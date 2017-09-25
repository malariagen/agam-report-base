# -*- coding: utf-8 -*-


# standard library imports
import os
import sys
import random
import inspect


# third party imports
import numpy as np
import pandas
import pyfasta
import allel
import h5py
import zarr


# internal imports
from .config import YAMLConfigFile, Key, Format
from .util import Logger, read_dataframe, guess_callset_format, open_callset, close_callset, \
    hash_params, get_genotype_array


class PopulationAnalysis(object):

    def __init__(self, path, overwrite=False, logfile=sys.stdout, cache_compressor=None):
        """TODO"""

        self.path = path

        # ensure directory exists
        if os.path.isfile(path):
            raise ValueError('Cannot create analysis directory, a file exists at path %r.' % path)
        os.makedirs(self.path, exist_ok=True)

        # setup config file
        config_path = os.path.join(self.path, 'config.yml')
        dump_kws = dict(default_flow_style=False)
        self.config = YAMLConfigFile(path=config_path, dump_kws=dump_kws, overwrite=overwrite)

        # setup logger
        self.log = Logger(logfile)

        # setup other instance attributes for lazy loading
        self._genome_assembly = None
        self._genome_annotation = None
        self._sample_table = None

        # setup cache compressor
        if cache_compressor is None:
            cache_compressor = zarr.Blosc(cname='zstd', clevel=1, shuffle=0)
        self.cache_compressor = cache_compressor

    def abs_path(self, path):
        # relative paths treated as relative to analysis directory
        if os.path.isabs(path):
            abs_path = path
        else:
            abs_path = os.path.abspath(os.path.join(self.path, path))
        return abs_path

    def set_genome_assembly(self,
                            path,
                            label=None,
                            description=None,
                            format=Format.FASTA):
        """Set the reference genome assembly.

        Parameters
        ----------
        path: str
            Path to assembly file. If given as a relative path, must be relative to the analysis
            directory.
        label: str
            Name of the assembly.
        description: str
            Description of the assembly.
        format: str
            File format. Only 'fasta' currently supported.

        """

        # check format
        format = str(format).lower()
        if format != Format.FASTA:
            raise NotImplementedError("Format %r not currently supported." % format)

        # check valid fasta
        try:
            abs_path = self.abs_path(path)
            assembly = pyfasta.Fasta(abs_path, key_fn=lambda v: v.split()[0])

        except Exception as e:
            raise ValueError('An error occurred opening the genome assembly file: {}\n'
                             'Does the path contain a valid fasta file?'
                             .format(e))

        # store for quick access
        self._genome_assembly = assembly

        # report some diagnostics
        n_chroms = len(assembly)
        n_bases = sum(len(assembly[chrom]) for chrom in assembly)
        self.log('Setting genome assembly; found {0} chromosomes, total size {1:,} bp.'
                 .format(n_chroms, n_bases))

        # store configuration
        config = {
            Key.PATH: path,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.FORMAT: format
        }
        self.config.set_genome_assembly(config)

        # TODO need to revalidate callsets, accessibility, ..., if change genome assembly?

    @property
    def genome_assembly(self):
        """TODO

        Returns
        -------
        pyfasta.Fasta

        """
        if self._genome_assembly is None:
            config = self.config.get_genome_assembly()
            path = self.abs_path(config[Key.PATH])
            self._genome_assembly = pyfasta.Fasta(path, key_fn=lambda v: v.split()[0])
        return self._genome_assembly

    def load_genome_assembly(self, chrom: str) -> np.ndarray:
        return np.asarray(self.genome_assembly[chrom])

    def set_genome_annotation(self,
                              path,
                              label=None,
                              description=None,
                              format=Format.GFF3):
        """TODO"""

        # check format
        format = str(format).lower()
        if format != Format.GFF3:
            raise NotImplementedError("Format %r not currently supported." % format)

        # check valid gff3
        try:
            abs_path = self.abs_path(path)
            annotation = allel.gff3_to_dataframe(abs_path, attributes=['ID', 'Parent'])

        except Exception as e:
            raise ValueError('An error occurred opening the genome assembly file: {}\n'
                             'Does the path contain a valid fasta file?'
                             .format(e))

        # store for quick access
        self._genome_annotation = annotation

        # report some diagnostics
        n_features = len(annotation)
        loc_gene = (annotation['type'] == 'gene')  # type: np.ndarray
        n_genes = np.count_nonzero(loc_gene)
        self.log('Setting genome annotation; found {0:,} features ({1:,} genes).'
                 .format(n_features, n_genes))

        # store configuration
        config = {
            Key.PATH: path,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.FORMAT: format,
        }
        self.config.set_genome_annotation(config)

    @property
    def genome_annotation(self):
        """TODO"""
        if self._genome_annotation is None:
            self._genome_annotation = self.load_genome_annotation(attributes=['ID', 'Parent'])
        return self._genome_annotation

    def load_genome_annotation(self, **kwargs):
        """TODO"""
        config = self.config.get_genome_annotation()
        path = self.abs_path(config[Key.PATH])
        kwargs.setdefault('attributes', ['ID', 'Parent'])
        annotation = allel.gff3_to_dataframe(path, **kwargs)
        return annotation

    def set_genome_accessibility(self, path, format=Format.HDF5, dataset_name='is_accessible',
                                 label=None, description=None):
        """TODO"""

        # check format
        format = str(format).lower()
        if format != Format.HDF5:
            raise NotImplementedError("Format %r not currently supported." % format)

        # check valid data
        abs_path = self.abs_path(path)
        try:
            accessibility = h5py.File(abs_path, mode='r')

        except Exception as e:
            raise ValueError('Error opening accessibility file: {}\n'
                             'Does the path contain a valid HDF5 file?'
                             .format(e))

        try:

            # check file validity
            n_chroms = 0
            n_bases_accessible = 0

            chroms = sorted(self.genome_assembly)
            for chrom in chroms:

                if chrom in accessibility:
                    group = accessibility[chrom]
                    n_chroms += 1
                else:
                    continue

                # check dataset exists
                if dataset_name not in group:
                    raise ValueError('Could not find accessibility dataset for chromosome '
                                     '{!r}.'.format(chrom))

                # check dataset type and shape
                is_accessible = group[dataset_name]
                if is_accessible.dtype != bool:
                    raise ValueError('Expected bool dtype, found {!r}'
                                     .format(is_accessible.dtype))
                if is_accessible.ndim != 1:
                    raise ValueError('Expected dataset with 1 dimension, found {}'
                                     .format(is_accessible.ndim))
                n_bases = len(self.genome_assembly[chrom])
                if is_accessible.shape[0] != n_bases:
                    raise ValueError('Expected dataset with length {}, found {}.'
                                     .format(n_bases, is_accessible.shape[0]))

                n_bases_accessible += np.count_nonzero(is_accessible[:])

            if n_chroms == 0:
                raise ValueError('No chromosome groups found in accessibility data.')

        finally:
            accessibility.close()

        # log some diagnostics
        self.log('Setting genome accessibility; found {:,} chromosomes, {:,} accessible bp.'
                 .format(n_chroms, n_bases_accessible))

        # store configuration
        config = {
            Key.PATH: path,
            Key.FORMAT: format,
            Key.DATASET_NAME: dataset_name,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
        }
        self.config.set_genome_accessibility(config)

    def load_genome_accessibility(self, chrom):
        """TODO"""
        conf = self.config.get_genome_accessibility()
        path = self.abs_path(conf[Key.PATH])
        dataset_name = conf[Key.DATASET_NAME]
        with h5py.File(path, mode='r') as accessibility:
            is_accessible = accessibility[chrom][dataset_name][:]
        return is_accessible

    def set_sample_table(self,
                         path,
                         index_col,
                         name=Key.MAIN,
                         format=Format.CSV,
                         read_kws=None,
                         label=None,
                         description=None):
        """TODO"""

        # check args
        if read_kws is None:
            read_kws = dict()
        else:
            read_kws = dict(read_kws)

        # force index_col into read keywords
        read_kws['index_col'] = index_col

        # help for tab-delimited formats
        format = str(format).lower()
        if format in {'tsv', 'tab', 'txt', 'text'}:
            format = Format.CSV
            read_kws.setdefault('sep', '\t')

        # check data can be read by pandas
        abs_path = self.abs_path(path)
        df = read_dataframe(path=abs_path, format=format, read_kws=read_kws)

        # log some diagnostics
        self.log('Setting sample table {!r}; found {:,} columns, {:,} rows.'
                 .format(name, len(df.columns), len(df)))

        # store as attribute for convenience
        if name == Key.MAIN:
            self._sample_table = df

        # store configuration
        config = {
            Key.PATH: path,
            Key.FORMAT: format,
            Key.READ_KWS: read_kws,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
        }
        self.config.set_sample_table(name, config)

    @property
    def sample_table(self):
        if self._sample_table is None:
            self._sample_table = self.load_sample_table()
        return self._sample_table

    def load_sample_table(self, name=Key.MAIN):
        """TODO"""
        config = self.config.get_sample_table(name)
        path = self.abs_path(config[Key.PATH])
        read_kws = config[Key.READ_KWS]
        format = config[Key.FORMAT]
        df = read_dataframe(path=path, format=format, read_kws=read_kws)
        return df

    def join_sample_tables(self, names=None, on=None, how='outer'):
        """TODO"""

        # find all sample table names
        if names:
            first = names[0]
            others = names[1:]
        else:
            names = sorted(self.config.get_sample_tables())
            if Key.MAIN in names:
                first = Key.MAIN
                others = [n for n in names if n != first]
            else:
                first = names[0]
                others = names[1:]

        # load and join
        df = self.load_sample_table(first)
        if others:
            dfs_others = [self.load_sample_table(n) for n in others]
            df = df.join(dfs_others, on=on, how=how)

        return df

    def set_callset(self, path, format=None, name=Key.MAIN, phased=False, label=None,
                    description=None, samples_path=None):
        """TODO"""

        # check format
        abs_path = self.abs_path(path)
        if format is None:
            format = guess_callset_format(abs_path)
        else:
            format = str(format).lower()
            if format not in {Format.HDF5, Format.ZARR}:
                raise NotImplementedError("Format %r not currently supported." % format)
        # TODO auto extract vcf

        # check valid callset
        try:
            root = open_callset(abs_path, format=format)

        except Exception as e:
            raise ValueError('An error occurred opening the callset: {}\n'
                             'Does the path contain a valid HDF5 file or Zarr directory?'
                             .format(e))

        try:

            n_chroms = 0
            n_variants = 0

            if samples_path is None and 'samples' in root:
                samples_path = 'samples'

            for chrom in sorted(self.genome_assembly):
                if chrom in root:
                    n_chroms += 1

                    # check variant positions
                    pos_path = 'variants/POS'
                    if pos_path not in root[chrom]:
                        raise ValueError("Callset is missing variants/POS for chromosome {!r}."
                                         .format(chrom))
                    pos = root[chrom][pos_path]
                    n_variants += pos.shape[0]

                    # check samples
                    if samples_path is None and 'samples' in root[chrom]:
                        samples_path = '/'.join([chrom, 'samples'])

            if n_chroms == 0:
                raise ValueError('No chromosome groups found in callset.')

            # check samples
            if samples_path:
                samples = root[samples_path][:].tolist()
                n_samples = len(samples)
            else:
                raise ValueError('No samples found in callset.')

        finally:
            close_callset(root)

        # report some diagnostics
        self.log('Setting callset {!r}; found {:,} chromosomes, {:,} samples, {:,} variants.'
                 .format(name, n_chroms, n_samples, n_variants))

        # store configuration
        config = {
            Key.PATH: path,
            Key.FORMAT: format,
            Key.PHASED: phased,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.SAMPLES_PATH: samples_path
        }
        self.config.set_callset(name, config)

    def set_population(self,
                       name,
                       query=None,
                       samples=None,
                       color=None,
                       label=None,
                       description=None):
        """TODO"""

        if query and samples:
            raise ValueError('Please provide either a query or a list of samples, not both.')

        config = {
            Key.COLOR: color,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
        }
        if samples:
            samples = list(samples)
            n_samples = len(samples)
            config[Key.SAMPLES] = samples

        else:
            df = self.join_sample_tables()
            df_pop = df.query(query)
            n_samples = len(df_pop)
            if n_samples == 0:
                raise ValueError('Query does not match any samples.')
            config[Key.QUERY] = query

        self.log('Setting population {!r}; found {:,} samples.'
                 .format(name, n_samples))
        self.config.set_population(name, config)

    def open_callset(self, callset=Key.MAIN):
        """TODO"""
        config = self.config.get_callset(callset)
        path = self.abs_path(config[Key.PATH])
        format = config[Key.FORMAT]
        if format == Format.HDF5:
            root = h5py.File(path, mode='r')
        elif format == Format.ZARR:
            root = zarr.open_group(path, mode='r')
        else:
            raise RuntimeError('Unexpected format: {!r}.'.format(format))
        return root

    def get_callset_samples(self, callset=Key.MAIN):
        """TODO"""
        config = self.config.get_callset(callset)
        samples_path = config[Key.SAMPLES_PATH]
        root = self.open_callset(callset)
        try:
            samples = root[samples_path][:]
            if samples.dtype.kind == 'S':
                samples = samples.astype('U')
            samples = samples.tolist()
        finally:
            close_callset(root)
        return samples

    def get_population_samples(self, pop):
        """TODO"""
        config = self.config.get_population(pop)
        query = config.get(Key.QUERY, None)
        if query:
            df = self.join_sample_tables()
            df_pop = df.query(query)
            samples = df_pop.index.values.tolist()
        else:
            samples = config[Key.SAMPLES]
        return samples

    def locate_samples(self, callset=Key.MAIN, pop=None, downsample=None, seed=None):
        """TODO"""

        callset_samples = self.get_callset_samples(callset)
        if pop:
            samples = self.get_population_samples(pop)
        else:
            samples = callset_samples

        # deal with downsampling
        if downsample:
            if downsample > len(samples):
                raise ValueError('Not enough samples ({}) to downsample to {}.'
                                 .format(len(samples), downsample))
            random.seed(seed)
            samples = sorted(random.sample(samples, downsample))

        # find sample indices
        if samples == callset_samples:
            # no selection, use None to indicate all samples
            samples = None
            sample_indices = None

        else:
            sample_indices = list()
            for s in samples:
                if s not in callset_samples:
                    raise ValueError('Sample {!r} not found in callset {!r}'
                                     .format(s, callset))
                sample_indices.append(callset_samples.index(s))

        return samples, sample_indices

    def cache_save(self, cache_group, cache_key, params_doc, data, names=None):
        # self.log('saving', cache_group, cache_key)
        cache = zarr.open_group(os.path.join(self.path, 'cache'), mode='a')
        group = cache.require_group(cache_group)
        if names is None:
            if cache_key in group:
                del group[cache_key]
            result = group.create_dataset(cache_key, dtype=data.dtype, shape=data.shape,
                                          compressor=self.cache_compressor)
            result[:] = data
        else:
            if cache_key in group:
                del group[cache_key]
            result = group.create_group(cache_key)
            for n, d in zip(names, data):
                ds = result.create_dataset(n, dtype=d.dtype, shape=d.shape,
                                           compressor=self.cache_compressor)
                ds[:] = d
        # mark success
        result.attrs['__params__'] = params_doc

    def cache_load(self, cache_group, cache_key, names=None):
        cache = zarr.open_group(os.path.join(self.path, 'cache'), mode='a')
        group = cache.require_group(cache_group)
        if cache_key not in group:
            self.log('computing', cache_group, cache_key)
            raise CacheMiss
        result = group[cache_key]
        if '__params__' not in result.attrs:
            self.log('computing', cache_group, cache_key)
            raise CacheMiss
        self.log('loading', cache_group, cache_key)
        if names is None:
            return result[:]
        else:
            return tuple(result[n][:] for n in names)

    def count_alleles(self, chrom, callset=Key.MAIN, pop=None, max_allele=3, dataset_name=None,
                      downsample=None, seed=None):
        """TODO"""

        # determine which samples to include
        samples, sample_indices = self.locate_samples(callset=callset, pop=pop,
                                                      downsample=downsample, seed=seed)

        # setup parameters
        params = dict(
            chrom=chrom,
            # N.B., include full callset config here in case it gets changed
            callset=self.config.get_callset(callset),
            samples=samples,
            max_allele=max_allele,
            dataset_name=dataset_name,
        )

        # check cache
        cache_group = inspect.currentframe().f_code.co_name
        params_doc, cache_key = hash_params(params)

        try:
            result = self.cache_load(cache_group, cache_key)

        except CacheMiss:

            root = self.open_callset(callset)

            try:
                # get genotype array
                gt = get_genotype_array(root, chrom, dataset_name)

                # run computation
                result = gt.count_alleles(max_allele=max_allele, subpop=sample_indices).compute()

            finally:
                close_callset(root)

            # save result
            self.cache_save(cache_group, cache_key, params_doc, result)

        return allel.AlleleCountsArray(result)

    # TODO __repr__

    def windowed_diversity(self, chrom, window_size, window_start=None, window_stop=None,
                           window_step=None, pop=None, callset=Key.MAIN, downsample=None,
                           seed=None, max_allele=3, dataset_name=None, equally_accessible=True):
        """TODO"""

        # determine which samples to include
        samples, sample_indices = self.locate_samples(callset=callset, pop=pop,
                                                      downsample=downsample, seed=seed)

        # setup parameters
        params = dict(
            chrom=chrom,
            window_size=window_size,
            window_start=window_start,
            window_stop=window_stop,
            window_step=window_step,
            # N.B., include full callset config here in case it gets changed
            callset=self.config.get_callset(callset),
            samples=samples,
            max_allele=max_allele,
            dataset_name=dataset_name,
            equally_accessible=equally_accessible,
        )

        # check cache
        cache_group = inspect.currentframe().f_code.co_name
        params_doc, cache_key = hash_params(params)

        try:
            result = self.cache_load(cache_group, cache_key, names=['windows', 'pi'])

        except CacheMiss:

            # load allele counts
            ac = self.count_alleles(chrom=chrom, callset=callset, pop=pop, downsample=downsample,
                                    seed=seed, max_allele=max_allele, dataset_name=dataset_name)

            # TODO load variant positions

            # TODO setup windows

            # TODO run windowed computation

            # TODO save result

        # TODO return result


class CacheMiss(Exception):
    pass
