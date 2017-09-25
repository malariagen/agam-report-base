# -*- coding: utf-8 -*-


# standard library imports
import os
import sys
from typing import Mapping, List
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
from .config import YAMLConfigFile
from .util import Logger, read_dataframe, guess_callset_format, open_callset, close_callset, \
    hash_params


class PopulationAnalysis(object):

    def __init__(self,
                 path: str,
                 overwrite: bool=False,
                 logfile=sys.stdout):
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
        self._genome_assembly: pyfasta.Fasta = None
        self._genome_annotation: pandas.DataFrame = None
        self._sample_table = None

    def get_abs_path(self, path):
        # relative paths treated as relative to analysis directory
        if os.path.isabs(path):
            abs_path = path
        else:
            abs_path = os.path.abspath(os.path.join(self.path, path))
        return abs_path

    def set_genome_assembly(self,
                            path: str,
                            label: str=None,
                            description: str=None,
                            format: str='fasta'):
        """Set the reference genome assembly.

        Parameters
        ----------
        path
            Path to assembly file. If given as a relative path, must be relative to the analysis
            directory.
        label
            Name of the assembly.
        description
            Description of the assembly.
        format
            File format. Only 'fasta' currently supported.

        """

        # check format
        format = format.lower()
        if format != 'fasta':
            raise NotImplementedError("Format %r not currently supported." % format)

        # check valid fasta
        try:
            abs_path = self.get_abs_path(path)
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
        config = dict(
            path=path, label=label, description=description, format=format
        )
        self.config.set_genome_assembly(config)

        # TODO need to revalidate callsets, accessibility, ..., if change genome assembly?

    @property
    def genome_assembly(self) -> pyfasta.Fasta:
        """TODO"""
        if self._genome_assembly is None:
            config = self.config.get_genome_assembly()
            path = self.get_abs_path(config['path'])
            self._genome_assembly = pyfasta.Fasta(path, key_fn=lambda v: v.split()[0])
        return self._genome_assembly

    def load_genome_assembly(self, chrom: str) -> np.ndarray:
        return np.asarray(self.genome_assembly[chrom])

    def set_genome_annotation(self,
                              path: str,
                              label: str=None,
                              description: str=None,
                              format: str='gff3'):
        """TODO"""

        # check format
        format = format.lower()
        if format != 'gff3':
            raise NotImplementedError("Format %r not currently supported." % format)

        # check valid gff3
        try:
            abs_path = self.get_abs_path(path)
            annotation = allel.gff3_to_dataframe(abs_path, attributes=['ID', 'Parent'])

        except Exception as e:
            raise ValueError('An error occurred opening the genome assembly file: {}\n'
                             'Does the path contain a valid fasta file?'
                             .format(e))

        # store for quick access
        self._genome_annotation = annotation

        # report some diagnostics
        n_features = len(annotation)
        n_genes = np.count_nonzero(annotation['type'] == 'gene')
        self.log('Setting genome annotation; found {0:,} features ({1:,} genes).'
                 .format(n_features, n_genes))

        # store configuration
        config = dict(
            path=path, label=label, description=description, format=format
        )
        self.config.set_genome_annotation(config)

    @property
    def genome_annotation(self):
        """TODO"""
        if self._genome_annotation is None:
            self._genome_annotation = self.load_genome_annotation(attributes=['ID', 'Parent'])
        return self._genome_annotation

    def load_genome_annotation(self, **kwargs) -> pandas.DataFrame:
        """TODO"""
        config = self.config.get_genome_annotation()
        path = self.get_abs_path(config['path'])
        annotation = allel.gff3_to_dataframe(path, **kwargs)
        return annotation

    def set_genome_accessibility(self,
                                 path: str,
                                 format: str='hdf5',
                                 dataset_name: str='is_accessible',
                                 label: str=None,
                                 description: str=None):
        """TODO"""

        # check format
        format = format.lower()
        if format != 'hdf5':
            raise NotImplementedError("Format %r not currently supported." % format)

        # check valid data
        abs_path = self.get_abs_path(path)
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

            for chrom in sorted(self.genome_assembly):

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
        config = dict(
            path=path, format=format, dataset_name=dataset_name, label=label,
            description=description
        )
        self.config.set_genome_accessibility(config)

    def load_genome_accessibility(self, chrom: str) -> np.ndarray:
        """TODO"""
        config = self.config.get_genome_accessibility()
        path = self.get_abs_path(config['path'])
        with h5py.File(path, mode='r') as accessibility:
            is_accessible = accessibility[chrom][config['dataset_name']][:]
        return is_accessible

    def set_sample_table(self,
                         path: str,
                         index_col: str,
                         name: str='main',
                         format: str='csv',
                         read_kws: dict=None,
                         label: str=None,
                         description: str=None):
        """TODO"""

        # check args
        if read_kws is None:
            read_kws = dict()
        else:
            read_kws = dict(read_kws)
        # force index_col into read keywords
        read_kws['index_col'] = index_col
        # help for tab-delimited formats
        format = format.lower()
        if format in {'tsv', 'tab', 'txt', 'text'}:
            format = 'csv'
            read_kws.setdefault('sep', '\t')

        # check data can be read by pandas
        abs_path = self.get_abs_path(path)
        df = read_dataframe(path=abs_path, format=format, read_kws=read_kws)

        # log some diagnostics
        self.log('Setting sample table {!r}; found {:,} rows.'
                 .format(name, len(df)))

        # store as attribute for convenience
        if name == 'main':
            self._sample_table = df

        # store configuration
        config = dict(
            path=path, format=format, read_kws=read_kws, label=label, description=description
        )
        self.config.set_sample_table(name, config)

    @property
    def sample_table(self):
        if self._sample_table is None:
            self._sample_table = self.load_sample_table()
        return self._sample_table

    def load_sample_table(self, name: str='main'):
        """TODO"""
        config = self.config.get_sample_table(name)
        path = self.get_abs_path(config['path'])
        read_kws = config['read_kws']
        format = config['format']
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
            if 'main' in names:
                first = 'main'
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

    def set_callset(self, path, format=None, name='main', phased=False, label=None,
                    description=None):
        """TODO"""

        # check format
        abs_path = self.get_abs_path(path)
        if format is None:
            format = guess_callset_format(abs_path)
        else:
            format = format.lower()
            if format not in {'hdf5', 'zarr'}:
                raise NotImplementedError("Format %r not currently supported." % format)
        # TODO auto extract vcf

        # check valid callset
        try:
            callset = open_callset(abs_path, format=format)

        except Exception as e:
            raise ValueError('An error occurred opening the callset: {}\n'
                             'Does the path contain a valid HDF5 file or Zarr directory?'
                             .format(e))

        try:

            n_chroms = 0
            n_variants = 0
            samples = None
            n_samples = 0
            if 'samples' in callset:
                samples = callset['samples'][:].tolist()
                n_samples = len(samples)

            for chrom in sorted(self.genome_assembly):
                if chrom in callset:
                    n_chroms += 1

                    # check variant positions
                    pos_path = 'variants/POS'
                    if pos_path not in callset[chrom]:
                        raise ValueError("Callset is missing POS for chromosome {!r}."
                                         .format(chrom))
                    pos = callset[chrom][pos_path]
                    n_variants += pos.shape[0]

                    # check samples
                    if 'samples' in callset[chrom] and samples is None:
                        samples = callset[chrom]['samples'][:].tolist()
                        n_samples = len(samples)

            if n_chroms == 0:
                raise ValueError('No chromosome groups found in callset.')

            if n_samples == 0:
                raise ValueError('No samples found in callset.')

        finally:
            close_callset(callset)

        # report some diagnostics
        self.log('Setting callset {!r}; found {:,} chromosomes, {:,} samples, {:,} variants.'
                 .format(name, n_chroms, n_samples, n_variants))

        # store configuration
        config = dict(
            path=path, format=format, phased=phased, label=label, description=description
        )
        self.config.set_callset(name, config)

    def set_population(self,
                       name: str,
                       query: str=None,
                       samples: List[str]=None,
                       color=None,
                       label: str=None,
                       description: str=None):
        """TODO"""

        if query and samples:
            raise ValueError('Please provide either a query or a list of samples, not both.')

        config = dict(color=color, label=label, description=description)
        if samples:
            samples = list(samples)
            n_samples = len(samples)
            config['samples'] = samples

        else:
            df = self.join_sample_tables()
            df_pop = df.query(query)
            n_samples = len(df_pop)
            if n_samples == 0:
                raise ValueError('Query does not match any samples.')
            config['query'] = query

        self.log('Setting population {!r}; found {:,} samples.'
                 .format(name, n_samples))
        self.config.set_population(name, config)

    def get_callset_samples(self, callset='main'):
        # TODO
        pass

    def get_population_samples(self, pop):
        # TODO
        pass

    def locate_samples(self, callset='main', pop=None, downsample=None, seed=None):

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

    def cache_save(self, group, key, result, names=None):
        # TODO
        pass

    def cache_load(self, group, key, names=None):
        # TODO
        pass

    def count_alleles(self, chrom, callset='main', pop=None, max_allele=3, dataset_name=None,
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
        cache_key = hash_params(params)

        try:
            result = self.cache_load(cache_group, cache_key)

        except CacheMiss:
            self.log('computing', cache_group, cache_key)

            root = self.open_callset(callset)

            try:
                # get genotype array
                gt = get_genotype_array(root, chrom, dataset_name)

                # run computation
                result = gt.count_alleles(max_allele=max_allele, subpop=sample_indices).compute()

            finally:
                close_callset(root)

            # save result
            self.cache_save(cache_group, cache_key, result)

        return allel.AlleleCountsArray(result)


class CacheMiss(Exception):
    pass
