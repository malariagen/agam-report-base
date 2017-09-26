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
from .config import YAMLConfigFile, Key, Format, Section
from .util import Logger, read_dataframe, guess_callset_format, open_callset, close_callset, \
    hash_params


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
            Human-readable label.
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
        conf = {
            Key.PATH: path,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.FORMAT: format
        }
        self.config.set_genome_assembly(conf)

        # TODO need to revalidate callsets, accessibility, ..., if change genome assembly?

    @property
    def genome_assembly(self):
        """TODO

        Returns
        -------
        pyfasta.Fasta

        """
        if self._genome_assembly is None:
            conf = self.config.get_genome_assembly()
            path = self.abs_path(conf[Key.PATH])
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
        conf = {
            Key.PATH: path,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.FORMAT: format,
        }
        self.config.set_genome_annotation(conf)

    @property
    def genome_annotation(self):
        """TODO"""
        if self._genome_annotation is None:
            self._genome_annotation = self.load_genome_annotation(attributes=['ID', 'Parent'])
        return self._genome_annotation

    def load_genome_annotation(self, **kwargs):
        """TODO"""
        conf = self.config.get_genome_annotation()
        path = self.abs_path(conf[Key.PATH])
        kwargs.setdefault('attributes', ['ID', 'Parent'])
        annotation = allel.gff3_to_dataframe(path, **kwargs)
        return annotation

    def set_genome_accessibility(self, path, format=Format.HDF5, dataset='is_accessible',
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
                if dataset not in group:
                    raise ValueError('Could not find accessibility dataset for chromosome '
                                     '{!r}.'.format(chrom))

                # check dataset type and shape
                is_accessible = group[dataset]
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
        conf = {
            Key.PATH: path,
            Key.FORMAT: format,
            Key.DATASET: dataset,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
        }
        self.config.set_genome_accessibility(conf)

    def load_genome_accessibility(self, chrom):
        """TODO"""
        conf = self.config.get_genome_accessibility()
        path = self.abs_path(conf[Key.PATH])
        dataset = conf[Key.DATASET]
        with h5py.File(path, mode='r') as accessibility:
            is_accessible = accessibility[chrom][dataset][:]
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

        # TODO use name as label if label is None

        # store configuration
        conf = {
            Key.PATH: path,
            Key.FORMAT: format,
            Key.READ_KWS: read_kws,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
        }
        self.config.set_sample_table(name, conf)

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
                    description=None, samples_path=None, genotype_dataset=None, max_allele=3):
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
            callset_root = open_callset(abs_path, format=format)

        except Exception as e:
            raise ValueError('An error occurred opening the callset: {}\n'
                             'Does the path contain a valid HDF5 file or Zarr directory?'
                             .format(e))

        try:

            n_chroms = 0
            n_variants = 0

            if samples_path is None and 'samples' in callset_root:
                samples_path = 'samples'

            for chrom in sorted(self.genome_assembly):
                if chrom in callset_root:
                    n_chroms += 1

                    # check variants and calldata groups
                    if 'variants' not in callset_root[chrom]:
                        raise ValueError("Callset is missing 'variants' group for chromosome {!r}."
                                         .format(chrom))
                    variants = callset_root[chrom]['variants']
                    if 'calldata' not in callset_root[chrom]:
                        raise ValueError("Callset is missing 'calldata' group for chromosome {!r}."
                                         .format(chrom))
                    calldata = callset_root[chrom]['calldata']

                    # check variant positions
                    pos_dataset = 'POS'
                    if pos_dataset not in variants:
                        raise ValueError("Callset is missing 'variants/POS' for chromosome {!r}."
                                         .format(chrom))
                    pos = variants[pos_dataset]
                    n_variants += pos.shape[0]

                    # check samples
                    if samples_path is None and 'samples' in callset_root[chrom]:
                        samples_path = '/'.join([chrom, 'samples'])

                    # check genotypes
                    if genotype_dataset is None:
                        for k in 'genotype', 'GT':
                            if k in calldata:
                                genotype_dataset = k
                                break
                    else:
                        if genotype_dataset not in calldata:
                            raise ValueError("Callset is missing 'calldata/{}' for chromosome {!r}."
                                             .format(genotype_dataset, chrom))

            if n_chroms == 0:
                raise ValueError('No chromosome groups found in callset.')

            # check samples
            if samples_path:
                samples = callset_root[samples_path][:].tolist()
                n_samples = len(samples)
            else:
                raise ValueError('No samples found in callset.')

            # check genotypes
            if genotype_dataset is None:
                raise ValueError('No genotypes found in callset.')

        finally:
            close_callset(callset_root)

        # report some diagnostics
        self.log('Setting callset {!r}; found {:,} chromosomes, {:,} samples, {:,} variants.'
                 .format(name, n_chroms, n_samples, n_variants))

        # TODO use name as label if label is None

        # store configuration
        conf = {
            Key.PATH: path,
            Key.FORMAT: format,
            Key.PHASED: phased,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.SAMPLES_PATH: samples_path,
            Key.GENOTYPE_DATASET: genotype_dataset,
            Key.MAX_ALLELE: max_allele,
        }
        self.config.set_callset(name, conf)

    def set_population(self,
                       name,
                       query=None,
                       samples=None,
                       color=None,
                       label=None,
                       description=None,
                       down_sample=None,
                       random_seed=None):
        """TODO"""

        if query and samples:
            raise ValueError('Please provide either a query or a list of samples, not both.')

        # TODO use name as label if label is None

        # common configuration
        conf = {
            Key.COLOR: color,
            Key.LABEL: label,
            Key.DESCRIPTION: description,
            Key.DOWN_SAMPLE: down_sample,
            Key.RANDOM_SEED: random_seed,
        }

        if samples:
            samples = list(samples)
            n_samples = len(samples)
            conf[Key.SAMPLES] = samples

        else:
            df = self.join_sample_tables()
            df_pop = df.query(query)
            n_samples = len(df_pop)
            if n_samples == 0:
                raise ValueError('Query does not match any samples.')
            conf[Key.QUERY] = query

        self.log('Setting population {!r}; found {:,} samples.'
                 .format(name, n_samples))
        self.config.set_population(name, conf)

    def get_population_samples(self, pop):
        """TODO"""

        # get population configuration
        conf = self.config.get_population(pop)
        query = conf.get(Key.QUERY, None)

        # handle query
        if query:
            df = self.join_sample_tables()
            df_pop = df.query(query)
            samples = df_pop.index.values.tolist()

        # handle explicit samples list
        else:
            samples = conf[Key.SAMPLES]

        # handle down-sampling
        down_sample = conf.get(Key.DOWN_SAMPLE, None)
        if down_sample:
            if down_sample > len(samples):
                raise ValueError('Not enough samples in population {!r} ({}) to downsample to {}.'
                                 .format(pop, len(samples), down_sample))
            random_seed = conf.get(Key.RANDOM_SEED, None)
            random.seed(random_seed)
            samples = sorted(random.sample(samples, down_sample))

        return samples

    def open_callset(self, callset=Key.MAIN):
        """TODO"""
        config = self.config.get_callset(callset)
        path = self.abs_path(config[Key.PATH])
        format = config[Key.FORMAT]
        if format == Format.HDF5:
            callset_root = h5py.File(path, mode='r')
        elif format == Format.ZARR:
            callset_root = zarr.open_group(path, mode='r')
        else:
            raise RuntimeError('Unexpected format: {!r}.'.format(format))
        return callset_root

    def get_callset_samples(self, callset=Key.MAIN):
        """TODO"""
        config = self.config.get_callset(callset)
        samples_path = config[Key.SAMPLES_PATH]
        callset_root = self.open_callset(callset)
        try:
            samples = callset_root[samples_path][:]
            if samples.dtype.kind == 'S':
                samples = samples.astype('U')
            samples = samples.tolist()
        finally:
            close_callset(callset_root)
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
        # self.log('saving {}/{}'.format(cache_group, cache_key))
        cache = zarr.open_group(os.path.join(self.path, 'cache'), mode='a')
        group = cache.require_group(cache_group)

        if names is None:
            # result is single array result
            if cache_key in group:
                del group[cache_key]
            result = group.create_dataset(cache_key, dtype=data.dtype, shape=data.shape,
                                          compressor=self.cache_compressor)
            result[:] = data

        else:
            # result is tuple of arrays
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

        # look for result in cache
        try:
            if cache_key not in group:
                raise CacheMiss
            result = group[cache_key]
            if '__params__' not in result.attrs:
                # incomplete save
                raise CacheMiss
        except CacheMiss:
            self.log('computing {}/{}'.format(cache_group, cache_key))
            raise

        # self.log('loading {}/{}'.format(cache_group, cache_key))
        if names is None:
            return result[:]
        else:
            return tuple(result[n][:] for n in names)

    def count_alleles(self, chrom, callset=Key.MAIN, pop=None):
        """TODO"""

        # determine which samples to include
        samples, sample_indices = self.locate_samples(callset=callset, pop=pop)

        # setup parameters
        callset_conf = self.config.get_callset(callset)
        params = dict(
            chrom=chrom,
            # N.B., include full callset config here in case it gets changed
            callset=callset_conf,
            samples=samples,
        )

        # check cache
        cache_group = inspect.currentframe().f_code.co_name
        params_doc, cache_key = hash_params(params)

        try:
            result = self.cache_load(cache_group, cache_key)

        except CacheMiss:

            genotype_dataset = callset_conf[Key.GENOTYPE_DATASET]
            max_allele = callset_conf[Key.MAX_ALLELE]
            callset_root = self.open_callset(callset)

            try:
                # get genotype array
                gt = allel.GenotypeDaskArray(callset_root[chrom]['calldata'][genotype_dataset])

                # run computation
                result = gt.count_alleles(max_allele=max_allele, subpop=sample_indices).compute()

            finally:
                close_callset(callset_root)

            # save result
            self.cache_save(cache_group, cache_key, params_doc, result)

        return allel.AlleleCountsArray(result)

    # TODO __repr__

    def load_variant_index(self, chrom, callset=Key.MAIN):
        """TODO"""
        callset_root = self.open_callset(callset)
        try:
            pos = callset_root[chrom]['variants']['POS'][:]
            pos = allel.SortedIndex(pos)
        finally:
            close_callset(callset_root)
        return pos

    def windowed_diversity(self, chrom, window_size,
                           # window_start=None, window_stop=None, window_step=None,
                           pop=None, callset=Key.MAIN, eqaccess=True):
        """TODO"""

        # handle multiple chromosomes
        if isinstance(chrom, (list, tuple)):
            results = [self.windowed_diversity(chrom=c, window_size=window_size, pop=pop,
                                               callset=callset, eqaccess=eqaccess)
                       for c in chrom]
            windows = np.concatenate([result[0] for result in results])
            pi = np.concatenate([result[1] for result in results])
            return windows, pi

        # determine which samples to include
        samples, sample_indices = self.locate_samples(callset=callset, pop=pop)

        # setup parameters
        params = dict(
            chrom=chrom,
            window_size=window_size,
            # window_start=window_start,
            # window_stop=window_stop,
            # window_step=window_step,
            # N.B., include full callset config here in case it gets changed
            callset=self.config.get_callset(callset),
            samples=samples,
            eqaccess=eqaccess,
        )

        # check cache
        cache_group = inspect.currentframe().f_code.co_name
        params_doc, cache_key = hash_params(params)
        result_names = 'windows', 'pi'

        try:

            # load from cache
            result = self.cache_load(cache_group, cache_key, names=result_names)

        except CacheMiss:

            # load allele counts
            ac = self.count_alleles(chrom=chrom, callset=callset, pop=pop)

            # load variant positions
            pos = self.load_variant_index(chrom=chrom, callset=callset)

            # load accessibility
            is_accessible = self.load_genome_accessibility(chrom)

            # setup windows
            # TODO add support for start, stop and step
            if eqaccess:
                windows = allel.equally_accessible_windows(is_accessible=is_accessible,
                                                           size=window_size)

            else:
                chrom_size = len(self.genome_assembly[chrom])
                windows = allel.stats.window.position_windows(pos, size=window_size, start=1,
                                                              stop=chrom_size)

            # run windowed computation
            pi, _, _, _ = allel.windowed_diversity(pos, ac, windows=windows,
                                                   is_accessible=is_accessible)

            # save result
            result = windows, pi
            self.cache_save(cache_group=cache_group, cache_key=cache_key, data=result,
                            names=result_names, params_doc=params_doc)

        return result

    def windowed_diversity_distplot(self, chrom, window_size,
                                    # window_start=None, window_stop=None, window_step=None,
                                    pop=None, callset=Key.MAIN, eqaccess=True,
                                    ax=None, plot_kws=None, average=np.median, ci_kws=None,
                                    xlim=None):

        # load data
        _, pi = self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop,
                                        callset=callset, eqaccess=eqaccess)
        if isinstance(chrom, (list, tuple)):
            chrom_label = 'Chromosomes: ' + ', '.join(chrom)
        else:
            chrom_label = 'Chromosome: ' + chrom

        # plot as percent
        pi = pi * 100

        # main plot
        import seaborn as sns
        if plot_kws is None:
            plot_kws = dict()
        plot_kws.setdefault('ax', ax)
        ax = sns.distplot(pi, **plot_kws)

        # tidy up
        ax.set_xlim(left=0)
        ax.set_xlabel(r'$\theta_{\pi}$ (%)')
        ax.set_ylabel('Frequency')
        pop_label = self.get_population_label(pop)
        ax.set_title('Population: {}; {}'
                     .format(pop_label, chrom_label))
        if xlim:
            ax.set_xlim(*xlim)

        # annotate average
        if average is not None:
            import scikits.bootstrap as bootstrap
            if ci_kws is None:
                ci_kws = dict()
            ci_kws['statfunction'] = average
            ci = bootstrap.ci(pi, **ci_kws)
            m = average(pi)
            ax.text(0.02, 0.98,
                    '{}={:.3f}%\n95% CI [{:.3f}-{:.3f}]'
                    .format(average.__name__, m, ci[0], ci[1]),
                    ha='left', va='top', transform=ax.transAxes)

    def get_population_label(self, pop):
        """TODO"""
        label = self.config.get_population(pop)[Key.LABEL]
        if label is None:
            label = pop
        n_samples = len(self.get_population_samples(pop))
        label += ' (n={})'.format(n_samples)
        return label

    def windowed_diversity_violinplot(self, chrom, window_size, pops=None,
                                      # window_start=None, window_stop=None, window_step=None,
                                      callset=Key.MAIN, eqaccess=True,
                                      ax=None, plot_kws=None, ylim=None):
        """TODO"""

        # load data
        if isinstance(pops, str):
            pops = [pops]
        results = [
            self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop, callset=callset,
                                    eqaccess=eqaccess)
            for pop in pops
        ]

        # plot as percent
        data = [pi * 100 for _, pi in results]

        # handle multiple chroms
        if isinstance(chrom, (list, tuple)):
            chrom_label = 'Chromosomes: ' + ', '.join(chrom)
        else:
            chrom_label = 'Chromosome: ' + chrom

        # main plot
        import seaborn as sns
        if plot_kws is None:
            plot_kws = dict()
        plot_kws['ax'] = ax
        ax = sns.violinplot(data=data, **plot_kws)

        # tidy up
        pop_labels = [self.get_population_label(pop) for pop in pops]
        ax.set_xticklabels(pop_labels)
        ax.set_xlabel('Population')
        ax.set_ylabel(r'$\theta_{\pi}$ (%)')
        if ylim:
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(bottom=0)
        ax.set_title(chrom_label)

    def windowed_diversity_average_ci(self, chrom, window_size, pop=None,
                                      # window_start=None, window_stop=None, window_step=None,
                                      callset=Key.MAIN, eqaccess=True,
                                      average=np.median, ci_kws=None):
        """TODO"""

        # load data
        _, pi = self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop,
                                        callset=callset, eqaccess=eqaccess)

        # compute CI
        import scikits.bootstrap as bootstrap
        if ci_kws is None:
            ci_kws = dict()
        ci_kws['statfunction'] = average
        ci = bootstrap.ci(pi, **ci_kws)
        m = average(pi)

        return m, list(ci)

    def windowed_diversity_compare(self, chrom, window_size, pops,
                                   # window_start=None, window_stop=None, window_step=None,
                                   callset=Key.MAIN, eqaccess=True,
                                   average=np.median, ci_kws=None):

        # load data
        results = {pop: self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop,
                                                callset=callset,
                                                eqaccess=eqaccess)
                   for pop in pops}

        # setup printing
        pop_labels = [self.get_population_label(pop) for pop in pops]
        label_len = max(len(l) for l in pop_labels)

        # handle multiple chroms
        if isinstance(chrom, (list, tuple)):
            chrom_label = 'chromosomes ' + ', '.join(chrom)
        else:
            chrom_label = 'Chromosome ' + chrom

        title = 'Nucleotide diversity; {}\n'.format(chrom_label)
        title += '-' * (len(title) - 1)
        self.log(title)

        # compute averages and CIs
        import scikits.bootstrap as bootstrap
        if ci_kws is None:
            ci_kws = dict()
        ci_kws['statfunction'] = average
        for pop in pops:
            # compute as percent
            result = results[pop]
            pi = result[1] * 100
            ci = bootstrap.ci(pi, **ci_kws)
            m = average(pi)
            label = self.get_population_label(pop)
            self.log('{} : {}={:.3f}%; 95% CI [{:.3f}-{:.3f}]'
                     .format(label.ljust(label_len), average.__name__, m, ci[0], ci[1]))

        # compare pairwise
        import itertools
        import scipy.stats
        for pop1, pop2 in itertools.combinations(pops, 2):
            pi1 = results[pop1][1]
            pi2 = results[pop2][1]
            lbl1 = self.get_population_label(pop1)
            lbl2 = self.get_population_label(pop2)
            t, p = scipy.stats.wilcoxon(pi1, pi2)
            p_fmt = '{:.2e}' if p < 1e-3 else '{:.3f}'
            p_str = p_fmt.format(p)
            self.log('{} versus {} : Wilcoxon signed rank test P={}; statistic={:.1f}'
                     .format(lbl1, lbl2, p_str, t))

    def windowed_diversity_chromplot(self, chrom, window_size, pop=None,
                                     callset=Key.MAIN, eqaccess=True,
                                     ax=None, plot_kws=None, legend=True, legend_kws=None):

        # load data
        windows, pi = self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop,
                                              callset=callset,
                                              eqaccess=eqaccess)
        # plot as percent
        pi = pi * 100

        # main plot
        import matplotlib.pyplot as plt
        import seaborn as sns
        if plot_kws is None:
            plot_kws = dict()
        plot_kws.setdefault('linewidth', 1)
        pop_label = self.get_population_label(pop)
        plot_kws.setdefault('label', pop_label)
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 2))
        x = windows.mean(axis=1)
        ax.plot(x, pi, **plot_kws)

        # tidy
        chrom_size = len(self.genome_assembly[chrom])
        ax.set_xlim(0, chrom_size)
        ax.set_ylim(bottom=0)
        ax.set_ylabel(r'$\theta_{\pi}$ (%)')
        xticks = ax.get_xticks()
        if chrom_size > 1e6:
            xticklabels = xticks / 1e6
            units = 'Mbp'
        elif chrom_size > 1e5:
            xticklabels = xticks / 1e5
            units = 'kbp'
        else:
            xticklabels = xticks
            units = 'bp'
        ax.set_xticklabels(xticklabels)
        ax.set_title(chrom)
        ax.set_xlabel('Position ({})'.format(units))
        if legend:
            if legend_kws is None:
                legend_kws = dict()
            legend_kws.setdefault('loc', 'upper left')
            legend_kws.setdefault('bbox_to_anchor', (1, 1))
            legend_kws.setdefault('title', 'Population')
            ax.legend(**legend_kws)

    def windowed_diversity_delta_chromplot(self, chrom, window_size, pop1, pop2,
                                           callset=Key.MAIN, eqaccess=True,
                                           ax=None, plot_kws=None, legend=True, legend_kws=None):

        # load data
        windows, pi1 = self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop1,
                                               callset=callset,
                                               eqaccess=eqaccess)
        _, pi2 = self.windowed_diversity(chrom=chrom, window_size=window_size, pop=pop2,
                                         callset=callset,
                                         eqaccess=eqaccess)

        # plot as percent
        delta = (pi1 - pi2) * 100

        # main plot
        import matplotlib.pyplot as plt
        import seaborn as sns
        if plot_kws is None:
            plot_kws = dict()
        plot_kws.setdefault('linewidth', 0)
        pop1_label = self.get_population_label(pop1)
        pop2_label = self.get_population_label(pop2)
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 2))
        x = windows.mean(axis=1)
        ax.fill_between(x, delta, where=delta>0, label=pop1_label)
        ax.fill_between(x, delta, where=delta<0, label=pop2_label)

        # tidy
        chrom_size = len(self.genome_assembly[chrom])
        ax.set_xlim(0, chrom_size)
        ax.set_ylabel(r'$\Delta \theta_{\pi}$ (%)')
        xticks = ax.get_xticks()
        if chrom_size > 1e6:
            xticklabels = xticks / 1e6
            units = 'Mbp'
        elif chrom_size > 1e5:
            xticklabels = xticks / 1e5
            units = 'kbp'
        else:
            xticklabels = xticks
            units = 'bp'
        ax.set_xticklabels(xticklabels)
        ax.set_title(chrom)
        ax.set_xlabel('Position ({})'.format(units))
        if legend:
            if legend_kws is None:
                legend_kws = dict()
            legend_kws.setdefault('loc', 'upper left')
            legend_kws.setdefault('bbox_to_anchor', (1, 1))
            legend_kws.setdefault('title', 'Population')
            ax.legend(**legend_kws)

    def genomeplot(self, plotf, chroms, fig=None, gridspec_kws=None, **kwargs):
        """TODO"""

        from matplotlib.gridspec import GridSpec
        import matplotlib.pyplot as plt
        import seaborn as sns
        if fig is None:
            fig = plt.figure(figsize=(8, 2))

        # setup subplots
        if fig.axes:
            # assume axes already set up
            axs = fig.axes
            if len(axs) != len(chroms):
                raise ValueError('Figure has wrong number of axes; expected {}, found {}.'
                                 .format(len(chroms), len(axs)))
        else:
            # setup gridspec
            chrom_sizes = [len(self.genome_assembly[chrom]) for chrom in chroms]
            if gridspec_kws is None:
                gridspec_kws = dict()
            gridspec_kws.setdefault('width_ratios', chrom_sizes)
            gs = GridSpec(nrows=1, ncols=len(chroms), **gridspec_kws)
            axs = []
            for i in range(len(chroms)):
                if i == 0:
                    ax = fig.add_subplot(gs[i])
                else:
                    ax = fig.add_subplot(gs[i], sharey=axs[0])
                axs.append(ax)

        # do plotting
        # see https://stackoverflow.com/questions/22511550/gridspec-with-shared-axes-in-python
        legend = kwargs.pop('legend', True)
        legend_kws = kwargs.pop('legend_kws', dict())
        for i, (ax, chrom) in enumerate(zip(axs, chroms)):
            plotf(chrom=chrom, ax=ax, legend=False, **kwargs)
            if i > 0:
                ax.set_xlabel('')
                ax.set_ylabel('')
                plt.setp(ax.get_yticklabels(), visible=False)

        # legend
        if legend:
            legend_kws.setdefault('loc', 'upper left')
            legend_kws.setdefault('bbox_to_anchor', (1, 1))
            axs[-1].legend(**legend_kws)

        # tidy up
        fig.tight_layout(w_pad=0)

    def windowed_diversity_genomeplot(self, chroms, fig=None, gridspec_kws=None, **kwargs):
        """TODO"""
        kwargs.setdefault('legend_kws', dict(title='Population'))
        self.genomeplot(self.windowed_diversity_chromplot, chroms=chroms, fig=fig,
                        gridspec_kws=gridspec_kws, **kwargs)

    def windowed_diversity_delta_genomeplot(self, chroms, fig=None, gridspec_kws=None, **kwargs):
        """TODO"""
        kwargs.setdefault('legend_kws', dict(title='Population'))
        self.genomeplot(self.windowed_diversity_delta_chromplot, chroms=chroms, fig=fig,
                        gridspec_kws=gridspec_kws, **kwargs)

    # TODO windowed_diversity_delta_genome_plot

    def __repr__(self):
        r = '<PopulationAnalysis at {!r}>'.format(self.path)
        return r


class CacheMiss(Exception):
    pass
