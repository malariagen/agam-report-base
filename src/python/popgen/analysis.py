# -*- coding: utf-8 -*-


# standard library imports
import os
import sys
from typing import Mapping


# third party imports
import numpy as np
import pandas
import pyfasta
import allel
import h5py
import zarr


# internal imports
from .config import YAMLConfigFile
from .util import Logger, read_dataframe, guess_callset_format, open_callset, close_callset


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
                            title: str=None,
                            description: str=None,
                            format: str='fasta'):
        """Set the reference genome assembly.

        Parameters
        ----------
        path
            Path to assembly file. If given as a relative path, must be relative to the analysis
            directory.
        title
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
        self.log('Setting genome assembly; found {0} sequences, total size {1:,} bp.'
                 .format(n_chroms, n_bases))

        # store configuration
        config = dict(
            path=path, title=title, description=description, format=format
        )
        self.config.set_genome_assembly(config)

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
                              title: str=None,
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
            path=path, title=title, description=description, format=format
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
                                 title: str=None,
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
            n_bases_total = 0
            n_bases_total_accessible = 0

            for chrom in sorted(self.genome_assembly):
                dataset_path = '/'.join([chrom, dataset_name])

                # check dataset exists
                if dataset_path not in accessibility:
                    raise ValueError('Could not find accessibility dataset for chromosome '
                                     '{!r}.'.format(chrom))

                # check dataset type and shape
                is_accessible = accessibility[dataset_path]
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

                n_bases_total += n_bases
                n_bases_total_accessible += np.count_nonzero(is_accessible[:])

        finally:
            accessibility.close()

        # log some diagnostics
        pc_genome_accessible = n_bases_total_accessible * 100 / n_bases_total
        self.log('Setting genome accessibility; found {:,} accessible bases '
                 '({:.1f}% of genome).'
                 .format(n_bases_total_accessible, pc_genome_accessible))

        # store configuration
        config = dict(
            path=path, format=format, dataset_name=dataset_name, title=title,
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
                         title: str=None,
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
            path=path, format=format, read_kws=read_kws, title=title, description=description
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

    def set_callset(self, path, format=None, name='main', phased=False, title=None,
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
            path=path, format=format, phased=phased, title=title, description=description
        )
        self.config.set_callset(name, config)
