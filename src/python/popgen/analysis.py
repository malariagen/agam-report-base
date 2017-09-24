# -*- coding: utf-8 -*-


# standard library imports
import os
import sys


# third party imports
import numpy as np
import pyfasta
import allel


# internal imports
from .config import YAMLConfigFile
from .util import Logger


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
        self._genome_assembly = None
        self._genome_annotation = None

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
            raise NotImplementedError("Only 'fasta' format currently supported.")

        # check valid fasta
        try:
            abs_path = self.get_abs_path(path)
            assembly = pyfasta.Fasta(abs_path, key_fn=lambda v: v.split()[0])

        except Exception as e:
            raise ValueError('An error occurred opening the genome assembly file: {}\n'
                             'Does the path contain a valid fasta file?'
                             .format(e))

        else:

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
    def genome_assembly(self):
        """TODO"""
        if self._genome_assembly is None:
            config = self.config.get_genome_assembly()
            path = self.get_abs_path(config['path'])
            self._genome_assembly = pyfasta.Fasta(path, key_fn=lambda v: v.split()[0])
        return self._genome_assembly

    def load_reference_sequence(self, chrom):
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
            raise NotImplementedError("Only 'gff3' format currently supported.")

        # check valid gff3
        try:
            abs_path = self.get_abs_path(path)
            annotation = allel.gff3_to_dataframe(abs_path, attributes=['ID', 'Parent'])

        except Exception as e:
            raise ValueError('An error occurred opening the genome assembly file: {}\n'
                             'Does the path contain a valid fasta file?'
                             .format(e))

        else:

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

    def load_genome_annotation(self, **kwargs):
        """TODO"""
        config = self.config.get_genome_annotation()
        path = self.get_abs_path(config['path'])
        annotation = allel.gff3_to_dataframe(path, **kwargs)
        return annotation
