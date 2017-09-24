# -*- coding: utf-8 -*-


# standard library imports
import os
import sys


# third party imports
import pyfasta


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

    def set_genome_assembly(self,
                            path: str,
                            title: str=None,
                            description: str=None):
        """Set the reference genome assembly.

        Parameters
        ----------
        path
            Path to FASTA file. If given as a relative path, must be relative to the analysis
            directory.
        title
            Name of the assembly.
        description
            Description of the assembly.

        """

        # normalise path - relative paths treated as relative to analysis directory
        if os.path.isabs(path):
            abs = path
        else:
            abs = os.path.abspath(os.path.join(self.path, path))

        # check valid FASTA
        try:
            assembly = pyfasta.Fasta(abs, key_fn=lambda v: v.split()[0])

        except Exception as e:
            raise ValueError('An error occurred opening the genome assembly file: {}\n'
                             'Is the file a valid FASTA file?'
                             .format(e))

        else:
            n_chroms = len(assembly)
            n_bases = sum(len(assembly[chrom]) for chrom in assembly)
            self.log('Setting genome assembly with {0} chromosomes, total size {1:,} bp.'
                     .format(n_chroms, n_bases))
            config = dict(
                path=path, title=title, description=description, format='fasta'
            )
            self.config.set_genome_assembly(config)
