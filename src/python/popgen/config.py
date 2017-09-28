# -*- coding: utf-8 -*-
from collections import MutableMapping
import os
from enum import Enum

import yaml
import json


class Section(object):
    GENOME = 'genome'
    SAMPLE_TABLES = 'sample_tables'
    POPULATIONS = 'populations'
    CALLSETS = 'callsets'


class Key(object):
    ASSEMBLY = 'assembly'
    ANNOTATION = 'annotation'
    ACCESSIBILITY = 'accessibility'
    MAIN = 'main'
    PATH = 'path'
    FORMAT = 'format'
    DATASET = 'dataset'
    SAMPLES = 'samples'
    SAMPLES_PATH = 'samples_path'
    PHASED = 'phased'
    READ_KWS = 'read_kws'
    LABEL = 'label'
    DESCRIPTION = 'description'
    QUERY = 'query'
    COLOR = 'color'
    MARKER = 'marker'
    DOWN_SAMPLE = 'down_sample'
    RANDOM_SEED = 'random_seed'
    GENOTYPE_DATASET = 'genotype_dataset'
    MAX_ALLELE = 'max_allele'


class Format(object):
    HDF5 = 'hdf5'
    ZARR = 'zarr'
    FASTA = 'fasta'
    GFF3 = 'gff3'
    CSV = 'csv'


class Statistic(object):
    DIVERSITY = 'diversity'
    WATTERSON_THETA = 'watterson_theta'
    TAJIMA_D = 'tajima_d'


class ConfigFile(object):

    def __init__(self,
                 path,
                 overwrite=False,
                 load_kws=None,
                 dump_kws=None):
        self.path = path
        if load_kws is None:
            load_kws = dict()
        self.load_kws = load_kws
        if dump_kws is None:
            dump_kws = dict()
        self.dump_kws = dump_kws
        if overwrite or not os.path.exists(path):
            self.dump(dict())

    def __repr__(self):
        r = '<{0} at {1!r}>\n'.format(type(self).__name__, self.path)
        with open(self.path, mode='r', encoding='utf8') as f:
            r += f.read()
        return r

    def __str__(self):
        with open(self.path, mode='r', encoding='utf8') as f:
            return f.read()

    def dump(self, d):
        """Write the contents of the file."""
        raise NotImplementedError

    def load(self):
        """Read the contents of the file."""
        raise NotImplementedError

    def set(self,
            section,
            key,
            value):
        """Set the value of a configuration option.

        Parameters
        ----------
        section
            Name of configuration section.
        key
            Name of configuration option.
        value
            Configuration value.

        """
        config = self.load()
        if section not in config:
            config[section] = dict()
        config[section][key] = value
        self.dump(config)

    def get(self,
            section,
            key=None):
        """Get a configuration section or option.

        Parameters
        ----------
        section
            Name of configuration section.
        key
            Name of configuration option.

        Returns
        -------
        value
            If key is None, returns the section, otherwise returns the value.

        """
        config = self.load()
        group = config.get(section, dict())
        if key is None:
            return group
        elif key not in group:
            raise KeyError('Item {!r} has not been set in configuration section {!r}.'
                           .format(key, section))
        else:
            return group[key]

    def set_genome_assembly(self, value):
        self.set(Section.GENOME, Key.ASSEMBLY, value)

    def get_genome_assembly(self):
        return self.get(Section.GENOME, Key.ASSEMBLY)

    def set_genome_annotation(self, value):
        self.set(Section.GENOME, Key.ANNOTATION, value)

    def get_genome_annotation(self):
        return self.get(Section.GENOME, Key.ANNOTATION)

    def set_genome_accessibility(self, value):
        self.set(Section.GENOME, Key.ACCESSIBILITY, value)

    def get_genome_accessibility(self):
        return self.get(Section.GENOME, Key.ACCESSIBILITY)

    def set_sample_table(self, key, value):
        self.set(Section.SAMPLE_TABLES, key, value)

    def get_sample_table(self, key):
        return self.get(Section.SAMPLE_TABLES, key)

    def get_sample_tables(self):
        return self.get(Section.SAMPLE_TABLES)

    def set_population(self, key, value):
        self.set(Section.POPULATIONS, key, value)

    def get_population(self, key):
        return self.get(Section.POPULATIONS, key)

    def get_populations(self):
        return self.get(Section.POPULATIONS)

    def set_callset(self, key, value):
        self.set(Section.CALLSETS, key, value)

    def get_callset(self, key):
        return self.get(Section.CALLSETS, key)

    def get_callsets(self):
        return self.get(Section.CALLSETS)


class YAMLConfigFile(ConfigFile):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def load(self):
        with open(self.path, mode='r', encoding='utf8') as f:
            return yaml.load(f)

    def dump(self, d):
        with open(self.path, mode='w', encoding='utf8') as f:
            yaml.dump(d, f, **self.dump_kws)


class JSONConfigFile(ConfigFile):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def load(self):
        with open(self.path, mode='r', encoding='utf8') as f:
            return json.load(f, **self.load_kws)

    def dump(self, d):
        with open(self.path, mode='w', encoding='utf8') as f:
            json.dump(d, f, **self.dump_kws)
