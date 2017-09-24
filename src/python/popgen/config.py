# -*- coding: utf-8 -*-
from collections import MutableMapping
import os
import yaml
import json


class ConfigFile(object):

    def __init__(self,
                 path: str,
                 overwrite: bool=False,
                 load_kws: dict=None,
                 dump_kws: dict=None):
        self.path = path
        if load_kws is None:
            load_kws = dict()
        self.load_kws: dict = load_kws
        if dump_kws is None:
            dump_kws = dict()
        self.dump_kws: dict = dump_kws
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

    def dump(self, d: dict):
        """Write the contents of the file."""
        raise NotImplementedError

    def load(self):
        """Read the contents of the file."""
        raise NotImplementedError

    def set(self,
            section: str,
            key: str,
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
            section: str,
            key: str=None):
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
            raise KeyError('Item {0!r} has not been set in configuration section {1!r}.'
                           .format(section, key))
        else:
            return group[key]

    def set_genome_assembly(self, value):
        self.set('genome', 'assembly', value)

    def get_genome_assembly(self):
        return self.get('genome', 'assembly')

    def set_genome_annotation(self, value):
        self.set('genome', 'annotation', value)

    def get_genome_annotation(self):
        return self.get('genome', 'annotation')

    def set_genome_accessibility(self, value):
        self.set('genome', 'accessibility', value)

    def get_genome_accessibility(self):
        return self.get('genome', 'accessibility')

    def set_sample_table(self, key: str, value):
        self.set('sample_tables', key, value)

    def get_sample_table(self, key: str):
        return self.get('sample_tables', key)

    def get_sample_tables(self):
        return self.get('sample_tables')

    def set_population(self, key: str, value):
        self.set('populations', key, value)

    def get_population(self, key: str):
        return self.get('populations', key)

    def get_populations(self):
        return self.get('populations')

    def set_callset(self, key: str, value):
        self.set('callsets', key, value)

    def get_callset(self, key: str):
        return self.get('callsets', key)

    def get_callsets(self):
        return self.get('callsets')


class YAMLConfigFile(ConfigFile):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def load(self):
        with open(self.path, mode='r', encoding='utf8') as f:
            return yaml.load(f)

    def dump(self, d: dict):
        with open(self.path, mode='w', encoding='utf8') as f:
            yaml.dump(d, f, **self.dump_kws)


class JSONConfigFile(ConfigFile):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def load(self):
        with open(self.path, mode='r', encoding='utf8') as f:
            return json.load(f, **self.load_kws)

    def dump(self, d: dict):
        with open(self.path, mode='w', encoding='utf8') as f:
            json.dump(d, f, **self.dump_kws)
