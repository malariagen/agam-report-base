# -*- coding: utf-8 -*-
from collections import MutableMapping
import json
import os


class JSONFile(MutableMapping):

    def __init__(self,
                 path: str,
                 overwrite: bool=False,
                 load_kws: dict=None,
                 dump_kws: dict=None):
        """
        A MutableMapping backed by a JSON file.

        Parameters
        ----------
        path
            File path.
        overwrite
            If True, replace existing file.
        load_kws
            Default keyword arguments when loading.
        dump_kws
            Default keyword arguments when dumping.

        """
        self.path = path
        if load_kws is None:
            load_kws = dict()
        self.load_kws = load_kws
        if dump_kws is None:
            dump_kws = dict()
        self.dump_kws = dump_kws
        if overwrite or not os.path.exists(self.path):
            self.dump(dict())

    def load(self, **kwargs):
        load_kws = self.load_kws.copy()
        load_kws.update(kwargs)
        with open(self.path, mode='r', encoding='utf8') as fp:
            obj = json.load(fp, **load_kws)
        return obj

    def dump(self, obj, **kwargs):
        dump_kws = self.dump_kws.copy()
        dump_kws.update(kwargs)
        with open(self.path, mode='w', encoding='utf8') as fp:
            json.dump(obj, fp, **dump_kws)

    def __getitem__(self, key):
        d = self.load()
        return d[key]

    def __iter__(self):
        d = self.load()
        return iter(d)

    def __len__(self):
        d = self.load()
        return len(d)

    def __setitem__(self, key, value):
        d = self.load()
        d[key] = value
        self.dump(d)

    def __delitem__(self, key):
        d = self.load()
        del d[key]
        self.dump(d)

    def __repr__(self):
        with open(self.path, mode='r', encoding='utf8') as fp:
            return fp.read()
