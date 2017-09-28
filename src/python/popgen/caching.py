# -*- coding: utf-8 -*-
from collections import deque, defaultdict
import sys
import zarr


class CacheMiss(Exception):

    def __init__(self, section, key):
        self.section = section
        self.key = key


class MemoryCache(object):

    def __init__(self, capacity):
        self.capacity = capacity
        self.size = 0
        self.queue = deque()
        self.index = defaultdict(dict)

    def get(self, section, key):
        try:
            result = self.index[section][key]
        except KeyError:
            # print('MemoryCache.get', section, key, 'miss')
            raise CacheMiss(section, key)
        else:
            # print('MemoryCache.get', section, key, 'hit')
            return result

    def put(self, section, key, result):
        # print('MemoryCache.put', section, key)
        # FIFO

        # determine size of result in bytes
        n = sizeof(result)

        # if too large by itself, no point in continuing
        if n > self.capacity:
            return

        # make space
        while self.queue and (self.size + n) > self.capacity:
            s, k, x = self.queue.pop()
            del self.index[s][k]
            self.size -= sizeof(x)

        # store
        self.queue.appendleft((section, key, result))
        self.index[section][key] = result
        self.size += n

    def clear(self, section):
        # print('MemoryCache.clear', section)
        for key in list(self.index[section].keys()):
            # remove from index
            result = self.index[section][key]
            n = sizeof(result)
            del self.index[section][key]
            # remove from queue
            self.queue.remove((section, key, result))
            # size accounting
            self.size -= n


class PersistentCache(object):

    def __init__(self, path, compressor=None):
        self.root = zarr.open_group(path, mode='a')
        self.compressor = compressor

    def get(self, section, key):
        group = self.root.require_group(section)

        # find result
        try:
            result = group[key]
        except KeyError:
            # print('PersistentCache.get', section, key, 'miss')
            raise CacheMiss(section, key)

        # check complete save
        if '__params__' not in result.attrs:
            # print('PersistentCache.get', section, key, 'miss (incomplete)')
            raise CacheMiss(section, key)

        # load
        # print('PersistentCache.get', section, key, 'hit')
        if isinstance(result, zarr.Array):
            return result[:]
        else:
            names = sorted(result.keys())
            return tuple(result[n][:] for n in names)

    def put(self, section, key, result, params):
        # print('PersistentCache.put', section, key)

        # ensure group
        group = self.root.require_group(section)

        # clear previous
        if key in group:
            del group[key]

        # handle multiple result
        if isinstance(result, (list, tuple)):

            # make group to hold results
            cached_result = group.create_group(key)

            # determine names for results
            names = ['f{:02d}'.format(i) for i in range(len(result))]

            # create datasets to store results
            for n, d in zip(names, result):
                cached_result.create_dataset(n, data=d, compressor=self.compressor)

        else:

            # create dataset to store result
            cached_result = group.create_dataset(key, data=result, compressor=self.compressor)

        # mark success
        cached_result.attrs['__params__'] = params

    def clear(self, section):
        # print('PersistentCache.clear', section)

        if section in self.root:
            del self.root[section]


def sizeof(result):
    if hasattr(result, 'nbytes'):
        return result.nbytes
    if isinstance(result, (list, tuple)):
        return sum(sizeof(x) for x in result)
    return sys.getsizeof(result)
