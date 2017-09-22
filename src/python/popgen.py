# -*- coding: utf-8 -*-
import os
from collections import MutableMapping
import json
import sys
from contextlib import contextmanager
import hashlib
import inspect
import allel
import zarr


class JSONFile(MutableMapping):

    def __init__(self,
                 path: str,
                 overwrite: bool=False,
                 load_kws: dict=None,
                 dump_kws: dict=None):
        """
        A MutableMapping interface to a JSON file.

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


def _check_path_exists(path):
    if not os.path.exists(path):
        raise ValueError('path does not exist: %r' % path)


def _check_isfile(path):
    if not os.path.exists(path):
        raise ValueError('file does not exist: %r' % path)
    if not os.path.isfile(path):
        raise ValueError('path is not a file: %r' % path)


def _check_pandas_reader(reader):
    import pandas
    if not reader.startswith('read_') or not hasattr(pandas, reader):
        raise ValueError('invalid reader %r; must be the name of a pandas function that can read '
                         'data into a data frame' % reader)


def _read_dataframe(path, reader, read_kws):
    _check_isfile(path)
    _check_pandas_reader(reader)
    import pandas
    read = getattr(pandas, reader)
    try:
        df = read(path, **read_kws)
    except Exception as e:
        raise ValueError('data could not be read by pandas: %s' % e)
    return df


def _open_callset(path, mode='r'):

    if (os.path.isfile(path) and
            (path.endswith('.hdf5') or path.endswith('.h5'))):
        import h5py
        return h5py.File(path, mode=mode)

    elif os.path.isdir(path):
        return zarr.open_group(path, mode=mode)

    # TODO zarr in zip file

    else:
        raise ValueError('cound not open callset at path %r' % path)


@contextmanager
def _open_callset_context(path, mode='r'):
    callset = _open_callset(path, mode=mode)
    if hasattr(callset, 'close'):
        try:
            yield callset
        finally:
            callset.close()
    else:
        yield callset


class PopulationAnalysis(object):

    def __init__(self,
                 path: str,
                 overwrite: bool=False,
                 logfile=sys.stdout,
                 cache_compressor=None):
        self.path = path
        # ensure directory exists
        os.makedirs(self.path, exist_ok=True)
        # setup config file
        config_path = os.path.join(self.path, 'config.json')
        dump_kws = dict(indent=4, sort_keys=True)
        self.config = JSONFile(path=config_path, dump_kws=dump_kws, overwrite=overwrite)
        # store reference to logfile
        self.logfile = logfile
        # setup cache compressor
        if cache_compressor is None:
            cache_compressor = zarr.Blosc(cname='zstd', clevel=1, shuffle=0)
        self.cache_compressor = cache_compressor

    def _get_config(self, key):
        return self.config.get(key, dict())

    def _set_config_field(self, key: str, field: str, obj: dict):
        config = self._get_config(key)
        config[field] = obj
        self.config[key] = config

    def _get_sample_data_sources(self):
        return self._get_config('sample_data')

    def _set_sample_data_source(self, name: str, obj: dict):
        self._set_config_field('sample_data', name, obj)

    def _get_populations(self):
        return self._get_config('populations')

    def _set_population(self, pop: str, obj: dict):
        self._set_config_field('populations', pop, obj)

    def _get_population_query(self, pop: str):
        pops = self._get_populations()
        if pop not in pops:
            raise ValueError('population %r has not been set' % pop)
        return pops[pop]['query']

    def _get_callsets(self):
        return self._get_config('callsets')

    def _get_callset(self, name):
        callsets = self._get_callsets()
        if name not in callsets:
            raise ValueError('callset %r has not been set' % name)
        return callsets[name]

    def _get_callset_path(self, name):
        return self._get_callset(name)['path']

    def _set_callset(self, name: str, obj: dict):
        self._set_config_field('callsets', name, obj)

    def set_sample_data(self,
                        path: str,
                        index_col: str,
                        name: str='main',
                        reader: str='read_csv',
                        read_kws: dict=None):
        """Add or set a sample data file.

        Parameters
        ----------
        path
            File path.
        index_col
            Index column, i.e., name of column containing unique sample identifier.
        name
            Name to reference the sample data file.
        reader
            Name of pandas function to use to read the data into a data frame.
        read_kws
            Passed through to pandas reader function.

        """

        # check args
        if read_kws is None:
            read_kws = dict()
        # force index_col into read keywords
        read_kws['index_col'] = index_col
        # check data can be read by pandas
        df = _read_dataframe(path, reader, read_kws)

        # print some diagnostics
        self.log('setting sample data %r with %s rows' % (name, len(df)))

        # store in configuration
        params = dict(path=path, reader=reader, read_kws=read_kws)
        self._set_sample_data_source(name, params)

    def load_sample_data(self, name: str='main', pop=None):

        # obtain read params
        config = self._get_sample_data_sources()
        params = config.get(name, None)
        if params is None:
            raise KeyError('sample data %r has not been set' % name)

        # read into dataframe
        df = _read_dataframe(params['path'], params['reader'], params['read_kws'])

        # subset to population
        if pop:
            query = self._get_population_query(pop)
            df = df.query(query)

        return df

    def load_sample_data_joined(self, on=None, how='outer', pop=None):

        # find all sample data configurations
        sample_data_config = self._get_sample_data_sources()
        names = sorted(sample_data_config.keys())
        if len(names) == 0:
            raise ValueError('no sample data has been set')

        # determine which to treat as the first dataframe in the join
        if 'main' in names:
            first = 'main'
        else:
            first = names[0]
        others = [n for n in names if n != first]

        # load and join
        df = self.load_sample_data(first)
        if others:
            df_others = [self.load_sample_data(n) for n in others]
            df = df.join(df_others, on=on, how=how)

        # subset to population
        if pop:
            query = self._get_population_query(pop)
            df = df.query(query)

        return df

    def set_population(self, name, query):

        # check query
        df_samples = self.load_sample_data_joined()
        df_pop = df_samples.query(query)
        if len(df_pop) == 0:
            raise ValueError('query does not match any samples')
        else:
            self.log('setting population %r with %s samples' % (name, len(df_pop)))

        # store config
        params = dict(query=query)
        self._set_population(name, params)

    def log(self, *msg):
        if isinstance(self.logfile, str):
            with open(self.logfile, mode='a') as f:
                print(*msg, file=f)
        else:
            print(*msg, file=self.logfile)

    def set_callset(self, path, name='main'):
        """TODO"""

        # check args
        _check_path_exists(path)
        # check can open
        with _open_callset_context(path) as callset:

            # print some diagnostics
            keys = list(callset.keys())
            if 'samples' not in keys:
                raise ValueError("callset does not contain required field 'samples'")
            contigs = [k for k in keys if k != 'samples']
            samples = callset['samples'][:]
            self.log('setting callset %r with %s samples and %s contigs' %
                     (name, len(samples), len(contigs)))

        # store in config
        self._set_callset(name, dict(path=path))

    def open_callset(self, name='main'):

        # obtain params
        config = self._get_callsets()
        params = config.get(name, None)
        if params is None:
            raise KeyError('callset %r has not been set' % name)

        # open
        path = params['path']
        callset = _open_callset(path)
        return callset

    def __repr__(self):
        r = ""
        r += "Population Analysis\n"
        r += "-------------------\n"
        r += "Path: %s\n" % self.path

        try:

            # setup for information about samples and populations
            df_samples = self.load_sample_data_joined()

        except ValueError:
            pass

        else:
            r += "Total no. samples: %s\n" % len(df_samples)

            # population sample sizes
            pops = self._get_populations()
            if pops:
                r += "Populations:\n"
                for population, params in sorted(pops.items()):
                    query = params['query']
                    r += "  - %s (n=%s) %r\n" % (population, len(df_samples.query(query)), query)

        return r

    def allele_counts(self, contig, pop=None, callset='main', max_allele=3, gt_name=None):

        # setup parameters
        params = dict(
            contig=contig,
            pop=pop,
            callset=self._get_callset(callset),
            max_allele=max_allele,
            gt_name=gt_name,
        )

        # hash and check cache
        cache_group = inspect.currentframe().f_code.co_name
        cache_key = _hash_params(params)
        try:
            ac = self.cache_load(cache_group, cache_key)

        except CacheMiss:
            self.log('computing', cache_group, cache_key)

            callset_path = self._get_callset_path(callset)
            with _open_callset_context(callset_path) as callset:

                # get genotype data
                gt = _get_genotype_dataset(callset, contig, gt_name)

                # deal with pop
                sample_indices = None
                if pop:
                    sample_indices = self._get_sample_indices(callset, pop)

                # run computation
                ac = gt.count_alleles(max_allele=max_allele, subpop=sample_indices).compute()

                # save in cache
                self.cache_save(cache_group, cache_key, ac)

        return allel.AlleleCountsArray(ac)

    def _get_sample_indices(self, callset, pop):
        samples = _get_samples_list(callset)
        df_pop = self.load_sample_data_joined(pop=pop)
        pop_samples = df_pop.index.values
        sample_indices = [samples.index(s) for s in pop_samples]
        return sample_indices

    def cache_save(self, cache_group, cache_key, data, names=None):
        self.log('saving', cache_group, cache_key)
        cache = zarr.open_group(os.path.join(self.path, 'cache'), mode='a')
        group = cache.require_group(cache_group)
        if names is None:
            result = group.require_dataset(cache_key, dtype=data.dtype, shape=data.shape,
                                           compressor=self.cache_compressor)
            result[:] = data
        else:
            result = group.require_group(cache_key)
            for n, d in zip(names, data):
                ds = result.require_dataset(n, dtype=d.dtype, shape=d.shape,
                                            compressor=self.cache_compressor)
                ds[:] = d
        # mark success
        result.attrs['__ok__'] = True

    def cache_load(self, cache_group, cache_key, names=None):
        cache = zarr.open_group(os.path.join(self.path, 'cache'), mode='a')
        group = cache.require_group(cache_group)
        if cache_key not in group:
            raise CacheMiss
        result = group[cache_key]
        if '__ok__' not in result.attrs:
            raise CacheMiss
        self.log('loading', cache_group, cache_key)
        if names is None:
            return result[:]
        else:
            return tuple(result[n][:] for n in names)


class CacheMiss(Exception):
    pass


def _get_genotype_dataset(callset, contig, dataset_name):
    # acommodate different names for the genotype dataset
    gt = None
    if dataset_name:
        gt_path = '/'.join([contig, 'calldata', dataset_name])
        if gt_path in callset:
            gt = callset[gt_path]
    else:
        for dataset_name in 'genotype', 'GT':
            gt_path = '/'.join([contig, 'calldata', dataset_name])
            if gt_path in callset:
                gt = callset[gt_path]
    if gt is None:
        raise RuntimeError('could not find genotype dataset for contig %r' % contig)
    else:
        return allel.GenotypeDaskArray(gt)


def _get_samples_list(callset):
    samples = callset['samples'][:]
    if samples.dtype.kind == 'S':
        samples = samples.astype('U')
    return samples.tolist()


def _hash_params(params):
    s = json.dumps(params, indent=4, sort_keys=True, ensure_ascii=True, skipkeys=False)
    k = hashlib.md5(s.encode('ascii')).hexdigest()
    return k
