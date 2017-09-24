# -*- coding: utf-8 -*-
import os
import sys


from .jsonfile import JSONFile
from . import util


class AnalysisSetup(object):

    def __init__(self,
                 path: str,
                 overwrite: bool=False,
                 logfile=sys.stdout):
        self.path = path
        store_path = os.path.join(self.path, 'setup.json')
        self.store = JSONFile(store_path, overwrite=overwrite)
        self.log = util.Logger(logfile)

    def _get_dict(self, field):
        return self.store.get(field, dict())

    def _set_dict_item(self, field, key, value):
        d = self._get_dict(field)
        d[key] = value
        self.store[field] = d

    def get_sample_data_setups(self):
        return self._get_dict('sample_data')

    def get_sample_data_setup(self, name):
        setups = self.get_sample_data_setups()
        setup = setups.get(name, None)
        if setup is None:
            raise KeyError('sample data %r has not been set' % name)
        return setup

    def _set_sample_data_setup(self, name, setup):
        self._set_dict_item('sample_data', name, setup)

    def set_sample_data(self,
                        path: str,
                        index_col: str,
                        name: str='main',
                        reader: str='read_csv',
                        read_kws: dict=None):
        """Add or set a sample data source.

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
        df = util.read_dataframe(path, reader, read_kws)

        # diagnostics
        self.log('setting sample data %r with %s rows' % (name, len(df)))

        # store in configuration
        setup = dict(path=path, reader=reader, read_kws=read_kws)
        self._set_sample_data_setup(name, setup)

    def get_populations_setups(self):
        return self._get_dict('populations')

    def get_populations_setup(self, pop):
        setups = self.get_populations_setups()
        if pop not in setups:
            raise ValueError('population %r has not been set' % pop)
        return setups[pop]

    def get_population_query(self, pop):
        setup = self.get_populations_setup(pop)
        return setup['query']

    def _set_populations_setup(self, pop, setup):
        self._set_dict_item('populations', pop, setup)

    def set_population(self, name, query):
        """TODO"""

        # check query
        df_samples = self.load_sample_data_joined()
        df_pop = df_samples.query(query)
        if len(df_pop) == 0:
            raise ValueError('query does not match any samples')
        else:
            self.log('setting population %r with %s samples' % (name, len(df_pop)))

        # store config
        params = dict(query=query)
        self._set_populations_setup(name, params)

    def load_sample_data(self, name: str='main', pop=None):
        """TODO"""

        # obtain read params
        setup = self.get_sample_data_setup(name)

        # read into dataframe
        df = util.read_dataframe(setup['path'], setup['reader'], setup['read_kws'])

        # subset to population
        if pop:
            query = self.get_population_query(pop)
            df = df.query(query)

        return df

    def load_sample_data_joined(self, on=None, how='outer', pop=None):
        """TODO"""

        # find all sample data configurations
        sample_data_config = self.get_sample_data_setups()
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
            query = self.get_population_query(pop)
            df = df.query(query)

        return df




