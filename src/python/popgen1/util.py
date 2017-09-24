# -*- coding: utf-8 -*-
import os
import pandas


def check_path_exists(path):
    if not os.path.exists(path):
        raise ValueError('path does not exist: %r' % path)


def check_isfile(path):
    if not os.path.exists(path):
        raise ValueError('file does not exist: %r' % path)
    if not os.path.isfile(path):
        raise ValueError('path is not a file: %r' % path)


def check_pandas_reader(reader):
    if not reader.startswith('read_') or not hasattr(pandas, reader):
        raise ValueError('invalid reader %r; must be the name of a pandas function that can read '
                         'data into a data frame' % reader)


def read_dataframe(path, reader, read_kws):
    check_isfile(path)
    check_pandas_reader(reader)
    read = getattr(pandas, reader)
    try:
        df = read(path, **read_kws)
    except Exception as e:
        raise ValueError('data could not be read by pandas: %s' % e)
    return df


class Logger(object):

    def __init__(self, logfile):
        self.logfile = logfile

    def __call__(self, *args):
        if self.logfile is not None:
            if isinstance(self.logfile, str):
                with open(self.logfile, mode='a') as f:
                    print(*args, file=f)
            else:
                print(*args, file=self.logfile)
