# -*- coding: utf-8 -*-
import os
import pandas
import h5py
import zarr


def check_path_exists(path):
    if not os.path.exists(path):
        raise ValueError('Path does not exist: {!r}.'.format(path))


def check_isfile(path):
    if not os.path.exists(path):
        raise ValueError('File does not exist: {!r}'.format(path))
    if not os.path.isfile(path):
        raise ValueError('Path is not a file: {!r}'.format(path))


def get_pandas_reader(format):
    reader = 'read_' + format
    if not hasattr(pandas, reader):
        raise ValueError('Format not supported: {!r}.'.format(format))
    read = getattr(pandas, reader)
    return read


def read_dataframe(path, format, read_kws):
    check_isfile(path)
    read = get_pandas_reader(format)
    try:
        df = read(path, **read_kws)
    except Exception as e:
        raise ValueError('Data could not be read by pandas: %s.' % e)
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


def guess_callset_format(path):
    if not os.path.exists(path):
        raise ValueError('Path does not exist: {!r}.'.format(path))
    if os.path.isfile(path):
        try:
            with h5py.File(path, mode='r') as h5f:
                return 'hdf5'
        except:
            pass
        # TODO Zarr zip format
    elif os.path.isdir(path):
        try:
            zarr.open_group(path, mode='r')
            return 'zarr'
        except:
            pass
    raise NotImplementedError('Unsupported format for callset at path {!r}'.format(path))


def open_callset(path, format, mode='r'):
    if format == 'hdf5':
        return h5py.File(path, mode=mode)
    if format == 'zarr':
        return zarr.open_group(path, mode=mode)
    raise NotImplementedError('Format {!r} not currently supported.'.format(format))


def close_callset(callset):
    if hasattr(callset, 'close'):
        callset.close()
