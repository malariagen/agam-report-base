# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import os


import re
import zarr
import petl as etl
import pandas
import h5py


title = 'Phase 1 selection data release'


# noinspection PyGlobalUndefined
def init(release_dir):
    """Initialise data resources.

    Parameters
    ----------
    release_dir : string
        Local filesystem path where data from the release are stored.

    """

    # windowed summary data
    ###########

    global ihs_windowed, xpehh_windowed, xpclr_windowed, h12_windowed, dTjD_windowed 
    ws_dir = os.path.join(release_dir, 'windowed')
    
    for metric in ("ihs", "xpehh", "xpclr", "hstats", "dTjD"):

        fn = os.path.join(ws_dir, '{metric}.txt.gz').format(metric=metric)

        if os.path.exists(fn):
            exec("{metric}_windowed=pandas.read_csv('{path}', sep='\t')".format(metric=metric, path=fn), globals())


    # raw data 
    ##########
    table_stems = {"hstats": "hstats_table_([A-Za-z0-9]+)_(.+)\.txt\.gz",
                   "xpclr": "[A-Za-z0-9]+_(.+vs+.)\.xpclr\.txt\.gz"}

    global ihs_raw, xpehh_raw, xpclr_raw, hstats_raw, dTjD_raw 
    raw_dir = os.path.join(release_dir, 'raw')
    
    for metric in ("ihs", "xpehh", "xpclr", "hstats", "dTjD"):

        fn = os.path.join(raw_dir, '{metric}/output.zarr').format(metric=metric)

        if os.path.exists(fn):
            exec("{metric}_raw=zarr.open_group('{path}', 'r')".format(metric=metric, path=fn), globals())
            continue

        ## otherwise we can load tables?
        #temp = {}
        #output_dir = os.path.join(raw_dir, metric, 'output')
        #files = os.listdir(output_dir)
        #for f in files:
        #    mm = re.search(table_stems[metric], f)
        #    if mm is not None:
        #        chrom, comp = mm.groups()
        #        if chrom not in temp:
        #            temp[chrom] = {}
        #        temp[chrom][comp] = lambda: pandas.read_csv(os.path.join(output_dir, f), sep="\t")
        #
        #print("exec", metric, len(files), list(temp.keys())) 
        #exec("{metric}_raw=temp".format(metric=metric))
