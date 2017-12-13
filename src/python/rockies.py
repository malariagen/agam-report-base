# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import bisect
import collections


import numpy as np
import matplotlib.pyplot as plt


def extract_signal(df, col, chrom, recmap, spacing=20000):
    """Extract a genome-wide selection scan signal from a dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        A dataframe with windowed data. Expected to contain columns "chrom",
        "start", "stop", as well as `col`. Expect that "start" and "stop" are
        GFF-style 1-based stop-inclusive coords.
    col : str
        Name of column containing signal.
    chrom : str
        Chromosome to extract. May be a tuple, e.g., ('2R', '2L'), in which
        case the data for each chromosome will be concatenated.
    recmap : dict [str -> array]
        Recombination map. A dictionary mapping chromosome names onto arrays,
        where each array holds the absolute recombination rate for each base
        in the chromosome.
    spacing : int, optional
        Amount of physical distance to insert when concatenating data from
        multiple chromosome arms.

    Returns
    -------
    ppos
    gpos
    signal

    """

    # handle multiple chromosomes
    if isinstance(chrom, (tuple, list)):
        assert len(chrom) == 2, 'can only concatenate two chromosomes'
        (ppos1, gpos1, signal1), (ppos2, gpos2, signal2) = \
            [extract_signal(df, col, chrom=c, recmap=recmap) for c in chrom]
        ppos = np.concatenate([ppos1, ppos2 + np.max(ppos1) + spacing])
        gpos = np.concatenate([gpos1, (gpos2 +
                                       np.max(gpos1) +
                                       recmap[chrom[1]][0] * spacing)])
        signal = np.concatenate([signal1, signal2])
        return ppos, gpos, signal

    # extract data for single chromosome arm
    df = df.reset_index().set_index('chrom')
    df_chrom = df.loc[chrom]

    # compute the physical position for each window
    ppos = (df_chrom.start + df_chrom.stop) / 2

    # compute genetic length of each window
    glen = df_chrom.apply(
        lambda row: np.sum(recmap[chrom][int(row.start-1):int(row.stop)]),
        axis=1
    )

    # compute the genetic length position for each window
    gpos = np.cumsum(glen)

    # extract the signal column
    signal = df_chrom[col]

    return ppos, gpos, signal


def plot_signal(df, col, chrom, recmap, spacing=20000, start=None, stop=None,
                figsize=(16, 3), distance='physical', ax=None):
    """Convenience function to plot a windowed selection signal."""

    # extract data
    ppos, gpos, signal = extract_signal(df, col, chrom, recmap, spacing=spacing)

    # determine x coord
    if distance == 'genetic':
        x = gpos
    else:
        x = ppos

    if start is not None or stop is not None:
        flt = np.ones(x.shape, dtype=bool)
        if start is not None:
            flt = flt & (x >= start)
        if stop is not None:
            flt = flt & (x <= stop)
        x = x[flt]
        signal = signal[flt]

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.plot(x, signal, linestyle=' ', marker='o', mfc='none', markersize=3)
    ax.set_title(col)
    ax.set_xlabel('Chromosome {} ({} distance)'.format(chrom, distance))
    if fig is not None:
        fig.tight_layout()
    return ax


Fit = collections.namedtuple('Fit',
                             'delta_aic peak_result_l peak_result_r '
                             'null_result_l null_result_r x y xl yl xr yr')


def fit_peak(gpos, signal, idx_center, flank_size, peak_model, peak_params,
             null_model, null_params):
    """Fit peak and null models to the data at window index `i`.

    Parameters
    ----------
    gpos : ndarray
        Genetic length position.
    signal : ndarray
        Selection signal.
    idx_center : int
        Index of window to fit peak on.
    flank_size : float
        Size of flanks to include when fitting.
    peak_model : Model
    peak_params : Parameters
    null_model : Model
    null_params : Parameters

    Returns
    -------
    TODO

    """

    # central position
    gpos_center = gpos[idx_center]

    # obtain data for region to fit on
    idx_left = bisect.bisect_left(gpos, gpos_center - flank_size)
    idx_right = bisect.bisect_right(gpos, gpos_center + flank_size)
    # include center in both flanks
    loc_left = slice(idx_left, idx_center + 1)
    loc_right = slice(idx_center, idx_right)
    loc = slice(idx_left, idx_right)
    xl = gpos_center - gpos[loc_left]
    xr = gpos[loc_right] - gpos_center
    x = gpos[loc] - gpos_center
    yl = signal[loc_left]
    yr = signal[loc_right]
    y = signal[loc]

    # strip nans
    missing = np.isnan(y)
    x = x[~missing]
    y = y[~missing]
    missingl = np.isnan(yl)
    xl = xl[~missingl]
    yl = yl[~missingl]
    missingr = np.isnan(yr)
    xr = xr[~missingr]
    yr = yr[~missingr]

    # setup output variables
    delta_aic = 0
    peak_result_l = peak_result_r = None
    null_result_l = null_result_r = None

    # check we have some data to fit
    if xl.shape[0] > 2 and xr.shape[0] > 2:

        # fit null model - allow for one outlier, remove maximum
        no_outlier_l = xl < xl.max()
        no_outlier_r = xr < xr.max()
        null_result_l = null_model.fit(yl[no_outlier_l], x=xl[no_outlier_l],
                                       params=null_params)
        null_result_r = null_model.fit(yr[no_outlier_r], x=xr[no_outlier_r],
                                       params=null_params)

        # fit peak model
        peak_result_l = peak_model.fit(yl, x=xl, params=peak_params)
        peak_result_r = peak_model.fit(yr, x=xr, params=peak_params)

        # compute delta AIC as minimum of the flanks, requires both flanks to
        # have a good fit
        delta_aic = min(null_result_l.aic - peak_result_l.aic,
                        null_result_r.aic - peak_result_r.aic)

    return Fit(delta_aic, peak_result_l, peak_result_r, null_result_l,
               null_result_r, x, y, xl, yl, xr, yr)


def plot_peak(gpos, signal, idx_center, flank_size, peak_model, peak_params,
              null_model, null_params, plot=True, report=True,
              figsize=(12, 4)):

    fit = fit_peak(gpos, signal, idx_center, flank_size, peak_model,
                   peak_params, null_model, null_params)

    if report:
        print(fit.peak_result_l.fit_report())
        print(fit.peak_result_r.fit_report())
        print(fit.null_result_l.fit_report())
        print(fit.null_result_r.fit_report())

    # plot the fit
    if plot:
        fig, ax = plt.subplots(figsize=figsize)
        peak_kws = dict(color='k', linestyle='--', lw=2)
        null_kws = dict(color='r', linestyle='--', lw=2)
        signal_kws = dict(marker='o', linestyle=' ', mfc='none', mew=2)
        ax.plot(-fit.xl, fit.peak_result_l.best_fit, **peak_kws)
        ax.plot(fit.xr, fit.peak_result_r.best_fit, **peak_kws)
        ax.plot(-fit.xl, [fit.peak_result_l.params['c']] * fit.xl.shape[0],
                **peak_kws)
        ax.plot(fit.xr, [fit.peak_result_r.params['c']] * fit.xr.shape[0],
                **peak_kws)
        ax.plot(-fit.xl, [fit.null_result_l.best_fit] * fit.xl.shape[0],
                **null_kws)
        ax.plot(fit.xr, [fit.null_result_r.best_fit] * fit.xr.shape[0],
                **null_kws)
        ax.plot(fit.x, fit.y, **signal_kws)
        ax.set_ylim(bottom=0)
        ax.text(.98, .98,
                ('[[Left Flank]]\n'
                 '[[AIC]] Peak: {:.1f}; Null: {:.1f}; Delta: {:.1f}\n'
                 '[[BIC]] Peak: {:.1f}; Null: {:.1f}; Delta: {:.1f}\n'
                 .format(fit.peak_result_r.aic, fit.null_result_r.aic,
                         fit.null_result_r.aic - fit.peak_result_r.aic,
                         fit.peak_result_r.bic, fit.null_result_r.bic,
                         fit.null_result_r.bic - fit.peak_result_r.bic)),
                ha='right', va='top', transform=ax.transAxes)
        ax.text(.02, .98,
                ('[[Right Flank]]\n'
                 '[[AIC]] Peak: {:.1f}; Null: {:.1f}; Delta: {:.1f}\n'
                 '[[BIC]] Peak: {:.1f}; Null: {:.1f}; Delta: {:.1f}\n'
                 .format(fit.peak_result_l.aic, fit.null_result_l.aic,
                         fit.null_result_l.aic - fit.peak_result_l.aic,
                         fit.peak_result_l.bic, fit.null_result_l.bic,
                         fit.null_result_l.bic - fit.peak_result_l.bic)),
                ha='left', va='top', transform=ax.transAxes)
        fig.tight_layout()
