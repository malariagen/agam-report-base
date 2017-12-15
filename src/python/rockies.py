# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from bisect import bisect_left, bisect_right
import collections


import numpy as np
import matplotlib.pyplot as plt
import lmfit
import seaborn as sns
palette = sns.color_palette()


def extract_signal(df, col, chrom, recmap, spacing=0):
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
        Chromosome to extract. May also be a tuple, e.g., ('2R', '2L'), in which
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
    starts
    stops
    gpos
    signal

    """

    # handle multiple chromosomes
    if isinstance(chrom, (tuple, list)):
        assert len(chrom) == 2, 'can only concatenate two chromosomes'
        (starts1, stops1, gpos1, signal1), (starts2, stops2, gpos2, signal2) = \
            [extract_signal(df, col, chrom=c, recmap=recmap) for c in chrom]
        starts = np.concatenate([starts1, starts2 + np.max(stops1) + spacing])
        stops = np.concatenate([stops1, stops2 + np.max(stops1) + spacing])
        gpos = np.concatenate([gpos1, (gpos2 +
                                       np.max(gpos1) +
                                       recmap[chrom[1]][0] * spacing)])
        signal = np.concatenate([signal1, signal2])
        return starts, stops, gpos, signal

    # extract data for single chromosome arm
    df = df.reset_index().set_index('chrom')
    df_chrom = df.loc[chrom]

    # extract window starts and stops
    starts = np.asarray(df_chrom.start)
    stops = np.asarray(df_chrom.stop)

    # compute genetic length of each window
    glen = np.array([
        np.sum(recmap[chrom][int(start-1):int(stop)])
        for start, stop in zip(starts, stops)
    ])

    # compute the genetic length position for each window
    gpos = np.cumsum(glen)

    # extract the signal column
    signal = np.asarray(df_chrom[col])

    return starts, stops, gpos, signal


def plot_signal(df, col, chrom, recmap, spacing=20000, start=None, stop=None,
                figsize=(16, 3), distance='physical', ax=None):
    """Convenience function to plot a windowed selection signal."""

    # extract data
    starts, stops, gpos, signal = extract_signal(df, col, chrom, recmap,
                                                 spacing=spacing)

    # determine x coord
    if distance == 'genetic':
        x = gpos
    else:
        x = (starts + stops) / 2

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

    ax.plot(x, signal, linestyle=' ', marker='o', mfc='none', markersize=3,
            color=palette[0])
    ax.set_title(col)
    ax.set_xlabel('Chromosome {} ({} distance)'.format(chrom, distance))
    if fig is not None:
        fig.tight_layout()
    return ax


def exponential(x, amplitude, decay, c, cap):
    """Exponential decay function.

    Parameters
    ----------
    x : ndarray
        Independent variable.
    amplitude : float
        Amplitude parameter.
    decay : float
        Decay parameter.
    c : float
        Constant baseline.
    cap : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # compute exponential
    y = c + amplitude * np.exp(-x / decay)

    # apply cap
    y = y.clip(None, cap)

    return y


def symmetric_exponential_peak(x, center, amplitude, decay, c, cap):
    """Symmetric exponential decay peak function.

    Parameters
    ----------
    x : ndarray
        Independent variable.
    center : int or float
        The center of the peak.
    amplitude : float
        Amplitude parameter.
    decay : float
        Decay parameter.
    c : float
        Constant baseline.
    cap : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # locate the center
    ix_cen = bisect_left(x, center)

    # compute left flank
    xl = center - x[:ix_cen]
    yl = c + amplitude * np.exp(-xl / decay)

    # compute right flank
    xr = x[ix_cen:] - center
    yr = c + amplitude * np.exp(-xr / decay)

    # prepare output
    y = np.concatenate([yl, yr])

    # apply cap
    y = y.clip(None, cap)

    return y


def asymmetric_exponential_peak(x, center, amplitude, left_decay, right_decay,
                                c, cap):
    """Asymmetric exponential decay peak function.

    Parameters
    ----------
    x : ndarray
        Independent variable.
    center : int or float
        The center of the peak.
    amplitude : float
        Amplitude parameter.
    left_decay : float
        Decay for left-hand flank.
    right_decay : float
        Decay for right-hand flank.
    c : float
        Constant baseline.
    cap : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # locate the center
    ix_cen = bisect_left(x, center)

    # compute left flank
    xl = center - x[:ix_cen]
    yl = c + amplitude * np.exp(-xl / left_decay)

    # compute right flank
    xr = x[ix_cen:] - center
    yr = c + amplitude * np.exp(-xr / right_decay)

    # prepare output
    y = np.concatenate([yl, yr])

    # apply cap
    y = y.clip(None, cap)

    return y


FitResult = collections.namedtuple(
    'FitResult',
    'delta_aic peak_result null_result loc xx yy best_fit peak residual '
    'peak_start_ix peak_stop_ix peak_start_x peak_stop_x baseline '
    'baseline_stderr'
)


def find_peak_limits(best_fit, baseline, stderr):
    ix_peak_start = ix_peak_stop = None
    for i in range(best_fit.shape[0]):
        v = best_fit[i]
        if ix_peak_start is None:
            if v > baseline + 3 * stderr:
                ix_peak_start = i
        elif ix_peak_stop is None:
            if v < baseline + 3 * stderr:
                ix_peak_stop = i
                break
    return ix_peak_start, ix_peak_stop


# noinspection PyUnresolvedReferences
class PeakFitter(object):
    """Abstract base class for peak fitters."""

    def fit(self, x, y, center, flank):

        # slice out the region of data to fit against
        ix_left = bisect_left(x, center - flank)
        ix_right = bisect_right(x, center + flank)
        loc = slice(ix_left, ix_right)
        xx = x[loc] - center  # make relative to center
        yy = y[loc]

        # fit the null model - remove one outlier
        no_outlier = yy < yy.max()
        null_result = self.null_model.fit(yy[no_outlier], x=xx[no_outlier],
                                          params=self.null_params)

        # fit the peak model
        peak_result = self.peak_model.fit(yy, x=xx, params=self.peak_params)

        # obtain difference in AIC
        delta_aic = null_result.aic - peak_result.aic

        # obtain best fit for peak data for subtracting from signal
        baseline = peak_result.params['c'].value
        baseline_stderr = peak_result.params['c'].stderr
        best_fit = peak_result.best_fit
        peak = best_fit - baseline
        residual = yy - peak

        # figure out the width of the peak
        peak_start_ix = peak_stop_ix = peak_start_x = peak_stop_x = None
        rix_peak_start, rix_peak_stop = find_peak_limits(best_fit, baseline,
                                                         baseline_stderr)
        if rix_peak_start is not None and rix_peak_stop is not None:
            peak_start_ix = ix_left + rix_peak_start
            peak_stop_ix = ix_left + rix_peak_stop
            peak_start_x = xx[rix_peak_start]
            peak_stop_x = xx[rix_peak_stop]

        return FitResult(delta_aic, peak_result, null_result, loc, xx, yy,
                         best_fit, peak, residual, peak_start_ix, peak_stop_ix,
                         peak_start_x, peak_stop_x, baseline, baseline_stderr)


class GaussianPeakFitter(PeakFitter):

    def __init__(self, amplitude, sigma, c):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = c
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = lmfit.models.GaussianModel() + lmfit.models.ConstantModel()
        peak_params = lmfit.Parameters()
        peak_params['center'] = lmfit.Parameter(value=0, vary=False)
        peak_params['amplitude'] = amplitude
        peak_params['sigma'] = sigma
        peak_params['c'] = c
        self.peak_model = peak_model
        self.peak_params = peak_params


class LorentzianPeakFitter(PeakFitter):

    def __init__(self, amplitude, sigma, c):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = c
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = (lmfit.models.LorentzianModel() +
                      lmfit.models.ConstantModel())
        peak_params = lmfit.Parameters()
        peak_params['center'] = lmfit.Parameter(value=0, vary=False)
        peak_params['amplitude'] = amplitude
        peak_params['sigma'] = sigma
        peak_params['c'] = c
        self.peak_model = peak_model
        self.peak_params = peak_params


class SymmetricExponentialPeakFitter(PeakFitter):

    def __init__(self, amplitude, decay, c, cap):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = c
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = lmfit.Model(symmetric_exponential_peak)
        peak_params = lmfit.Parameters()
        peak_params['center'] = lmfit.Parameter(value=0, vary=False)
        peak_params['amplitude'] = amplitude
        peak_params['decay'] = decay
        peak_params['c'] = c
        peak_params['cap'] = cap
        self.peak_model = peak_model
        self.peak_params = peak_params


class AsymmetricExponentialPeakFitter(PeakFitter):

    def __init__(self, amplitude, left_decay, right_decay, c, cap):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = c
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = lmfit.Model(asymmetric_exponential_peak)
        peak_params = lmfit.Parameters()
        peak_params['center'] = lmfit.Parameter(value=0, vary=False)
        peak_params['amplitude'] = amplitude
        peak_params['left_decay'] = left_decay
        peak_params['right_decay'] = right_decay
        peak_params['c'] = c
        peak_params['cap'] = cap
        self.peak_model = peak_model
        self.peak_params = peak_params


class PairExponentialPeakFitter(PeakFitter):

    def __init__(self, amplitude, decay, c, cap):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = c
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = lmfit.Model(exponential)
        peak_params = lmfit.Parameters()
        peak_params['amplitude'] = amplitude
        peak_params['decay'] = decay
        peak_params['c'] = c
        peak_params['cap'] = cap
        self.peak_model = peak_model
        self.peak_params = peak_params

    def fit(self, x, y, center, flank):

        # slice out the region of data to fit against
        ix_left = bisect_left(x, center - flank)
        ix_right = bisect_right(x, center + flank)
        loc = slice(ix_left, ix_right)
        xx = x[loc] - center
        yy = y[loc]

        # split into left and right flanks
        ix_center = bisect_right(xx, 0)
        xl = -xx[:ix_center]
        yl = yy[:ix_center]
        xr = xx[ix_center:]
        yr = yy[ix_center:]

        # check there's some data on both flanks
        if xl.shape[0] > 3 and xr.shape[0] > 3:
            # fit each flank separately

            # find outliers
            no_outlier_l = yl < yl.max()
            no_outlier_r = yr < yr.max()

            # fit the null model - allow one outlier
            null_result_l = self.null_model.fit(yl[no_outlier_l],
                                                x=xl[no_outlier_l],
                                                params=self.null_params)
            null_result_r = self.null_model.fit(yr[no_outlier_r],
                                                x=xr[no_outlier_r],
                                                params=self.null_params)

            # fit the peak model
            peak_result_l = self.peak_model.fit(yl, x=xl,
                                                params=self.peak_params)
            peak_result_r = self.peak_model.fit(yr, x=xr,
                                                params=self.peak_params)

            # obtain difference in AIC
            delta_aic_l = null_result_l.aic - peak_result_l.aic
            delta_aic_r = null_result_r.aic - peak_result_r.aic
            delta_aic = min([delta_aic_l, delta_aic_r])

            # determine baseline
            baseline_l = peak_result_l.params['c'].value
            baseline_r = peak_result_r.params['c'].value
            baseline_stderr_l = peak_result_l.params['c'].stderr
            baseline_stderr_r = peak_result_r.params['c'].stderr
            baseline = max([baseline_l, baseline_r])
            baseline_stderr = max([baseline_stderr_l, baseline_stderr_r])

            # obtain best fit for peak data for subtracting from signal
            best_fit_l = peak_result_l.best_fit
            peak_l = best_fit_l - baseline
            residual_l = yl - peak_l
            best_fit_r = peak_result_r.best_fit
            peak_r = best_fit_r - baseline
            residual_r = yr - peak_r

            # prepare output
            min_result_ix = np.argmin([delta_aic_l, delta_aic_r])
            peak_result = [peak_result_l, peak_result_r][min_result_ix]
            null_result = [null_result_l, null_result_r][min_result_ix]
            best_fit = np.concatenate([best_fit_l, best_fit_r])
            peak = np.concatenate([peak_l, peak_r])
            residual = np.concatenate([residual_l, residual_r])

            # figure out the width of the peak
            peak_start_ix = peak_stop_ix = peak_start_x = peak_stop_x = None
            rix_peak_start, rix_peak_stop = find_peak_limits(best_fit, baseline,
                                                             baseline_stderr)
            if rix_peak_start is not None and rix_peak_stop is not None:
                peak_start_ix = ix_left + rix_peak_start
                peak_stop_ix = ix_left + rix_peak_stop
                peak_start_x = xx[rix_peak_start]
                peak_stop_x = xx[rix_peak_stop]

        else:
            delta_aic = peak_result = null_result = best_fit = peak = \
                residual = peak_start_ix = peak_stop_ix = peak_start_x = \
                peak_stop_x = baseline = baseline_stderr = None

        return FitResult(delta_aic, peak_result, null_result, loc, xx, yy,
                         best_fit, peak, residual, peak_start_ix,
                         peak_stop_ix, peak_start_x, peak_stop_x, baseline,
                         baseline_stderr)


def plot_peak(fit, figsize=(12, 3.5)):
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=figsize)

    ax = axs[0]
    # plot width of the peak
    if fit.peak_start_x and fit.peak_stop_x:
        ax.axvspan(fit.peak_start_x, fit.peak_stop_x, facecolor=palette[0],
                   alpha=.2)
    ax.axvline(0, color='k', lw=1, linestyle='--')
    # plot the fit
    ax.plot(fit.xx, fit.best_fit, lw=1, linestyle='--', color='k')
    # # plot the baseline
    # if fit.baseline is not None:
    #     ax.axhline(fit.baseline, lw=1, linestyle='--', color='k')
    # plot the data
    ax.plot(fit.xx, fit.yy, marker='o', linestyle=' ', markersize=3,
            mfc='none', color=palette[0], mew=1)
    ax.text(.02, .98, r'$\Delta_{i}$ : %.1f' % fit.delta_aic,
            transform=ax.transAxes, ha='left', va='top')
    ax.set_title('Best Fit')
    ax.set_ylabel('Selection statistic')
    ax.set_xlabel('Genetic distance (relative to center)')
    ax.set_ylim(bottom=0)

    ax = axs[1]
    ax.plot(fit.xx, fit.residual, marker='o', linestyle=' ', markersize=3,
            mfc='none', color=palette[0], mew=1)
    if fit.baseline is not None:
        ax.axhline(fit.baseline, lw=1, linestyle='--', color='k')
    ax.set_ylim(axs[0].get_ylim())
    ax.set_title('Residual')
    ax.set_ylabel('Selection statistic')
    ax.set_xlabel('Genetic distance (relative to center)')
    ax.set_ylim(bottom=0)

    ax = axs[2]
    ax.hist(fit.residual.clip(0, None), bins=30)
    if fit.baseline is not None:
        ax.axvline(fit.baseline, lw=1, linestyle='--', color='k')
    ax.set_xlabel('Selection statistic')
    ax.set_ylabel('Frequency')
    ax.set_title('Residual')

    fig.tight_layout()


def scan_fit(x, y, flank, fitter, centers, delta_aics, fits,
             start_index=None, stop_index=None, debug=False):

    # N.B., there may be more centers than data points, because nans must
    # have been removed from data, but we will fit at all windows
    assert not np.any(np.isnan(y))
    assert x.shape == y.shape
    assert (centers.shape == delta_aics.shape == fits.shape)

    # determine region to scan over
    n = centers.shape[0]
    if start_index is None:
        start_index = 0
    if stop_index is None:
        stop_index = n
    assert start_index >= 0
    assert stop_index <= n

    # iterate and fit
    for i in range(start_index, stop_index):
        if debug and i % 100 == 0:
            print('scan progress', i, centers[i])

        # central position to fit at
        center = centers[i]

        # fit the peak
        fit = fitter.fit(x, y, center, flank)

        # store the results
        fits[i] = fit
        if fit.delta_aic is not None:
            delta_aics[i] = fit.delta_aic


def find_peaks(window_starts, window_stops, gpos, signal, flank, fitter,
               min_delta_aic=50, max_iter=20, extend_delta_aic=10,
               max_param_stderr=2, debug=False):
    """TODO"""

    def log(*args):
        if debug:
            print(*args)

    window_starts = np.asarray(window_starts)
    window_stops = np.asarray(window_stops)
    gpos = np.asarray(gpos)
    signal = np.asarray(signal)
    assert (window_starts.shape == window_stops.shape == gpos.shape ==
            signal.shape)
    n = gpos.shape[0]

    # setup working data structures
    delta_aics = np.zeros(n, dtype='f8')
    fits = np.empty(n, dtype=object)

    # strip out missing data
    missing = np.isnan(signal)
    x = gpos[~missing]
    y = signal[~missing]
    starts_nomiss = window_starts[~missing]
    stops_nomiss = window_stops[~missing]
    ppos_nomiss = (starts_nomiss + stops_nomiss) / 2

    # first pass model fits
    scan_fit(x, y, flank=flank, fitter=fitter, centers=gpos,
             delta_aics=delta_aics, fits=fits, debug=debug)

    # keep track of which iteration we're on
    iteration = 0

    # find the first peak
    best_ix = np.argmax(delta_aics)
    best_fit = fits[best_ix]
    best_delta_aic = delta_aics[best_ix]
    log('first pass', best_ix, best_delta_aic)

    while best_delta_aic > min_delta_aic and iteration < max_iter:

        if debug:
            print('*' * 80)
            print('Iteration', iteration)
            print('Peak index:', best_ix)
            print('Delta AIC:', best_delta_aic)
            print('Window:', window_starts[best_ix],
                  window_stops[best_ix])
            print(best_fit.peak_result.fit_report())
            print(best_fit.null_result.fit_report())
            print('*' * 80)

            # plot landscape of delta AIC
            # noinspection PyTypeChecker
            fig, axs = plt.subplots(nrows=2, figsize=(12, 4), sharex=True)
            ax = axs[0]
            ax.set_ylim(bottom=0)
            ax.plot(x, y, marker='o', linestyle=' ', markersize=2, mfc='none')
            ax = axs[1]
            ax.plot(gpos, delta_aics, lw=1)
            ax.axvline(gpos[best_ix], color='k', linestyle='--', lw=2)
            ax.set_ylim(bottom=0)
            ax.set_ylabel(r'$\Delta_{i}$')
            ax.set_xlabel('Position')

            # plot the best peak fit
            plot_peak(best_fit)
            plt.show()

        log('find extent of region under selection')
        hit_start = window_starts[best_ix]
        hit_stop = window_stops[best_ix]
        # search left
        i = best_ix - 1
        while 0 <= i < n:
            if best_delta_aic - fits[i].delta_aic < extend_delta_aic:
                if debug:
                    print('extend hit left', fits[i].delta_aic)
                hit_start = window_starts[i]
                i -= 1
            else:
                break
        # search right
        i = best_ix + 1
        while 0 <= i < n:
            if best_delta_aic - fits[i].delta_aic < extend_delta_aic:
                log('extend hit right', fits[i].delta_aic)
                hit_stop = window_stops[i]
                i += 1
            else:
                break

        if debug:
            print('find flanking region')
        peak_start = peak_stop = None
        if best_fit.peak_start_ix is not None:
            peak_start = starts_nomiss[best_fit.peak_start_ix]
        if best_fit.peak_stop_ix is not None:
            peak_stop = stops_nomiss[best_fit.peak_stop_ix]

        if debug:
            fig, ax = plt.subplots(figsize=(12, 5))
            # plot region under peak
            if (best_fit.peak_start_ix is not None and
                    best_fit.peak_stop_ix is not None):
                ax.axvspan(starts_nomiss[best_fit.peak_start_ix],
                           hit_start, color=palette[0], alpha=.2)
                ax.axvspan(hit_stop,
                           stops_nomiss[best_fit.peak_stop_ix],
                           color=palette[0], alpha=.2)
            # plot hit region
            ax.axvspan(hit_start, hit_stop, color=palette[3], alpha=.5)
            # plot best window
            ax.axvspan(window_starts[best_ix], window_stops[best_ix],
                       color=palette[3], alpha=1)
            # plot fit
            ax.plot(ppos_nomiss[best_fit.loc], best_fit.best_fit,
                    linestyle='--', color='k', lw=1)
            # # plot the baseline
            # if best_fit.baseline is not None:
            #     ax.axhline(best_fit.baseline, lw=1, linestyle='--', color='k')
            # plot data
            ax.plot(ppos_nomiss[best_fit.loc], best_fit.yy, marker='o',
                    linestyle=' ', color=palette[0], mew=1, mfc='none',
                    markersize=3)
            ax.set_xlabel('Physical position')
            ax.set_ylabel('Selection statistic')
            plt.show()

        # filter out hits with too much parameter uncertainty
        params_check = True
        for k in best_fit.peak_result.params:
            p = best_fit.peak_result.params[k]
            if (p.stderr / p.value) > max_param_stderr:
                params_check = False
                log('failed param check: ', p)

        if params_check:
            yield (best_fit, best_delta_aic, best_ix,
                   window_starts[best_ix], window_stops[best_ix],
                   hit_start, hit_stop, peak_start, peak_stop)

        # subtract peak from signal
        y[best_fit.loc] -= best_fit.peak

        # rescan region of the peak
        center = gpos[best_ix]
        ix_rescan_left = bisect_left(gpos, center - (flank * 2))
        ix_rescan_right = bisect_right(gpos, center + (flank * 2))
        log('rescan', ix_rescan_left, ix_rescan_right)
        scan_fit(x, y, flank=flank, fitter=fitter, centers=gpos,
                 delta_aics=delta_aics, fits=fits,
                 start_index=ix_rescan_left, stop_index=ix_rescan_right,
                 debug=debug)

        # find next peak
        best_ix = np.argmax(delta_aics)
        best_fit = fits[best_ix]
        best_delta_aic = delta_aics[best_ix]

        iteration += 1
