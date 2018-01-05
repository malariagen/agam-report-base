# -*- coding: utf-8 -*-
from bisect import bisect_left, bisect_right
import collections
import numpy
import matplotlib.pyplot as plt
import lmfit
import scipy
import seaborn as sns
palette = sns.color_palette()


# workaround PyCharm issue with finding numpy functions
class NumpyModule(object):

    def __getattr__(self, item):
        return getattr(numpy, item)


np = NumpyModule()


def load_values(df, seqid, recmap, genome, spacing=0, values_col='value',
                seqid_col='seqid', start_col='start', end_col='end'):
    """Extract genome-wide selection scan values from a dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        A dataframe with windowed data. Expected to contain columns "chrom" (sequence
        identifier), "start", "stop", as well as `col`. Expect that "start" and "stop" are
        GFF-style 1-based stop-inclusive coords.
    values_col : str
        Name of column containing statistic values.
    seqid : str
        Sequence identifier. May also be a tuple, e.g., ('2R', '2L'), in which
        case the data for each sequence will be concatenated.
    recmap : dict [str -> array]
        Recombination map. A dictionary mapping sequence identifiers onto arrays,
        where each array holds the absolute recombination rate for each base
        in the sequence.
    genome : dict-like [str -> array]
        Genome sequences.
    spacing : int, optional
        Amount of physical distance to insert when concatenating data from
        multiple chromosome arms.
    seqid_col, start_col, end_col : str
        Column names for window sequence identifier, start and end coordinates.

    Returns
    -------
    starts
        Window physical start positions.
    ends
        Window physical end positions.
    gstarts
        Window genetic start positions.
    gends
        Window genetic end positions.
    values
        Statistic values.

    """

    # handle multiple sequences
    if isinstance(seqid, (tuple, list)):
        assert len(seqid) == 2, 'can only concatenate two sequences'
        ((starts1, ends1, gstarts1, gends1, values1, percentiles1),
         (starts2, ends2, gstarts2, gends2, values2, percentiles2)) = \
            [load_values(df, values_col=values_col, seqid=c, recmap=recmap,
                         genome=genome, seqid_col=seqid_col, start_col=start_col,
                         end_col=end_col)
             for c in seqid]
        seq1_plen = len(genome[seqid[0]])
        seq1_glen = np.sum(recmap[seqid[0]])
        starts = np.concatenate([starts1, starts2 + seq1_plen + spacing])
        ends = np.concatenate([ends1, ends2 + seq1_plen + spacing])
        gspacing = recmap[seqid[1]][0] * spacing
        gstarts = np.concatenate([gstarts1, gstarts2 + seq1_glen + gspacing])
        gends = np.concatenate([gends1, gends2 + seq1_glen + gspacing])
        values = np.concatenate([values1, values2])
        percentiles = np.concatenate([percentiles1, percentiles2])
        return starts, ends, gstarts, gends, values, percentiles

    # extract data for single seqid
    df = df.reset_index().set_index(seqid_col)
    df['percentile'] = scipy.stats.rankdata(df[values_col]) / len(df)

    df_seq = df.loc[seqid]

    # extract window starts and ends
    starts = np.asarray(df_seq[start_col])
    ends = np.asarray(df_seq[end_col])

    # compute genetic length of each window
    glen = np.array([
        np.sum(recmap[seqid][int(start - 1):int(end)])
        for start, end in zip(starts, ends)
    ])

    # compute the genetic length position for each window
    gbins = np.append([0], np.cumsum(glen))
    gstarts = gbins[:-1]
    gends = gbins[1:]

    # extract the values column
    values = np.asarray(df_seq[values_col])
    percentiles = np.asarray(df_seq['percentile'])

    return starts, ends, gstarts, gends, values, percentiles


def plot_values(starts, ends, values, figsize=(16, 3), xlabel='Physical position (bp)',
                ax=None):
    """Convenience function to plot values from a windowed selection scan."""

    # determine x coord
    x = (starts + ends) / 2

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.plot(x, values, linestyle=' ', marker='o', mfc='none', markersize=3,
            color=palette[0], mew=1)
    ax.set_xlabel(xlabel)
    ax.set_xlim(starts[0], ends[-1])
    if fig is not None:
        fig.tight_layout()
    return ax


def exponential(x, amplitude, decay, baseline, floor, ceiling):
    """Exponential decay function.

    Parameters
    ----------
    x : ndarray
        Independent variable.
    amplitude : float
        Amplitude parameter.
    decay : float
        Decay parameter.
    baseline : float
        Baseline parameter.
    floor : float
        Minimum value that the result can take.
    ceiling : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # compute exponential
    y = baseline + amplitude * np.exp(-x / decay)

    # apply cap
    y = y.clip(floor, ceiling)

    return y


def symmetric_exponential_peak(x, center, amplitude, decay, baseline, floor, ceiling):
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
    baseline : float
        Baseline parameter.
    floor : float
        Minimum value that the result can take.
    ceiling : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # locate the center
    ix_cen = bisect_right(x, center)

    # compute left flank
    xl = center - x[:ix_cen]
    yl = baseline + amplitude * np.exp(-xl / decay)

    # compute right flank
    xr = x[ix_cen:] - center
    yr = baseline + amplitude * np.exp(-xr / decay)

    # prepare output
    y = np.concatenate([yl, yr])

    # apply cap
    y = y.clip(floor, ceiling)

    return y


def asymmetric_exponential_peak(x, center, amplitude, decay_left, decay_right,
                                baseline, floor, ceiling):
    """Asymmetric exponential decay peak function.

    Parameters
    ----------
    x : ndarray
        Independent variable.
    center : int or float
        The center of the peak.
    amplitude : float
        Amplitude parameter.
    decay_left : float
        Decay parameter (left flank).
    decay_right : float
        Decay parameter (right flank).
    baseline : float
        Baseline parameter.
    floor : float
        Minimum value that the result can take.
    ceiling : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # locate the center
    ix_cen = bisect_right(x, center)

    # compute left flank
    xl = center - x[:ix_cen]
    yl = baseline + amplitude * np.exp(-xl / decay_left)

    # compute right flank
    xr = x[ix_cen:] - center
    yr = baseline + amplitude * np.exp(-xr / decay_right)

    # prepare output
    y = np.concatenate([yl, yr])

    # apply cap
    y = y.clip(floor, ceiling)

    return y


def uncoupled_exponential_peak(x, center, amplitude_left, amplitude_right, decay_left,
                               decay_right, baseline_left, baseline_right, floor_left,
                               floor_right, ceiling_left, ceiling_right):
    """Exponential decay peak function where flanks are fitted completely
    independently, ie.., no shared parameters.

    Parameters
    ----------
    x : ndarray
        Independent variable.
    center : int or float
        The center of the peak.
    amplitude_left : float
        Amplitude parameter (left flank).
    amplitude_right : float
        Amplitude parameter (right flank).
    decay_left : float
        Decay parameter (left flank).
    decay_right : float
        Decay parameter (right flank).
    baseline_left : float
        Baseline parameter (left flank).
    baseline_right : float
        Baseline parameter (right flank).
    floor_left : float
        Minimum value that the result can take.
    floor_right : float
        Minimum value that the result can take.
    ceiling_left : float
        Maximum value that the result can take.
    ceiling_right : float
        Maximum value that the result can take.

    Returns
    -------
    y : ndarray

    """

    # locate the center
    ix_cen = bisect_right(x, center)

    # compute left flank
    xl = center - x[:ix_cen]
    yl = baseline_left + amplitude_left * np.exp(-xl / decay_left)
    yl = yl.clip(floor_left, ceiling_left)

    # compute right flank
    xr = x[ix_cen:] - center
    yr = baseline_right + amplitude_right * np.exp(-xr / decay_right)
    yr = yr.clip(floor_right, ceiling_right)

    # prepare output
    y = np.concatenate([yl, yr])

    return y


PeakFitResult = collections.namedtuple(
    'PeakFitResult',
    'null_result peak_result null_result_left null_result_right peak_result_left '
    'peak_result_right xx yy loc best_fit peak residual peak_start_ix peak_end_ix '
    'peak_start_x peak_end_x baseline baseline_stderr delta_aic delta_aic_left '
    'delta_aic_right'
)


def find_peak_limits(best_fit, threshold):
    ix_peak_start = ix_peak_end = None
    # work forward to find peak start
    for i in range(best_fit.shape[0]):
        v = best_fit[i]
        if ix_peak_start is None:
            if v > threshold:
                ix_peak_start = i
    # work backwards to find peak end
    for i in range(best_fit.shape[0] - 1, -1, -1):
        v = best_fit[i]
        if ix_peak_end is None:
            if v > threshold:
                ix_peak_end = i
                break
    assert ix_peak_start is not None
    assert ix_peak_end is not None
    return ix_peak_start, ix_peak_end


class SymmetricExponentialPeakFitter(object):

    def __init__(self, amplitude, decay, baseline, floor, ceiling, null):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = null
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = lmfit.Model(symmetric_exponential_peak)
        peak_params = lmfit.Parameters()
        peak_params['center'] = lmfit.Parameter(value=0, vary=False)
        peak_params['amplitude'] = amplitude
        peak_params['decay'] = decay
        peak_params['baseline'] = baseline
        peak_params['ceiling'] = ceiling
        peak_params['floor'] = floor
        self.peak_model = peak_model
        self.peak_params = peak_params

        # initialise flank model
        self.flank_model = lmfit.Model(exponential)

    def fit(self, x, y, center, flank):

        # slice out the region of data to fit against
        ix_left = bisect_left(x, center - flank)
        ix_right = bisect_right(x, center + flank)
        loc = slice(ix_left, ix_right)
        xx = x[loc] - center  # make relative to center
        yy = y[loc]

        # fit the null model - allow one outlier
        no_outlier = yy < yy.max()
        null_result = self.null_model.fit(yy[no_outlier], x=xx[no_outlier],
                                          params=self.null_params)

        # fit the peak model
        peak_result = self.peak_model.fit(yy, x=xx, params=self.peak_params)

        # obtain best fit for peak data for subtracting from values
        baseline = peak_result.params['baseline'].value
        baseline_stderr = peak_result.params['baseline'].stderr
        best_fit = peak_result.best_fit
        peak = best_fit - baseline
        residual = yy - peak

        # figure out the width of the peak
        rix_peak_start, rix_peak_end = find_peak_limits(best_fit,
                                                        (baseline + 3 * baseline_stderr))
        peak_start_ix = ix_left + rix_peak_start
        peak_end_ix = ix_left + rix_peak_end
        peak_start_x = xx[rix_peak_start]
        peak_end_x = xx[rix_peak_end]

        # split data into two flanks
        ix_cen = bisect_right(xx, 0)
        xx_left = -xx[:ix_cen]
        yy_left = yy[:ix_cen]
        xx_right = xx[ix_cen:]
        yy_right = yy[ix_cen:]

        # obtain fit for left and right flanks
        flank_params = lmfit.Parameters()
        amplitude = peak_result.params['amplitude'].value
        decay = peak_result.params['decay'].value
        baseline = peak_result.params['baseline'].value
        floor = peak_result.params['floor'].value
        ceiling = peak_result.params['ceiling'].value
        flank_params['amplitude'] = lmfit.Parameter(value=amplitude, vary=False)
        flank_params['decay'] = lmfit.Parameter(value=decay, vary=False)
        flank_params['baseline'] = lmfit.Parameter(value=baseline, vary=False)
        flank_params['floor'] = lmfit.Parameter(value=floor, vary=False)
        flank_params['ceiling'] = lmfit.Parameter(value=ceiling, vary=False)
        peak_result_left = self.flank_model.fit(yy_left, x=xx_left, params=flank_params)
        null_result_left = self.null_model.fit(yy_left, x=xx_left, params=self.null_params)
        peak_result_right = self.flank_model.fit(yy_right, x=xx_right,
                                                 params=flank_params)
        null_result_right = self.null_model.fit(yy_right, x=xx_right,
                                                params=self.null_params)

        # obtain differences in AIC
        delta_aic = null_result.aic - peak_result.aic
        delta_aic_left = null_result_left.aic - peak_result_left.aic
        delta_aic_right = null_result_right.aic - peak_result_right.aic

        return PeakFitResult(null_result=null_result,
                             peak_result=peak_result,
                             null_result_left=null_result_left,
                             null_result_right=null_result_right,
                             peak_result_left=peak_result_left,
                             peak_result_right=peak_result_right,
                             xx=xx, yy=yy, loc=loc, best_fit=best_fit,
                             baseline=baseline, baseline_stderr=baseline_stderr,
                             peak=peak, peak_start_ix=peak_start_ix,
                             peak_end_ix=peak_end_ix, peak_start_x=peak_start_x,
                             peak_end_x=peak_end_x, residual=residual,
                             delta_aic=delta_aic, delta_aic_left=delta_aic_left,
                             delta_aic_right=delta_aic_right)


class AsymmetricExponentialPeakFitter(object):

    def __init__(self, amplitude, decay_left, decay_right, baseline, floor, ceiling,
                 null):

        # initialise null model
        null_model = lmfit.models.ConstantModel()
        null_params = lmfit.Parameters()
        null_params['c'] = null
        self.null_model = null_model
        self.null_params = null_params

        # initialise peak model
        peak_model = lmfit.Model(asymmetric_exponential_peak)
        peak_params = lmfit.Parameters()
        peak_params['center'] = lmfit.Parameter(value=0, vary=False)
        peak_params['amplitude'] = amplitude
        peak_params['decay_left'] = decay_left
        peak_params['decay_right'] = decay_right
        peak_params['baseline'] = baseline
        peak_params['ceiling'] = ceiling
        peak_params['floor'] = floor
        self.peak_model = peak_model
        self.peak_params = peak_params

        # initialise flank model
        self.flank_model = lmfit.Model(exponential)

    def fit(self, x, y, center, flank):

        # slice out the region of data to fit against
        ix_left = bisect_left(x, center - flank)
        ix_right = bisect_right(x, center + flank)
        loc = slice(ix_left, ix_right)
        xx = x[loc] - center  # make relative to center
        yy = y[loc]

        # fit the null model - allow one outlier
        no_outlier = yy < yy.max()
        null_result = self.null_model.fit(yy[no_outlier], x=xx[no_outlier],
                                          params=self.null_params)

        # fit the peak model
        peak_result = self.peak_model.fit(yy, x=xx, params=self.peak_params)

        # obtain best fit for peak data for subtracting from values
        baseline = peak_result.params['baseline'].value
        baseline_stderr = peak_result.params['baseline'].stderr
        best_fit = peak_result.best_fit
        peak = best_fit - baseline
        residual = yy - peak

        # figure out the width of the peak
        rix_peak_start, rix_peak_end = find_peak_limits(best_fit,
                                                        (baseline + 3 * baseline_stderr))
        peak_start_ix = ix_left + rix_peak_start
        peak_end_ix = ix_left + rix_peak_end
        peak_start_x = xx[rix_peak_start]
        peak_end_x = xx[rix_peak_end]

        # split data into two flanks
        ix_cen = bisect_right(xx, 0)
        xx_left = -xx[:ix_cen]
        yy_left = yy[:ix_cen]
        xx_right = xx[ix_cen:]
        yy_right = yy[ix_cen:]

        # obtain fit for left and right flanks
        amplitude = peak_result.params['amplitude'].value
        decay_left = peak_result.params['decay_left'].value
        decay_right = peak_result.params['decay_right'].value
        baseline = peak_result.params['baseline'].value
        floor = peak_result.params['floor'].value
        ceiling = peak_result.params['ceiling'].value

        flank_params = lmfit.Parameters()
        flank_params['amplitude'] = lmfit.Parameter(value=amplitude, vary=False)
        flank_params['decay'] = lmfit.Parameter(value=decay_left, vary=False)
        flank_params['baseline'] = lmfit.Parameter(value=baseline, vary=False)
        flank_params['floor'] = lmfit.Parameter(value=floor, vary=False)
        flank_params['ceiling'] = lmfit.Parameter(value=ceiling, vary=False)
        peak_result_left = self.flank_model.fit(yy_left, x=xx_left, params=flank_params)
        null_result_left = self.null_model.fit(yy_left, x=xx_left, params=self.null_params)
        flank_params = lmfit.Parameters()
        flank_params['amplitude'] = lmfit.Parameter(value=amplitude, vary=False)
        flank_params['decay'] = lmfit.Parameter(value=decay_right, vary=False)
        flank_params['baseline'] = lmfit.Parameter(value=baseline, vary=False)
        flank_params['floor'] = lmfit.Parameter(value=floor, vary=False)
        flank_params['ceiling'] = lmfit.Parameter(value=ceiling, vary=False)
        peak_result_right = self.flank_model.fit(yy_right, x=xx_right,
                                                 params=flank_params)
        null_result_right = self.null_model.fit(yy_right, x=xx_right,
                                                params=self.null_params)

        # obtain differences in AIC
        delta_aic = null_result.aic - peak_result.aic
        delta_aic_left = null_result_left.aic - peak_result_left.aic
        delta_aic_right = null_result_right.aic - peak_result_right.aic

        return PeakFitResult(null_result=null_result,
                             peak_result=peak_result,
                             null_result_left=null_result_left,
                             null_result_right=null_result_right,
                             peak_result_left=peak_result_left,
                             peak_result_right=peak_result_right,
                             xx=xx, yy=yy, loc=loc, best_fit=best_fit,
                             baseline=baseline, baseline_stderr=baseline_stderr,
                             peak=peak, peak_start_ix=peak_start_ix,
                             peak_end_ix=peak_end_ix, peak_start_x=peak_start_x,
                             peak_end_x=peak_end_x, residual=residual,
                             delta_aic=delta_aic, delta_aic_left=delta_aic_left,
                             delta_aic_right=delta_aic_right)


def plot_peak_fit(fit, figsize=(8, 2.5), iter_out_dir=None,
                  xlabel='Genetic distance (cM)', ylabel='Selection statistic', dpi=None):
    # noinspection PyTypeChecker
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=figsize, facecolor='w', dpi=dpi)

    ax = axs[0]
    # plot width of the peak
    if fit.peak_start_x and fit.peak_end_x:
        ax.axvspan(fit.peak_start_x, fit.peak_end_x, facecolor=palette[0],
                   alpha=.2)
    # ax.axvline(0, color='k', lw=.5, linestyle='--')
    # plot the fit
    ax.plot(fit.xx, fit.best_fit, lw=.5, linestyle='--', color='k')
    # # plot the baseline
    # if fit.baseline is not None:
    #     ax.axhline(fit.baseline, lw=1, linestyle='--', color='k')
    # plot the data
    ax.plot(fit.xx, fit.yy, marker='o', linestyle=' ', markersize=3,
            mfc='none', color=palette[0], mew=.5)
    ax.text(.02, .98, r'$\Delta_{i}$ : %.1f' % fit.delta_aic_left,
            transform=ax.transAxes, ha='left', va='top')
    ax.text(.98, .98, r'$\Delta_{i}$ : %.1f' % fit.delta_aic_right,
            transform=ax.transAxes, ha='right', va='top')
    ax.text(.5, 1.02, r'$\Delta_{i}$ : %.1f' % fit.delta_aic,
            transform=ax.transAxes, ha='center', va='bottom')
    ax.set_title('Peak fit', loc='left')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_ylim(bottom=0)

    ax = axs[1]
    if fit.baseline is not None:
        ax.axhline(fit.baseline, lw=.5, linestyle='--', color='k')
    ax.plot(fit.xx, fit.residual, marker='o', linestyle=' ', markersize=3,
            mfc='none', color=palette[0], mew=.5)
    ax.set_ylim(axs[0].get_ylim())
    ax.set_title('Residual', loc='left')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_ylim(bottom=0)

    ax = axs[2]
    ax.hist(fit.residual, bins=30)
    # if fit.baseline is not None:
    #     ax.axvline(fit.baseline, lw=.5, linestyle='--', color='k')
    ax.set_xlabel(ylabel)
    ax.set_ylabel('Frequency')
    ax.set_title('Residual', loc='left')

    fig.tight_layout()

    if iter_out_dir:
        fig.savefig(os.path.join(iter_out_dir, 'peak_fit.png'),
                    bbox_inches='tight', dpi=150, facecolor='w')
        plt.close()
