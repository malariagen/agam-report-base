# -*- coding: utf-8 -*-
# cython: wraparound=False
# cython: boundscheck=False
from __future__ import absolute_import, print_function, division


cimport numpy as cnp
import numpy as np


def count_gametes(cnp.int8_t[:, :] h):
    """Count the number of each gametic type observed for each pair of variants.
    Observation of all four gametic types for any pair of variants is evidence for
    recombination."""

    cdef:
        Py_ssize_t n, m, i, j, k
        cnp.uint8_t[:, :] d
        cnp.uint32_t[:, :, :, :] count

    n = h.shape[0]
    m = h.shape[1]
    count = np.zeros((n, n, 2, 2), dtype='u4')
    for i in range(n):
        for j in range(i+1, n):
            for k in range(m):
                count[i, j, h[i, k], h[j, k]] += 1

    return np.asarray(count)


import random


cdef tuple _locate_alleles(cnp.int8_t[:, :] h,
                           Py_ssize_t i,
                           cnp.uint8_t[:] cc0,
                           cnp.uint8_t[:] cc1):
    cdef:
        Py_ssize_t nc0, nc1
        cnp.int8_t a
    nc0 = nc1 = 0
    for k in range(h.shape[1]):
        a = h[i, k]
        if a == 0:
            cc0[k] = 1
            nc0 += 1
        elif a == 1:
            cc1[k] = 1
            nc1 += 1
    return nc0, nc1


cdef tuple _count_gametes(c0, c1, p0, p1):
    g00 = c0 & p0
    g01 = c0 & p1
    g10 = c1 & p0
    g11 = c1 & p1
    n_gametes = np.array([np.sum(g00), np.sum(g01), np.sum(g10), np.sum(g11)], dtype='i4')
    return g00, g01, g10, g11, n_gametes


cdef object _locate_recombinants(g00, g01, g10, g11, n_gametes, randomize):
    smallest_cls = np.nonzero(n_gametes == np.min(n_gametes))[0]
    if randomize and len(smallest_cls) > 1:
        smallest_cls = random.choice(smallest_cls)
    else:
        smallest_cls = smallest_cls[0]
    if smallest_cls == 0:
        loc_recombinants = g00
    elif smallest_cls == 1:
        loc_recombinants = g01
    elif smallest_cls == 2:
        loc_recombinants = g10
    else:
        loc_recombinants = g11
    return loc_recombinants


def locate_breakpoints_by_4gametes_opt(haps,
                                       randomize=True,
                                       seed=None):

    cdef:
        Py_ssize_t i, j, k, n, m, nc0, nc1, np0, np1
        cnp.int8_t[:, :] h
        # cnp.int32_t[:] breaks
        cnp.int32_t[:] gametes
        cnp.int32_t[:] n_gametes
        cnp.uint8_t[:] cc0, cc1

    # take a copy of the haplotype data so we can modify
    h = np.array(haps, copy=True, dtype='i1')

    # set up randomization
    if randomize and seed is not None:
        random.seed(seed)

    # setup counting variables
    n = h.shape[0]
    m = h.shape[1]

    # setup breakpoints
    breaks = np.full(m, fill_value=n, dtype='i4')

    # remember which variants are informative for recombination
    unique_splits = set()
    splits = list()

    for i in range(n):

        # find occurrences of reference and alternate alleles
        c0 = np.zeros(m, dtype='bool')
        c1 = np.zeros(m, dtype='bool')
        nc0, nc1 = _locate_alleles(h, i, c0.view('u1'), c1.view('u1'))

        if nc0 + nc1 >= 4:

            if nc0 >= 2 and nc1 >= 2:
                # at least 2 occurrences of both ref and alt alleles, possible to observe 4 gametes

                # search all previous splits
                for p0, p1, np0, np1 in splits:

                    g00, g01, g10, g11, n_gametes = _count_gametes(c0, c1, p0, p1)

                    # apply 4 gamete test
                    if (n_gametes[3] == 0 or
                            n_gametes[2] == 0 or
                            n_gametes[1] == 0 or
                            n_gametes[0] == 0):
                        pass

                    else:

                        # locate recombinants
                        loc_recombinants = _locate_recombinants(g00, g01, g10, g11, n_gametes,
                                                                randomize)

                        # set breakpoints
                        breaks[loc_recombinants] = i

                        # remove recombinants from further analysis
                        np.asarray(h)[i:, loc_recombinants] = -1

                # add to list of splits
                split = tuple(c0), tuple(c1)
                if split not in unique_splits:
                    unique_splits.add(split)
                    splits.append((c0, c1, nc0, nc1))

        else:

            # less than 4 haplotypes remaining, cannot observe 4 gametes, special case

            if nc0 + nc1 == 1:

                # only one haplotype left, break now
                loc_recombinant = c0 | c1
                breaks[loc_recombinant] = i
                break

            elif nc0 == 1 and nc1 == 1:

                # break both in the pair
                loc_recombinants = c0 | c1
                breaks[loc_recombinants] = i
                break

            elif nc0 >= 1 and nc1 >= 1:

                # break haplotype carrying minor allele
                loc_recombinants = c0 if nc0 < nc1 else c1
                if np.any(loc_recombinants):
                    breaks[loc_recombinants] = i
                    np.asarray(h)[i:, loc_recombinants] = -1

    print('foo')
    return np.asarray(breaks)
