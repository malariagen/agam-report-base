# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import numpy as np
from numpy.testing import assert_array_equal


import hapclust


def test_locate_breakpoints_by_4gametes():

    for reverse in False, True:

        def f(haps, **kwargs):
            kwargs.setdefault('randomize', False)
            return hapclust.locate_breakpoints_by_4gametes(haps, reverse=reverse, **kwargs)

        h = [[0, 0, 1, 1],
             [0, 1, 0, 1]]
        ex = [1, 2, 2, 2]
        assert_array_equal(ex, f(h))

        h = [[0, 0, 1, 1],
             [0, 0, 0, 0],
             [0, 1, 0, 1]]
        ex = [2, 3, 3, 3]
        assert_array_equal(ex, f(h))

        h = [[0, 0, 0, 0],
             [0, 0, 1, 1],
             [0, 1, 0, 1]]
        ex = [2, 3, 3, 3]
        assert_array_equal(ex, f(h))

        # handle less than 4 haplotypes left
        h = [[0, 0, 1, 1],
             [0, 1, 0, 1],
             [0, 0, 0, 0],
             [1, 0, 0, 0],
             [0, 1, 0, 0],
             [0, 1, 0, 0],
             [0, 0, 1, 0],
             [0, 0, 0, 0]]
        ex = [1, 4, 6, 6]
        assert_array_equal(ex, f(h))

        # test with p_parsimony
        h = [[0, 0, 1, 1],
             [0, 1, 0, 1],
             [0, 1, 0, 1]]
        ex = [1, 3, 2, 3]
        assert_array_equal(ex, f(h, p_parsimony=[1, 1, 1]))
        ex = [3, 3, 3, 3]
        assert_array_equal(ex, f(h, p_parsimony=[0, 0, 0]))
        ex = [2, 3, 3, 3]
        assert_array_equal(ex, f(h, p_parsimony=[1, 0, 1]))

#     h = [[0, 0, 0]]
#     ex = [0, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0],
#          [0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[1, 1, 1],
#          [1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0],
#          [1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1],
#          [0, 0, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1],
#          [1, 1, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 1, 1],
#          [1, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 1, 1],
#          [0, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 1, 1],
#          [0, 0, 1]]
#     assert_array_equal(ex, f(h))
#
#     # N.B., singletons can always be rearranged and forwarded to the end of the chain, and so do
#     # not imply any informative split.
#     h = [[1, 0, 0],
#          [1, 1, 0]]
#     assert_array_equal(ex, f(h))
#
#
# def test_locate_breakpoints_by_parsimony_3haps_recom():
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     # There is evidence for recombination in all of the following cases.
#
#     # Derived allele (1) at second variant occurs with both ancestral and derived allele at first
#     # variant, implies recombination in either first or second haplotype.
#     h = [[0, 1, 1],
#          [1, 1, 0]]
#     #     R
#     ex = [1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0],
#          [0, 1, 1],
#          [1, 1, 0]]
#     #     R
#     ex = [2, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0],
#          [0, 1, 1],
#          [1, 1, 0]]
#     #     R
#     ex = [2, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0],
#          [0, 1, 1],
#          [1, 1, 1],
#          [1, 1, 0]]
#     #     R
#     ex = [3, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 1, 1],
#          [0, 0, 0],
#          [1, 1, 1],
#          [1, 1, 0]]
#     #     R
#     ex = [3, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     # change order of haplotypes carrying the mutation - here second or third haplotypes could be
#     #  the recombinant
#     h = [[1, 1, 0],
#          [0, 1, 1]]
#     #        R
#     ex = [0, 1, 0]
#     assert_array_equal(ex, f(h))
#
#     # # break down all the way
#     # TODO figure out what happens when only two singleton clusters remain
#     # h = [[1, 1, 0],
#     #      [0, 1, 1],
#     #      #   R
#     #      [1, 0, 1],
#     #      #R
#     #
#     #      ]
#
#
# def test_locate_breakpoints_by_parsimony_4haps_norec():
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     # There is no evidence for recombination in any of the following cases. Expected breakpoints,
#     # (0) implies no breakpoint found:
#     ex = [0, 0, 0, 0]
#
#     h = [[0, 0, 0, 0],
#          [1, 1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1, 1],
#          [0, 0, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1, 1],
#          [1, 1, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[1, 0, 0, 0],
#          [1, 1, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[1, 0, 0, 0],
#          [1, 1, 1, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0],
#          [0, 0, 0, 0],
#          [1, 1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1, 1],
#          [1, 1, 0, 0],
#          [0, 0, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1],
#          [0, 0, 1, 0],
#          [0, 1, 0, 0],
#          [0, 1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     # break down all the way
#     h = [[0, 0, 1, 1],
#          [0, 1, 1, 0],
#          #   R
#          [1, 0, 1, 1],  # fixed in all remaining lineages, => no recombination
#          [1, 0, 0, 1],  # see three gametes => recombination
#          #R
#          [0, 0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#
# def test_locate_breakpoints_by_parsimony_4haps_recom():
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     # There is evidence for recombination in all of the following cases.
#
#     h = [[0, 0, 1, 1],
#          [0, 1, 1, 0]]
#     #        R
#     ex = [0, 1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[1, 1, 0, 0],
#          [0, 1, 1, 0]]
#     #        R
#     ex = [0, 1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1, 1],
#          [0, 1, 1, 1]]
#     #        R
#     ex = [0, 1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     # N.B., smallest number of recombinations
#     h = [[0, 0, 1, 1],
#          [1, 1, 1, 0]]
#     #           R
#     ex = [0, 0, 1, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 1, 1],
#          [0, 0, 0, 0],
#          [0, 1, 1, 0]]
#     #        R
#     ex = [0, 2, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     # N.B., smallest number of recombinations
#     h = [[1, 1, 0, 0],
#          [0, 0, 1, 1],
#          [1, 1, 1, 0]]
#     #           R
#     ex = [0, 0, 2, 0]
#     assert_array_equal(ex, f(h))
#
#     # N.B., smallest number of recombinations
#     h = [[0, 0, 1, 1],
#          [1, 1, 0, 0],
#          [1, 1, 1, 0]]
#     #           R
#     ex = [0, 0, 2, 0]
#     assert_array_equal(ex, f(h))
#
#
# def test_locate_breakpoints_by_parsimony_5haps_norec():
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     # There is no evidence for recombination in any of the following cases. Expected breakpoints,
#     # (0) implies no breakpoint found:
#     ex = [0, 0, 0, 0, 0]
#
#     h = [[0, 0, 0, 0, 0],
#          [1, 1, 1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1],
#          [0, 0, 0, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1],
#          [1, 1, 0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     # N.B., singletons can be forwarded
#     h = [[0, 0, 0, 1, 1],
#          [1, 0, 0, 0, 0],
#          [1, 1, 1, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     # N.B., singletons can be forwarded
#     h = [[1, 1, 1, 0, 0],
#          [1, 0, 0, 0, 0],
#          [1, 1, 0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#
# def test_locate_breakpoints_by_parsimony_5haps_recom():  # flake8: noqa
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     # There is evidence for recombination in all of the following cases.
#
#     h = [[0, 0, 0, 1, 1],
#          [0, 0, 1, 1, 0]]
#     #           R
#     ex = [0, 0, 1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1],
#          [0, 0, 0, 0, 0],
#          [0, 0, 1, 1, 0]]
#     #           R
#     ex = [0, 0, 2, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1],
#          [0, 1, 1, 0, 0],
#          [0, 0, 1, 1, 0]]
#     #           R
#     ex = [0, 0, 2, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     # N.B., smallest number of recombinations
#     h = [[0, 0, 0, 1, 1],
#          [0, 1, 1, 1, 0]]
#     #              R
#     ex = [0, 0, 0, 1, 0]
#     assert_array_equal(ex, f(h))
#
#     # N.B., smallest number of recombinations
#     h = [[0, 0, 0, 1, 1],
#          [1, 1, 1, 1, 1],
#          [0, 1, 1, 1, 0]]
#     #              R
#     ex = [0, 0, 0, 2, 0]
#     assert_array_equal(ex, f(h))
#
#     # Break down all clusters. N.B., when there is only a single cluster with a single haplotype
#     # remaining, have to terminate, otherwise would continue on forever.
#     h = [[0, 0, 0, 1, 1],
#          [0, 1, 1, 1, 0],
#          #         R
#          [1, 1, 1, 0, 0],
#          #R
#          [0, 0, 0, 0, 1],
#          [0, 0, 1, 0, 1],
#          #      R
#          [1, 0, 1, 0, 0],
#          [0, 1, 0, 0, 1],
#          #   R        R
#          ]
#     ex = [2, 6, 4, 1, 6]
#     assert_array_equal(ex, f(h))
#
#
# def test_locate_breakpoints_by_parsimony_7haps_norec():  # flake8: noqa
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     # There is no evidence for recombination in any of the following cases. Expected breakpoints,
#     # (0) implies no breakpoint found:
#     ex = [0, 0, 0, 0, 0]
#
#     h = [[0, 0, 0, 0, 1, 1, 1],
#          [0, 0, 0, 0, 1, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0, 1, 1, 1],
#          [0, 0, 0, 0, 0, 1, 1]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0, 1, 1, 1],
#          [1, 1, 0, 0, 0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0, 0, 1, 1],
#          [1, 1, 1, 0, 0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0, 1, 1, 1],
#          [1, 1, 1, 0, 0, 0, 0]]
#     assert_array_equal(ex, f(h))
#
#
# def test_locate_breakpoints_by_parsimony_7haps_recom():  # flake8: noqa
#
#     def f(haps):
#         # To get predictable results for testing, treat "1" as the derived allele, and alyaws
#         # pick the first cluster to recombine if two clusters carry equal numbers of mutations.
#         return locate_breakpoints_by_parsimony(haps, polarized=True, randomize=False)
#
#     h = [[0, 0, 0, 0, 0, 1, 1],
#          [0, 0, 0, 0, 1, 1, 1]]
#     #                 R
#     ex = [0, 0, 0, 0, 1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0, 0, 1, 1],
#          [0, 0, 0, 1, 1, 1, 0]]
#     #                    R
#     ex = [0, 0, 0, 0, 0, 1, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 0, 1, 1, 1],
#          [0, 0, 0, 1, 1, 0, 0]]
#     #              R
#     ex = [0, 0, 0, 1, 0, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1, 1, 1],
#          [0, 1, 1, 1, 1, 0, 0]]
#     #        R  R
#     ex = [0, 1, 1, 0, 0, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1, 1, 1],
#          [0, 1, 1, 1, 1, 1, 0]]
#     #        R  R
#     ex = [0, 1, 1, 0, 0, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1, 1, 1],
#          [1, 1, 1, 1, 1, 0, 0]]
#     #              R  R
#     ex = [0, 0, 0, 1, 1, 0, 0]
#     assert_array_equal(ex, f(h))
#
#     h = [[0, 0, 0, 1, 1, 1, 1],
#          [1, 1, 1, 1, 1, 1, 0]]
#     #     R  R  R
#     ex = [1, 1, 1, 0, 0, 0, 0]
#     assert_array_equal(ex, f(h))
