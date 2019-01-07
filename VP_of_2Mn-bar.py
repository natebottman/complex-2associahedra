# =============================================================================
# A somewhat-brute-force approach to computing the
# virtual Poincare polynomial of 2Mn-bar.
#
# Because of the fiber-product aspect, what we
# compute in general is the VP of a fiber product
# of 2Mn-bars, for different n's of the same length.
#
# Now it's memoized!
# =============================================================================

import sympy as sp
from itertools import product
import copy
import warnings as w

x = sp.symbols('x')


# =============================================================================
# These 6 functions allow us to produce partitions and 2-partitions
# of various kinds.
# =============================================================================


def binary_seqs(s):
    # Returns list of binary sequences of length s, for use in parts

    if s == 0:
        return [[]]
    elif s >= 1:
        lists = [[[i] + seq for seq in binary_seqs(s-1)] for i in [0, 1]]
        return lists[0]+lists[1]


def parts(nums):
    # Returns all partitions of a list

    if len(nums) == 0:
        return [[]]
    elif len(nums) >= 1:
        new_parts = []
        for part in parts(nums[:-1]):
            new_parts += [part + [[nums[-1]]]]
            temp_part = copy.deepcopy(part)
            for i in range(len(temp_part)):
                new_parts += [[group if j != i else (group+[nums[-1]])
                               for (j, group) in enumerate(temp_part)]]
        return new_parts


def parts_rf(nums):

    return [tuple(tuple(gp) for gp in part) for part in parts(nums)]


def twoparts_fused(twonums):
    # Idea: consider a single witch ball, and suppose that all seams fuse.
    # How can marked points distribute themselves?
    #
    # Takes in a dictionary of the form eg {1:[a,b],2:[c]}, where the
    # keys are the seam labels and the values are the marked points.
    # Returns a list of lists of dicts, eg [{1:[a],2:[]},{1:[b],2:[c]}].

    nums = twonums.keys()
    num_pts = sum([len(val) for val in twonums.values()])
    if num_pts == 0:
        return [[]]
    elif num_pts >= 1:
        i = [j for j in nums if len(twonums[j]) >= 1][0]
        pt = twonums[i][-1]
        smaller = {k: (v[:-1] if k == i else v)
                   for k, v in twonums.items()}
        small_twoparts = twoparts_fused(smaller)
        return [[{k: (v+[pt] if (k, bub) == (i, bub0) else v)
                  for k, v in bub.items()}
                 for bub in twopart]
                for twopart in small_twoparts for bub0 in twopart] \
            + [twopart + [{k: ([pt] if k == i else [])
                          for k in nums}]
               for twopart in small_twoparts]


def twoparts_fused_rf(twonums):

    nums = twonums.keys()
    return [{tuple(nums): twopart}
            for twopart in twoparts_fused(twonums)]


def twoparts_fixed_part(twonums, part):

    if set(twonums.keys()) != set(elt for gp in part for elt in gp):
        raise Exception('second argument is not a partition of the seams')

    twoparts = list(product(*[twoparts_fused_rf({k: twonums[k] for k in gp})
                              for gp in part]))

    return [{k: v for d in twopart for k, v in d.items()}
            for twopart in twoparts]


# =============================================================================
# These 2 functions compute the VP's of (a) ordered configurations of points
# in k-times-punctured C, and (b) the interior of a fiber product of 2Mn-bars.
# =============================================================================


def VP_OC_C(k, m=0):
    # Returns the VP of the ordered configuration space of k points in C
    # minus m points.

    return sp.expand(sp.prod([x**2-i for i in range(m, k+m)]))


def VP_fiber_prod_interior(m):
    # Returns the VP of the fiber product of 2-associahedra corresponding
    # to m = [m1,...,ms]. All mi's must be of the same length.

    s, k = len(m[0]), len(m)

    if s == 1:
        return sp.expand(sp.prod([VP_OC_C(m[i][0]-2, 2) for i in range(k)]))

    elif s >= 2:
        a = [sorted(mi)[::-1] for mi in m]
        poly = VP_OC_C(s-2, 2) \
            * sp.prod([VP_OC_C(a[i][0]-1, 1) for i in range(k)]) \
            * sp.prod([sp.prod([VP_OC_C(a[i][j], 0) for j in range(1, s)])
                       for i in range(k)])
        return sp.expand(poly)


# =============================================================================
# These 3 functions are helpers for the Repo class.
# =============================================================================


def pts_from_nums(m, s):

    if any([len(mi) != s for mi in m]):
        raise Exception("not all mi's have length s")

    nums = list(range(1, s+1))
    twonums_seq = [{i+1: [(a+1, i+1, j+1) for j in range(m[a][i])]
                    for i in range(s)} for a in range(len(m))]
    return [twonums_seq, nums]


def all_n(dim):
    # Recursively generates all lists such that the corresponding
    # 2-associahedron has dimension dim. We return only
    # weakly-increasing sequences.

    if dim == 0:
        return [[2], [1, 0]]

    elif dim >= 1:
        old_inc = [[n[:k] + [n[k]+1] + n[k+1:] for k in range(0, len(n))]
                   for n in all_n(dim-1)]
        old_inc = [n for l in old_inc for n in l]
        only_one_pt = [[0]*k + [1] + [0]*(dim+1-k) for k in range(0, dim+2)]
        with_dup = old_inc + only_one_pt
        return [sorted(list(n)) for n in set([tuple(n) for n in with_dup])]


def reformat_n(n):

    if len(n) == 0:
        return ()
    elif len(n) == 1:
        return (tuple(sorted(list(n[0]))),)
    elif len(n[0]) == 2:
        return tuple(sorted(list(tuple(sorted(list(ni))) for ni in n)))
    else:
        return tuple(sorted(list(n)))


# =============================================================================
# The Repo class is where we finally compute VP's of fiber product of 2Mn-bars.
# =============================================================================


class Repo:

    def __init__(self):
        self.ass = {1: 1}
        self.twoass = {}
        self.twoass_roots_unfused = {}

    def return_all(self):
        return [self.ass,
                self.twoass,
                self.twoass_roots_unfused]

    def get_ass(self, s):
        if s not in self.ass:
            self.update_ass(s)
        return self.ass[s]

    def update_ass(self, s):
        if s in self.ass:
            w.warn("an unnecessary computation is being done")
        ps = parts(list(range(1, s+1)))[:-1]
        poly = sum([VP_OC_C(len(p)-2, 2) *
                    sp.prod([self.get_ass(len(c)) for c in p]) for p in ps])
        self.ass[s] = sp.expand(poly)

    def get_2ass(self, n, s):

        m = reformat_n(n)
        if (m, s) not in self.twoass:
            self.update_2ass(m, s)
        return self.twoass[(m, s)]

    def update_2ass(self, n, s):

        m = reformat_n(n)
        if any([len(mi) != s for mi in m]):
            raise Exception('one of the mi\'s does not have length s')
        if (m, s) in self.twoass.keys():
            w.warn("an unnecessary computation is being done")

        if len(m) == 0:
            self.twoass[((), s)] = self.get_ass(s)
        elif s == 1:
            self.twoass[(m, 1)] = sp.prod([self.get_ass(mi[0]) for mi in m])
        else:
            twonums_seq, nums = pts_from_nums(m, s)
            P = 0
            for twopart_seq in product(*[twoparts_fused(twonums)
                                         for twonums in twonums_seq]):
                P += sp.prod([self.get_ass(len(twopart))
                              for twopart in twopart_seq]) \
                 * self.get_2ass_roots_unfused(tuple(tuple(len(bub[s])
                                                           for s in bub.keys())
                                                     for twopart in twopart_seq
                                                     for bub in twopart))
            self.twoass[(m, s)] = sp.expand(P)

    def get_2ass_roots_unfused(self, n):

        m = reformat_n(n)
        if m not in self.twoass_roots_unfused:
            self.update_2ass_roots_unfused(m)
        return self.twoass_roots_unfused[m]

    def update_2ass_roots_unfused(self, n):

        m = reformat_n(n)
        s = len(m[0])
        if any([len(mi) != s for mi in m]):
            raise Exception('the mi\'s are not all of the same length')
        if not (s >= 2):
            raise Exception('should be at least two seams')
        if m in self.twoass_roots_unfused:
            w.warn("an unnecessary computation is being done")

        P = 0
        twonums_seq, nums = pts_from_nums(m, s)
        for part in parts_rf(nums)[:-1]:
            for twopart_seq in product(*[twoparts_fixed_part(twonums, part)
                                         for twonums in twonums_seq]):
                P += VP_fiber_prod_interior([[len(twopart[gp]) for gp in part]
                                             for twopart in twopart_seq]) \
                 * sp.prod([self.get_2ass(tuple(tuple(len(bubble[s])
                                                      for s in bubble.keys())
                                                for twopart in twopart_seq
                                                for bubble in twopart[gp]),
                                          len(gp))
                            for gp in part])
        self.twoass_roots_unfused[m] = sp.expand(P)

    def get_single_2ass(self, m):
        return self.get_2ass((m,), len(m))

    def add_specific_dim(self, d):
        for n in all_n(d):
            if ((tuple(n),), len(n)) not in self.twoass:
                self.update_2ass((tuple(n),), len(n))
