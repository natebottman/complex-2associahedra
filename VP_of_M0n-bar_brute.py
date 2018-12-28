# a somewhat brute force approach to computing the
# virtual Poincare polynomial of M0n. A less-brute-
# force approach is implemented in another file.

from sympy import *
from math import factorial as fac
import copy

x = symbols('x')

def binary_seqs(s):
    # returns list of binary sequences of length s, for use in parts
    if s==0:
        return [[]]
    elif s>=1:
        lists = [[[i] + seq for seq in binary_seqs(s-1)] for i in [0,1]]
        return lists[0]+lists[1]

def parts(nums):
    # Returns all partitions of a list. Note that the first partition
    # is all singletons, and the last partition has everything together.
    if len(nums) in [0,1]:
        return [] if nums==[] else [[nums]]
    elif len(nums)>=2:
        new_parts = []
        for part in parts(nums[:-1]):
            new_parts += [part + [[nums[-1]]]]
            temp_part = copy.deepcopy(part)
            for i in range(len(temp_part)):
                new_parts += [[group if j!=i else (group+[nums[-1]]) \
                             for (j,group) in enumerate(temp_part)]]
        return new_parts
    
def VP_of_C_minus_d(d):
    return prod([x**2-(2+i) for i in range(0,d)])

class Repo:
    # This is a repository of VP's of M_{0,r}-bar.
    # Setting it up this way allows us to compute more and more
    # VP's, in a memoized fashion.

    def __init__(self):
        self.polys = { 1 : 1, 2 : 1 }
        
    def return_all(self):
        return self.polys
    
    def update(self,r):
        m = max(self.polys.keys())
        for s in range(m+1,r+1):
            ps = parts(list(range(1,s+1)))[:-1]
            poly = sum([VP_of_C_minus_d(len(p)-2)* \
                        prod([self.polys[len(c)] for c in p]) for p in ps])
            self.polys[s] = expand(poly)