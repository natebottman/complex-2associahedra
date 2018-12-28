from sympy import *
from math import factorial as fac

x = symbols('x')

def multinom(a,nums):
    return fac(a)//prod([fac(el) for el in nums])

def sum_lists(p1,p2):
    return [p1[i]+p2[i] for i in range(len(p1))]

def add_zeros(p1,k):
    return p1+[0]*k

def part_sum(p1,p2):
    # Sum of two partitions, not necessarily of
    # the same length.
    if len(p1) == len(p2):
        return sum_lists(p1,p2)
    elif len(p1)>len(p2):
        return sum_lists(p1,add_zeros(p2,len(p1)-len(p2)))
    else:
        return part_sum(p2,p1)

def parts(a,b):
    # returns list of all partitions of a
    # whose largest piece(s) has size at most b.
    L = [[a]]
    if a==0:
        return [[]]
    elif b==1:
        return L
    else:
        for i in range(2,b+1):
            for j in range(1,a//i+1):
                j_chunk = [0]*(i-1)+[j]
                L.extend([part_sum(part,j_chunk) \
                    for part in parts(a-j*i,i-1)])
        return L
    
def coeff(r,p):
    return multinom(r,p)//prod([fac(i+1)**k for (i,k) in enumerate(p)])

def VP_of_C_minus_d(d):
    return prod([x**2-(2+i) for i in range(0,d)])
    
def VP(n):
    # returns list of virtual Poincare polynomials of M_{0,k+1}-bar,
    # for k between 1 and n
    L = [1,1]
    for r in range(3,n+1):
        L.append(expand(sum([coeff(r,p) * \
                             VP_of_C_minus_d(sum(p)-2)* \
                             prod([L[i]**k for (i,k) in enumerate(p)]) \
                                 for p in parts(r,r-1)])))
    return L

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
            poly = sum([coeff(s,p) * VP_of_C_minus_d(sum(p)-2)* \
                        prod([self.polys[i+1]**k for (i,k) in enumerate(p)]) \
                        for p in parts(s,s-1)])
            self.polys[s] = expand(poly)