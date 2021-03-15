from copy import deepcopy
from itertools import product

def parts(nums):
    # Returns all partitions of a list

    if len(nums) == 1:
        return [[nums]]
    elif len(nums) >= 2:
        new_parts = []
        for part in parts(nums[:-1]):
            # we take all partitions of the first n-1 elements
            
            new_parts += [part + [[nums[-1]]]]
            # for every such partition, we include the new partition
            # where the n-th element is on its own
            
            temp_part = deepcopy(part)
            temp_part[-1] += [nums[-1]]
            new_parts += [temp_part]
        return new_parts
    
def cleaned_parts(nums):
    l = parts(nums)
    l.remove([nums])
    return l   

def nicer_prod(*args):
    return list(list(t) for t in product(*args))

def flatten(l):
    return [subitem for item in l for subitem in item]

def K_list_of_lists(nums):
    if len(nums) == 1:
        return [[[nums[0]]]]
    elif len(nums) == 2:
        return [[[nums[0]],[nums[1]]]]
    elif len(nums) >= 3:
        bs_list = []
        for outer_bs in cleaned_parts(nums):
            sub_bs_list = nicer_prod(*[K_list_of_lists(part) \
                                       for part in outer_bs])
            sub_bs_list = [flatten(item) for item in sub_bs_list]
            bs_list += [sub_bs + [outer_b for outer_b in outer_bs \
                                  if (outer_b not in sub_bs)] \
                            for sub_bs in sub_bs_list]
        return bs_list

# def K_list_of_sets(nums):
#     return [set(tuple(b) for b in bs) for bs in K_list_of_lists(nums)]
    
# def K_poset_helper(bs):
#     # returns all the bracketings that are just finer than the given one
    
#     bs_list = []
#     for b in bs:
#         for part in cleaned_parts(b):
#             bs_list += [bs.union(set([tuple(b) for b in part]))]
#     return bs_list

# def K_poset(nums):
#     if len(nums) == 1:
#         return {nums[0]:[]}
#     elif len(nums) >= 2:
        