'''
    A function to output k-mer distribution histogram
    w.r.t different k (k-mer size) and n (sketch size)
'''

import os
import sys
from collections import Counter
import numpy as np
import pandas as pd
import screed
# import MinHash
try:
    from CMash import MinHash as MH
except ImportError:
    try:
        import MinHash as MH
    except ImportError:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from CMash import MinHash as MH


def k_mer_sketch_histogram(n, k, genome, histogram_make=True, histogram_name=None, rev_comp=False):
    MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
    dist = Counter(MHS._counts)
    # save histogram
    if histogram_make:
        if histogram_name:
            figure_name = histogram_name
        else:
            species_name = MHS.input_file_name.split('/')[-1].split('.')[0]
            figure_name = species_name + '-k' + str(k) + '-n' + str(n)
        # TODO: histogram function not implemented yet
    return dist  # np.array(list(dist))


def k_mer_global_histogram(k, genome, histogram_make=True, histogram_name=None, rev_comp=False):
    k_mer_dict = {}
    # TODO: use KMC instead to make it faster
    def add_k_mer_dict(self, kmer, rev_comp=False):
        try:
            k_mer_dict[kmer] += 1
        except KeyError:
            k_mer_dict[kmer] = 1
        pass
    MH.CountEstimator.add = add_k_mer_dict
    MH.CountEstimator(n=1, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
    dist = Counter(k_mer_dict.values())
    return dist  # np.array(list(dist))



if __name__ == "__main__":
    n_list = [10000]  # n = number of hash functions
    k_list = [10]  # k = k-mer size
    file = "/storage/home/xbz5174/scratch/short_term_work_Feb/tests/Organisms/taxid_1307_414_genomic.fna.gz"
    print('Current working dir:'); os.system('pwd')
    # max_prime=9999999999971
    # get sketch k-mer distributions
    for k in k_list:
        # dist_global = k_mer_global_histogram(k=k, genome=file, histogram_make=True)
        # print('k=%d k-mer global dist:' % k, dist_global)
        for n in n_list:
            dist = k_mer_sketch_histogram(n=n, k=k, genome=file, histogram_make=True)
            print('k=%d n=%d k-mer dist:' %(k, n), dist)
    # get global sketch k-mer distributions
    # Must not get global dist before sketch dist as global MH funcitons are overide
    for k in k_list:
        dist_global = k_mer_global_histogram(k=k, genome=file, histogram_make=True)
        print('k=%d k-mer global dist:' %k, dist_global)

