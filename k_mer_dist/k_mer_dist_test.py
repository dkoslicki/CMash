'''
    A function to output k-mer distribution histogram
    w.r.t different k (k-mer size) and n (sketch size)
'''

import os
import sys
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance
from scipy.spatial import distance
# import KMC python API (compilied)
sys.path.insert(1, '/storage/home/xbz5174/work/tools/KMC-3.1.1/')
import py_kmc_api as kmc


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

def k_mer_global_histogram_KMC(k, genome, histogram_make=True, histogram_name=None):
    kmc_file = genome
    outname = genome.split('/')[-1]+'.res'
    # creat KMC database
    # TODO: not exactly sure but KMC seems doesn't accept sequences spanning more than 4 lines
    os.system('kmc -m24 -fa -ci0 -k%d %f %f ./kmc_global_count/' %(k, genome, outname))
    outname = './kmc_global_count/' + outname
    # read KMC data base to get count values
    kmer_data_base = kmc.KMCFile()
    kmer_data_base.OpenForListing(outname)
    kmer_object = kmc.KmerAPI(kmer_data_base.Info().kmer_length)
    counter = kmc.Count()
    counter_list = []
    while kmer_data_base.ReadNextKmer(kmer_object, counter):
        counter_list.append(int(counter.value))
    dist = Counter(counter_list)
    return dist


def total_variation_Metric(histo1, histo2):
    # histo1, histo2 are normalized distributions
    return L1_metric(histo1, histo2)/2


def wasserstein_metric(histo1, histo2):
    # histo1, histo2 are normalized distributions
    return wasserstein_distance(histo1, histo2)


def L1_metric (histo1, histo2):
    # histo1, histo2 are normalized distributions
    return distance.cityblock(histo1, histo2)



if __name__ == "__main__":
    n_list = [10000]  # n = number of hash functions
    k_list = [18]  # k = k-mer size
    file = "/storage/home/xbz5174/scratch/short_term_work_Feb/tests/Organisms/taxid_1307_414_genomic.fna.gz"
    print('Current working dir:'); os.system('pwd')
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

