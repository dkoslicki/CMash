'''
    A function to output k-mer distribution histogram
    w.r.t different k (k-mer size) and n (sketch size)
'''
import os
import sys
import numpy as np
from scipy.stats import wasserstein_distance
from scipy.spatial import distance
# import KMC python API (compiled)
sys.path.insert(1, '/storage/home/xbz5174/work/tools/KMC-3.1.1/')
sys.path.insert(1, '/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/')
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


def k_mer_sketch_histogram(n, k, genome, true_histogram=True, histogram_name=None, rev_comp=False):
    # input: n - sketch size (# Hash function), k - k-mer size, genome - fasta(.gz)
    # return np.array of distribution and histogram
    MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
    # turn array of counts of k-mers into occurence of k-mers with the counts
    counts = MHS._counts
    dist = np.zeors(max(counts))
    for _c in counts:
        dist[_c] = dist[_c] + 1
    dist_norm = dist / np.sum(dist)
    if true_histogram:
        histogram_draw(dist_norm, genome, histogram_name)
    return dist, dist_norm  # np.array(list(dist))


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


def k_mer_global_histogram_KMC(k, genome, true_histogram=True, histogram_name=None, runKMC = False):
    # TODO: not exactly sure but KMC seems doesn't accept sequences spanning more than 4 lines
    # create KMC database
    KMC_outname = genome.split('/')[-1] + '.res'
    outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
    if runKMC:
        # -ci2 - exclude k-mers occurring less than 2 times
        if '.fastq' in genome or '.fq' in genome:
            os.system('/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/kmc -fq -ci2 -k%d %s %s %s'
                      %(k, genome, outpath + KMC_outname, outpath))
        elif '.fasta' in genome or '.fa' in genome or '.fna' in genome:
            os.system('/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/kmc -fa -ci2 -k%d %s %s %s'
                      %(k, genome, outpath + KMC_outname, outpath))
        else:
            print("Is file fa/fq? Check file and its name!", file=sys.stderr)
            exit(2)
    # read KMC database to get count values
    kmer_data_base = kmc.KMCFile()
    kmer_data_base.OpenForListing(outpath + KMC_outname)
    kmer_object = kmc.KmerAPI(kmer_data_base.Info().kmer_length)
    counter = kmc.Count()
    counter_dict = {}
    while kmer_data_base.ReadNextKmer(kmer_object, counter):
        try:
            counter_dict[int(counter.value)] = counter_dict[int(counter.value)] + 1
        except KeyError:
            counter_dict[int(counter.value)] = 1
    # get distribution
    dist = np.zeors(max(counter_dict.values()))
    for _k, _v in counter_dict.items():
        dist[_k] = _v
    dist_norm = dist / np.sum(dist)
    return dist, dist_norm


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
        dist_global = k_mer_global_histogram_KMC(k=k, genome=file, histogram_make=True, runKMC=True)
        print('k=%d k-mer global dist:' % k, dist_global)