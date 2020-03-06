'''
    A function to output k-mer distribution histogram
    w.r.t different k (k-mer size) and n (sketch size)
'''
import os
import sys
import numpy as np
import pickle
from scipy.stats import wasserstein_distance
from scipy.spatial import distance
import matplotlib.pyplot as plt; plt.rcdefaults()
import seaborn as sns
import matplotlib.pyplot as plt
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


def k_mer_sketch_histogram(n, k, genome, rev_comp=False):
    n = int(n); k = int(k)
    # input: n - sketch size (# Hash function), k - k-mer size, genome - fasta(.gz)
    # return np.array of distribution and histogram
    KMC_outname = genome.split('/')[-1] + '.ksize' + str(k) + '.res'
    outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
    # if the value not stored, compute it, else load it
    if not os.path.isfile(outpath + KMC_outname + '.sketch' + str(n) + '.pickle'):
        # if MinHash Estimator with larger sketch size doesn't exists, compute it with current sketch size
        MHS_filenames = os.listdir(outpath + 'MH_counts/')
        if MHS_filenames:
            try:
                # get min sketch sizes of existing MinHash Estimator which is greater than n
                sketch_size_existing = [int(_.split('.sketch')[-1].split('.MHScounts.pickle')[0]) for _ in MHS_filenames
                                        if (_.endswith('.MHScounts.pickle') and '.ksize' + str(k) + '.' in _
                                            and KMC_outname in _)]
                sketch_size_existing_greater_than_n = min([_ for _ in sketch_size_existing if _ >= n])
                MHS_count_name = outpath + 'MH_counts/' + KMC_outname + '.sketch' + str(sketch_size_existing_greater_than_n) + '.MHScounts.pickle'
                with open(MHS_count_name, 'rb') as MHS_sketch_count_file:
                    MHS_count = pickle.load(MHS_sketch_count_file)
                    counts = MHS_count[:n]
            # sketch_size_existing_greater_than_n is empty
            except (ValueError, FileNotFoundError):
                MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='n', input_file_name=genome, rev_comp=rev_comp)
                counts = MHS._counts
        else:
            MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='n', input_file_name=genome, rev_comp=rev_comp)
            counts = MHS._counts
        # check if MHS counts with k & n is saved nor not
        MHS_count_name = outpath + 'MH_counts/' + KMC_outname + '.sketch' + str(n) + '.MHScounts.pickle'
        if not os.path.isfile(MHS_count_name):
            with open(MHS_count_name, 'wb') as MHS_sketch_count_file:
                pickle.dump(counts, MHS_sketch_count_file)
        # turn array of counts of k-mers into occurrence of k-mers with the counts
        dist = np.zeros(max(counts))
        for _c in counts:
            dist[_c - 1] = dist[_c - 1] + 1
        dist_norm = dist / np.sum(dist)
        with open(outpath + KMC_outname + '.sketch' + str(n) + '.pickle', 'wb') as config_sketch_file:
            pickle.dump([dist, dist_norm], config_sketch_file)
    else:
        with open(outpath + KMC_outname + '.sketch' + str(n) + '.pickle', 'rb') as config_sketch_file:
            dist, dist_norm = pickle.load(config_sketch_file)
    return dist, dist_norm  # np.array(list(dist))


def k_mer_global_histogram_KMC(k, genome, runKMC=False):
    # create KMC database
    KMC_outname = genome.split('/')[-1] + '.ksize' + str(k) + '.res'
    outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
    # if the value not stored, compute it, else load it
    if not os.path.isfile(outpath + KMC_outname + '.global.pickle'):
        # check if KMC database exists
        if runKMC or not os.path.isfile(outpath + KMC_outname + '.kmc_pre'):
            # -ci2 - exclude k-mers occurring less than 2 times
            if '.fastq' in genome or '.fq' in genome:
                os.system('/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/kmc -fq -v -ci0 -b -cs3000 -k%d %s %s %s'
                          % (k, genome, outpath + KMC_outname, outpath))
            elif '.fasta' in genome or '.fa' in genome or '.fna' in genome:
                os.system('/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/kmc -fm -v -ci0 -b -cs3000 -k%d %s %s %s'
                          % (k, genome, outpath + KMC_outname, outpath))
            else:
                print("Is file fa/fq? Check file and its name!", file=sys.stderr)
                exit(2)
        # read KMC database to get count values
        # TODO: KMC doesnot count k-mers more than 255 times
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
        dist = np.zeros(max(counter_dict.keys()))
        for _k, _v in counter_dict.items():
            dist[_k - 1] = _v
        dist_norm = dist / np.sum(dist)
        with open(outpath + KMC_outname + '.global.pickle', 'wb') as config_global_file:
            pickle.dump([dist, dist_norm], config_global_file)
    else:
        with open(outpath + KMC_outname + '.global.pickle', 'rb') as config_global_file:
            dist, dist_norm = pickle.load(config_global_file)
    return dist, dist_norm


def histogram_draw(dist, genome, histogram_name=None, k=0, n=0, rm_occur_leq=0, to_normalize=False, to_log=False):
    if histogram_name:
        figure_name = histogram_name
    else:
        species_name = genome.split('/')[-1].split('.')[0]
        figure_name = 'Dist_k-mer_%s_ksize%d_sketch%d.png' % (species_name, k, n)

    dist = np.array(dist)
    # remove counts if the occurrence is smaller than or equal to rm_occur_leq
    dist = dist[rm_occur_leq:]
    if to_normalize:
        dist = dist/sum(dist)
    if to_log:
        # dist = np.nan_to_num(np.log(dist), nan=0.0, neginf=0)
        dist = np.log(dist)
    # index + 1 is the occurrence of k-mer

    plt.clf()
    plt.rcdefaults()
    plt.grid(axis='y', alpha=0.75)
    y_pos = np.arange(rm_occur_leq + 1, len(dist) + rm_occur_leq + 1)
    if to_log:
        plt.ylabel('log # such k-mers')
    else:
        plt.ylabel('# such k-mers')
    plt.xlabel('# occurrences of k-mer')
    plt.title(figure_name)
    ax = sns.scatterplot(x=y_pos, y=dist, marker='+', s=1)
    plt.savefig(os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/' + figure_name)


def histogram_draw_multi(dist_norm_1, dist_norm_2, genome, histogram_name=None, k=0, n=0):
    # TODO: put two histo in one figure
    # bar chart with 2 bars
    pass


def total_variation_Metric(histo1, histo2, occur_at_least=1, to_normalize=True):
    # histo1, histo2 are normalized distributions
    return L1_metric(histo1, histo2, occur_at_least, to_normalize)/2


def wasserstein_metric(histo1, histo2, occur_at_least=1, to_normalize=True):
    # histo1, histo2 are normalized distributions
    # make two distribution have the same dim
    histo1[:occur_at_least] = np.zeros(occur_at_least)
    histo2[:occur_at_least] = np.zeros(occur_at_least)
    if to_normalize:
        histo1 = histo1/sum(histo1)
        histo2 = histo2/sum(histo2)
    if len(histo1) > len(histo2):
        _padding = np.zeros(len(histo1)-len(histo2))
        histo2 = np.concatenate((histo2, _padding))
    elif len(histo1) < len(histo2):
        _padding = np.zeros(len(histo2) - len(histo1))
        histo1 = np.concatenate((histo1, _padding))
    return wasserstein_distance(histo1, histo2)


def L1_metric (histo1, histo2, occur_at_least=1, to_normalize=True):
    # histo1, histo2 are normalized distributions
    # make two distribution have the same dim
    histo1[:occur_at_least] = np.zeros(occur_at_least)
    histo2[:occur_at_least] = np.zeros(occur_at_least)
    if to_normalize:
        histo1 = histo1/sum(histo1)
        histo2 = histo2/sum(histo2)
    if len(histo1) > len(histo2):
        _padding = np.zeros(len(histo1)-len(histo2))
        histo2 = np.concatenate((histo2, _padding))
    elif len(histo1) < len(histo2):
        _padding = np.zeros(len(histo2) - len(histo1))
        histo1 = np.concatenate((histo1, _padding))
    return distance.cityblock(histo1, histo2)


if __name__ == "__main__":
    # n = number of hash functions
    n_list = [10]
    # k = k-mer size
    k_list = [1]
    exit('Use k_mer_dist_test_run.py to run this script.')
