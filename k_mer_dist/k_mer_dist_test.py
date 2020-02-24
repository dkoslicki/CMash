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


def k_mer_sketch_histogram(n, k, genome, true_histogram=True, histogram_name=None, rev_comp=False):
    # input: n - sketch size (# Hash function), k - k-mer size, genome - fasta(.gz)
    # return np.array of distribution and histogram
    KMC_outname = genome.split('/')[-1] + '.ksize' + str(k) + '.res'
    outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
    # if the value not stored, compute it, else load it
    if not os.path.isfile(outpath + KMC_outname + '.sketch' + str(n) + '.pickle'):
        # if MinHash Estimator with larger sketch size doesn't exists, compute it with current sketch size
        MHS_filenames = os.listdir(outpath + '/MH_counts/')
        if MHS_filenames:
            try:
                # get min sketch sizes of existing MinHash Estimator which is greater than n
                sketch_size_existing = [int(_.split('.sketch')[-1].split('.MHScounts.pickle')[0]) for _ in MHS_filenames
                                        if _.endswith('.MHScounts.pickle') and '.ksize' + str(k) + '.' in _]
                sketch_size_existing_greater_than_n = min([_ for _ in sketch_size_existing if _ >= n])
                MHS_count_name = outpath + '/MH_counts/' + KMC_outname + '.sketch' + str(sketch_size_existing_greater_than_n) + '.MHScounts.pickle'
                with open(MHS_count_name, 'rb') as MHS_sketch_count_file:
                    MHS_count = pickle.load(MHS_sketch_count_file)
                    counts = MHS_count[:n]
            # sketch_size_existing_greater_than_n is empty
            except ValueError:
                MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
                counts = MHS._counts
        else:
            MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
            counts = MHS._counts
        # check if MHS counts with k & n is saved nor not
        MHS_count_name = outpath + '/MH_counts/' + KMC_outname + '.sketch' + str(n) + '.MHScounts.pickle'
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
    if true_histogram:
        histogram_draw(dist_norm, genome, histogram_name, k, n)
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
    from collections import Counter
    dist = Counter(k_mer_dict.values())
    return dist  # np.array(list(dist))


def k_mer_global_histogram_KMC(k, genome, true_histogram=True, histogram_name=None, runKMC = False):
    # TODO: not exactly sure but KMC seems doesn't accept sequences spanning more than 4 lines
    # create KMC database
    KMC_outname = genome.split('/')[-1] + '.ksize' + str(k) + '.res'
    outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
    # if the value not stored, compute it, else load it
    if not os.path.isfile(outpath + KMC_outname + '.global.pickle'):
        # check if KMC database exists
        if runKMC or not os.path.isfile(outpath + KMC_outname + '.kmc_pre'):
            # -ci2 - exclude k-mers occurring less than 2 times
            if '.fastq' in genome or '.fq' in genome:
                os.system('/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/kmc -fq -ci2 -k%d %s %s %s'
                          % (k, genome, outpath + KMC_outname, outpath))
            elif '.fasta' in genome or '.fa' in genome or '.fna' in genome:
                os.system('/storage/home/xbz5174/work/tools/KMC-3.1.1/bin/kmc -fa -ci2 -k%d %s %s %s'
                          % (k, genome, outpath + KMC_outname, outpath))
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
        dist = np.zeros(max(counter_dict.values()))
        for _k, _v in counter_dict.items():
            dist[_k - 1] = _v
        dist_norm = dist / np.sum(dist)
        with open(outpath + KMC_outname + '.global.pickle', 'wb') as config_global_file:
            pickle.dump([dist, dist_norm], config_global_file)
    else:
        with open(outpath + KMC_outname + '.global.pickle', 'rb') as config_global_file:
            dist, dist_norm = pickle.load(config_global_file)
    if true_histogram:
        histogram_draw(dist_norm, genome, histogram_name, k)
    return dist, dist_norm


def histogram_draw(dist_norm, genome, histogram_name=None, k=0, n=0, rm_occur_leq=0):
    if histogram_name:
        figure_name = histogram_name
    else:
        species_name = genome.split('/')[-1].split('.')[0]
        figure_name = species_name + '-k' + str(k) + '-n' + str(n)
    # remove counts if the occurrence is smaller than or equal to rm_occur_leq
    dist_norm = dist_norm[rm_occur_leq:]
    y_pos = np.arange(rm_occur_leq, len(dist_norm))
    plt.bar(y_pos, dist_norm, align='center', alpha=0.5)
    plt.xticks(y_pos, y_pos)
    plt.ylabel('# counts of k-mers occurring such times')
    plt.xlabel('# k-mer occur')
    plt.title('Dist k-mer %f ksize%d sketch%d' %(figure_name, k, n))
    plt.savefig(os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'+'Dist k-mer %f ksize%d sketch%d.png'
                %(figure_name, k, n))
    plt.clf()
    plt.close()
    return 0

def histogram_draw_two(dist_norm_1, dist_norm_2, genome, histogram_name=None, k=0, n=0):
    # TODO: put two histo in one figure
    # bar chart with 2 bars
    pass


def total_variation_Metric(histo1, histo2):
    # histo1, histo2 are normalized distributions
    return L1_metric(histo1, histo2)/2


def wasserstein_metric(histo1, histo2):
    # histo1, histo2 are normalized distributions
    # make two distribution have the same dim
    if len(histo1) > len(histo2):
        _padding = np.zeros(len(histo1)-len(histo2))
        histo2 = np.concatenate((histo2, _padding))
    elif len(histo1) < len(histo2):
        _padding = np.zeros(len(histo2) - len(histo1))
        histo1 = np.concatenate((histo1, _padding))
    return wasserstein_distance(histo1, histo2)


def L1_metric (histo1, histo2):
    # histo1, histo2 are normalized distributions
    # make two distribution have the same dim
    if len(histo1) > len(histo2):
        _padding = np.zeros(len(histo1)-len(histo2))
        histo2 = np.concatenate((histo2, _padding))
    elif len(histo1) < len(histo2):
        _padding = np.zeros(len(histo2) - len(histo1))
        histo1 = np.concatenate((histo1, _padding))
    return distance.cityblock(histo1, histo2)


if __name__ == "__main__":

    # n = number of hash functions
    n_list = [int(1e2), int(1e4), int(1e6), int(1e8)]
    # k = k-mer size
    k_list = [21, 60, 120]
    # Use abs path to avoid confusion
    genome_mid_1 = '/storage/home/xbz5174/scratch/short_term_work_Feb/data_real/SRR072232.fastq.gz'
    genome_mid_2 = '/storage/home/xbz5174/scratch/short_term_work_Feb/data_real/SRR172902.fastq.gz'
    genome = genome_mid_1

    print('Current working dir:')
    os.system('pwd')
    # get sketch k-mer distributions
    for k in k_list:
        KMC_outname = genome.split('/')[-1] + '.ksize' + str(k) + '.res'
        outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
        # if the value not stored, compute it, else load it
        if not os.path.isfile(outpath + KMC_outname + '.global.pickle'):
            dist_global, dist_global_norm = k_mer_global_histogram_KMC(k, genome, true_histogram=True, runKMC=False)
            with open(outpath + KMC_outname + '.global.pickle', 'wb') as config_global_file:
                pickle.dump([dist_global, dist_global_norm], config_global_file)
        else:
            with open(outpath + KMC_outname + '.global.pickle', 'rb') as config_global_file:
                dist_global, dist_global_norm = pickle.load(config_global_file)
        print('ksize:',k,' - global dist', dist_global_norm)
        for n in n_list:
            # if the value not stored, compute it, else load it
            if not os.path.isfile(outpath + KMC_outname + '.sketch' + str(n) +'.pickle'):
                dist_sketch, dist_sketch_norm = k_mer_sketch_histogram(n, k, genome, true_histogram=True, rev_comp=False)
                with open(outpath + KMC_outname + '.sketch' + str(n) +'.pickle', 'wb') as config_sketch_file:
                    pickle.dump([dist_sketch, dist_sketch_norm], config_sketch_file)
            else:
                with open(outpath + KMC_outname +  '.sketch' + str(n) +'.pickle', 'rb') as config_sketch_file:
                    dist_sketch, dist_sketch_norm = pickle.load(config_sketch_file)
            print('ksize, n:', k, n, ' - sketch dist', dist_sketch_norm)
            print('distances',
                  total_variation_Metric(dist_global_norm, dist_sketch_norm),
                  wasserstein_metric(dist_global_norm, dist_sketch_norm),
                  L1_metric(dist_global_norm, dist_sketch_norm)
                  )