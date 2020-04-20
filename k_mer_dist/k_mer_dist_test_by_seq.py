'''
    A function to output k-mer distribution histogram
    w.r.t different k (k-mer size) and n (sketch size)
    using CMash which counts by k-mer sequences
'''
import os
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt; plt.rcdefaults()
# import MinHash_by_seq
try:
    from CMash import MinHash_by_seq as MH
except ImportError:
    try:
        import MinHash_by_seq as MH
    except ImportError:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from CMash import MinHash_by_seq as MH


def k_mer_sketch_histogram(n, k, genome, rev_comp=False):
    n = int(n); k = int(k)
    # input: n - sketch size (# Hash function), k - k-mer size, genome - fasta(.gz)
    # return np.array of abundance and normalized abundance distribution
    KMC_outname = genome.split('/')[-1] + '.ksize' + str(k) + '.res'
    outpath = os.path.dirname(os.path.realpath(__file__)) + '/kmc_global_count/'
    # if the value not stored, compute it, else load it
    if not os.path.isfile(outpath + KMC_outname + '.sketch' + str(n) + '.byseq.pickle'):
        # if MinHash Estimator with larger sketch size doesn't exists, compute it with current sketch size
        MHS_filenames = os.listdir(outpath + 'MH_counts/')
        if MHS_filenames:
            try:
                # get min sketch sizes of existing MinHash Estimator which is greater than n
                sketch_size_existing = [int(_.split('.sketch')[-1].split('.MHScounts.byseq.pickle')[0]) for _ in MHS_filenames
                                        if (_.endswith('.MHScounts.byseq.pickle') and '.ksize' + str(k) + '.' in _
                                            and KMC_outname in _)]
                sketch_size_existing_greater_than_n = min([_ for _ in sketch_size_existing if _ >= n])
                MHS_count_name = outpath + 'MH_counts/' + KMC_outname + '.sketch' + str(sketch_size_existing_greater_than_n) + '.MHScounts.byseq.pickle'
                with open(MHS_count_name, 'rb') as MHS_sketch_count_file:
                    MHS_count = pickle.load(MHS_sketch_count_file)
                    counts = MHS_count[:n]
            # sketch_size_existing_greater_than_n is empty
            except (ValueError, FileNotFoundError):
                MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
                counts = MHS._counts
        else:
            MHS = MH.CountEstimator(n=n, ksize=k, save_kmers='y', input_file_name=genome, rev_comp=rev_comp)
            counts = MHS._counts
        # check if MHS counts with k & n is saved nor not
        MHS_count_name = outpath + 'MH_counts/' + KMC_outname + '.sketch' + str(n) + '.MHScounts.byseq.pickle'
        if not os.path.isfile(MHS_count_name):
            with open(MHS_count_name, 'wb') as MHS_sketch_count_file:
                pickle.dump(counts, MHS_sketch_count_file)
        # turn array of counts of k-mers into occurrence of k-mers with the counts
        dist = np.zeros(max(counts))
        for _c in counts:
            dist[_c - 1] = dist[_c - 1] + 1
        dist_norm = dist / np.sum(dist)
        with open(outpath + KMC_outname + '.sketch' + str(n) + '.byseq.pickle', 'wb') as config_sketch_file:
            pickle.dump([dist, dist_norm], config_sketch_file)
    else:
        with open(outpath + KMC_outname + '.sketch' + str(n) + '.byseq.pickle', 'rb') as config_sketch_file:
            dist, dist_norm = pickle.load(config_sketch_file)
    return dist, dist_norm  # np.array(list(dist))


if __name__ == "__main__":
    exit('Use k_mer_dist_test_run.py to run this script.')