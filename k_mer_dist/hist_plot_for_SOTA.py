'''
    A helper script to plot and compute distance matrix of the results of KmerEstimate, ntCard
'''
import seaborn as sns
import os
import sys
from collections import Counter
import numpy as np
from scipy.stats import wasserstein_distance
from scipy.spatial import distance
import matplotlib.pyplot as plt; plt.rcdefaults()
# use pickle to store/load ground truth
import pickle
# import k_mer_dist_test functions
try:
    from k_mer_dist import k_mer_dist_test as kmd
except ModuleNotFoundError:
    try:
        import k_mer_dist_test as kmd
    except ModuleNotFoundError:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from k_mer_dist import k_mer_dist_test as kmd


if __name__ == "__main__":
    genome_path = '/gpfs/scratch/xbz5174/short_term_work_Feb/data_real/SRR072232.fastq.gz'
    genome = 'SRR072232'
    k_sizes = [21, 60]
    # read and compute for kmerEst
    for k in k_sizes:
        dist, _ = kmd.k_mer_global_histogram_KMC(k, genome_path, runKMC=False)
        print(dist)

        with open('sota_res/%s.k%d.kmerEst.hist' % (genome, k), 'r') as k_est_f:
            # remove top 2 lines (F1, F0)
            dist_est = np.array([int(line.split()[-1]) for line in k_est_f.read().splitlines()[2:]])
            print('KmerEstimate k=', k, '\n', dist_est)
            print(kmd.total_variation_Metric(dist, dist_est, occur_at_least=0),
                  kmd.wasserstein_metric(dist, dist_est, occur_at_least=0),
                  kmd.L1_metric(dist, dist_est, occur_at_least=0),
                  kmd.weighted_L1_matrix(dist, dist_est, occur_at_least=0),
                  sep='\t')
            kmd.histogram_draw(dist_est, genome, k=k, rm_occur_leq=0, to_normalize=False, to_log=True,
                               histogram_name='kmerEst.%s.k%g.png' % (genome, k))
            print(Counter(dist_est))

        with open('sota_res/%s.k%d.ntcard_k%d.hist' % (genome, k, k), 'r') as nt_card_f:
            # remove top 2 lines (F1, F0)
            dist_est = np.array([int(line.split()[-1]) for line in nt_card_f.read().splitlines()[2:]])
            print('ntCard k=', k, '\n', dist_est)
            print(kmd.total_variation_Metric(dist, dist_est, occur_at_least=0),
                  kmd.wasserstein_metric(dist, dist_est, occur_at_least=0),
                  kmd.L1_metric(dist, dist_est, occur_at_least=0),
                  kmd.weighted_L1_matrix(dist, dist_est, occur_at_least=0),
                  sep='\t')
            kmd.histogram_draw(dist_est, genome, k=k, rm_occur_leq=0, to_normalize=False, to_log=True, histogram_name='ntCard.%s.k%g.png' % (genome, k))
            print(Counter(dist_est))