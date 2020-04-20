'''
    A helper script to run/debug/plot k_mer_dist_test.py
'''
import seaborn as sns
import os
import sys
import numpy as np
import matplotlib.pyplot as plt; plt.rcdefaults()
# import k_mer_dist_test functions
try:
    from k_mer_dist import k_mer_dist_test as kmd
except ModuleNotFoundError:
    try:
        import k_mer_dist_test as kmd
    except ModuleNotFoundError:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from k_mer_dist import k_mer_dist_test as kmd
# import k_mer_dist_test_by_seq
# this is another version of CMash using sequence of kmers to count (instead of their hash values)
try:
    from k_mer_dist import k_mer_dist_test_by_seq as kmdq
except ModuleNotFoundError:
    try:
        import k_mer_dist_test_by_seq as kmdq
    except ModuleNotFoundError:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from k_mer_dist import k_mer_dist_test_by_seq as kmdq


def plot(h1, h2, h3, h4, x_low=8, x_high=64):
    # plot KMC, CMash, KmerEst and ntCard
    # all distributions normalized
    h1 = h1/sum(h1)
    h2 = h2/sum(h2)
    h3 = h3/sum(h3)
    h4 = h4/sum(h4)
    # make them the same length
    max_len = max(len(h1), len(h2), len(h3), len(h4))
    h1 = np.concatenate((h1, np.zeros(max_len - len(h1))))
    h2 = np.concatenate((h2, np.zeros(max_len - len(h2))))
    h3 = np.concatenate((h3, np.zeros(max_len - len(h3))))
    h4 = np.concatenate((h4, np.zeros(max_len - len(h4))))
    # concatenate y, x, hue
    y = np.concatenate((h1[x_low-1:x_high], h2[x_low-1:x_high], h3[x_low-1:x_high], h4[x_low-1:x_high]))
    x = np.arange(x_low, x_high+1)
    x = np.concatenate((x, x, x, x))
    hue = ['KMC'] * (x_high-x_low+1) + ['CMash'] * (x_high-x_low+1) + ['KmerEst'] * (x_high-x_low+1) + ['ntCard'] * (x_high-x_low+1)
    # make plot, set y-axis in log scale
    plt.figure(figsize=(9, 3))
    ax = sns.barplot(x, y, hue)
    ax.set_yscale("log")
    ax.set(ylim=(0.000005, 0.01))
    plt.xticks(np.arange(x_low-x_low, x_high-x_low, 10), np.arange(x_low, x_high+1, 10))
    plt.text(x_high-x_low-15, 0.002, "k="+str(k)+"\nn="+str(int(n)))
    plt.xlabel('Occurrences of k-mers')
    plt.ylabel('Proportion of such k-mer')
    # plt.show()
    plt.savefig(genome+'.k'+str(k)+'n'+str(n)+'.png')


if __name__ == "__main__":
    g1 = '/gpfs/scratch/xbz5174/short_term_work_Feb/data_real/SRR072232.fastq.gz'
    g2 = '/gpfs/scratch/xbz5174/short_term_work_Feb/data_real/SRR172902.fastq.gz'
    n_list = [int(1e4), int(1e3)]  # sketch size

    # make plot
    for k in [120]: #[21, 60, 120]:
        for g in [g1]: #[g1, g2]
            if g == g1 or k !=120:
                # get global distribution
                dist_global, _ = kmd.k_mer_global_histogram_KMC(k, g, runKMC=False)
                for n in n_list:
                    dist_sk, _ = kmd.k_mer_sketch_histogram(n, k, g, rev_comp=False)
                    # dist_skq, _ = kmdq.k_mer_sketch_histogram(n, k, g, rev_comp=False)  # k-mer abundance by CMash with sequences
                    # get name of genome
                    genome = g.split('/')[-1].split('.')[0]
                    # read KmerEst output
                    with open('sota_res/%s.k%d.kmerEst.hist' % (genome, k), 'r') as k_est_f:
                        # remove top 2 lines (F1, F0)
                        dist_ke = np.array([int(line.split()[-1]) for line in k_est_f.read().splitlines()[2:]])
                    # read ntCard output
                    with open('sota_res/%s.k%d.ntcard_k%d.hist' % (genome, k, k), 'r') as nt_card_f:
                        # remove top 2 lines (F1, F0)
                        dist_nt = np.array([int(line.split()[-1]) for line in nt_card_f.read().splitlines()[2:]])
                    plot(dist_global, dist_sk, dist_ke, dist_nt, 4, 48)

    # distrance matrices
    for k in [21, 60, 120]:
        for g in [g1, g2]:
            if g == g1 or k !=120:
                # get global distribution
                dist_global, _ = kmd.k_mer_global_histogram_KMC(k, g, runKMC=False)
                genome = g.split('/')[-1].split('.')[0]
                # read KmerEst output
                with open('sota_res/%s.k%d.kmerEst.hist' % (genome, k), 'r') as k_est_f:
                    # remove top 2 lines (F1, F0)
                    dist_ke = np.array([int(line.split()[-1]) for line in k_est_f.read().splitlines()[2:]])
                # read ntCard output
                with open('sota_res/%s.k%d.ntcard_k%d.hist' % (genome, k, k), 'r') as nt_card_f:
                    # remove top 2 lines (F1, F0)
                    dist_nt = np.array([int(line.split()[-1]) for line in nt_card_f.read().splitlines()[2:]])
                for n in [1e3, 1e4, 1e5, 1e6]:
                    dist_sk, _ = kmd.k_mer_sketch_histogram(n, k, g, rev_comp=False)
                    print(genome, k, n,
                          kmd.L1_metric(dist_global, dist_sk),
                          kmd.L1_metric(dist_global, dist_ke),
                          kmd.L1_metric(dist_global, dist_nt),
                          sep='\t')

