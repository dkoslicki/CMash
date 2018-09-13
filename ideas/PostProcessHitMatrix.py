# find the unique k-mers
# TP's should include at least some unique k-mer hits
# possibly include some with overwhelming amount of non-unique hits?
# or sets of ones where there is no unique hits among them
# scipy.sparse.save_npz

from CMash import MinHash as MH
import numpy as np
import scipy as sp
from scipy.io import loadmat
import pandas as pd
import os


# First, read in the sketches of just the appropriate critters
#cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.csv'
cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/test.csv'
training_base_name = '/nfs1/Koslicki_Lab/koslickd/RepoPhlAn-7-24-18/out/microbes_24072018/fna/'
#training_hdf_file = '/home/dkoslicki/Data/MiCOPMinHash/AllBacteria.hd5'
training_hdf_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/Test.h5'
hit_matrices_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/test_hit_matrix.npz'
coverage_threshold = 0
sort_key = 'k=5'
location_of_thresh = -1

# read in the file and sort as needed
df = pd.read_csv(cmash_out_file, index_col=0)

# Need to do this after the fact, otherwise my basis could get screwed up
#max_key = df.keys()[-1]
#if not sort_key:
#	sort_key = max_key
#cmash_thresh = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)

names_passed_thresh = list(df.index)
names_passed_thresh_with_path = []
for name in names_passed_thresh:
	names_passed_thresh_with_path.append(training_base_name + name)
CEs = MH.import_multiple_from_single_hdf5(training_hdf_file, import_list=names_passed_thresh_with_path)
training_file_names = [c.input_file_name for c in CEs]

# import the hit matrices
hit_matrices_dict = loadmat(hit_matrices_file)

# now, for each one of the sketches, look for unique k-mer in it, set non-unique to zero
k_range = sorted([int(i.split('=')[1]) for i in df.keys()])
num_unique = dict()
for i in range(len(CEs)):
	for k_size in k_range:
		current_kmers = [k[:k_size] for k in CEs[i]._kmers]
		current_kmers_set = set([k[:k_size] for k in CEs[i]._kmers])
		other_kmers_set = set()
		for j in range(len(CEs)):
			if j != i:
				other_kmers_set.update([k[:k_size] for k in CEs[j]._kmers])
		non_unique = current_kmers_set.intersection(other_kmers_set)
		to_zero_indicies = [current_kmers.index(kmer) for kmer in non_unique]
		hit_matrices_dict['k=%d' % k_size][i, to_zero_indicies] = 0  # set these to zero since they show up in other sketches (so not informative)
		num_unique[i, k_range.index(k_size)] = len(current_kmers_set) - len(non_unique)  # keep track of the size of the unique k-mers

containment_indices = np.zeros((len(CEs), len(k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
for k_size_loc in range(len(k_range)):
	k_size = k_range[k_size_loc]
	containment_indices[:, k_size_loc] = (hit_matrices_dict['k=%d' % k_size].sum(axis=1).ravel()) #/float(num_hashes))

for k_size_loc in range(len(k_range)):
	k_size = k_range[k_size_loc]
	for hash_loc in np.where(containment_indices[:, k_size_loc])[0]:  # find the genomes with non-zero containment
		unique_kmers = set()
		for kmer in CEs[hash_loc]._kmers:
			unique_kmers.add(kmer[:k_size])  # find the unique k-mers
		containment_indices[hash_loc, k_size_loc] /= float(len(unique_kmers))  # TODO: this doesn't seem like the right way to normalize, but the below gives numbers > 1
		#containment_indices[hash_loc, k_size_loc] /= float(num_unique[hash_loc, k_size_loc])  # divide by the unique num of k-mers

results = dict()
for k_size_loc in range(len(k_range)):
	ksize = k_range[k_size_loc]
	key = 'k=%d' % ksize
	results[key] = containment_indices[:, k_size_loc]
df = pd.DataFrame(results, map(os.path.basename, training_file_names))
df = df.reindex(labels=['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
sort_key = 'k=%d' % k_range[location_of_thresh]
max_key = 'k=%d' % k_range[-1]
filtered_results = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)  # only select those where the highest k-mer size's count is above the threshold
#filtered_results.to_csv(results_file, index=True, encoding='utf-8')
