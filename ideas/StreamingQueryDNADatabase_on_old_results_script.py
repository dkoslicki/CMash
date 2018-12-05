#! /usr/bin/env python
import khmer
import marisa_trie as mt
import numpy as np
import os
import sys
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
except ImportError:
	try:
		import MinHash as MH
	except ImportError:
		sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
		from CMash import MinHash as MH
import multiprocessing
import pandas as pd
import argparse
from argparse import ArgumentTypeError
import re
import matplotlib.pyplot as plt
from hydra import WritingBloomFilter, ReadingBloomFilter
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse import save_npz
from scipy.io import savemat
import timeit
from itertools import islice

data_dir = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/"
output_dir = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/SaveIntermediateFiles/RunAgain/"

# read in the arguments
k_range = [30, 40, 50, 60]
training_data = os.path.join(data_dir, "micopdb_n_1000_k_60.h5")
query_file = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/CAMITest/RL_S001__insert_270.fq.gz"
results_file = os.path.join(output_dir, "output.csv")
npz_file = os.path.splitext(results_file)[0] + "_hit_matrix.npz"
num_threads = 48
location_of_thresh = -1
coverage_threshold = 0
streaming_database_file = os.path.splitext(training_data)[0] + ".tst"  # name of the tst training file
streaming_database_file = os.path.abspath(streaming_database_file)
hydra_file = os.path.join(data_dir, "micopdb_n_1000_k_60_prefilter_30_60_10.bf")
verbose = True
num_reads_per_core = 100000
sensitive = False
match_tuples_file = os.path.join(data_dir, "SaveIntermediateFiles/output_match_tuples.txt")

# Training data
if verbose:
	print("Reading in sketches")
	t0 = timeit.default_timer()
sketches = MH.import_multiple_from_single_hdf5(training_data)
if sketches[0]._kmers is None:
	raise Exception(
		"For some reason, the k-mers were not saved when the database was created. Try running MakeStreamingDNADatabase.py again.")
num_hashes = len(sketches[0]._kmers)  # note: this is relying on the fact that the sketches were properly constructed
max_ksize = sketches[0].ksize

# adjust the k-range if necessary
k_range = [val for val in k_range if val <= max_ksize]

# adjust location of thresh if necessary
if location_of_thresh:
	if location_of_thresh >= len(k_range):
		print("Warning, k_range is of length %d, reducing location of threshold from %d to %d" % (len(k_range), location_of_thresh, len(k_range)))
		location_of_thresh = len(k_range) - 1

# Get names of training files for use as rows in returned tabular data
training_file_names = []
for i in range(len(sketches)):
	training_file_names.append(sketches[i].input_file_name)

if verbose:
	print("Finished reading in sketches")
	t1 = timeit.default_timer()
	print("Time: %f" % (t1 - t0))

if verbose:
	print("Reading in/creating ternary search tree")
	t0 = timeit.default_timer()
# Make the Marissa tree
if streaming_database_file is None:
	streaming_database_file = os.path.splitext(training_data)[0] + ".tst"
	streaming_database_file = os.path.abspath(streaming_database_file)
	print("It appears a tst training file has not been created (did you remember to use MakeStreamingDNADatabase.py?).")
	print("I'm creating one anyway at: %s" % streaming_database_file)
	print("This may take a while...")
	to_insert = set()
	for i in range(len(sketches)):
		for kmer_index in range(len(sketches[i]._kmers)):
			kmer = sketches[i]._kmers[kmer_index]
			to_insert.add(kmer + 'x' + str(i) + 'x' + str(kmer_index))  # format here is kmer+x+hash_index+kmer_index
	tree = mt.Trie(to_insert)
	tree.save(streaming_database_file)
else:
	tree = mt.Trie()
	tree.load(streaming_database_file)

# all the k-mers of interest in a set (as a pre-filter)
if not hydra_file:  # create one
	try:
		all_kmers_bf = WritingBloomFilter(len(sketches)*len(k_range)*num_hashes*2, 0.01)
		for sketch in sketches:
			for kmer in sketch._kmers:
				for ksize in k_range:
					all_kmers_bf.add(kmer[0:ksize])  # put all the k-mers and the appropriate suffixes in
					all_kmers_bf.add(khmer.reverse_complement(kmer[0:ksize]))  # also add the reverse complement
	except IOError:
		print("No such file or directory/error opening file: %s" % hydra_file)
		sys.exit(1)
else:  # otherwise read it in
	try:
		all_kmers_bf = ReadingBloomFilter(hydra_file)
	except IOError:
		print("No such file or directory/error opening file: %s" % hydra_file)
		sys.exit(1)
if verbose:
	print("Finished reading in/creating ternary search tree")
	t1 = timeit.default_timer()
	print("Time: %f" % (t1 - t0))
# Seen k-mers (set of k-mers that already hit the trie, so don't need to check again)
seen_kmers = set()



if verbose:
	print("Forming hit matrix")
	t0 = timeit.default_timer()
#print("Len matches: %d" % len(match_tuples))
# create k_range spare matrices. Rows index by genomes (sketch/hash index), columns index by k_mer_loc
row_ind_dict = dict()
col_ind_dict = dict()
value_dict = dict()
unique_kmers = dict()  # this will keep track of the unique k-mers seen in each genome (sketch/hash loc)
for k_size in k_range:
	row_ind_dict[k_size] = []
	col_ind_dict[k_size] = []
	value_dict[k_size] = []

match_tuples = []
with open(match_tuples_file, "r") as fid:
	for line in fid.readlines():
		line_split = line.strip().split()
		hash_loc = int(line_split[0])
		k_size_loc = int(line_split[1])
		kmer_loc = int(line_split[2])
		match_tuples.append((hash_loc, k_size_loc, kmer_loc))

for hash_loc, k_size_loc, kmer_loc in match_tuples:
	if hash_loc not in unique_kmers:
		unique_kmers[hash_loc] = set()
	k_size = k_range[k_size_loc]
	kmer = sketches[hash_loc]._kmers[kmer_loc][:k_size]
	if kmer not in unique_kmers[hash_loc]:  # if you've seen this k-mer before, don't add it. NOTE: this makes sure we don't over count
		row_ind_dict[k_size].append(hash_loc)
		col_ind_dict[k_size].append(kmer_loc)
		value_dict[k_size].append(1)
		unique_kmers[hash_loc].add(kmer)

hit_matrices = []

for k_size in k_range:
	mat = csc_matrix((value_dict[k_size], (row_ind_dict[k_size], col_ind_dict[k_size])), shape=(len(sketches), num_hashes))
	hit_matrices.append(mat)
	############################################################
	save_npz(os.path.splitext(results_file)[0] + "_first_hit_matrix_%d.npz" % k_size, mat)
	############################################################

if verbose:
	print("Finished forming hit matrix")
	t1 = timeit.default_timer()
	print("Time: %f" % (t1 - t0))

if verbose:
	print("Computing containment indicies")
	t0 = timeit.default_timer()
containment_indices = np.zeros((len(sketches), len(k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
for k_size_loc in range(len(k_range)):
	containment_indices[:, k_size_loc] = (hit_matrices[k_size_loc].sum(axis=1).ravel()) #/float(num_hashes))

for k_size_loc in range(len(k_range)):
	k_size = k_range[k_size_loc]
	for hash_loc in np.where(containment_indices[:, k_size_loc])[0]:  # find the genomes with non-zero containment
		unique_kmers = set()
		for kmer in sketches[hash_loc]._kmers:
			unique_kmers.add(kmer[:k_size])  # find the unique k-mers
		containment_indices[hash_loc, k_size_loc] /= float(len(unique_kmers))  # divide by the unique num of k-mers
if verbose:
	print("Finished computing containment indicies")
	t1 = timeit.default_timer()
	print("Time: %f" % (t1 - t0))
############################################################
savemat(os.path.splitext(results_file)[0] + "_first_containment_indices.mat", {'containment_indices':containment_indices})
############################################################


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

if True:
	if verbose:
		print("Exporting results")
		t0 = timeit.default_timer()
	filtered_results.to_csv(results_file+"non_filtered.csv", index=True, encoding='utf-8')

# TODO: may not have to do this if I do the pos-processing directly in here
# export the reduced hit matrices
# first, get the basis of the reduced data frame
to_select_names = list(filtered_results.index)
all_names = map(os.path.basename, training_file_names)
rows_to_select = []
for name in to_select_names:
	rows_to_select.append(all_names.index(name))
hit_matrices_dict = dict()
# the reduce the hit matrix to this basis
for i in range(len(k_range)):
	k_size = k_range[i]
	hit_matrices_dict['k=%d' % k_size] = hit_matrices[i][rows_to_select, :]
# then export  # TODO: not necessary if I do the post-processing right here
savemat(npz_file+"reduced_hit_matrix.npz", hit_matrices_dict, appendmat=False, do_compression=True)


###############################################################################
if verbose and not sensitive:
	print("Starting the post-processing")
	t0 = timeit.default_timer()

if not sensitive:
	# Do the post-processing
	# Make the hit matrices dense
	hit_matrices_dense_dict = dict()
	for k_size in k_range:
		hit_matrices_dense_dict['k=%d' % k_size] = hit_matrices_dict['k=%d' % k_size].todense()

	hit_matrices_dict = hit_matrices_dense_dict

	# get the count estimators of just the organisms of interest
	CEs = MH.import_multiple_from_single_hdf5(training_data, import_list=to_select_names)  # TODO: could make it a tad more memory efficient by sub-selecting the 'sketches'

	all_kmers_with_counts = dict()
	is_unique_kmer = set()
	is_unique_kmer_per_ksize = dict()
	for k_size in k_range:
		is_unique_kmer_per_ksize[k_size] = set()
		for i in range(len(CEs)):
			for big_kmer in CEs[i]._kmers:
				kmer = big_kmer[:k_size]
				if kmer in all_kmers_with_counts:
					all_kmers_with_counts[kmer] += 1
				else:
					all_kmers_with_counts[kmer] = 1

	for kmer in all_kmers_with_counts.keys():
		if all_kmers_with_counts[kmer] == 1:
			k_size = len(kmer)
			is_unique_kmer_per_ksize[k_size].add(kmer)
			is_unique_kmer.add(kmer)

	num_unique = dict()
	for i in range(len(CEs)):
		for k_size in k_range:
			current_kmers = [k[:k_size] for k in CEs[i]._kmers]
			current_kmers_set = set(current_kmers)
			non_unique = set()

			for kmer in current_kmers:
				if kmer not in is_unique_kmer_per_ksize[k_size]:
					non_unique.add(kmer)

			to_zero_indicies = [ind for ind, kmer in enumerate(current_kmers) if kmer in non_unique]
			hit_matrices_dict['k=%d' % k_size][i, to_zero_indicies] = 0  # set these to zero since they show up in other sketches (so not informative)
			num_unique[i, k_range.index(k_size)] = len(current_kmers_set) - len(non_unique)  # keep track of the size of the unique k-mers

	# sum the modified hit matrices to get the size of the intersection
	containment_indices = np.zeros((len(to_select_names), len(k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
	for k_size_loc in range(len(k_range)):
		k_size = k_range[k_size_loc]
		containment_indices[:, k_size_loc] = (
			hit_matrices_dict['k=%d' % k_size].sum(axis=1).ravel())  # /float(num_hashes))

	# then normalize by the number of unique k-mers (to get the containment index)
	for k_size_loc in range(len(k_range)):
		k_size = k_range[k_size_loc]
		for hash_loc in np.where(containment_indices[:, k_size_loc])[0]:  # find the genomes with non-zero containment
			unique_kmers = set()
			for kmer in CEs[hash_loc]._kmers:
				unique_kmers.add(kmer[:k_size])  # find the unique k-mers
			containment_indices[hash_loc, k_size_loc] /= float(
				len(unique_kmers))  # TODO: this doesn't seem like the right way to normalize, but apparently it is!
	# containment_indices[hash_loc, k_size_loc] /= float(num_unique[hash_loc, k_size_loc])  # divide by the unique num of k-mers

	results2 = dict()
	for k_size_loc in range(len(k_range)):
		ksize = k_range[k_size_loc]
		key = 'k=%d' % ksize
		results2[key] = containment_indices[:, k_size_loc]
	df2 = pd.DataFrame(results2, map(os.path.basename, to_select_names))
	df2 = df2.reindex(labels=['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
	sort_key = 'k=%d' % k_range[location_of_thresh]
	max_key = 'k=%d' % k_range[-1]
	filtered_results = df2[df2[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)  # only select those where the highest k-mer size's count is above the threshold

	filtered_results.to_csv(results_file, index=True, encoding='utf-8')
	if verbose:
		t1 = timeit.default_timer()
		print("Finished thresholding. Time: %f" % (t1 - t0))

