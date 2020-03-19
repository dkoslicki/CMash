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


class Create:
	def __init__(self, training_database_file: str, bloom_filter_file: str, TST_file: str, k_range: list):
		self.bloom_filter_file = bloom_filter_file
		self.TST_file = TST_file
		self.k_range = k_range
		self.training_database = training_database_file
		self.seen_kmers = set() # Seen k-mers (set of k-mers that already hit the trie, so don't need to check again)

	def import_TST(self) -> None:
		# no more safety net for those that didn't create a TST properly with the CreateStreamingQueryDNADatabase.py
		self.tree = mt.Trie()
		self.tree.load(self.TST_file)

	def create_BF_prefilter(self) -> None:
		tree = self.tree
		k_range = self.k_range
		if not self.bloom_filter_file:  # create one
			try:
				# Get all the k-mers in the TST, put them in a bloom filter
				# all_kmers_bf = WritingBloomFilter(len(sketches) * len(k_range) * num_hashes * 20, 0.01)
				self.all_kmers_bf = WritingBloomFilter(len(tree.keys()) * len(k_range) * 5, 0.01, ignore_case=True)  # fudge factor of 5 will make the BF larger, but also slightly faster
				for kmer_info in tree.keys():
					kmer = kmer_info.split('x')[0]  # remove the location information and just get the kmer
					for ksize in k_range:
						self.all_kmers_bf.add(kmer[0:ksize])
						self.all_kmers_bf.add(khmer.reverse_complement(kmer[0:ksize]))

			except IOError:
				print("No such file or directory/error opening file: %s" % self.bloom_filter_file)
				sys.exit(1)
		else:  # otherwise read it in
			try:
				self.all_kmers_bf = ReadingBloomFilter(self.bloom_filter_file)
			except IOError:
				print("No such file or directory/error opening file: %s" % self.bloom_filter_file)
				sys.exit(1)


# shared object that will update the intersection counts
class Counters:
	def __init__(self, tree: mt.Trie, k_range: list, seen_kmers: set, all_kmers_bf: WritingBloomFilter):
		self.tree = tree
		self.k_range =k_range
		self.seen_kmers = seen_kmers
		self.all_kmers_bf = all_kmers_bf

	# This class is basically an array of counters (on the same basis as the sketches)
	# it's used to keep track (in a parallel friendly way) of which streamed k-mers went into the training file sketches
	#def __init__(self, training_database_file=None, bloom_filter_file=None, TST_file=None, k_range=None):
	#	super().__init__(training_database_file, bloom_filter_file, TST_file, k_range)

	def return_matches(self, input_kmer: str, k_size_loc: int) -> tuple:
		""" Get all the matches in the trie with the kmer prefix"""
		match_info = set()
		to_return = []
		saw_match = False
		tree = self.tree

		# look for matches to both the kmer and its reverse complement in the TST as we can't assume
		# directionality of reads (and training database is constructed without reverse complements)
		for kmer in [input_kmer, khmer.reverse_complement(input_kmer)]:
			prefix_matches = tree.keys(kmer)  # get all the k-mers whose prefix matches
			# get the location of the found kmers in the counters
			for item in prefix_matches:
				split_string = item.split('x')  # first is the hash location, second is which k-mer
				hash_loc = int(split_string[1])
				kmer_loc = int(split_string[2])
				match_info.add((hash_loc, k_size_loc, kmer_loc))
			saw_match = False
			if match_info:
				saw_match = True
				for tup in match_info:
					to_return.append(tup)
			if saw_match:  # Only need to see a match to the original kmer or the reverse complement, don't return both otherwise you over-count
				break
		return to_return, saw_match

	def process_seq(self, seq: str) -> list:
		k_range = self.k_range
		seen_kmers = self.seen_kmers
		all_kmers_bf = self.all_kmers_bf
		#  start with small kmer size, if see match, then continue looking for longer k-mer sizes, otherwise move on
		small_k_size = k_range[0]  # start with the small k-size
		to_return = []
		for i in range(len(seq) - small_k_size + 1):  # look at all k-mers
			kmer = seq[i:i + small_k_size]
			possible_match = False
			if kmer not in seen_kmers:  # if we should process it
				if kmer in all_kmers_bf:  # if we should process it
					match_list, saw_match = self.return_matches(kmer, 0)
					if saw_match:
						seen_kmers.add(kmer)
						seen_kmers.add(khmer.reverse_complement(kmer))
						to_return.extend(match_list)
					possible_match = True
			# TODO: note: I could (since it'd only be for a single kmer size, keep a set of *all* small_kmers I've tried and use this as another pre-filter
			else:
				possible_match = True  # FIXME: bug introduced here in cf64b7aace5eadf738b920109d6419c9d930a1dc

			# start looking at the other k_sizes, don't overhang len(seq)
			if possible_match:
				for other_k_size in [x for x in k_range[1:] if i + x <= len(seq)]:
					kmer = seq[i:i + other_k_size]
					if kmer in all_kmers_bf:
						# if True:
						k_size_loc = k_range.index(other_k_size)
						match_list, saw_match = self.return_matches(kmer, k_size_loc)
						if saw_match:
							to_return.extend(match_list)
					else:
						pass  # if you didn't see a match at a smaller k-length, you won't at a larger one
		return to_return

# class to take the processed data and turn it into the containment indicies matrices
class Containment:
	def __init__(self, k_range: list, match_tuples: list, sketches: list, num_hashes: int):  # TODO: would like to indicate that sketches should be a list of CEs from MH
		self.k_range = k_range
		self.match_tuples = match_tuples
		self.sketches = sketches
		self.num_hashes = num_hashes

	def create_to_hit_matrices(self) -> None:
		k_range = self.k_range
		match_tuples = self.match_tuples
		sketches = self.sketches
		num_hashes = self.num_hashes

		# create k_range spare matrices. Rows index by genomes (sketch/hash index), columns index by k_mer_loc
		row_ind_dict = dict()
		col_ind_dict = dict()
		value_dict = dict()
		unique_kmers = dict()  # this will keep track of the unique k-mers seen in each genome (sketch/hash loc)
		for k_size in k_range:
			row_ind_dict[k_size] = []
			col_ind_dict[k_size] = []
			value_dict[k_size] = []

		match_tuples = set(match_tuples)  # uniquify, so we don't make the row/col ind dicts too large

		# convert the match tuples to the necessary format to be turned into a matrix/tensor
		for hash_loc, k_size_loc, kmer_loc in match_tuples:
			if hash_loc not in unique_kmers:
				unique_kmers[hash_loc] = set()
			k_size = k_range[k_size_loc]
			kmer = sketches[hash_loc]._kmers[kmer_loc][:k_size]
			if kmer not in unique_kmers[
				hash_loc]:  # if you've seen this k-mer before, don't add it. NOTE: this makes sure we don't over count
				row_ind_dict[k_size].append(hash_loc)
				col_ind_dict[k_size].append(kmer_loc)
				value_dict[k_size].append(1)  # only counting presence/absence, so just a 1 for the value
				unique_kmers[hash_loc].add(kmer)

		# list of matrices that contain the hits: len(hit_matrices) == k_sizes
		# each hit_matrices[i] has rows indexed by which genome/sketch they belong to
		# columns indexed by where the k-mer appeared in the sketch/hash list
		self.hit_matrices = []

		for k_size in k_range:
			# convert to matrices
			mat = csc_matrix((value_dict[k_size], (row_ind_dict[k_size], col_ind_dict[k_size])),
							 shape=(len(sketches), num_hashes))
			self.hit_matrices.append(mat)

	def create_containment_indicies(self) -> None:
		sketches = self.sketches
		k_range = self.k_range
		hit_matrices = self.hit_matrices
		# prep the containment indicies matrix: rows are genome/sketch, one column for each k-mer size in k_range
		self.containment_indices = np.zeros((len(sketches), len(
			k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
		for k_size_loc in range(len(k_range)):
			# sum over the columns: i.e. total up the number of matches in the sketch/hash list
			self.containment_indices[:, k_size_loc] = (hit_matrices[k_size_loc].sum(axis=1).ravel())  # /float(num_hashes))

		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			for hash_loc in np.where(self.containment_indices[:, k_size_loc])[
				0]:  # find the genomes with non-zero containment
				unique_kmers = set()
				for kmer in sketches[hash_loc]._kmers:
					# find the unique k-mers: for smaller k-mer truncation sizes, the length of the sketch might have
					# been reduced due to duplicates, so adjust for this factor to get the correct denominator
					unique_kmers.add(kmer[:k_size])
				self.containment_indices[hash_loc, k_size_loc] /= float(
					len(unique_kmers))  # divide by the unique num of k-mers

	def create_data_frame(self, training_file_names: list, location_of_thresh: int, coverage_threshold: float) -> None:
		k_range = self.k_range
		containment_indices = self.containment_indices
		results = dict()
		for k_size_loc in range(len(k_range)):
			ksize = k_range[k_size_loc]
			key = 'k=%d' % ksize
			results[key] = containment_indices[:, k_size_loc]
		df = pd.DataFrame(results, map(os.path.basename, training_file_names))
		df = df.reindex(labels=['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
		sort_key = 'k=%d' % k_range[location_of_thresh]
		max_key = 'k=%d' % k_range[-1]
		self.filtered_results = df[df[sort_key] > coverage_threshold].sort_values(max_key,
																			 ascending=False)  # only select those where the highest k-mer size's count is above the threshold

class PostProcess:
	def __init__(self, fil):
		pass

	def prepare_post_process(self):
		to_select_names = list(filtered_results.index)
		all_names = list(map(os.path.basename, training_file_names))
		rows_to_select = []
		for name in to_select_names:
			rows_to_select.append(all_names.index(name))
		hit_matrices_dict = dict()
		# the reduce the hit matrix to this basis
		for i in range(len(k_range)):
			k_size = k_range[i]
			hit_matrices_dict['k=%d' % k_size] = hit_matrices[i][rows_to_select, :]

		# Make the hit matrices dense
		hit_matrices_dense_dict = dict()
		for k_size in k_range:
			hit_matrices_dense_dict['k=%d' % k_size] = hit_matrices_dict['k=%d' % k_size].todense()

		hit_matrices_dict = hit_matrices_dense_dict

	def find_kmers_in_filtered_results(self):
		# get the count estimators of just the organisms of interest
		CEs = MH.import_multiple_from_single_hdf5(training_database_file,
												  import_list=to_select_names)  # TODO: could make it a tad more memory efficient by sub-selecting the 'sketches'

		# get all the kmers (for each kmer size) and form their counts in the subset of predicted sketches to be in the sample
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

	def find_unique_kmers(self):
		# Use this to identify which k-mers are unique (i.e. show up in exactly one sketch)
		for kmer in all_kmers_with_counts.keys():
			if all_kmers_with_counts[kmer] == 1:
				k_size = len(kmer)
				is_unique_kmer_per_ksize[k_size].add(kmer)
				is_unique_kmer.add(kmer)

	def find_non_unique_kmers(self):
		# Also keep track of which kmers appear in more than one sketch (not unique)
		num_unique = dict()
		for i in range(len(CEs)):
			for k_size in k_range:
				current_kmers = [k[:k_size] for k in CEs[i]._kmers]
				current_kmers_set = set(current_kmers)
				non_unique = set()

				for kmer in current_kmers:
					if kmer not in is_unique_kmer_per_ksize[k_size]:
						non_unique.add(kmer)
				# reduce the hit matrices by removing the hits corresponding to non-unique k-mers
				to_zero_indicies = [ind for ind, kmer in enumerate(current_kmers) if kmer in non_unique]
				# if you use a really small initial kmer size, some of the hit matrices may be empty
				# due to all k-mers being shared in common
				if hit_matrices_dict['k=%d' % k_size].size > 0:
					hit_matrices_dict['k=%d' % k_size][
						i, to_zero_indicies] = 0  # set these to zero since they show up in other sketches (so not informative)
				num_unique[i, k_range.index(k_size)] = len(current_kmers_set) - len(
					non_unique)  # keep track of the size of the unique k-mers

	def create_post_containment_indicies(self):
		# sum the modified hit matrices to get the size of the intersection
		containment_indices = np.zeros((len(to_select_names), len(
			k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			containment_indices[:, k_size_loc] = (
				hit_matrices_dict['k=%d' % k_size].sum(axis=1).ravel())  # /float(num_hashes))

		# then normalize by the number of unique k-mers (to get the containment index)
		# In essence, this is the containment index, restricted to unique k-mers. This effectively increases the specificity,
		# but also increases the variance/confidence interval, since this decreases the size of the sketch.
		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			for hash_loc in np.where(containment_indices[:, k_size_loc])[
				0]:  # find the genomes with non-zero containment
				unique_kmers = set()
				for kmer in CEs[hash_loc]._kmers:
					unique_kmers.add(kmer[:k_size])  # find the unique k-mers
				containment_indices[hash_loc, k_size_loc] /= float(len(
					unique_kmers))  # FIXME: this doesn't seem like the right way to normalize, but apparently it is!
				# containment_indices[hash_loc, k_size_loc] /= float(num_unique[hash_loc, k_size_loc])  # FIXME: in small tests, this seems to give better results. To be revisted.

	def create_data_frame(self):
		# spit out these results
		results = dict()
		for k_size_loc in range(len(k_range)):
			ksize = k_range[k_size_loc]
			key = 'k=%d' % ksize
			results[key] = containment_indices[:, k_size_loc]
		df = pd.DataFrame(results, map(os.path.basename, to_select_names))
		df = df.reindex(labels=['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
		sort_key = 'k=%d' % k_range[location_of_thresh]
		max_key = 'k=%d' % k_range[-1]
		# TODO: might not want to filter at this point again, or adjust the coverage_threshold since now it's only based on unique k-mers
		filtered_results = df[df[sort_key] > coverage_threshold].sort_values(max_key,
																			 ascending=False)  # only select those where the highest k-mer size's count is above the threshold



def main():
	"""
	Basically a bunch of simple command line tests
	"""
	top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	TST_file = os.path.join(top_dir, 'tests/TrainingDatabase.tst')
	training_database_file = os.path.join(top_dir, 'tests/TrainingDatabase.h5')
	k_range = [10, 12, 14, 16, 18, 20]

	C = Create(training_database_file=training_database_file, bloom_filter_file="", TST_file=TST_file, k_range=k_range)

	# test import of TST
	C.import_TST()
	print(f"number of keys in tree: {len(C.tree.keys())}")

	# test creation of BF
	C.create_BF_prefilter()
	print(f"Number of buckets in BF: {C.all_kmers_bf.buckets()}")

	# test import of Counters class

	counters = Counters(tree=C.tree, k_range=C.k_range, seen_kmers=C.seen_kmers, all_kmers_bf=C.all_kmers_bf)
	print(f"Processed sequence with counter: {counters.process_seq('AGTCCGCGCCACTGGCAGTGACCATCGACACGCAGACGGAGATTAACAACATTGTACTGGTCAATGATACCGGTATGCCG')}")





# simple way to do testing
if __name__ == "__main__":
	main()
