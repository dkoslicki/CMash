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


class Query:
	def __init__(self, training_database_file=None, bloom_filter_file=None, TST_file=None, k_range=None):
		self.bloom_filter_file = bloom_filter_file
		self.TST_file = TST_file
		self.k_range = k_range
		self.training_database = training_database_file
		self.seen_kmers = set() # Seen k-mers (set of k-mers that already hit the trie, so don't need to check again)
		pass  # TBD what needs to be passed


	def import_TST(self):  # no more safety net for those that didn't create a TST properly with the CreateStreamingQueryDNADatabase.py
		self.tree = mt.Trie()
		self.tree.load(self.TST_file)

	def create_BF_prefilter(self):
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
	def __init__(self, tree=None, k_range=None, seen_kmers=None, all_kmers_bf=None):
		self.tree = tree
		self.k_range =k_range
		self.seen_kmers = seen_kmers
		self.all_kmers_bf = all_kmers_bf
		pass
	# This class is basically an array of counters (on the same basis as the sketches)
	# it's used to keep track (in a parallel friendly way) of which streamed k-mers went into the training file sketches
	#def __init__(self, training_database_file=None, bloom_filter_file=None, TST_file=None, k_range=None):
	#	super().__init__(training_database_file, bloom_filter_file, TST_file, k_range)

	def return_matches(self, input_kmer, k_size_loc):
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

	def process_seq(self, seq):
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



def main():
	"""
	Basically a bunch of simple command line tests
	"""
	top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	TST_file = os.path.join(top_dir, 'tests/TrainingDatabase.tst')
	training_database_file = os.path.join(top_dir, 'tests/TrainingDatabase.h5')
	k_range = [10, 12, 14, 16, 18, 20]

	Q = Query(training_database_file=training_database_file, bloom_filter_file=None, TST_file=TST_file, k_range=k_range)

	# test import of TST
	Q.import_TST()
	print(f"number of keys in tree: {len(Q.tree.keys())}")

	# test creation of BF
	Q.create_BF_prefilter()
	print(f"Number of buckets in BF: {Q.all_kmers_bf.buckets()}")

	# test import of Counters class

	C = Counters(tree=Q.tree, k_range=Q.k_range, seen_kmers=Q.seen_kmers, all_kmers_bf=Q.all_kmers_bf)
	print(f"Processed sequence with counter: {C.process_seq('AGTCCGCGCCACTGGCAGTGACCATCGACACGCAGACGGAGATTAACAACATTGTACTGGTCAATGATACCGGTATGCCG')}")



# simple way to do testing
if __name__ == "__main__":
	main()
