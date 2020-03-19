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
	def __init__(self, training_database=None, bloom_filter_file=None, TST_file=None, k_range=None):
		self.bloom_filter_file = bloom_filter_file
		self.TST_file = TST_file
		self.k_range = k_range
		self.training_database = training_database
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



def main():
	"""
	Basically a bunch of simple command line tests
	"""
	top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	TST_file = os.path.join(top_dir, 'tests/TrainingDatabase.tst')
	training_database = os.path.join(top_dir, 'tests/TrainingDatabase.h5')
	k_range = [10, 12, 14, 16, 18, 20]

	Q = Query(training_database=training_database, bloom_filter_file=None, TST_file=TST_file, k_range=k_range)
	Q.import_TST()
	print(f"number of keys in tree: {len(Q.tree.keys())}")

# simple way to do testing
if __name__ == "__main__":
	main()
