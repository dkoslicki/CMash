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
	def __init__(self, bloom_filter_file=None, TST_tree=None, k_range=None):
		self.bloom_filter_file = bloom_filter_file
		self.tree = TST_tree
		self.k_range = k_range
		pass  # TBD what needs to be passed


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

