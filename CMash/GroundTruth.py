import khmer
import marisa_trie as mt
import numpy as np
import os
import sys
import pandas as pd
from hydra import WritingBloomFilter, ReadingBloomFilter
from scipy.sparse import csc_matrix
import re
import screed
from argparse import ArgumentTypeError
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
	from CMash import Query
except ImportError:
	try:
		import MinHash as MH
		import Query
	except ImportError:
		sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
		from CMash import MinHash as MH
		from CMash import Query


# FIXME: Make sure I am *identifying* the rc k-mers with the k-mers and not counting them as distinct entities


notACTG = re.compile('[^ACTG]')  # look for any not ACTG

class TrueContainment:
	"""
	This class has functionality to compute the ground truth containment indicies and return them in the same format
	as the scripts (to ease future testing). It is **only** intended for:
	1. Small training databases
	2. Databases that were formed using genomes that you have direct access to (i.e. live on your file system)
	"""
	def __init__(self, training_database_file: str, query_file: str, k_sizes: str):
		self.training_database_file = training_database_file
		self.query_file = query_file
		self.CEs = []
		self.k_sizes = self.parseNumList(k_sizes)
		self.query_kmers = dict()

	def import_database(self):
		self.CEs = MH.import_multiple_from_single_hdf5(self.training_database_file)

	def parseNumList(self, k_sizes_str: str) -> list:
		"""
		Parses a string like 10-21-1 and turn it into a list like [10, 11, 12,...,21]
		:param k_sizes_str: the <start>-<end>-<increment> string
		:type k_sizes_str: str
		:return: list of k-mer sizes
		:rtype: list
		"""
		m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', k_sizes_str)
		# ^ (or use .split('-'). anyway you like.)
		if not m:
			raise ArgumentTypeError(
				"'" + k_sizes_str + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
		start = int(m.group(1))
		end = int(m.group(2))
		if m.group(3):
			increment = int(m.group(3))
		else:
			increment = 1
		return list(range(start, end + 1, increment))

	@staticmethod
	def kmers(seq, ksize):
		"""yield all k-mers of len ksize from seq.
		Returns an iterable object
		"""
		for i in range(len(seq) - ksize + 1):
			yield seq[i:i + ksize]

	def return_kmers(self, input_file):
		k_sizes = self.k_sizes
		k_size_to_kmers = dict()
		# initialize all the k-mer sizes for the query file
		for k_size in k_sizes:
			k_size_to_kmers[k_size] = set()

		# iterate over the file only once
		for record in screed.open(input_file):
			seq = record.sequence  # get the sequence
			seq = seq.upper()  # convert to upper-case
			seq_split_onlyACTG = notACTG.split(seq)  # split it on non-ACTG sequences
			# if there are no non-ACTG's we get only a single one
			if len(seq_split_onlyACTG) == 1:
				for k_size in k_sizes:
					for kmer in self.kmers(seq, k_size):
						if kmer:
							k_size_to_kmers[k_size].add(kmer)  # add the kmer
							k_size_to_kmers[k_size].add(khmer.reverse_complement(kmer))  # add the reverse complement
			# otherwise, we need to do the same thing for each of the subsequences
			else:
				for sub_seq in seq_split_onlyACTG:
					if sub_seq:
						for k_size in k_sizes:
							for kmer in self.kmers(seq, k_size):
								if kmer:
									k_size_to_kmers[k_size].add(kmer)  # add the kmer
									k_size_to_kmers[k_size].add(khmer.reverse_complement(kmer))  # add the reverse complement
		return k_size_to_kmers





def main():
	training_database_file = "/home/dkoslicki/Desktop/CMash/tests/script_tests/TrainingDatabase.h5"
	query_file = "/home/dkoslicki/Desktop/CMash/tests/Organisms/taxid_1192839_4_genomic.fna.gz"

if __name__ == "__main__":
	main()