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
import multiprocessing

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

# Note: here, I am using canonical k-mers: i.e. to disambiguate between a k-mer and it's reverse complement, I
# simply take as the representative whichever one is lexicographically smaller.


notACTG = re.compile('[^ACTG]')  # look for any not ACTG


class TrueContainment:
	"""
	This class has functionality to compute the ground truth containment indicies and return them in the same format
	as the scripts (to ease future testing). It is **only** intended for:
	1. Small training databases
	2. Databases that were formed using genomes that you have direct access to (i.e. live on your file system)
	"""

	def __init__(self, training_database_file: str, k_sizes: str):
		self.training_database_file = training_database_file
		self.k_sizes = self.__parseNumList(k_sizes)
		self.CEs = self.__import_database()
		self.training_file_names = self.__return_file_names()
		self.training_file_to_ksize_to_kmers = self.__compute_all_training_kmers()

	def __import_database(self) -> list:
		CEs = MH.import_multiple_from_single_hdf5(self.training_database_file)
		return CEs

	def __return_file_names(self):
		training_file_names = list(map(lambda x: x.input_file_name.decode('utf-8'), self.CEs))
		return training_file_names

	def __parseNumList(self, k_sizes_str: str) -> list:
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
	def __kmers(seq, ksize):
		"""yield all k-mers of len ksize from seq.
		Returns an iterable object
		"""
		for i in range(len(seq) - ksize + 1):
			yield seq[i:i + ksize]

	def _return_ksize_to_kmers(self, input_file: str) -> dict:
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
					for kmer in self.__kmers(seq, k_size):
						if kmer:
							# Use canonical k-mers
							temp_kmer = kmer
							temp_kmer_rc = khmer.reverse_complement(kmer)
							if temp_kmer < temp_kmer_rc:
								k_size_to_kmers[k_size].add(temp_kmer)  # add the kmer
							else:
								k_size_to_kmers[k_size].add(temp_kmer_rc)  # add the reverse complement
			# otherwise, we need to do the same thing for each of the subsequences
			else:
				for sub_seq in seq_split_onlyACTG:
					if sub_seq:
						for k_size in k_sizes:
							for kmer in self.__kmers(seq, k_size):
								if kmer:
									# Use canonical k-mers
									temp_kmer = kmer
									temp_kmer_rc = khmer.reverse_complement(kmer)
									if temp_kmer < temp_kmer_rc:
										k_size_to_kmers[k_size].add(temp_kmer)  # add the kmer
									else:
										k_size_to_kmers[k_size].add(temp_kmer_rc)  # add the reverse complement
		return k_size_to_kmers

	@staticmethod
	def __return_containment_index(set1: set, set2: set):
		return len(set1.intersection(set2)) / float(len(set1))

	def __compute_all_training_kmers(self):
		training_file_to_ksize_to_kmers = dict()
		num_threads = multiprocessing.cpu_count()
		pool = multiprocessing.Pool(processes=num_threads)
		# res is returned in the same order as self.training_file_names according to the docs
		res = pool.map(self._return_ksize_to_kmers, self.training_file_names)
		for (item, file_name) in zip(res, self.training_file_names):
			training_file_to_ksize_to_kmers[file_name] = item
		pool.close()
		return training_file_to_ksize_to_kmers

	def __return_containment_indicies(self, query_file: str) -> np.ndarray:
		training_file_names = self.training_file_names
		k_sizes = self.k_sizes
		training_file_to_ksize_to_kmers = self.training_file_to_ksize_to_kmers
		num_files = len(training_file_names)
		# rows are the files, columns are the k-mer sizes
		containment_indicies = np.zeros((num_files, len(k_sizes)))
		# if the query file is part of the training files, then nothing extra to do
		if query_file in training_file_names:
			for (j, k_size) in enumerate(k_sizes):
				query_kmers = training_file_to_ksize_to_kmers[query_file][k_size]
				for (i, file_name) in enumerate(training_file_names):
					training_kmers = training_file_to_ksize_to_kmers[file_name][k_size]
					# | train \cap query| / | train |
					containment_indicies[i, j] = self.__return_containment_index(training_kmers, query_kmers)
		else:
			# need to compute the k-mers in the query file
			query_file_to_ksize_to_kmers = self._return_ksize_to_kmers(query_file)
			for (j, k_size) in enumerate(k_sizes):
				query_kmers = query_file_to_ksize_to_kmers[query_file][k_size]
				for (i, file_name) in enumerate(training_file_names):
					training_kmers = training_file_to_ksize_to_kmers[file_name][k_size]
					# | train \cap query| / | train |
					containment_indicies[i, j] = self.__return_containment_index(training_kmers, query_kmers)
		return containment_indicies

	def return_containment_data_frame(self, query_file: str, location_of_thresh: int, coverage_threshold: float) -> pd.DataFrame:
		k_range = self.k_sizes
		training_file_names = self.training_file_names
		containment_indices = self.__return_containment_indicies(query_file)
		df = Query.return_data_frame(training_file_names=training_file_names,
		                             k_range=k_range,
		                             location_of_thresh=location_of_thresh,
		                             containment_indices=containment_indices,
		                             coverage_threshold=coverage_threshold)
		return df


def main():
	training_database_file = "/home/dkoslicki/Desktop/CMash/tests/script_tests/TrainingDatabase.h5"
	query_file = "/home/dkoslicki/Desktop/CMash/tests/Organisms/taxid_1192839_4_genomic.fna.gz"
	k_range = "4-5-1"  # small test
	#k_range = "10-21-2"  # same as in tests/script_tests/run_small_tests.sh
	g = TrueContainment(training_database_file, k_range)
	print(g.return_containment_data_frame(query_file, -1, .1))  # defaults



if __name__ == "__main__":
	main()
