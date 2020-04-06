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
import argparse

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
		"""
		Private function that imports the HDF5 training file.
		:return: a list of CountEstimators
		:rtype: MinHash.CountEstimator
		"""
		CEs = MH.import_multiple_from_single_hdf5(self.training_database_file)
		return CEs

	def __return_file_names(self):
		"""
		Private function that gets all the files names contained in the training data.
		:return: a list of file names
		:rtype: list
		"""
		training_file_names = list(map(lambda x: x.input_file_name.decode('utf-8'), self.CEs))
		return training_file_names

	@staticmethod
	def __parseNumList(k_sizes_str: str) -> list:
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
		"""
		yield all k-mers of len ksize from seq.
		Returns an iterable object
		:param seq: a DNA sequence
		:type seq: str
		:param ksize: a k-mer size
		:type ksize: int
		"""
		for i in range(len(seq) - ksize + 1):
			yield seq[i:i + ksize]

	def _return_ksize_to_kmers(self, input_file: str) -> dict:
		"""
		Enumerates all the k-mers specified by self.k_sizes in the genome/metagenome specified by input_file.
		:param input_file: a file path to a fna/fq/gzipped file containing DNA sequences
		:type input_file: str
		:return: a dictionary with keys corresponding to k-mer sizes in self.k_sizes, and values dictionaries containing canonical k-mers
		:rtype: dict
		"""
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
		"""
		Computes the containment index
		:param set1: a set of k-mers
		:type set1: set
		:param set2: another set of k-mers
		:type set2: set
		:return: containment index  |set1 \cap set2| / |set 1|
		:rtype: float
		"""
		return len(set1.intersection(set2)) / float(len(set1))

	def __compute_all_training_kmers(self):
		"""
		In a parallelized fashion, enumerate all the k-mers for the given self.k_sizes in the input training genomes.
		:return: a dictionary with keys given by self.training_file_names, values are dictionaries: keys are k_sizes, values are sets of canonical k-mers
		:rtype: dict
		"""
		training_file_to_ksize_to_kmers = dict()
		num_threads = multiprocessing.cpu_count()
		pool = multiprocessing.Pool(processes=int(min(num_threads, len(self.training_file_names))))
		# res is returned in the same order as self.training_file_names according to the docs
		res = pool.map(self._return_ksize_to_kmers, self.training_file_names)
		for (item, file_name) in zip(res, self.training_file_names):
			training_file_to_ksize_to_kmers[file_name] = item
		pool.close()
		return training_file_to_ksize_to_kmers

	def __return_containment_indicies(self, query_file: str) -> np.ndarray:
		"""
		Creates a matrix of containment indicies:
			for each i in self.training_file_names:
				for each k in k_sizes:
					containment_indicies[i ,k] = |query_file_k-mers \cap training_file_i_k-mers| / |training_file_i_k-mers|
		:param query_file: a file pointing to a fasta/q (maybe compressed) file
		:type query_file: str
		:return: a numpy matrix of containment indicies: containment_indicies[i ,k] = |query_file_k-mers \cap training_file_i_k-mers| / |training_file_i_k-mers|
		:rtype: np.ndarray
		"""
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
				query_kmers = query_file_to_ksize_to_kmers[k_size]
				for (i, file_name) in enumerate(training_file_names):
					training_kmers = training_file_to_ksize_to_kmers[file_name][k_size]
					# | train \cap query| / | train |
					containment_indicies[i, j] = self.__return_containment_index(training_kmers, query_kmers)
		return containment_indicies

	def return_containment_data_frame(self, query_file: str, location_of_thresh: int, coverage_threshold: float) -> pd.DataFrame:
		"""
		Returns a Pandas Data frame with rows indexed by training file names, columns indicated by k-mer sizes, and entries the
		containment indicies for the give query_file. Same exact format as CMash/Query.py and scripts/StreamingQueryDNADatabase.py
		:param query_file: a file pointing to a fasta/q (maybe compressed) file. Need not be in the training data
		:type query_file: str
		:param location_of_thresh: where in self.k_sizes the thresholding should take place (-1 means the last one)
		:type location_of_thresh: int
		:param coverage_threshold: filter out those results that have containment indicies strictly below this threshold
		:type coverage_threshold: float
		:return: Returns a Pandas Data frame with rows indexed by training file names, columns indicated by k-mer sizes, and entries the
		containment indicies for the give query_file.
		:rtype: pandas.DataFrame
		"""
		k_range = self.k_sizes
		training_file_names = self.training_file_names
		containment_indices = self.__return_containment_indicies(query_file)
		df = Query.return_data_frame(training_file_names=training_file_names,
		                             k_range=k_range,
		                             location_of_thresh=location_of_thresh,
		                             containment_indices=containment_indices,
		                             coverage_threshold=coverage_threshold)
		return df


class TrueContainmentKMC:
	"""
	This class has functionality to compute the ground truth containment indicies and return them in the same format
	as the scripts (to ease future testing). It is only intended for:
	1. Small-ish training databases
	2. Databases that were formed using genomes that you have direct access to (i.e. live on your file system)
	"""

	def __init__(self, training_database_file: str, k_sizes: str):
		self.training_database_file = training_database_file
		self.k_sizes = self.__parseNumList(k_sizes)
		self.CEs = self.__import_database()
		self.training_file_names = self.__return_file_names()
		self.training_file_to_ksize_to_kmers = self.__compute_all_training_kmers()

	def __import_database(self) -> list:
		"""
		Private function that imports the HDF5 training file.
		:return: a list of CountEstimators
		:rtype: MinHash.CountEstimator
		"""
		CEs = MH.import_multiple_from_single_hdf5(self.training_database_file)
		return CEs

	def __return_file_names(self):
		"""
		Private function that gets all the files names contained in the training data.
		:return: a list of file names
		:rtype: list
		"""
		training_file_names = list(map(lambda x: x.input_file_name.decode('utf-8'), self.CEs))
		return training_file_names

	@staticmethod
	def __parseNumList(k_sizes_str: str) -> list:
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
		"""
		yield all k-mers of len ksize from seq.
		Returns an iterable object
		:param seq: a DNA sequence
		:type seq: str
		:param ksize: a k-mer size
		:type ksize: int
		"""
		for i in range(len(seq) - ksize + 1):
			yield seq[i:i + ksize]

	def _return_ksize_to_kmers(self, input_file: str) -> dict:
		"""
		Enumerates all the k-mers specified by self.k_sizes in the genome/metagenome specified by input_file.
		:param input_file: a file path to a fna/fq/gzipped file containing DNA sequences
		:type input_file: str
		:return: a dictionary with keys corresponding to k-mer sizes in self.k_sizes, and values dictionaries containing canonical k-mers
		:rtype: dict
		"""
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
		"""
		Computes the containment index
		:param set1: a set of k-mers
		:type set1: set
		:param set2: another set of k-mers
		:type set2: set
		:return: containment index  |set1 \cap set2| / |set 1|
		:rtype: float
		"""
		return len(set1.intersection(set2)) / float(len(set1))

	def __compute_all_training_kmers(self):
		"""
		In a parallelized fashion, enumerate all the k-mers for the given self.k_sizes in the input training genomes.
		:return: a dictionary with keys given by self.training_file_names, values are dictionaries: keys are k_sizes, values are sets of canonical k-mers
		:rtype: dict
		"""
		training_file_to_ksize_to_kmers = dict()
		num_threads = multiprocessing.cpu_count()
		pool = multiprocessing.Pool(processes=int(min(num_threads, len(self.training_file_names))))
		# res is returned in the same order as self.training_file_names according to the docs
		res = pool.map(self._return_ksize_to_kmers, self.training_file_names)
		for (item, file_name) in zip(res, self.training_file_names):
			training_file_to_ksize_to_kmers[file_name] = item
		pool.close()
		return training_file_to_ksize_to_kmers

	def __return_containment_indicies(self, query_file: str) -> np.ndarray:
		"""
		Creates a matrix of containment indicies:
			for each i in self.training_file_names:
				for each k in k_sizes:
					containment_indicies[i ,k] = |query_file_k-mers \cap training_file_i_k-mers| / |training_file_i_k-mers|
		:param query_file: a file pointing to a fasta/q (maybe compressed) file
		:type query_file: str
		:return: a numpy matrix of containment indicies: containment_indicies[i ,k] = |query_file_k-mers \cap training_file_i_k-mers| / |training_file_i_k-mers|
		:rtype: np.ndarray
		"""
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
				query_kmers = query_file_to_ksize_to_kmers[k_size]
				for (i, file_name) in enumerate(training_file_names):
					training_kmers = training_file_to_ksize_to_kmers[file_name][k_size]
					# | train \cap query| / | train |
					containment_indicies[i, j] = self.__return_containment_index(training_kmers, query_kmers)
		return containment_indicies

	def return_containment_data_frame(self, query_file: str, location_of_thresh: int, coverage_threshold: float) -> pd.DataFrame:
		"""
		Returns a Pandas Data frame with rows indexed by training file names, columns indicated by k-mer sizes, and entries the
		containment indicies for the give query_file. Same exact format as CMash/Query.py and scripts/StreamingQueryDNADatabase.py
		:param query_file: a file pointing to a fasta/q (maybe compressed) file. Need not be in the training data
		:type query_file: str
		:param location_of_thresh: where in self.k_sizes the thresholding should take place (-1 means the last one)
		:type location_of_thresh: int
		:param coverage_threshold: filter out those results that have containment indicies strictly below this threshold
		:type coverage_threshold: float
		:return: Returns a Pandas Data frame with rows indexed by training file names, columns indicated by k-mer sizes, and entries the
		containment indicies for the give query_file.
		:rtype: pandas.DataFrame
		"""
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
	parser = argparse.ArgumentParser(
		description="This script calculates the *ground truth* containment indicies for each of the training/reference sketches"
		            " via brute force enumeration of all the (canonical) k-mers", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-c', '--containment_threshold', type=float,
	                    help="Only return results with containment index above this "
	                         "threshold at the maximum k-mer size.", default=0.1)
	parser.add_argument('-l', '--location_of_thresh', type=int,
	                    help="Location in range to apply the threshold passed by the -c flag. -l 2 -c 5-50-10 means the"
	                         " threshold will be applied at k-size 25. Default is largest size.", default=-1)
	parser.add_argument('in_file', help="Input file: FASTA/Q file to be processes")
	parser.add_argument('reference_file',
	                    help='Training database/reference file (in HDF5 format). Created with MakeStreamingDNADatabase.py')
	parser.add_argument('out_file', help='Output csv file with the containment indices.')
	parser.add_argument('k_range', type=str,
	                    help="Range of k-mer sizes in the formate <start>-<end>-<increment>."
	                         " So 5-10-2 means [5, 7, 9]. If <end> is larger than the k-mer size"
	                         "of the training data, this will automatically be reduced.")

	# get the args
	args = parser.parse_args()
	k_sizes = args.k_range
	training_database_file = args.reference_file
	query_file = args.in_file
	out_file = args.out_file
	location_of_thresh = args.location_of_thresh
	coverage_threshold = args.containment_threshold

	# pre-compute the kmers in the training database
	g = TrueContainment(training_database_file=training_database_file, k_sizes=k_sizes)

	# compute the containment indicies
	df = g.return_containment_data_frame(query_file=query_file, location_of_thresh=location_of_thresh, coverage_threshold=coverage_threshold)

	# save them
	df.to_csv(out_file, index=True, encoding='utf-8')


if __name__ == "__main__":
	main()
