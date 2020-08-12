import khmer
import numpy as np
import os
import sys
import pandas as pd
import re
import screed
from argparse import ArgumentTypeError
import multiprocessing
import argparse
import subprocess
import json
from itertools import starmap
import tempfile

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
	as the scripts (to ease future testing). It is only intended for:
	1. Small-ish training databases
	2. Databases that were formed using genomes that you have direct access to (i.e. live on your file system)
	"""

	def __init__(self, training_database_file: str, k_sizes: str, temp_dir: str):
		self.training_database_file = training_database_file
		self.k_sizes = self.__parseNumList(k_sizes)
		self.CEs = self.__import_database()
		self.training_file_names = self.__return_file_names()
		self.temp_dir = temp_dir
		if not os.path.exists(temp_dir):
			os.mkdir(temp_dir)
		# compute all the training k-mers up front
		self._compute_all_training_kmers()

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
	def _kmc_count(input_file_name: str, output_file_name: str, kmer_size: int, threads=1) -> None:
		"""
		Calls KMC to compute the k-mers for a given input file name
		:param input_file_name:
		:type input_file_name:
		:param output_file_name:
		:type output_file_name:
		:param kmer_size:
		:type kmer_size:
		"""
		input_types = ['-fm', '-fq', '-fa', '-fbam']
		success = False
		for input_type in input_types:
			res = subprocess.run(f"kmc -k{kmer_size} {input_type} -r -t{threads} -ci0 -cs3 -j{output_file_name}.log {input_file_name} {output_file_name} .", shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
			if res.returncode == 0:
				success = True
				break
		if not success:
			raise Exception(f"Unknown sequence format: must be one of multifasta, fastq, fasta, or BAM (gzipped or uncompressed). Culprit file is {input_file_name}. Command was {res.args}")

	@staticmethod
	def _kmc_return_distinct_kmers(kmc_log_file: str) -> int:
		"""
	Parses the KMC log file to return the number of distinct k-mers
		:param kmc_log_file:
		:type kmc_log_file:
		:return:
		:rtype:
		"""
		with open(kmc_log_file, 'r') as fid:
			res = json.load(fid)
			return res['Stats']['#Unique_k-mers']

	@staticmethod
	def _kmc_return_intersection_count(kmc_input_file1: str, kmc_input_file2: str) -> int:
		"""
		Takes two kmc counted files, returns the number of k-mers in their intersection
		:param kmc_input_file1:
		:type kmc_input_file1:
		:param kmc_input_file2:
		:type kmc_input_file2:
		"""
		dir_name = os.path.dirname(kmc_input_file1)
		intersect_file = os.path.join(dir_name, f"{os.path.basename(kmc_input_file1)}_intersect_{os.path.basename(kmc_input_file2)}")
		dump_file = os.path.join(dir_name, f"{os.path.basename(kmc_input_file1)}_intersect_{os.path.basename(kmc_input_file2)}_dump")
		res = subprocess.run(f"kmc_tools simple {kmc_input_file1} -ci1 {kmc_input_file2} -ci1 intersect {intersect_file}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")
		res = subprocess.run(f"kmc_dump {intersect_file} {dump_file}; cat {dump_file} | wc -l", shell=True, capture_output=True)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")

		intersect_count = int(res.stdout)

		# then clean up the mess
		os.remove(f"{intersect_file}.kmc_pre")
		os.remove(f"{intersect_file}.kmc_suf")
		os.remove(dump_file)

		return intersect_count

	def __kmc_output_name_converter(self, input_file: str, k_size: str) -> str:
		temp_dir = self.temp_dir
		return f"{os.path.join(temp_dir, os.path.basename(input_file))}_k_{k_size}"

	def _compute_all_training_kmers(self):
		num_threads = 48
		to_compute = []
		# create the tuples to be computed on: (input file, ouput_kmc_file, k_kmer_size)
		for training_file in self.training_file_names:
			for k_size in self.k_sizes:
				output_file = self.__kmc_output_name_converter(training_file, k_size)
				to_compute.append((training_file, output_file, k_size))
		pool = multiprocessing.Pool(processes=int(min(num_threads, len(self.training_file_names))))
		res = pool.starmap(self._kmc_count, to_compute)
		# consume everything so we know the process has completed
		for it in res:
			pass
		pool.close()

	def _return_containment_index(self, query_file: str, i: int, j: int) -> tuple:
		k_size = self.k_sizes[j]
		training_file = self.training_file_names[i]
		training_kmc_output = self.__kmc_output_name_converter(training_file, k_size)
		query_kmc_output = self.__kmc_output_name_converter(query_file, k_size)
		numerator = self._kmc_return_intersection_count(query_kmc_output, training_kmc_output)
		denomenator = self._kmc_return_distinct_kmers(f"{training_kmc_output}.log")
		return (i, j, numerator / float(denomenator))  # | train \cap query| / | train |

	def _return_containment_indicies(self, query_file: str) -> np.ndarray:
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
		# compute the k-mers in the query file
		for k_size in k_sizes:
			# store the query file kmc outputs to a dict for future convenience
			self._kmc_count(query_file, self.__kmc_output_name_converter(query_file, k_size), k_size, threads=48)

		# compute the containment indicies
		# rows are the files, columns are the k-mer sizes
		containment_indicies = np.zeros((len(training_file_names), len(k_sizes)))

		# serial version
		#for (j, k_size) in enumerate(k_sizes):
		#	query_kmc_output = self.__kmc_output_name_converter(query_file, k_size)
		#	for (i, training_file) in enumerate(training_file_names):
		#		training_kmc_output = self.__kmc_output_name_converter(training_file, k_size)
		#		numerator = self._kmc_return_intersection_count(query_kmc_output, training_kmc_output)
		#		denomenator = self._kmc_return_distinct_kmers(f"{training_kmc_output}.log")
		#		containment_indicies[i, j] = numerator / float(denomenator)  # | train \cap query| / | train |

		# parallel version
		to_compute = []
		for i in range(len(training_file_names)):
			for j in range(len(k_sizes)):
				to_compute.append((query_file, i, j))
		pool = multiprocessing.Pool(processes=int(min(multiprocessing.cpu_count()/float(4), len(self.training_file_names))))
		res = pool.starmap(self._return_containment_index, to_compute)
		for (i, j, ci) in res:
			containment_indicies[i, j] = ci
		pool.close()

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
		containment_indices = self._return_containment_indicies(query_file)
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
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
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
	temp_dir = tempfile.TemporaryDirectory()
	# pre-compute the kmers in the training database
	g = TrueContainment(training_database_file=training_database_file, k_sizes=k_sizes, temp_dir=temp_dir.name)

	# compute the containment indicies
	df = g.return_containment_data_frame(query_file=query_file, location_of_thresh=location_of_thresh, coverage_threshold=coverage_threshold)

	# save them
	df.to_csv(out_file, index=True, encoding='utf-8')


if __name__ == "__main__":
	main()
