import khmer
import marisa_trie as mt
import numpy as np
import os
import sys
import pandas as pd
from hydra import WritingBloomFilter, ReadingBloomFilter
from scipy.sparse import csc_matrix
from h5py import File
import subprocess
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
except ImportError:
	try:
		import MinHash as MH
	except ImportError:
		sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
		from CMash import MinHash as MH


def return_data_frame(training_file_names: list, k_range: list, location_of_thresh: int, containment_indices: np.ndarray, coverage_threshold: float) -> pd.DataFrame:
	"""
	Creates a nicely formatted Pandas data frame from the self.containment_indicies.
	:param training_file_names: the file names that were used to create the training database
	:type training_file_names: list
	:param k_range: a range of k-mer sizes
	:type k_range: list
	:param location_of_thresh: where in self.k_range the thresholding should take place (-1 means the last one)
	:type location_of_thresh: int
	:param containment_indices: the containment indicies matrix you wish to convert to a pandas data frame
	:type containment_indices: numpy.ndarray
	:param coverage_threshold: filter out those results that have containment indicies strictly below this threshold
	:type coverage_threshold: float
	"""
	results = dict()
	for k_size_loc in range(len(k_range)):
		ksize = k_range[k_size_loc]
		key = 'k=%d' % ksize
		results[key] = containment_indices[:, k_size_loc]
	df = pd.DataFrame(results, map(os.path.basename, training_file_names))
	df = df.reindex(labels=['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
	sort_key = 'k=%d' % k_range[location_of_thresh]
	max_key = 'k=%d' % k_range[-1]
	# only select those where the highest k-mer size's count is above the threshold
	filtered_results = df[df[sort_key] >= coverage_threshold].sort_values(max_key, ascending=False)
	return filtered_results


class Create:
	"""
	This class has functionality to:
	1. Import the ternary search tree created in the training step
	2. Create or import the bloom filter pre-filter file
	"""
	def __init__(self, training_database_file: str, bloom_filter_file: str, TST_file: str, k_range: list):
		"""
		Initializes the class
		:param training_database_file: file pointing to the HDF5 training database created with MakeStreamingDNADatabase.py
		:param bloom_filter_file: (optional) file pointing to file created with MakeStreamingPrePfilter.py. If empty string, one will be created
		:param TST_file: file pointing to the TST file (ternary search tree) that was created with MakeStreamingDNADatabase.py
		:param k_range: range of k-mer sizes. eg [10, 20, 30]
		"""
		self.bloom_filter_file = bloom_filter_file
		self.TST_file = TST_file
		self.k_range = k_range
		self.training_database = training_database_file
		self.tree = None  # populated by import_TST
		self.all_kmers_bf = None  # populated by create_BF_prefilter

	def import_TST(self) -> None:
		"""
		Imports the ternary search tree
		"""
		# no more safety net for those that didn't create a TST properly with the CreateStreamingQueryDNADatabase.py
		self.tree = mt.Trie()
		self.tree.load(self.TST_file)

	def create_BF_prefilter(self, result_file=None) -> None:
		"""
		Imports or creates the pre-filter Bloom filter
		:param result_file: (optional) if you'd like to export the bloom filter, populate that here
		:type result_file: str
		"""
		tree = self.tree
		k_range = self.k_range
		if not self.bloom_filter_file:  # create one
			try:
				# Get all the k-mers in the TST, put them in a bloom filter
				# all_kmers_bf = WritingBloomFilter(len(sketches) * len(k_range) * num_hashes * 20, 0.01)
				if result_file:
					# save it to the file
					self.all_kmers_bf = WritingBloomFilter(len(tree.keys()) * len(k_range) * 5, 0.01, ignore_case=True, filename=result_file)  # fudge factor of 5 will make the BF larger, but also slightly faster
				else:
					# keep it in memory
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


class Intersect:
	"""
	This class computes the intersection of the kmers in the input sample
	and the kmers in the training sample prior to testing each training sample
	for membership in the input sample. This eliminates the need to pre-screen kmers
	for presence in the TST via the prefilter since every kmer in the input must be
	in both the input and the training samples.
	"""
	def __init__(self, reads_path, training_path, input_type='fasta', threads=8, temp_dir=None, verbose=False):
		# handle file paths
		#reads
		self.reads_path = os.path.abspath(reads_path)
		if not os.path.exists(self.reads_path):
			raise Exception(f"It appears that the input reads {reads_path} does not exist.")

		# training database
		if not os.path.exists(training_path):
			raise Exception(f"It appears that the training database {training_path} does not exist.")
		# get the full path to the training database
		self.cmashDatabase = os.path.abspath(training_path)
		# get just the name of the training database
		self.cmashBaseName = os.path.splitext(os.path.abspath(training_path))[0]
		if self.cmashBaseName == "" or self.cmashBaseName == None:
			raise Exception("Invalid training database file name/path")
		# these are names we are setting ourselves. Might want to consider picking a name that won't cause conflicts
		# if running in parallel
		self.cmashDump = self.cmashBaseName + "_dump.fa"
		self.db_kmers_loc = self.cmashBaseName + '_dump'

		# if no temp directory was provided, place it in the same place the database was found
		if not temp_dir:
			self.temp_dir = os.path.dirname(self.cmashDatabase)
		else:
			# otherwise, just use what was provided
			self.temp_dir = temp_dir

		# handle the kmc commands
		self.kmc_dir = ""
		self.kmc = self.kmc_dir + "kmc"
		self.kmc_tools = self.kmc_dir + "kmc_tools"
		self.kmc_dump = self.kmc_dir + "kmc_dump"
		# test if KMC is installed correctly
		res = subprocess.run("kmc -h", shell=True, stdout=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception("It appears KMC is not installed. Please try '$ conda install -c bioconda KMC' and try again.")

		# set remaining properties
		self.threads = threads
		self.ksize = self.get_kmer_size()
		self.verbose = verbose
		self.input_type = input_type
		if input_type != "fasta" and input_type != "fastq":
			raise Exception("only allowable input_types are fasta and fastq")

		# read counts file
		self.reads_kmc_out_file = os.path.join(self.temp_dir, 'reads_' + str(self.ksize) + 'mers_dump')
		# intersection dump file
		self.intersection_kmc_out_file = os.path.join(self.temp_dir, str(self.ksize) + 'mers_intersection')
		self.intersection_kmc_dump_file = os.path.join(self.temp_dir, str(self.ksize) + 'mers_intersection_dump')

	def get_kmer_size(self):
		"""Reads the training database (in HDF5 format)
		and returns the uniform size of the k-mers.
		The ksize attribute is stored in the first
		subgroup of the CountEstimators group.
		"""
		db = File(self.cmashDatabase, "r")
		ce = db["CountEstimators"]
		ksize = None
		for it in ce.values():
			ksize = it.attrs["ksize"]
			break
		if not ksize:
			raise Exception(f"Failed to automatically detect k-mer size in {self.cmashDatabase}")
		return ksize


	def dump_training_kmers(self):
		"""
		This code is based on dump_kmers.py in the Metalign repository.
		Dumps the kmers from the input's CountEstimator to a fasta file.
		"""
		if self.verbose:
			print("dumping training k-mers")

		# delete the dump if it already exists
		subprocess.run(f"rm {self.cmashDump}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		training_database = self.cmashDatabase  # first input is the training file name
		dump_file = self.cmashDump  # second input is the desired output dump file
		CEs = MH.import_multiple_from_single_hdf5(training_database)
		with open(dump_file, 'w') as fid:
			i = 0
			for CE in CEs:
				for kmer in CE._kmers:
					fid.write('>seq%d\n' % i)
					fid.write('%s\n' % kmer)
					i += 1

	def count_training_kmers(self):
		"""
		Opens a KMC process to count the dumped training kmers.
		-ci1 excludes all kmers which appear less than one time (excludes no kmers).
		"""
		if self.verbose:
			print("counting training k-mers")

		out_path = self.db_kmers_loc
		# remove the old dumped database if it exists
		subprocess.run(f"rm {out_path}.kmc_pre", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		subprocess.run(f"rm {out_path}.kmc_suf", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		# count the training k-mers
		if self.verbose:
			res = subprocess.run(f"{self.kmc} -v -k{self.ksize} -fm -ci0 -cs3 -t{self.threads} {self.cmashDump} {out_path} .", shell=True)
		else:
			res = subprocess.run(f"{self.kmc} -v -k{self.ksize} -fm -ci0 -cs3 -t{self.threads} {self.cmashDump} {out_path} .", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")

	def count_input_kmers(self):
		"""
		Opens a KMC process to count the dumped kmers from the input sample.
		-ci1 excludes all kmers which appear less than one time (excludes no kmers).
		"""
		if self.verbose:
			print("counting input k-mers")

		if self.input_type == 'fastq':
				type_arg = '-fq'
		elif self.input_type == 'fasta':
				type_arg = '-fa'
		else:
			raise Exception("Only allowable input_types are 'fasta' and 'fastq'.")

		out_path = self.reads_kmc_out_file

		#count kmers in the input sample (not the training db)
		# TODO: will need to try with/without -fm if type_arg == '-fa'
		if self.verbose:
			res = subprocess.run(f"{self.kmc} -v -k{self.ksize} {type_arg} -ci0 -cs3 -t{self.threads} -fm {self.reads_path} {out_path} {self.temp_dir}", shell=True)
		else:
			res = subprocess.run(f"{self.kmc} -v -k{self.ksize} {type_arg} -ci0 -cs3 -t{self.threads} -fm {self.reads_path} {out_path} {self.temp_dir}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")


	def intersect(self):
		"""
		This source is from the Metalign repository.
		Opens a new KMC process to intersect the kmers
		from the input sample and the training database.
		The kmers in the intersection are written to a file as-is. Then,
		the contents of the file are read and rewritten in FASTA format.
		"""

		db_kmers_loc = self.db_kmers_loc
		in_path = self.reads_kmc_out_file
		out_path = self.intersection_kmc_out_file
		dump_path = self.intersection_kmc_dump_file

		if self.verbose:
			print("Intersecting input sample & training sample...")

		# intersect kmers
		if self.verbose:
			res = subprocess.run(f"{self.kmc_tools} simple {db_kmers_loc} -ci0 {in_path} -ci0 intersect {out_path} -ci0 ", shell=True)
		else:
			res = subprocess.run(f"{self.kmc_tools} simple {db_kmers_loc} -ci0 {in_path} -ci0 intersect {out_path} -ci0 ", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")

		#dump intersection
		if self.verbose:
			print("dumping intersection to FASTA file")
		if self.verbose:
			res = subprocess.run(f"{self.kmc_dump} -ci0 {out_path} {dump_path}", shell=True)
		else:
			res = subprocess.run(f"{self.kmc_dump} -ci0 {out_path} {dump_path}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		if res.returncode != 0:
			raise Exception(f"The command {res.args} failed to run and returned the returncode={res.returncode}")

		#read intersection & rewrite in fasta format
		with(open(dump_path, 'r')) as infile:
				with(open(dump_path + '.fa', 'w')) as fasta:
						for line in infile:
								seq = line.split()[0]
								fasta.write('>seq' + '\n' + seq + '\n')
		self.out_file = dump_path + '.fa'

	def compute_intersection(self):
		self.dump_training_kmers()
		self.count_training_kmers()
		self.count_input_kmers()
		self.intersect()


# shared object that will update the intersection counts
class Counters:
	"""
	This class is mainly to facilitate the parallelization of enumerating input query k-mer hits to the ternary
	search tree
	"""
	def __init__(self, tree: mt.Trie, k_range: list, all_kmers_bf: WritingBloomFilter):
		"""
		Initialize the class
		:param tree: The actual ternary search tree (TST) that k-mers and their pre-fixes will be searched against
		:type tree: marisa_trie.Tree
		:param k_range: a range of k-mer sizes
		:type k_range: list
		:param all_kmers_bf: a pre-filter that checks in O(1) time if an input k-mer may be in the TST
		:type all_kmers_bf: hydra.WritingBloomFilter
		"""
		self.tree = tree
		self.k_range =k_range
		self.seen_kmers = set()
		self.all_kmers_bf = all_kmers_bf

	# This class is basically an array of counters (on the same basis as the sketches)
	# it's used to keep track (in a parallel friendly way) of which streamed k-mers went into the training file sketches
	#def __init__(self, training_database_file=None, bloom_filter_file=None, TST_file=None, k_range=None):
	#	super().__init__(training_database_file, bloom_filter_file, TST_file, k_range)

	def return_matches(self, input_kmer: str, k_size_loc: int) -> tuple:
		"""
		Get all the matches in the TST with the kmer prefix
		:param input_kmer: an input k-mer
		:type input_kmer: str
		:param k_size_loc: where in self.k_range this k-mer (via it's length) belongs
		:type k_size_loc: inst
		:return: a tuple: first of which is a list of strings (all the matches in the TST), and the second is a Boolean indicating if you saw a match
		:rtype: tuple
		"""
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
		"""
		Takes an input sequence, breaks it into its k-mers (for every size self.k_range), and after some filtering and
		checking, sends it to return_matches to query the TST
		:param seq: an input DNA sequence
		:type seq: string
		:return: a list of keys indicating all the TST hits for all the k-mers in seq
		:rtype: list
		"""
		k_range = self.k_range
		seen_kmers = self.seen_kmers
		all_kmers_bf = self.all_kmers_bf
		#  start with small kmer size, if see match, then continue looking for longer k-mer sizes, otherwise move on
		small_k_size = k_range[0]  # start with the small k-size
		to_return = []
		seq = seq.upper()
		# TODO: could, for efficiency, also remove non-ACTG, but those won't match anyways since they aren't in the TST
		#  might not actually be more efficient to search for non-ACTG too
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
				possible_match = True  # FIXME: bug introduced here in cf64b7aace5eadf738b920109d6419c9d930a1dc, make sure it didn't happen again

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
	"""
	This class handles all the conversion from raw TST hits, to containment index computations
	"""
	# TODO: would like to indicate that sketches should be a list of CEs from MH
	def __init__(self, k_range: list, match_tuples: list, sketches: list, num_hashes: int):
		"""
		Initialize the class
		:param k_range: a range of k-mer sizes
		:type k_range: list
		:param match_tuples: all the tuples returned by Query.return_matches after processing each sequence
		:type match_tuples:  list
		:param sketches: the MinHash sketches
		:type sketches: list[CMash.CountEstimator]
		:param num_hashes: number of hashes that each count estimator has
		:type num_hashes: int
		"""
		self.k_range = k_range
		self.match_tuples = match_tuples
		self.sketches = sketches
		self.num_hashes = num_hashes
		self.hit_matrices = []
		# TODO: could make this thing sparse, or do the filtering for above threshold here
		# as it's probably a memory hog
		self.containment_indices = np.zeros((len(sketches), len(k_range)))
		self.filtered_results = None  # will be filled in when I do the exporting

	def create_to_hit_matrices(self) -> None:
		"""
		Converts the match tuples into a list of matrices, one for each k-mer size in self.k_range
		"""
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
		for k_size in k_range:
			# convert to matrices
			mat = csc_matrix((value_dict[k_size], (row_ind_dict[k_size], col_ind_dict[k_size])),
							 shape=(len(sketches), num_hashes))
			self.hit_matrices.append(mat)

	def create_containment_indicies(self) -> None:
		"""
		Utilizes the self.hit_matrices to compute the actual containment indicies in self.containment_indicies
		"""
		sketches = self.sketches
		k_range = self.k_range
		hit_matrices = self.hit_matrices
		# prep the containment indicies matrix: rows are genome/sketch, one column for each k-mer size in k_range
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
					if kmer[:k_size]:
						unique_kmers.add(kmer[:k_size])
				self.containment_indices[hash_loc, k_size_loc] /= float(
					len(unique_kmers))  # divide by the unique num of k-mers

	def create_data_frame(self, training_file_names: list, location_of_thresh: int, coverage_threshold: int) -> None:
		"""
		Creates a nicely formatted Pandas data frame from the self.containment_indicies.
		"""
		self.filtered_results = return_data_frame(training_file_names, self.k_range, location_of_thresh, self.containment_indices, coverage_threshold)


class PostProcess:
	"""
	A class to perform the post-processing for more specific (less sensitive) results.
	Main idea here is to only concentrate on the unique k-mers: those that don't show up in more than one genome
	as they are more specific to the presence of that genome being present in the sample
	"""
	def __init__(self, filtered_results: pd.DataFrame, training_file_names: list, k_range: list, hit_matrices: list):
		"""
		Initialize the class
		:param filtered_results: the Containment.filtered_results pandas data frame
		:type filtered_results: pandas.DataFrame
		:param training_file_names: the file names that were used to create the training database
		:type training_file_names: list
		:param k_range: a range of k-mer sizes
		:type k_range: list
		:param hit_matrices: the list of hit matrices from Containment.hit_matrices
		:type hit_matrices: list
		"""
		self.filtered_results = filtered_results
		self.training_file_names = training_file_names
		self.k_range = k_range
		self.hit_matrices = hit_matrices
		self.to_select_names = None
		self.all_kmers_with_counts = None
		self.is_unique_kmer_per_ksize = dict()
		self.CEs = None
		self.hit_matrices_dict = None
		self.containment_indices = None
		self.num_unique_dict = dict()

	def prepare_post_process(self) -> None:
		"""
		Converts the filtered results to a dictionary mapping k-mer size to dense matrices. Stores in self.hit_matrices_dict
		"""
		filtered_results = self.filtered_results
		training_file_names = self.training_file_names
		k_range = self.k_range
		hit_matrices = self.hit_matrices
		self.to_select_names = list(filtered_results.index)
		all_names = list(map(os.path.basename, training_file_names))
		rows_to_select = []
		for name in self.to_select_names:
			rows_to_select.append(all_names.index(name))
		hit_matrices_dict_temp = dict()

		# TODO: could make this much cleaner by not using the intermediate hit_matrices_dict_temp and just converting
		# TODO: to dense immediately then reduce the hit matrix to this basis
		for i in range(len(k_range)):
			k_size = k_range[i]
			hit_matrices_dict_temp['k=%d' % k_size] = hit_matrices[i][rows_to_select, :]

		# Make the hit matrices dense
		hit_matrices_dense_dict = dict()
		for k_size in k_range:
			hit_matrices_dense_dict['k=%d' % k_size] = hit_matrices_dict_temp['k=%d' % k_size].todense()

		self.hit_matrices_dict = hit_matrices_dense_dict

	def find_kmers_in_filtered_results(self, training_database_file: str) -> None:
		"""
		For each of the genomes that showed up in self.filtered_results, collects all their k-mers and counts
		and puts it in self.all_kmers_with_counts.
		:param training_database_file: file pointing to the HDF5 training database created with MakeStreamingDNADatabase.py
		:type training_database_file: string
		"""
		to_select_names = self.to_select_names
		k_range = self.k_range
		#is_unique_kmer_per_ksize = self.is_unique_kmer_per_ksize

		# get the count estimators of just the organisms of interest
		# TODO: could make it a LOT more memory efficient by sub-selecting the 'sketches'
		self.CEs = MH.import_multiple_from_single_hdf5(training_database_file, import_list=to_select_names)

		# get all the kmers (for each kmer size) and form their counts in the subset of predicted sketches to be in the sample
		self.all_kmers_with_counts = dict()
		for k_size in k_range:
			#self.is_unique_kmer_per_ksize[k_size] = set()
			for i in range(len(self.CEs)):
				for big_kmer in self.CEs[i]._kmers:
					kmer = big_kmer[:k_size]
					if kmer in self.all_kmers_with_counts:
						self.all_kmers_with_counts[kmer] += 1
					else:
						self.all_kmers_with_counts[kmer] = 1

	def find_unique_kmers(self) -> None:
		"""
		Finds the k-mers that showed up in exactly one sketch (i.e. uniquely identified the related genome)
		"""
		# Use this to identify which k-mers are unique (i.e. show up in exactly one sketch)
		k_range = self.k_range
		is_unique_kmer_per_ksize = self.is_unique_kmer_per_ksize
		for k_size in k_range:
			self.is_unique_kmer_per_ksize[k_size] = set()
		is_unique_kmer = set()  # TODO: might not actually need this
		all_kmers_with_counts = self.all_kmers_with_counts
		for kmer in all_kmers_with_counts.keys():
			if all_kmers_with_counts[kmer] == 1:
				k_size = len(kmer)
				is_unique_kmer_per_ksize[k_size].add(kmer)
				is_unique_kmer.add(kmer)

	def find_non_unique_kmers_reduce_hit_matrices(self) -> None:
		"""
		Finds the k-mers that showed up in more than one sketch/genome, use these to reduce the hit matrices to set
		the corresponding entries to zero.
		"""
		# Also keep track of which kmers appear in more than one sketch (not unique)
		CEs = self.CEs
		k_range = self.k_range
		is_unique_kmer_per_ksize = self.is_unique_kmer_per_ksize
		hit_matrices_dict = self.hit_matrices_dict
		num_unique_dict = self.num_unique_dict
		for i in range(len(CEs)):
			for k_size in k_range:
				current_kmers = [k[:k_size] for k in CEs[i]._kmers]
				current_kmers_set = set(current_kmers)
				non_unique_set = set()

				for kmer in current_kmers:
					if kmer not in is_unique_kmer_per_ksize[k_size]:
						non_unique_set.add(kmer)
				# reduce the hit matrices by removing the hits corresponding to non-unique k-mers
				to_zero_indicies = [ind for ind, kmer in enumerate(current_kmers) if kmer in non_unique_set]
				# if you use a really small initial kmer size, some of the hit matrices may be empty
				# due to all k-mers being shared in common
				# set these to zero since they show up in other sketches (so not informative)
				if hit_matrices_dict['k=%d' % k_size].size > 0:
					hit_matrices_dict['k=%d' % k_size][i, to_zero_indicies] = 0
				# keep track of the size of the unique k-mers
				num_unique_dict[i, k_range.index(k_size)] = len(current_kmers_set) - len(non_unique_set)

	def create_post_containment_indicies(self) -> None:
		"""
		Goes through the hit_matrices_dict and sets to 0 those that showed up in more than one sketch
		(i.e. did not uniquely identify the associated genome)
		"""
		# sum the modified hit matrices to get the size of the intersection
		to_select_names = self.to_select_names
		k_range = self.k_range
		hit_matrices_dict = self.hit_matrices_dict
		CEs = self.CEs

		# TODO: could make this thing sparse, or do the filtering for above threshold here
		self.containment_indices = np.zeros((len(to_select_names), len(k_range)))
		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			self.containment_indices[:, k_size_loc] = (
				hit_matrices_dict['k=%d' % k_size].sum(axis=1).ravel())  # /float(num_hashes))

		# then normalize by the number of unique k-mers (to get the containment index)
		# In essence, this is the containment index, restricted to unique k-mers. This effectively increases the specificity,
		# but also increases the variance/confidence interval, since this decreases the size of the sketch.
		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			for hash_loc in np.where(self.containment_indices[:, k_size_loc])[
				0]:  # find the genomes with non-zero containment
				unique_kmers = set()
				for kmer in CEs[hash_loc]._kmers:
					unique_kmers.add(kmer[:k_size])  # find the unique k-mers
				# FIXME: this doesn't seem like the right way to normalize, but apparently it is!
				self.containment_indices[hash_loc, k_size_loc] /= float(len(unique_kmers))
				# FIXME: in small tests, this seems to give better results. To be revisted.
				#self.containment_indices[hash_loc, k_size_loc] /= float(self.num_unique_dict[hash_loc, k_size_loc])

	def create_data_frame(self, training_file_names: list, location_of_thresh: int, coverage_threshold: int) -> None:
		"""
		Creates a nicely formatted Pandas data frame from the self.containment_indicies.
		:param training_file_names: the file names that were used to create the training database
		:type training_file_names: list
		:param location_of_thresh: here in self.k_range the thresholding should take place (-1 means the last one)
		:type location_of_thresh: int
		:param coverage_threshold: filter out those results that have containment indicies below this threshold
		:type coverage_threshold: float
		"""
		self.filtered_results = return_data_frame(training_file_names, self.k_range, location_of_thresh,
												  self.containment_indices, coverage_threshold)

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

	#counters = Counters(tree=C.tree, k_range=C.k_range, seen_kmers=C.seen_kmers, all_kmers_bf=C.all_kmers_bf)
	#print(f"Processed sequence with counter: {counters.process_seq('AGTCCGCGCCACTGGCAGTGACCATCGACACGCAGACGGAGATTAACAACATTGTACTGGTCAATGATACCGGTATGCCG')}")





# simple way to do testing
if __name__ == "__main__":
	main()
