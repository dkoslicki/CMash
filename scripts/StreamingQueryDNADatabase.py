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

# TODO: export hit matrices


def parseNumList(input):
	"""Thank you stack overflow"""
	m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', input)
	# ^ (or use .split('-'). anyway you like.)
	if not m:
		raise ArgumentTypeError("'" + input + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end+1, increment))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script calculates containment indicies for each of the training/reference sketches"
									" by streaming through the query file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
	parser.add_argument('-c', '--containment_threshold', type=float, help="Only return results with containment index above this "
															  "threshold at the maximum k-mer size.", default=0.1)
	parser.add_argument('-p', '--plot_file', action="store_true", help="Optional flag to specify that a plot of the "
																	   "k-mer curves should also be saved (same basename"
																	   "as the out_file).")
	parser.add_argument('-r', '--reads_per_core', type=int, help="Number of reads per core in each chunk of parallelization."
																 " Set as high as memory will allow (eg. 1M on 256GB, 48 core machine)", default=100000)
	parser.add_argument('-f', '--filter_file',
						help="Location of pre-filter bloom filter. Use only if you absolutely know what you're doing "
							"(hard to error check bloom filters).")
	parser.add_argument('-l', '--location_of_thresh', type=int,
						help="Location in range to apply the threshold passed by the -c flag. -l 2 -c 5-50-10 means the"
							" threshold will be applied at k-size 25. Default is largest size.", default=-1)
	parser.add_argument('-v', '--verbose', action="store_true", help="Print out progress report/timing information")
	parser.add_argument('in_file', help="Input file: FASTA/Q file to be processes")
	parser.add_argument('reference_file', help='Training database/reference file (in HDF5 format). Created with MakeStreamingDNADatabase.py')
	parser.add_argument('out_file', help='Output csv file with the containment indices.')
	parser.add_argument('range', type=parseNumList, help="Range of k-mer sizes in the formate <start>-<end>-<increment>."
														   " So 5-10-2 means [5, 7, 9]. If <end> is larger than the k-mer size"
														   "of the training data, this will automatically be reduced.")

	# read in the arguments
	args = parser.parse_args()
	k_range = args.range
	if k_range is None:
		raise Exception("The --range argument is required, no matter what the help menu says.")
	training_data = args.reference_file
	query_file = args.in_file
	results_file = args.out_file
	npz_file = os.path.splitext(results_file)[0] + "_hit_matrix.npz"
	num_threads = args.threads
	location_of_thresh = args.location_of_thresh
	coverage_threshold = args.containment_threshold
	streaming_database_file = os.path.splitext(training_data)[0] + ".tst"  # name of the tst training file
	streaming_database_file = os.path.abspath(streaming_database_file)
	hydra_file = args.filter_file
	verbose = args.verbose
	num_reads_per_core = args.reads_per_core
	if not os.path.exists(streaming_database_file):
		streaming_database_file = None
	if args.plot_file:
		plot_file = os.path.abspath(os.path.splitext(results_file)[0] + ".png")

	# Import data and error checking
	# Query file
	if not os.path.exists(query_file):
		raise Exception("Query file %s does not exist." % query_file)
	if not os.path.exists(training_data):
		raise Exception("Training/reference file %s does not exist." % training_data)

	# Training data
	if verbose:
		print("Reading in sketches")
		t0 = timeit.default_timer()
	sketches = MH.import_multiple_from_single_hdf5(training_data)
	if sketches[0]._kmers is None:
		raise Exception(
			"For some reason, the k-mers were not saved when the database was created. Try running MakeDNADatabase.py again.")
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
			all_kmers_bf = WritingBloomFilter(len(sketches)*len(k_range)*num_hashes, 0.01)
			for sketch in sketches:
				for kmer in sketch._kmers:
					for ksize in k_range:
						all_kmers_bf.add(kmer[0:ksize])  # put all the k-mers and the appropriate suffixes in
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

	# shared object that will update the intersection counts
	class Counters(object):
		# This class is basically an array of counters (on the same basis as the sketches)
		# it's used to keep track (in a parallel friendly way) of which streamed k-mers went into the training file sketches
		def __init__(self):
			pass

		def return_matches(self, kmer, k_size_loc):
			""" Get all the matches in the trie with the kmer prefix"""
			prefix_matches = tree.keys(kmer)  # get all the k-mers whose prefix matches
			match_info = set()
			# get the location of the found kmers in the counters
			for item in prefix_matches:
				split_string = item.split('x')  # first is the hash location, second is which k-mer
				hash_loc = int(split_string[1])
				kmer_loc = int(split_string[2])
				match_info.add((hash_loc, k_size_loc, kmer_loc))
			to_return = []
			saw_match = False
			if match_info:
				saw_match = True
				for tup in match_info:
					to_return.append(tup)
			return to_return, saw_match

		def process_seq(self, seq):
			#  start with small kmer size, if see match, then continue looking for longer k-mer sizes, otherwise move on
			small_k_size = k_range[0]  # start with the small k-size
			to_return = []
			for i in range(len(seq) - small_k_size + 1):  # look at all k-mers
				kmer = seq[i:i + small_k_size]
				possible_match = False
				if kmer not in seen_kmers:  # if we should process it
					if kmer in all_kmers_bf:  # if we should process it
						match_list, saw_match = self.return_matches(kmer, 0)
						if saw_match:  #  TODO: note, I *could* add all the trie matches and their sub-kmers to the seen_kmers
							seen_kmers.add(kmer)
							to_return.extend(match_list)
						possible_match = True
					# TODO: note: I could (since it'd only be for a single kmer size, keep a set of *all* small_kmers I've tried and use this as another pre-filter
				else:
					possible_match = True
				# start looking at the other k_sizes, don't overhang len(seq)
				if possible_match:
					for other_k_size in [x for x in k_range[1:] if i+x <= len(seq)]:
						kmer = seq[i:i + other_k_size]
						if kmer in all_kmers_bf:
							k_size_loc = k_range.index(other_k_size)
							match_list, saw_match = self.return_matches(kmer, k_size_loc)
							if saw_match:
								to_return.extend(match_list)
						else:
							break
			return to_return


	# Initialize the counters
	# TODO: note, I could be doing a partial dedup here, just to reduce the memory usage...
	counter = Counters()
	def map_func(sequence):
		return counter.process_seq(sequence)

	pool = multiprocessing.Pool(processes=num_threads)

	if verbose:
		print("Start streaming")
		t0 = timeit.default_timer()
	# populate the queue
	fid = khmer.ReadParser(query_file)  # This is faster than screed
	match_tuples = []
	#num_reads_per_core = 100000
	num_reads_per_chunk = num_reads_per_core * num_threads
	to_proc = [record.sequence for record in islice(fid, num_reads_per_chunk)]
	i = 0
	while to_proc:
			i += len(to_proc)
			if verbose:
				print("Read in %d sequences" % i)
			res = pool.map(map_func, to_proc, chunksize=max(1, min(num_reads_per_core, len(to_proc)/num_threads)))
			flattened_res = [item for sublist in res if sublist for item in sublist]
			flattened_res = list(set(flattened_res))  # dedup it
			match_tuples.extend(flattened_res)
			to_proc = [record.sequence for record in islice(fid, num_reads_per_chunk)]
	fid.close()
	#print(match_tuples)
	if verbose:
		print("Finished streaming")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))

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

	match_tuples = set(match_tuples)  # uniquify, so we don't make the row/col ind dicts too large

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

	if verbose:
		print("Exporting results")
		t0 = timeit.default_timer()
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
	filtered_results.to_csv(results_file, index=True, encoding='utf-8')

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
	# then export
	savemat(npz_file, hit_matrices_dict, appendmat=False, do_compression=True)


	# If requested, plot the results
	if args.plot_file:
		df = pd.read_csv(results_file)  # annoyingly, I have to read it back in to get the format into something I can work with
		dft = df.transpose()
		fig = plt.figure()
		for key in dft.keys():
			plt.plot(k_range, dft[key].values[1:])  # [1:] since the first entry is the name of the file
		plt.legend(dft.values[0])
		plt.xlabel('K-mer size')
		plt.ylabel('Containment Index')
		fig.savefig(plot_file)

	if verbose:
		print("Finished exporting results")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))
