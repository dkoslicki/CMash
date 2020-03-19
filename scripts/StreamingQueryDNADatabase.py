#! /usr/bin/env python
import khmer
import numpy as np
import os
import sys
import multiprocessing
import pandas as pd
import argparse
from argparse import ArgumentTypeError
import re
import matplotlib.pyplot as plt
import timeit
from itertools import islice

# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
	from Query import Create
	from Query import Counters
	from Query import Containment
except ImportError:
	try:
		import MinHash as MH
		import Create
		import Counters
		import Containment
	except ImportError:
		try:
			sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
			from CMash import MinHash as MH
			from CMash.Query import Create  # fix relative imports
			from CMash.Query import Counters
			from CMash.Query import Containment
		except ImportError:
			raise Exception("Unable to import necessary classes")


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
	parser = argparse.ArgumentParser(
		description="This script calculates containment indicies for each of the training/reference sketches"
					" by streaming through the query file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use",
						default=multiprocessing.cpu_count())
	parser.add_argument('-c', '--containment_threshold', type=float,
						help="Only return results with containment index above this "
							 "threshold at the maximum k-mer size.", default=0.1)
	parser.add_argument('-p', '--plot_file', action="store_true", help="Optional flag to specify that a plot of the "
																	   "k-mer curves should also be saved (same basename"
																	   "as the out_file).")
	parser.add_argument('-r', '--reads_per_core', type=int,
						help="Number of reads per core in each chunk of parallelization."
							 " Set as high as memory will allow (eg. 1M on 256GB, 48 core machine)", default=100000)
	parser.add_argument('-f', '--filter_file',
						help="Location of pre-filter bloom filter. Use only if you absolutely know what you're doing "
							 "(hard to error check bloom filters).")
	parser.add_argument('-l', '--location_of_thresh', type=int,
						help="Location in range to apply the threshold passed by the -c flag. -l 2 -c 5-50-10 means the"
							 " threshold will be applied at k-size 25. Default is largest size.", default=-1)
	parser.add_argument('--sensitive', action="store_true",
						help="Operate in sensitive mode. Marginally more true positives with significantly more false positives. Use with caution.",
						default=False)
	parser.add_argument('-v', '--verbose', action="store_true", help="Print out progress report/timing information")
	parser.add_argument('in_file', help="Input file: FASTA/Q file to be processes")
	parser.add_argument('reference_file',
						help='Training database/reference file (in HDF5 format). Created with MakeStreamingDNADatabase.py')
	parser.add_argument('out_file', help='Output csv file with the containment indices.')
	parser.add_argument('range', type=parseNumList,
						help="Range of k-mer sizes in the formate <start>-<end>-<increment>."
							 " So 5-10-2 means [5, 7, 9]. If <end> is larger than the k-mer size"
							 "of the training data, this will automatically be reduced.")

	# read in the arguments
	args = parser.parse_args()
	k_range = args.range
	if k_range is None:
		raise Exception("The --range argument is required, no matter what the help menu says.")
	training_database_file = args.reference_file
	query_file = args.in_file
	results_file = args.out_file
	npz_file = os.path.splitext(results_file)[0] + "_hit_matrix.npz"
	num_threads = args.threads
	location_of_thresh = args.location_of_thresh
	coverage_threshold = args.containment_threshold
	TST_file = os.path.splitext(training_database_file)[0] + ".tst"  # name of the tst training file
	TST_file = os.path.abspath(TST_file)
	hydra_file = args.filter_file
	verbose = args.verbose
	num_reads_per_core = args.reads_per_core
	sensitive = args.sensitive
	if not os.path.exists(TST_file):
		TST_file = None
	if args.plot_file:
		plot_file = os.path.abspath(os.path.splitext(results_file)[0] + ".png")

	# Import data and error checking
	# Query file
	if not os.path.exists(query_file):
		raise Exception("Query file %s does not exist." % query_file)
	if not os.path.exists(training_database_file):
		raise Exception("Training/reference file %s does not exist." % training_database_file)

	# Training data
	if verbose:
		print("Reading in sketches")
		t0 = timeit.default_timer()
	sketches = MH.import_multiple_from_single_hdf5(training_database_file)
	if sketches[0]._kmers is None:
		raise Exception(
			"For some reason, the k-mers were not saved when the database was created. Try running MakeStreamingDNADatabase.py again.")
	num_hashes = len(sketches[0]._kmers)  # note: this is relying on the fact that the sketches were properly constructed
	max_ksize = sketches[0].ksize

	def keyfunction(item):
		return os.path.basename(item.input_file_name)

	sketches = sorted(sketches, key=keyfunction)  # sort the sketches by the basename of input file

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
		training_file_names.append(str(sketches[i].input_file_name.decode('utf-8')))

	training_file_names = sorted(training_file_names, key=os.path.basename)  # sort based on base name

	if verbose:
		print("Finished reading in sketches")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))

	if verbose:
		print("Reading in/creating ternary search tree")
		t0 = timeit.default_timer()

	# Import the query class to handle the creation of all the necessary data structures
	C = Create(training_database_file=training_database_file, bloom_filter_file=hydra_file, TST_file=TST_file, k_range=k_range)

	# Make the Marissa tree
	C.import_TST()

	# create the Bloom Filter prefilter
	C.create_BF_prefilter()

	if verbose:
		print("Finished reading in/creating ternary search tree")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))

	# Initialize the counters
	counter = Counters(tree=C.tree, k_range=C.k_range, seen_kmers=C.seen_kmers, all_kmers_bf=C.all_kmers_bf)

	def map_func(sequence):
		return counter.process_seq(sequence)

	pool = multiprocessing.Pool(processes=num_threads)

	if verbose:
		print("Start streaming")
		t0 = timeit.default_timer()
	# Open the file to prepare for processing
	fid = khmer.ReadParser(query_file)  # This is faster than screed
	match_tuples = []
	num_reads_per_chunk = num_reads_per_core * num_threads
	to_proc = [record.sequence for record in islice(fid, num_reads_per_chunk)]
	i = 0
	# process in chunks for less I/O time
	while to_proc:
			i += len(to_proc)
			if verbose:
				print("Read in %d sequences" % i)
			res = pool.map(map_func, to_proc, chunksize=int(max(1, min(num_reads_per_core, len(to_proc)/num_threads))))
			flattened_res = [item for sublist in res if sublist for item in sublist]
			flattened_res = list(set(flattened_res))  # dedup it
			match_tuples.extend(flattened_res)
			to_proc = [record.sequence for record in islice(fid, num_reads_per_chunk)]
	fid.close()
	pool.close()

	if verbose:
		print("Finished streaming")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))


	# initialize the containment class which will handle the conversion of the match tuples into hit matrices
	# then into a containment indicies, and finally create a pandas DF for export
	containment = Containment(k_range=k_range, match_tuples=match_tuples, sketches=sketches, num_hashes=num_hashes)

	if verbose:
		print("Forming hit matrix")
		t0 = timeit.default_timer()

	# create the hit matrices from the match_tuples
	containment.create_to_hit_matrices()

	if verbose:
		print("Finished forming hit matrix")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))

	if verbose:
		print("Computing containment indicies")
		t0 = timeit.default_timer()

	# Create the containment indicies
	containment.create_containment_indicies()

	if verbose:
		print("Finished computing containment indicies")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))

	# prepare a nice pandas data-frame for export
	containment.create_data_frame(training_file_names=training_file_names, location_of_thresh=location_of_thresh, coverage_threshold=coverage_threshold)


	if sensitive:
		if verbose:
			print("Exporting results")
			t0 = timeit.default_timer()
		# export results
		containment.filtered_results.to_csv(results_file, index=True, encoding='utf-8')

	# If requested, plot the results
	# FIXME: plot
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


	# FIXME: post-process class
	###############################################################################
	# Start the post-processing
	if verbose and not sensitive:
		print("Starting the post-processing")
		t0 = timeit.default_timer()

	if not sensitive:
		# Do the post-processing
		# Main idea here is to only concentrate on the unique k-mers: those that don't show up in more than one genome
		# as they are more specific to the presence of that genome being present in the sample

		# FIXME: prepare_post_process
		# from the non-filtered out genomes, create a dictionary that maps k-mer size to the properly reduced hit_matrix
		# first, get the basis of the reduced data frame
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

		# FIXME: find_kmers_in_filtered_results
		# get the count estimators of just the organisms of interest
		CEs = MH.import_multiple_from_single_hdf5(training_database_file, import_list=to_select_names)  # TODO: could make it a tad more memory efficient by sub-selecting the 'sketches'

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

		# FIXME: find_unique_kmers
		# Use this to identify which k-mers are unique (i.e. show up in exactly one sketch)
		for kmer in all_kmers_with_counts.keys():
			if all_kmers_with_counts[kmer] == 1:
				k_size = len(kmer)
				is_unique_kmer_per_ksize[k_size].add(kmer)
				is_unique_kmer.add(kmer)

		# FIXME: find_non_unique_kmers
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
					hit_matrices_dict['k=%d' % k_size][i, to_zero_indicies] = 0  # set these to zero since they show up in other sketches (so not informative)
				num_unique[i, k_range.index(k_size)] = len(current_kmers_set) - len(non_unique)  # keep track of the size of the unique k-mers

		# FIXME: create_post_containment_indicies
		# sum the modified hit matrices to get the size of the intersection
		containment_indices = np.zeros((len(to_select_names), len(k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			containment_indices[:, k_size_loc] = (
				hit_matrices_dict['k=%d' % k_size].sum(axis=1).ravel())  # /float(num_hashes))

		# then normalize by the number of unique k-mers (to get the containment index)
		# In essence, this is the containment index, restricted to unique k-mers. This effectively increases the specificity,
		# but also increases the variance/confidence interval, since this decreases the size of the sketch.
		for k_size_loc in range(len(k_range)):
			k_size = k_range[k_size_loc]
			for hash_loc in np.where(containment_indices[:, k_size_loc])[0]:  # find the genomes with non-zero containment
				unique_kmers = set()
				for kmer in CEs[hash_loc]._kmers:
					unique_kmers.add(kmer[:k_size])  # find the unique k-mers
				containment_indices[hash_loc, k_size_loc] /= float(len(unique_kmers))  # FIXME: this doesn't seem like the right way to normalize, but apparently it is!
				#containment_indices[hash_loc, k_size_loc] /= float(num_unique[hash_loc, k_size_loc])  # FIXME: in small tests, this seems to give better results. To be revisted.

		# FIXME: create_data_frame
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
		filtered_results = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)  # only select those where the highest k-mer size's count is above the threshold

		filtered_results.to_csv(results_file, index=True, encoding='utf-8')
		if verbose:
			t1 = timeit.default_timer()
			print("Finished thresholding. Time: %f" % (t1 - t0))

