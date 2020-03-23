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
	from Query import PostProcess
except ImportError:
	try:
		import MinHash as MH
		import Create
		import Counters
		import Containment
		import PostProcess
	except ImportError:
		try:
			sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
			from CMash import MinHash as MH
			from CMash.Query import Create  # fix relative imports
			from CMash.Query import Counters
			from CMash.Query import Containment
			from CMash.Query import PostProcess
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
	num_threads = args.threads
	location_of_thresh = args.location_of_thresh
	coverage_threshold = args.containment_threshold
	TST_file = os.path.splitext(training_database_file)[0] + ".tst"  # name of the tst training file
	TST_file = os.path.abspath(TST_file)
	bloom_filter_file = args.filter_file
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

	sketches = sorted(sketches, key=lambda x: os.path.basename(x.input_file_name))

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
	C = Create(training_database_file=training_database_file, bloom_filter_file=bloom_filter_file, TST_file=TST_file, k_range=k_range)

	# Make the Marissa tree
	C.import_TST()

	# create the Bloom Filter prefilter
	C.create_BF_prefilter()

	if verbose:
		print("Finished reading in/creating ternary search tree")
		t1 = timeit.default_timer()
		print("Time: %f" % (t1 - t0))

	# Initialize the counters
	counter = Counters(tree=C.tree, k_range=C.k_range, all_kmers_bf=C.all_kmers_bf)

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


	# Start the post-processing, if requested
	if verbose and not sensitive:
		print("Starting the post-processing")
		t0 = timeit.default_timer()

	if not sensitive:
		# Do the post-processing
		# Main idea here is to only concentrate on the unique k-mers: those that don't show up in more than one genome
		# as they are more specific to the presence of that genome being present in the sample
		post_process = PostProcess(filtered_results=containment.filtered_results,
								   training_file_names=training_file_names, k_range=k_range,
								   hit_matrices=containment.hit_matrices)

		# import the TST and import/create the BF pre-filter
		post_process.prepare_post_process()

		# find all the k-mers in the current results
		post_process.find_kmers_in_filtered_results(training_database_file=training_database_file)

		# find those k-mers that are unique to their genome/sketch
		post_process.find_unique_kmers()

		# find k-mers that are shared between more than one genome/sketch
		post_process.find_non_unique_kmers_reduce_hit_matrices()

		# reduce the hit_matrices/containment indicies by not considering those k-mers that are shared between more
		# than one genome/sketch
		post_process.create_post_containment_indicies()

		# create the data frame
		to_select_names = post_process.to_select_names
		post_process.create_data_frame(training_file_names=to_select_names, location_of_thresh=location_of_thresh, coverage_threshold=coverage_threshold)

		# and then export it
		post_process.filtered_results.to_csv(results_file, index=True, encoding='utf-8')
		if verbose:
			t1 = timeit.default_timer()
			print("Finished thresholding. Time: %f" % (t1 - t0))

