#! /usr/bin/env python
import khmer
import numpy as np
from khmer.khmer_args import optimal_size
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
import pandas as pd
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
import threading
from itertools import *
import argparse
import screed


# Helper function that uses equations (2.1) and (2.7) that tells you where you need
# to set the threshold to ensure (with confidence 1-t) that you got all organisms
# with coverage >= c
def threshold_calc(k, c, p, confidence):
	delta = c*(1-np.sqrt(-2*np.log(1 - confidence)*(c+p) / float(c**2 * k)))
	if delta < 0:
		delta = 0
	return delta


# This will calculate the similarity indicies between one sketch and the sample NodeGraph
def compute_indicies(sketch, num_sample_kmers, true_fprate):
	#global sample_kmers
	num_hashes = len(sketch._kmers)
	num_training_kmers = sketch._true_num_kmers
	count = 0
	adjust = 0
	for kmer in sketch._kmers:
		if kmer != '':
			count += sample_kmers.get(kmer)
		else:
			adjust += 1  # Hash wasn't full, adjust accordingly
	count -= int(np.round(true_fprate*num_hashes))  # adjust by the FP rate of the bloom filter
	intersection_cardinality = count
	containment_index = count / float(num_hashes - adjust)
	jaccard_index = num_training_kmers * containment_index / float(num_training_kmers + num_sample_kmers - num_training_kmers * containment_index)
	#print("Train %d sample %d" % (num_training_kmers, num_sample_kmers))
	# It can happen that the query file has less k-mers in it than the training file, so just round to nearest reasonable value
	if containment_index > 1:
		containment_index = 1
	elif containment_index < 0:
		containment_index = 0
	if jaccard_index > 1:
		jaccard_index = 1
	elif jaccard_index < 0:
		jaccard_index = 0
	return intersection_cardinality, containment_index, jaccard_index


def unwrap_compute_indicies(arg):
	return compute_indicies(*arg)


def restricted_float(x):
	x = float(x)
	if x < 0.0 or x >= 1.0:
		raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0)" % (x,))
	return x


# Can read in with: test=pd.read_csv(os.path.abspath('../data/results.csv'),index_col=0)
def main():
	parser = argparse.ArgumentParser(description="This script creates a CSV file of similarity indicies between the"
									" input file and each of the sketches in the training/reference file.",
									formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
	parser.add_argument('-f', '--force', action="store_true", help="Force creation of new NodeGraph.")
	parser.add_argument('-fp', '--fp_rate', type=restricted_float, help="False positive rate.", default=0.0001)
	parser.add_argument('-ct', '--containment_threshold', type=restricted_float,
						help="Only return results with containment index above this value", default=0.02)
	parser.add_argument('-c', '--confidence', type=restricted_float,
						help="Desired probability that all results were returned with containment index above threshold [-ct]",
						default=0.95)
	parser.add_argument('-ng', '--node_graph', help="NodeGraph/bloom filter location. Used if it exists; if not, one "
													"will be created and put in the same directory as the specified "
													"output CSV file.", default=None)
	parser.add_argument('-b', '--base_name', action="store_true",
						help="Flag to indicate that only the base names (not the full path) should be saved in the output CSV file")
	parser.add_argument('-i', '--intersect_nodegraph', action="store_true",
						help="Option to only insert query k-mers in bloom filter if they appear anywhere in the training"
							 " database. Note that the Jaccard estimates will now be "
							 "J(query intersect union_i training_i, training_i) instead of J(query, training_i), "
							 "but will use significantly less space.")
	parser.add_argument('in_file', help="Input file: FASTQ/A file (can be gzipped).")
	parser.add_argument('training_data', help="Training/reference data (HDF5 file created by MakeTrainingDatabase.py)")
	parser.add_argument('out_csv', help='Output CSV file')

	# Parse and check args
	args = parser.parse_args()
	base_name = args.base_name
	training_data = os.path.abspath(args.training_data)
	if not os.path.exists(training_data):
		raise Exception("Training/reference file %s does not exist." % training_data)
	# Let's get the k-mer sizes in the training database
	ksizes = set()
	# Import all the training data
	sketches = MH.import_multiple_from_single_hdf5(training_data)
	# Check for issues with the sketches (can also check if all the kmers make sense (i.e. no '' or non-ACTG characters))
	if sketches[0]._kmers is None:
		raise Exception("For some reason, the k-mers were not saved when the database was created. Try running MakeDNADatabase.py again.")
	num_hashes = len(sketches[0]._kmers)
	for i in range(len(sketches)):
		sketch = sketches[i]
		if sketch._kmers is None:
			raise Exception(
				"For some reason, the k-mers were not saved when the database was created. Try running MakeDNADatabase.py again.")
		if len(sketch._kmers) != num_hashes:
			raise Exception("Unequal number of hashes for sketch of %s" % sketch.input_file_name)
		ksizes.add(sketch.ksize)
		if len(ksizes) > 1:
			raise Exception("Training/reference data uses different k-mer sizes. Culprit was %s." % (sketch.input_file_name))
	# Get the appropriate k-mer size
	ksize = ksizes.pop()
	# Get number of threads to use
	num_threads = args.threads
	# Check and parse the query file
	query_file = os.path.abspath(args.in_file)
	if not os.path.exists(query_file):
		raise Exception("Query file %s does not exist." % query_file)
	# Node graph is stored in the output folder with name <InputFASTQ/A>.NodeGraph.K<k_size>
	if args.node_graph is None:  # If no node graph is specified, create one
		node_graph_out = os.path.join(os.path.dirname(os.path.abspath(args.out_csv)),
									os.path.basename(query_file) + ".NodeGraph.K" + str(ksize))
		if not os.path.exists(node_graph_out):  # Don't complain if the default location works
			print("Node graph not provided (via -ng). Creating one at: %s" % node_graph_out)
	elif os.path.exists(args.node_graph):  # If one is specified and it exists, use it
		node_graph_out = args.node_graph
	else:  # Otherwise, the specified one doesn't exist
		raise Exception("Provided NodeGraph %s does not exist." % args.node_graph)
	# import and check the intersect nodegraph
	if args.intersect_nodegraph is True:
		intersect_nodegraph_file = os.path.splitext(training_data)[0] + ".intersect.Nodegraph"
	else:
		intersect_nodegraph_file = None
	intersect_nodegraph = None
	if intersect_nodegraph_file is not None:
		if not os.path.exists(intersect_nodegraph_file):
			raise Exception("Intersection nodegraph does not exist. Please re-run MakeDNADatabase.py with the -i flag.")
		try:
			intersect_nodegraph = khmer.load_nodegraph(intersect_nodegraph_file)
			if intersect_nodegraph.ksize() != ksize:
				raise Exception("Given intersect nodegraph %s has K-mer size %d while the database K-mer size is %d"
								% (intersect_nodegraph_file, intersect_nodegraph.ksize(), ksize))
		except:
			raise Exception("Could not load given intersect nodegraph %s" % intersect_nodegraph_file)
	results_file = os.path.abspath(args.out_csv)
	force = args.force
	fprate = args.fp_rate
	coverage_threshold = args.containment_threshold  # desired coverage cutoff
	confidence = args.confidence  # desired confidence that you got all the organisms with coverage >= desired coverage

	# Get names of training files for use as rows in returned tabular data
	training_file_names = []
	for i in range(len(sketches)):
		training_file_names.append(sketches[i].input_file_name)

	# Only form the Nodegraph if we need to
	global sample_kmers
	if not os.path.exists(node_graph_out) or force is True:
		hll = khmer.HLLCounter(0.01, ksize)
		hll.consume_seqfile(query_file)
		full_kmer_count_estimate = hll.estimate_cardinality()
		res = optimal_size(full_kmer_count_estimate, fp_rate=fprate)
		if intersect_nodegraph is None:  # If no intersect list was given, just populate the bloom filter
			sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
			#sample_kmers.consume_seqfile(query_file)
			rparser = khmer.ReadParser(query_file)
			threads = []
			for _ in range(num_threads):
				cur_thrd = threading.Thread(target=sample_kmers.consume_seqfile_with_reads_parser, args=(rparser,))
				threads.append(cur_thrd)
				cur_thrd.start()
			for thread in threads:
				thread.join()
		else:  # Otherwise, only put a k-mer in the bloom filter if it's in the intersect list
			# (WARNING: this will cause the Jaccard index to be calculated in terms of J(query\intersect hash_list, training)
			#  instead of J(query, training)
			# (TODO: fix this after khmer is updated)
			#intersect_nodegraph_kmer_count = intersect_nodegraph.n_unique_kmers()  # Doesnt work due to khmer bug
			intersect_nodegraph_kmer_count = intersect_nodegraph.n_occupied()  # Not technically correct, but I need to wait until khmer is updated
			if intersect_nodegraph_kmer_count < full_kmer_count_estimate:  # At max, we have as many k-mers as in the union of the training database (But makes this always return 0)
				res = optimal_size(intersect_nodegraph_kmer_count, fp_rate=fprate)
				sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
			else:
				sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
			for record in screed.open(query_file):
				seq = record.sequence
				for i in range(len(seq) - ksize + 1):
					kmer = seq[i:i + ksize]
					if intersect_nodegraph.get(kmer) > 0:
						sample_kmers.add(kmer)
		# Save the sample_kmers
		sample_kmers.save(node_graph_out)
		true_fprate = khmer.calc_expected_collisions(sample_kmers, max_false_pos=0.99)
	else:
		sample_kmers = khmer.load_nodegraph(node_graph_out)
		node_ksize = sample_kmers.ksize()
		if node_ksize != ksize:
			raise Exception("Node graph %s has wrong k-mer size of %d (input was %d). Try --force or change -k." % (
			node_graph_out, node_ksize, ksize))
		true_fprate = khmer.calc_expected_collisions(sample_kmers, max_false_pos=0.99)

	#num_sample_kmers = sample_kmers.n_unique_kmers()  # For some reason this only works when creating a new node graph, use the following instead
	num_sample_kmers = sample_kmers.n_occupied()

	# Compute all the indicies for all the training data
	pool = Pool(processes=num_threads)
	res = pool.map(unwrap_compute_indicies, zip(sketches, repeat(num_sample_kmers), repeat(true_fprate)))

	# Gather up the results in a nice form
	intersection_cardinalities = np.zeros(len(sketches))
	containment_indexes = np.zeros(len(sketches))
	jaccard_indexes = np.zeros(len(sketches))
	for i in range(len(res)):
		(intersection_cardinality, containment_index, jaccard_index) = res[i]
		intersection_cardinalities[i] = intersection_cardinality
		containment_indexes[i] = containment_index
		jaccard_indexes[i] = jaccard_index

	d = {'intersection': intersection_cardinalities, 'containment index': containment_indexes, 'jaccard index': jaccard_indexes}
	# Use only the basenames to label the rows (if requested)
	if base_name is True:
		df = pd.DataFrame(d, map(os.path.basename, training_file_names))
	else:
		df = pd.DataFrame(d, training_file_names)

	# Only get the rows above a certain threshold
	if coverage_threshold <= 0:
		est_threshold = 0
	else:
		est_threshold = threshold_calc(num_hashes, coverage_threshold, fprate, confidence)
	filtered_results = df[df['containment index'] > est_threshold].sort_values('containment index', ascending=False)
	# Export the results
	filtered_results.to_csv(results_file, index=True, encoding='utf-8')


if __name__ == "__main__":
	main()
