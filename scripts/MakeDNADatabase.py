#! /usr/bin/env python
# This script will make a training database of hashes
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
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
from itertools import *
import argparse
import khmer
from khmer.khmer_args import optimal_size

# This function will make a single min hash sketch upon being given a file name, sketch size, prime, and k-mer size
def make_minhash(genome, max_h, prime, ksize):
	MHS = MH.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, save_kmers='y', input_file_name=genome)
	# Just use HLL to estimate the number of kmers, no need to get exact count
	hll = khmer.HLLCounter(0.01, ksize)
	hll.consume_seqfile(genome)
	MHS._true_num_kmers = hll.estimate_cardinality()
	MHS.input_file_name = genome
	return MHS


# Unwrap for Python2.7
def make_minhash_star(arg):
	return make_minhash(*arg)


def main():
	parser = argparse.ArgumentParser(description="This script creates training/reference sketches for each FASTA/Q file"
									" listed in the input file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
	parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=500)
	parser.add_argument('-k', '--k_size', type=int, help="K-mer size", default=21)
	parser.add_argument('-i', '--intersect_nodegraph', action="store_true",
						help="Optional flag to export Nodegraph file containing all k-mers in the training database. Saved in same location as out_file.")
	parser.add_argument('in_file', help="Input file: file containing (absolute) file names of training genomes.")
	parser.add_argument('out_file', help='Output training database/reference file (in HDF5 format)')
	args = parser.parse_args()
	num_threads = args.threads
	prime = args.prime  # taking hashes mod this prime
	ksize = args.k_size
	max_h = args.num_hashes
	input_file_names = os.path.abspath(args.in_file)
	if not os.path.exists(input_file_names):
		raise Exception("Input file %s does not exist." % input_file_names)
	out_file = os.path.abspath(args.out_file)
	if args.intersect_nodegraph is True:
		intersect_nodegraph_file = os.path.splitext(out_file)[0] + ".intersect.Nodegraph"
	else:
		intersect_nodegraph_file = None

	file_names = list()
	fid = open(input_file_names, 'r')
	for line in fid.readlines():
		line = line.strip()
		if not os.path.exists(line):
			raise Exception("Training genome %s does not exist." % line)
		file_names.append(line)
	fid.close()

	# Open the pool and make the sketches
	pool = Pool(processes=num_threads)
	genome_sketches = pool.map(make_minhash_star, zip(file_names, repeat(max_h), repeat(prime), repeat(ksize)))

	# Export all the sketches
	MH.export_multiple_to_single_hdf5(genome_sketches, out_file)

	# If requested, save all the k-mers into a big Nodegraph (unfortunately, need to pass through the data again since we
	# a-priori don't know how big of a table we need to make
	if intersect_nodegraph_file is not None:
		total_num_kmers = 0
		for sketch in genome_sketches:
			total_num_kmers += sketch._true_num_kmers
		res = optimal_size(total_num_kmers, fp_rate=0.001)
		intersect_nodegraph = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
		for file_name in file_names:
			intersect_nodegraph.consume_seqfile(file_name)
		intersect_nodegraph.save(intersect_nodegraph_file)

if __name__ == "__main__":
	main()
