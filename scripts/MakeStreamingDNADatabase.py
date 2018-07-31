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
import tst

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
	parser.add_argument('-k', '--k_size', type=int, help="k-mer size", default=21)
	parser.add_argument('in_file', help="Input file: file containing (absolute) file names of training genomes.")
	parser.add_argument('out_file', help='Output training database/reference file (in HDF5 format). An additional file '
										 '(ending in .tst) will also be created in the same directory with the same base name.')
	args = parser.parse_args()
	num_threads = args.threads
	prime = args.prime  # taking hashes mod this prime
	ksize = args.k_size
	max_h = args.num_hashes
	input_file_names = os.path.abspath(args.in_file)
	if not os.path.exists(input_file_names):
		raise Exception("Input file %s does not exist." % input_file_names)
	out_file = os.path.abspath(args.out_file)

	# check for and make filename for tst file
	streaming_database_file = os.path.splitext(out_file)[0] + ".tst"
	streaming_database_file = os.path.abspath(streaming_database_file)

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

	# Save the ternary search tree
	tree = tst.TST()  # tst array
	for i in range(len(genome_sketches)):
		for kmer_index in range(len(genome_sketches[i]._kmers)):
			kmer = genome_sketches[i]._kmers[kmer_index]
			tree[kmer + 'x' + str(i) + 'x' + str(kmer_index)] = True  # format here is kmer+x+hash_index+kmer_index
	tree.write_to_file(streaming_database_file)


if __name__ == "__main__":
	main()

