#! /usr/bin/env python
import os
import sys
import khmer
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
except ImportError:
	try:
		import MinHash as MH
	except ImportError:
		sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
		from CMash import MinHash as MH
import argparse
from argparse import ArgumentTypeError
import re
from hydra import WritingBloomFilter


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
	parser = argparse.ArgumentParser(description="This script pre-computes a prefilter for the StreamingQueryDNADatabase (for a fixed k-range).", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('reference_file', help='Training database/reference file (in HDF5 format). Created with MakeStreamingDNADatabase.py')
	parser.add_argument('out_file', help='Output prefilter file.')
	parser.add_argument('range', type=parseNumList, help="Range of k-mer sizes in the formate <start>-<end>-<increment>."
														   " So 5-10-2 means [5, 7, 9]. If <end> is larger than the k-mer size"
														   "of the training data, this will automatically be reduced.")

	# read in the arguments
	args = parser.parse_args()
	k_range = args.range
	if k_range is None:
		raise Exception("The --range argument is required, no matter what the help menu says.")
	training_data = args.reference_file
	results_file = args.out_file
	results_file = os.path.abspath(results_file)

	# Import data and error checking
	# Query file
	if not os.path.exists(training_data):
		raise Exception("Training/reference file %s does not exist." % training_data)

	# Training data
	sketches = MH.import_multiple_from_single_hdf5(training_data)
	if sketches[0]._kmers is None:
		raise Exception(
			"For some reason, the k-mers were not saved when the database was created. Try running MakeDNADatabase.py again.")
	num_hashes = len(sketches[0]._kmers)  # note: this is relying on the fact that the sketches were properly constructed
	max_ksize = sketches[0].ksize

	# adjust the k-range if necessary
	k_range = [val for val in k_range if val <= max_ksize]

	# all the k-mers of interest in a set (as a pre-filter)
	try:
		all_kmers_bf = WritingBloomFilter(len(sketches)*len(k_range)*num_hashes*2, 0.01, results_file)
		for sketch in sketches:
			for kmer in sketch._kmers:
				for ksize in k_range:
					kmer_str = kmer[0:ksize].decode('utf-8')
					all_kmers_bf.add(kmer_str)  # put all the k-mers and the appropriate suffixes in
					all_kmers_bf.add(khmer.reverse_complement(kmer_str))  # also add the reverse complement
	except IOError:
		print("No such file or directory/error opening file: %s" % results_file)
		sys.exit(1)
