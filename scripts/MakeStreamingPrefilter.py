#! /usr/bin/env python
import os
import sys
import khmer
# The following is for ease of development (so I don't need to keep re-installing the tool)
#try:
#	from CMash import MinHash as MH
#except ImportError:
#	try:
#		import MinHash as MH
#	except ImportError:
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from CMash import MinHash as MH
from CMash.Query import Create
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
	# Training file
	if not os.path.exists(training_data):
		raise Exception("Training/reference file %s does not exist." % training_data)

	# get the name of the TST
	TST_file = os.path.splitext(training_data)[0] + ".tst"
	if not os.path.exists(TST_file):
		raise Exception(f"Ternary search tree {TST_file} does not exist. Did you run CreateStreamingDNADatabase.py?")

	# TODO: adjust the k-range if necessary
	#  k_range = [val for val in k_range if val <= max_ksize]  # max_ksize is unknown without reading in the whole
	#  training database, which I don't want to do

	# create the class
	C = Create(training_database_file=training_data, bloom_filter_file="", TST_file=TST_file, k_range=k_range)
	# import the TST
	C.import_TST()
	# create and export the BF
	C.create_BF_prefilter(result_file=results_file)
