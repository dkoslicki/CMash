# This script will test the issue where the file name order when you form the training database
# does not match the file name order when you read the training database back in
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

num_train = 100

# Read in all the file names
all_file_names = []
with open("/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/absolute_file_names.txt", "r") as fid:
	for line in fid.readlines():
		all_file_names.append(line.strip())

# form the training database on a few of them
subset_file_names_file = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/TestFileNameOrder/FileNames.txt"
with open(subset_file_names_file, "w") as fid:
	for i in range(num_train):
		fid.write("%s\n" % all_file_names[i])

out_hdf5_file = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/TestFileNameOrder/TrainingData.h5"
python = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/CMashVE/bin/python "
script = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/CMash/scripts/MakeStreamingDNADatabase.py "
script_args = subset_file_names_file + " " + out_hdf5_file + " -n 1000 -k 60"
os.system(python + script + script_args)

# Import the HDF5 file
sketches = MH.import_multiple_from_single_hdf5(out_hdf5_file)
tree = mt.Trie()
tree.load(out_hdf5_file.split('.')[0] + ".tst")

for sketch_index in range(num_train):
	for kmer in sketches[sketch_index]._kmers:
		is_correct = False
		for hit in tree.keys(kmer):
			hit_split = hit.split('x')
			tree_sketch_index = int(hit_split[1])
			if tree_sketch_index == sketch_index:
				is_correct = True
				break
		if not is_correct:
			raise Exception("Mismatch: sketch index was %d while in the tree it's %d: %s" % (sketch_index, tree_sketch_index, hit))