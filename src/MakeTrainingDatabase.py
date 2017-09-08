# This script will make a training database of hashes
import os
import screed
import MinHash as MH
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
from itertools import *


num_threads = multiprocessing.cpu_count()
prime = 9999999999971  # taking hashes mod this prime
ksize = 21  # k-mer length
max_h = 500  # max number of hashes in sketch
input_file_names = os.path.abspath('../data/FileNames.txt')
out_file = os.path.abspath('../data/AllSketches.h5')

file_names = list()
fid = open(input_file_names, 'r')
for line in fid.readlines():
	file_names.append(line.strip())
fid.close()


# This function will make a single min hash sketch upon being given a file name, sketch size, prime, and k-mer size
def make_minhash(genome, max_h, prime, ksize):
	kmers = set()
	name = os.path.basename(genome)
	MHS = MH.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, save_kmers='y')
	for record in screed.open(genome):
		seq = record.sequence
		for i in range(len(seq) - ksize + 1):
			kmer = seq[i:i+ksize]
			#kmer_rev = khmer.reverse_complement(kmer)  # No need to care about revcomp since nodegraph ID's them
			#if kmer < kmer_rev:
			kmers.add(kmer)
			MHS.add(kmer)
			#else:
			#	kmers.add(kmer_rev)
			#	MHS.add(kmer_rev)
	MHS._true_num_kmers = len(kmers)
	MHS.input_file_name = genome
	# Export the hash k-mers
	#fid = open(os.path.abspath(os.path.join('../data/Viruses/', name + ".Hash21mers.fa")), 'w')
	#for kmer in MHS._kmers:
	#	fid.write(">\n%s\n" % kmer)
	#fid.close()
	return MHS


# Unwrap for Python2.7
def make_minhash_star(arg):
	return make_minhash(*arg)

pool = Pool(processes=num_threads)
genome_sketches = pool.map(make_minhash_star, zip(file_names, repeat(max_h), repeat(prime), repeat(ksize)))

# Export all the sketches
MH.export_multiple_to_single_hdf5(genome_sketches, out_file)
