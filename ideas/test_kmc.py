# This code will test out the idea of using kmc to
# 1. quickly enumerate the k-mers
# 2. intersect these with the training database, output as fasta
# 3. use that reduced fasta of intersecting kmers as the query to CMash

####################################################################
# First, I will need to dump the training database to a fasta file
from CMash import MinHash as MH
import os
import blist
import itertools

training_out_file = '/nfs1/Koslicki_Lab/koslickd/KMC_test/NathanRefSeqTraining60mers.fa'
training_data ='/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/micopdb_n_1000_k_60.h5'
training_file_names = "/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/absolute_file_names.txt"

print("reading file list")
file_names = []
with open(training_file_names, 'r') as fid:
	for line in fid.readlines():
		line = line.strip()
		file_names.append(os.path.basename(line))

print("importing kmers")
chunk_size = 1000
iter = 0
with open(training_out_file, 'w') as fid:
	for file_iter in xrange(0, len(file_names), chunk_size):
		print("on file: %d" % file_iter)
		file_names_subset = file_names[file_iter:file_iter+chunk_size]
		sketchs = MH.import_multiple_from_single_hdf5(training_data, import_list=file_names_subset)
		all_kmers = itertools.chain.from_iterable(sketch._kmers for sketch in sketchs)
		print("forming the set")
		all_kmers_set = set(all_kmers)
		to_write = ""
		for kmer in all_kmers_set:
			to_write += ">seq%d\n" % iter
			to_write += "%s\n" % kmer
			iter += 1
		fid.write(to_write)

