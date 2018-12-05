# This code will test out the idea of using kmc to
# 1. quickly enumerate the k-mers
# 2. intersect these with the training database, output as fasta
# 3. use that reduced fasta of intersecting kmers as the query to CMash

####################################################################
# First, I will need to dump the training database to a fasta file
from CMash import MinHash as MH
import os
import blist

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
all_kmers = blist.blist()
iter = 0
for file_name in file_names:
	iter += 1
	sketch = MH.import_multiple_from_single_hdf5(training_data, import_list=[file_name])[0]
	all_kmers += sketch._kmers
	if iter % 1000 == 0:
		print("Finished %d files" % iter)

print("forming the set")
all_kmers_set = set(all_kmers)

num_kmers = len(all_kmers_set)
print("exporting results")
with open(training_out_file, 'w') as fid:
	iter = 0
	for kmer in all_kmers_set:
		fid.write(">seq_%d\n" % iter)
		fid.write("%s\n" % kmer)
		iter += 1
		if iter % 1000000 == 0:
			print("On kmer %d of %d" % (iter, num_kmers))

##########################################################################

