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
chunk_size = 5000
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

###########################################################################################
# Next, run kmc on this thing
# with a bash script, use the following:
# /home/pi/koslickd/KMC/./kmc -v -k60 -m200 -sm -fm -ci0 -cs3 -t48 -jlog_train NathanRefSeqTraining60mers.fa NathanRefSeq60mers /scratch/kmc_temp/
#

############################################################################################
# Next, run kmc on the input file
# in bash:
# /home/pi/koslickd/KMC/./kmc -v -k60 -m200 -sm -fq -ci2 -cs3 -t48 -jlog_sample /nfs1/Koslicki_Lab/koslickd/MiCOPCMash/CAMITest/RL_S001__insert_270.fq.gz CAMILow60mers /scratch/kmc_temp/
# 7 min, 38 seconds

#############################################################################################
# Next, intersect them
# in bash:
# /home/pi/koslickd/KMC/./kmc_tools simple NathanRefSeq60mers CAMILow60mers intersect intersection

##############################################################################################
# Dump them and turn into a fasta file
# /home/pi/koslickd/KMC/./kmc_dump intersection intersection_dump
# cat intersection_dump | cut -f1 | sed 's/^/>seq\n/g' > intersection_dump.fa

#############################################################################################
# Then run CMash on intersection_dump.fa
#/nfs1/Koslicki_Lab/koslickd/MiCOPCMash/CMashVE/bin/python /nfs1/Koslicki_Lab/koslickd/MiCOPCMash/CMash/scripts/StreamingQueryDNADatabase.py /nfs1/Koslicki_Lab/koslickd/KMC_test/intersection_dump.fa /nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/micopdb_n_1000_k_60.h5 /nfs1/Koslicki_Lab/koslickd/KMC_test/RL_S001__insert_270_classified_kmc.csv 30-60-10 -c 0 -r 1000000 -v -f /nfs1/Koslicki_Lab/koslickd/MiCOPCMash/TrainingData/NathanRefSeq/micopdb_n_1000_k_60_prefilter_30_60_10.bf
