# here I will manually re-create the steps that kmc is supposed to be doing, and compare where the differences are taking place

import os
import sys
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
	from CMash.Make import MakeTSTNew
	from Query import Create
	from Query import Intersect
	from Query import Counters
	from Query import Containment
	from Query import PostProcess
except ImportError:
	try:
		import MinHash as MH
		import Create
		import Intersect
		import Counters
		import Containment
		import PostProcess
	except ImportError:
		try:
			sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
			from CMash import MinHash as MH
			from CMash.Make import MakeTSTNew
			from CMash.Query import Create  # fix relative imports
			from CMash.Query import Intersect
			from CMash.Query import Counters
			from CMash.Query import Containment
			from CMash.Query import PostProcess
		except:
			print("Stupid IDE relative imports...")
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
from itertools import *
import argparse
import khmer
import marisa_trie as mt
import subprocess
import sys
import re

# for IDE REPL testing
os.chdir("/home/dkoslicki/Desktop/CMash/tests/script_tests_debug15/python_output")

# Import the database and dump the Kmers
reads_path = "../../Organisms/taxid_1192839_4_genomic.fna.gz"
training_path = "../TrainingDatabase.h5"
input_type = 'fasta'
threads = 16
temp_dir = "."
verbose = True
I = Intersect(reads_path, training_path, input_type=input_type, threads=threads, temp_dir=temp_dir, verbose=verbose)


#####################################################
# dumping database k-mers: count_training_kmers
# this is KMC-free, so I can just call Isaac's code
I.cmashDump = "TrainingDatabase_dump.fa"
I.dump_training_kmers()

# dump the k-mers using KMC
I.db_kmers_loc = "TrainingDatabase_dump"
I.count_training_kmers()
# FIXME: problem is here: the output of KMC is:
#Stats:
#   No. of k-mers below min. threshold :            0
#   No. of k-mers above max. threshold :            0
#   No. of unique k-mers               :            3  # <-------
#   No. of unique counted k-mers       :            3  # <-------
#   Total no. of k-mers                :            3  # <-------
#   Total no. of reads                 :            1  # <-------
#   Total no. of super-k-mers          :            0
# and:
# $ kmc_dump TrainingDatabase_dump /dev/fd/1
# AAAATCGCTC      1
# AAGTACTGAA      1
# ATACATAGCA      1

# NOTE: after changing kmc from -fa to -fm (since Isaac is dumping in multifasta format, we get:
#Stats:
#   No. of k-mers below min. threshold :            0
#   No. of k-mers above max. threshold :            0
#   No. of unique k-mers               :          544
#   No. of unique counted k-mers       :          544
#   Total no. of k-mers                :         1000
#   Total no. of sequences             :         1000
#   Total no. of super-k-mers          :            0

# which looks much more correct

# check this in python
python_training_kmers = set()
CEs = MH.import_multiple_from_single_hdf5(training_path)
for CE in CEs:
	for kmer in CE._kmers:
		kmer_rc = khmer.reverse_complement(kmer)
		if kmer < kmer_rc:
			python_training_kmers.add(kmer)
		else:
			python_training_kmers.add(kmer_rc)

print(f"Python's count of number of k-mers: {len(python_training_kmers)}")

#with open("python_training_canonical_kmers.txt", 'w') as fid:
#	for kmer in true_canonical_kmers:
#		fid.write(f"{kmer}\n")

# dump kmc version of the canonical kmers
result = subprocess.run(f"kmc_dump TrainingDatabase_dump /dev/fd/1", capture_output=True, shell=True)
kmc_training_kmers = set(map(lambda x: x.split('\t')[0], result.stdout.decode('utf-8').split('\n')))
kmc_training_kmers.remove('')
print(f"KMC's count of number of k-mers: {len(kmc_training_kmers)}")

if sorted(list(python_training_kmers)) == sorted(list(kmc_training_kmers)):
	print("Yes! Python and KMC agree on the database dumped k-mers")
else:
	raise Exception("NO! Python and KMC DO NOT agree on the database dumped k-mers")


#####################################################
# test the counting of input k-mers count_input_kmers()
I.count_input_kmers()

# do it in python
python_read_kmers = set()
fid = khmer.ReadParser(reads_path)
for record in fid:
	seq = record.sequence
	for kmer in MH.kmers(seq, I.ksize):
		kmer_rc = khmer.reverse_complement(kmer)
		if kmer < kmer_rc:
			python_read_kmers.add(kmer)
		else:
			python_read_kmers.add(kmer_rc)

# do it with KMC
I.count_input_kmers()
result = subprocess.run(f"kmc_dump {I.reads_kmc_out_file} /dev/fd/1", capture_output=True, shell=True)
kmc_reads_kmers = set(map(lambda x: x.split('\t')[0], result.stdout.decode('utf-8').split('\n')))
kmc_reads_kmers.remove('')
if sorted(list(python_read_kmers)) == sorted(list(kmc_reads_kmers)):
	print("Yes! Python and KMC agree on the read dumped k-mers")
else:
	raise Exception("NO! Python and KMC DO NOT agree on the read dumped k-mers")


#####################################################
# test the intersection: intersect()

I.intersect()

# do the intersection in python
python_intersection_kmers = python_read_kmers.intersection(python_training_kmers)

# do the intersection with KMC
# FIXME: kmc_dump doesn't like output to stdout
#result = subprocess.run(f"kmc_dump {I.intersection_kmc_dump_file} /dev/fd/1", capture_output=True, shell=True)
result = subprocess.run(f"kmc_dump {I.intersection_kmc_out_file} /dev/fd/1", capture_output=True, shell=True)
kmc_intersection_kmers = set(map(lambda x: x.split('\t')[0], result.stdout.decode('utf-8').split('\n')))
kmc_intersection_kmers.remove('')
if sorted(list(python_intersection_kmers)) == sorted(list(kmc_intersection_kmers)):
	print("Yes! Python and KMC agree on the intersection k-mers")
else:
	raise Exception("NO! Python and KMC DO NOT agree on the intersection k-mers")

# Make sure they have been written in FASTA correctly
fasta_intersection_kmers = set()
fid = khmer.ReadParser(I.intersection_kmc_dump_file+'.fa')
for record in fid:
	seq = record.sequence
	fasta_intersection_kmers.add(seq)

if sorted(list(python_intersection_kmers)) == sorted(list(fasta_intersection_kmers)):
	print("Yes! Python and fasta dump agree on the intersection k-mers")
else:
	raise Exception("NO! Python and fasta dump DO NOT agree on the intersection k-mers")


#########################################################
# test the StreamingQueryDNADatabas.py code on the original file vs the kmc reduced file
def parseNumList(input):
	"""Thank you stack overflow"""
	m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', input)
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end+1, increment))

TST_file = "../TrainingDatabase.tst"
k_range_str = '5-15-2'
k_range = parseNumList(k_range_str)
C = Create(training_database_file=training_path, bloom_filter_file="", TST_file=TST_file, k_range=k_range)
C.import_TST()
C.create_BF_prefilter()

# on the original reads
counter = Counters(tree=C.tree, k_range=C.k_range, all_kmers_bf=C.all_kmers_bf)
fid = khmer.ReadParser(reads_path)
to_proc = [record.sequence for record in fid]
fid.close()
res = map(counter.process_seq, to_proc)
flattened_res = [item for sublist in res if sublist for item in sublist]
all_reads_match_tuples = list(set(flattened_res))

# on the kmc dumped reads
counter = Counters(tree=C.tree, k_range=C.k_range, all_kmers_bf=C.all_kmers_bf)
fid = khmer.ReadParser(I.final_out_file)
to_proc = [record.sequence for record in fid]
fid.close()
res = map(counter.process_seq, to_proc)
flattened_res = [item for sublist in res if sublist for item in sublist]
kmc_reads_match_tuples = list(set(flattened_res))

if sorted(all_reads_match_tuples) == sorted(kmc_reads_match_tuples):
	print("Yes! match_tuples are the same for original reads or kmc reads")
else:
	raise Exception("NO! match_tuples are NOT the same for original reads or kmc reads")

missing_matches = set(all_reads_match_tuples) - set(kmc_reads_match_tuples)
missing_match = list(missing_matches)[0]
ce_index = missing_match[0]
k_range_loc = missing_match[1]
sketch_index = missing_match[2]
print(f"A missing kmer: {CEs[ce_index]._kmers[sketch_index]}")
print(f"A missing kmer_rc: {khmer.reverse_complement(CEs[ce_index]._kmers[sketch_index])}")
print(f"at kmer size: {k_range[k_range_loc]}")
missing_kmer_trunc = CEs[ce_index]._kmers[sketch_index][0:k_range[k_range_loc]]
print(f"kmer trunc: {missing_kmer_trunc}")
print(f"kmer_rc trunc: {khmer.reverse_complement(CEs[ce_index]._kmers[sketch_index])[0:k_range[k_range_loc]]}")
print(f"trunc kmer_rc: {khmer.reverse_complement(CEs[ce_index]._kmers[sketch_index][0:k_range[k_range_loc]])}")

# Look for that k-mer in the query file:
result = subprocess.run(f"zgrep -c -i {missing_kmer_trunc} {reads_path}", capture_output=True, shell=True)
print(f"matches in reads to missing_kmer_trunc: {result.stdout}")

# look for that k-mer in the intersection file
result = subprocess.run(f"zgrep -c -i {missing_kmer_trunc} {I.final_out_file}", capture_output=True, shell=True)
print(f"matches in intersection to missing_kmer_trunc: {result.stdout}")

# look for the k-mer in the reads dump
result = subprocess.run(f"kmc_dump {I.reads_kmc_out_file} /dev/fd/1 | cut -f 1 | grep -c -i {missing_kmer_trunc}", capture_output=True, shell=True)
print(f"matches in kmc reads dump to missing_kmer_trunc: {result.stdout}")
# so indeed, the kmer_trunc is in the reads dump, but it only matches as prefixes to k-mers that aren't in the training database:
# kmc_dump reads_15mers_dump /dev/fd/1 | cut -f 1 | grep -i TGCCCTGTGGC
#ATGCCCTGTGGCTGG
#ATGCTGCCCTGTGGC
#CCCCTGCCCTGTGGC
#CTGCCCTGTGGCAGC
#TGCCCTGTGGCAGCA  #<---------------
# grep -c TGCCCTGTGGCAGCA TrainingDatabase_dump.fa
# 0

