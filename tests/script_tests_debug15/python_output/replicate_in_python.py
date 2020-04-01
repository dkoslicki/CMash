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
true_canonical_kmers = set()
CEs = MH.import_multiple_from_single_hdf5(training_path)
for CE in CEs:
	for kmer in CE._kmers:
		kmer_rc = khmer.reverse_complement(kmer)
		if kmer < kmer_rc:
			true_canonical_kmers.add(kmer)
		else:
			true_canonical_kmers.add(kmer_rc)

print(f"Python's count of number of k-mers: {len(true_canonical_kmers)}")
#with open("python_training_canonical_kmers.txt", 'w') as fid:
#	for kmer in true_canonical_kmers:
#		fid.write(f"{kmer}\n")

# dump kmc version of the canonical kmers
result = subprocess.run(f"kmc_dump TrainingDatabase_dump /dev/fd/1", capture_output=True, shell=True)
kmc_canonical_kmers = set(map(lambda x: x.split('\t')[0], result.stdout.decode('utf-8').split('\n')))
kmc_canonical_kmers.remove('')

if sorted(list(true_canonical_kmers)) == sorted(list(kmc_canonical_kmers)):
	print("Yes! Python and KMC agree on the database dumped k-mers")
