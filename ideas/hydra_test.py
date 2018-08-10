import hydra
from _hydra import BloomCalculations, BloomFilter, \
	UnsupportedOperationException
from hydra import WritingBloomFilter, murmur_hash, ReadingBloomFilter
import sys
sys.path.append('/home/dkoslicki/Desktop/CMash/CMash')
import MinHash as MH
import timeit
import random
import numpy as np
import os


#num_in_bloom = 1000*100
num_in_bloom = 10000000

# delete the old bloom filter so we don't try to stuff more things into it
try:
	os.remove("test.bloom")
	os.remove("test.bloom.desc")
except:
	pass

bloom = WritingBloomFilter(num_in_bloom, 0.01, "test.bloom")

# Read it back in
#bloom = ReadingBloomFilter("test.bloom")

# If you want to add the real k-mers
#CEs=MH.import_multiple_from_single_hdf5('/home/dkoslicki/Desktop/CMash/data/SmallGenomes.h5')
#all_kmers = set()
#for CE in CEs:
#	for kmer in CE._kmers:
#		all_kmers.add("%s" % kmer)
#for CE in CEs:
#	for kmer in CE._kmers:
# 		bloom.add("%s" % kmer)

# If you want to add random kmers
all_kmers = set()
for _ in range(num_in_bloom):
	kmer = "".join(np.random.choice(["A", "C", "T", "G"], 60))
	bloom.add(kmer)
	all_kmers.add(kmer)

# Test the timing for a bloom query
N = 10000
i = 0
t0 = timeit.default_timer()
for _ in range(N):
	#kmer = random.sample(all_kmers, 1)[0]
	kmer = "".join(np.random.choice(["A", "C", "T", "G"], 60))
	if kmer in bloom:
		i += 1
	else:
		pass
t1 = timeit.default_timer()
print(t1-t0)

# Test the timing for a set
i = 0
t0 = timeit.default_timer()
for _ in range(N):
	#kmer = random.sample(all_kmers, 1)[0]
	kmer = "".join(np.random.choice(["A", "C", "T", "G"], 60))
	if kmer in all_kmers:
		i += 1
	else:
		pass
t1 = timeit.default_timer()
print(t1-t0)
