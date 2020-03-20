import os
import sys
import tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from CMash.Query import *
from CMash import MinHash as MH
import khmer

# FIXME: could probably do all the data creation, module initialization, and method calling, and then have the tests
# FIXME: just test the data

# create some test data
# First, the TST
seq1 = "ATCGTATGAGTATCGTCGATGCATGCATCGATGCATGCTACGTATCGCATGCATG"
seq2 = "ATCTACTCAACATTAACTACTCATATTAACTCACATTCATATCCATACTACTCGT"
seq3 = "ACTCATGTTAGATCGATATTGACTGATGACTCGTTGCACTGCATGCTGCATGATGC"
seq4 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
CE1 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y')
CE2 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y')
CE3 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y')
CE4 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y')
CE1.add_sequence(seq1)
CE2.add_sequence(seq2)
CE3.add_sequence(seq3)
CE4.add_sequence(seq4)
# CE's must have input names
CE1.input_file_name = "seq1"
CE2.input_file_name = "seq2"
CE3.input_file_name = "seq3"
CE4.input_file_name = "seq4"
CEs = [CE1, CE2, CE3, CE4]
temp_database_file = tempfile.mktemp()
MH.export_multiple_to_single_hdf5(CEs, temp_database_file)

# And create the TST
to_insert = set()
# add both the original k-mer and the reverse complement, as the MinHashes were created without reverse complement
for i in range(len(CEs)):
	for kmer_index in range(len(CEs[i]._kmers)):
		# normal kmer
		kmer = CEs[i]._kmers[kmer_index]
		if kmer:
			to_insert.add(kmer + 'x' + str(i) + 'x' + str(kmer_index))  # format here is kmer+x+hash_index+kmer_index
			# rev-comp kmer
			kmer = khmer.reverse_complement(CEs[i]._kmers[kmer_index])
			to_insert.add(kmer + 'x' + str(i) + 'x' + str(kmer_index))  # format here is kmer+x+hash_index+kmer_index

# export the TST
tree = mt.Trie(to_insert)
temp_TST_file = tempfile.mktemp()
tree.save(temp_TST_file)

# Create module tests
def test_initialize_Create():
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=[1, 3, 5])
	assert C.k_range == [1, 3, 5]
	assert C.TST_file == temp_TST_file
	pass

def test_Create_import_TST():
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=[1, 3, 5])
	C.import_TST()
	# make sure the k-mers and their reverse complements have been added
	assert C.tree.keys("AAA")
	assert C.tree.keys("TTT")
	assert C.tree.keys("ACT")
	assert C.tree.keys("AGT")
	# Make sure the correct kmers are being identified with the correct k-mers
	assert C.tree.keys("AAA")[0] == "AAAx3x0"
	assert C.tree.keys("TTT")[0] == "TTTx3x0"
	pass

def test_Create_BF_prefilter():
	k_range = [1, 3, 5]
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	C.import_TST()
	C.create_BF_prefilter()

	# Make sure each TST kmers has been inserted into the bloom tree
	for kmer_with_info in C.tree.keys():
		kmer = kmer_with_info.split('x')[0]
		assert kmer in C.all_kmers_bf

	# Make sure all the reverse complements are in there too
	for kmer_with_info in C.tree.keys():
		kmer = kmer_with_info.split('x')[0]
		kmer = khmer.reverse_complement(kmer)
		assert kmer in C.all_kmers_bf

	# go through each individual sequence in the sketches and make sure them and their rev-comps are in the BF
	for CE in CEs:
		for kmer in CE._kmers:
			if kmer:
				for k_size in k_range:
					assert kmer[0:k_size] in C.all_kmers_bf
					assert khmer.reverse_complement(kmer[0:k_size]) in C.all_kmers_bf


# Counters module tests
def test_initialize_Counters():
	k_range = [1, 3, 5]
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	C.import_TST()
	C.create_BF_prefilter()
	counters = Counters(tree=C.tree, k_range=k_range, all_kmers_bf=C.all_kmers_bf)

def test_Counters_return_matches():
	k_range = [1, 3, 5]
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	C.import_TST()
	C.create_BF_prefilter()
	counters = Counters(tree=C.tree, k_range=k_range, all_kmers_bf=C.all_kmers_bf)


def test_Counters_process_seq():
	pass


# Containment module tests
def test_initialize_Containment():
	pass

def test_Containment_create_to_hit_matrices():
	pass

def test_Containment_create_containment_indicies():
	pass

def test_Containment_create_data_frame():
	pass


# PostProcess module tests
def test_initialize_PostProcess():
	pass

def test_PostProcess_prepare_post_process():
	pass

def test_PostProcess_find_kmers_in_filtered_results():
	pass

def test_PostProcess_find_unique_kmers():
	pass

def test_PostProcess_find_non_unique_kmers_reduce_hit_matrices():
	pass

def test_PostProcess_create_post_containment_indicies():
	pass

def test_PostProces_create_data_frame():
	print("You are here")
	pass


# and for the lonely function
def test_return_data_frame():
	pass


# run through all the tests
def main():
	"""
	Main testing script. Run through all the tests above that start with "test_" and execute them.
	Designed to be run via the command line with $python test_Query.py
	:return:
	:rtype:
	"""
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	import test_Query  # curiously enough, this works!

	# get all the test functions
	test_functions = []
	for key in test_Query.__dict__:
		if key[0:5] == "test_":
			test_functions.append(key)

	# run through them. And they're in order of how they were written!
	for func in test_functions:
		getattr(test_Query, func)()


if __name__ == "__main__":
	main()
	print("All tests passed successfully!")
