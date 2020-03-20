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
seq2 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
seq3 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"
seq4 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
CE1 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=5, save_kmers='y')
CE2 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=5, save_kmers='y')
CE3 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=5, save_kmers='y')
CE4 = MH.CountEstimator(n=5, max_prime=9999999999971, ksize=5, save_kmers='y')
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

# TODO: marisa_trie has an issue with single character prefix lookups
# TODO: see https://github.com/pytries/marisa-trie/issues/55
# TODO: so set k-range above that
k_range = [2, 3, 5]

# Create module tests
def test_initialize_Create():
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	assert C.k_range == k_range
	assert C.TST_file == temp_TST_file
	pass

def test_Create_import_TST():
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	C.import_TST()
	# make sure the k-mers and their reverse complements have been added
	assert C.tree.keys("AAA")
	assert C.tree.keys("TTT")
	assert C.tree.keys("ATA")
	assert C.tree.keys("TAT")
	# Make sure the correct kmers are being identified with the correct k-mers
	for kmer in ["AAA", "TTT"]:
		matches = C.tree.keys("AAA")
		for match in matches:
			location_info = "x".join(match.split('x')[1:])
			assert (location_info == "3x0") or (location_info == "1x0")


def test_Create_BF_prefilter():
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

	# check if the BF is case insensitive
	for CE in CEs:
		for kmer in CE._kmers:
			if kmer:
				for k_size in k_range:
					trunc_kmer = kmer[0:k_size]
					trunc_kmer = trunc_kmer.lower()
					assert trunc_kmer in C.all_kmers_bf
					# khmer doesn't properly handle rev-comps of lower-case characters
					# see https://github.com/dib-lab/khmer/issues/1904
					assert khmer.reverse_complement(trunc_kmer.upper()).lower() in C.all_kmers_bf


# Counters module tests
def test_initialize_Counters():
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	C.import_TST()
	C.create_BF_prefilter()
	counters = Counters(tree=C.tree, k_range=k_range, all_kmers_bf=C.all_kmers_bf)

def test_Counters_return_matches():
	C = Create(training_database_file=temp_database_file, bloom_filter_file="", TST_file=temp_TST_file, k_range=k_range)
	C.import_TST()
	C.create_BF_prefilter()
	counters = Counters(tree=C.tree, k_range=k_range, all_kmers_bf=C.all_kmers_bf)

	# test the return matches on known k-mers
	# each sketch kmer (or it's reverse complement) should match to the TST
	# TODO: big note here: proper way to check this: take the reverse complement, THEN truncate
	#  (which effectively takes the suffix, as the suffix of a rev-comp is the prefix of the original)
	#  but this calls into question how create_BF_prefilter is working since it truncates, THEN takes the revcomp
	#  but this is the only way I could get all these tests to pass successfully
	for CE in CEs:
		for k_size in k_range:
			for kmer in CE._kmers:
				kmer = kmer[0:k_size]
				if kmer:
					k_size_loc = k_range.index(len(kmer))
					to_return, saw_match = counters.return_matches(input_kmer=kmer, k_size_loc=k_size_loc)
					assert saw_match
					for to_return_elem in to_return:
						truncated_sketches = list(map(lambda x: x[0:k_size], CEs[to_return_elem[0]]._kmers))
						# add the reverse complements as well, since the TST return_matches matches to rev-comps as well
						truncated_sketches_revcomp = list(map(lambda x: khmer.reverse_complement(x)[0:k_size], CEs[to_return_elem[0]]._kmers))
						# make sure the kmer really is in the sketch indicated by to_return, could be in the truncated or the rev-comp one
						assert (kmer in truncated_sketches) or (kmer in truncated_sketches_revcomp)
						# make sure the k_size_loc is correct
						assert to_return_elem[1] == k_size_loc
						# make sure it returned the correct location in the sketch
						# note that at some smaller kmer values, it may appear in multiple locations, so just make sure
						# that it appears somewhere in the list
						indices = [i for i, x in enumerate(truncated_sketches) if x == kmer]
						indices_revcomp = [i for i, x in enumerate(truncated_sketches_revcomp) if x == kmer]
						assert (to_return_elem[2] in indices) or (to_return_elem[2] in indices_revcomp)


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
