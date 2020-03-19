import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from CMash.Query import *

# Create module tests
def test_initialize_Create():
	# include error handling tests?
	pass

def test_Create_import_TST():
	pass

def test_Create_BF_prefilter():
	pass

# Counters module tests
def test_initialize_Counters():
	pass

def test_Counters_return_matches():
	pass

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

def test_PostProcess_find_non_unique_kmers():
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

	# run through them
	for func in test_functions:
		getattr(test_Query, func)()


if __name__ == "__main__":
	main()
