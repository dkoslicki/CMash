from __future__ import print_function
import khmer
import tst  # via pip install git+https://github.com/dkoslicki/pytst2.git#egg=pytst-1.18
import numpy as np
import os
import sys
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
	from CMash import MinHash as MH
except ImportError:
	try:
		import MinHash as MH
	except ImportError:
		sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
		from CMash import MinHash as MH
import time
from multiprocessing import Process, Value, Lock
import multiprocessing
import ctypes
import pandas as pd
import argparse
from argparse import ArgumentParser, ArgumentTypeError
import re
import matplotlib.pyplot as plt

def parseNumList(input):
	m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', input)
	# ^ (or use .split('-'). anyway you like.)
	if not m:
		raise ArgumentTypeError("'" + input + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end+1, increment))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script calculates containment indicies for each of the training/reference sketches"
									" by streaming through the query file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
	parser.add_argument('-c', '--containment_threshold', type=float, help="Only return results with containment index above this "
															  "threshold at the maximum k-mer size.", default=0.1)
	parser.add_argument('-p', '--plot_file', action="store_true", help="Optional flag to specify that a plot of the "
																	   "k-mer curves should also be saved (same basename"
																	   "as the out_file).")
	parser.add_argument('in_file', help="Input file: FASTA/Q file to be processes")
	parser.add_argument('reference_file', help='Training database/reference file (in HDF5 format). Created with MakeStreamingDNADatabase.py')
	parser.add_argument('out_file', help='Output csv file with the containment indices.')
	parser.add_argument('range', type=parseNumList, help="Range of k-mer sizes in the formate <start>-<end>-<increment>."
														   " So 5-10-2 means [5, 7, 9]. If <end> is larger than the k-mer size"
														   "of the training data, this will automatically be reduced.")

	# read in the arguments
	args = parser.parse_args()
	k_range = args.range
	if k_range is None:
		raise Exception("The --range argument is required, no matter what the help menu says.")
	training_data = args.reference_file
	query_file = args.in_file
	results_file = args.out_file
	num_threads = args.threads
	coverage_threshold = args.containment_threshold
	streaming_database_file = os.path.splitext(training_data)[0] + ".tst"  # name of the tst training file
	streaming_database_file = os.path.abspath(streaming_database_file)
	if not os.path.exists(streaming_database_file):
		streaming_database_file = None
	if args.plot_file:
		plot_file = os.path.abspath(os.path.splitext(results_file)[0] + ".png")

	# Import data and error checking
	# Query file
	if not os.path.exists(query_file):
		raise Exception("Query file %s does not exist." % query_file)
	if not os.path.exists(training_data):
		raise Exception("Training/reference file %s does not exist." % training_data)

	# Training data
	sketches = MH.import_multiple_from_single_hdf5(training_data)
	ksizes = set()
	if sketches[0]._kmers is None:
		raise Exception(
			"For some reason, the k-mers were not saved when the database was created. Try running MakeDNADatabase.py again.")
	num_hashes = len(sketches[0]._kmers)
	for i in range(len(sketches)):
		sketch = sketches[i]
		if sketch._kmers is None:
			raise Exception(
				"For some reason, the k-mers were not saved when the database was created. Try running MakeDNADatabase.py again.")
		if len(sketch._kmers) != num_hashes:
			raise Exception("Unequal number of hashes for sketch of %s" % sketch.input_file_name)
		ksizes.add(sketch.ksize)
		if len(ksizes) > 1:
			raise Exception(
				"Training/reference data uses different k-mer sizes. Culprit was %s." % (sketch.input_file_name))

	max_ksize = ksizes.pop()
	# adjust the k-range if necessary
	k_range = [val for val in k_range if val <= max_ksize]

	# Get names of training files for use as rows in returned tabular data
	training_file_names = []
	for i in range(len(sketches)):
		training_file_names.append(sketches[i].input_file_name)

	# Make the tst tree
	if streaming_database_file is None:
		streaming_database_file = os.path.splitext(training_data)[0] + ".tst"
		streaming_database_file = os.path.abspath(streaming_database_file)
		print("It appears a tst training file has not been created (did you remember to use MakeStreamingDNADatabase.py?).")
		print("I'm creating one anyway at: %s" % streaming_database_file)
		print("This may take a while...")
		tree = tst.TST()  # tst array
		for i in range(len(sketches)):
			for kmer_index in range(len(sketches[i]._kmers)):
				kmer = sketches[i]._kmers[kmer_index]
				tree[kmer + 'x' + str(i) + 'x' + str(kmer_index)] = True  # format here is kmer+x+hash_index+kmer_index
		tree.write_to_file(streaming_database_file)
	else:
		tree = tst.TST()
		tree.read_from_file(streaming_database_file)

	# all the k-mers of interest in a set (as a pre-filter)
	all_kmers_set = set()
	for sketch in sketches:
		for kmer in sketch._kmers:
			for ksize in k_range:
				all_kmers_set.add(kmer[0:ksize])  # put all the k-mers and the appropriate suffixes in

	# shared object that will update the intersection counts
	class Counters(object):
		# This class is basically an array of counters (on the same basis as the sketches)
		# it's used to keep track (in a parallel friendly way) of which streamed k-mers went into the training file sketches
		def __init__(self):
			self.ksize = max_ksize
			N = len(sketches)
			to_share = multiprocessing.Array(ctypes.c_bool, N * num_hashes * len(k_range))  # vector
			to_share_np = np.frombuffer(to_share.get_obj(), dtype=ctypes.c_bool)  # get it
			self.vals = to_share_np.reshape(N, num_hashes, len(k_range))
			self.lock = Lock()  # don't need since I don't care about collisions (once matched, always matched)

		def increment(self, hash_index, kmer_index, k_size_loc):
			self.vals[hash_index, kmer_index, k_size_loc] = True  # array

		def value(self, hash_index, k_size_loc):
			return np.sum(self.vals[hash_index, :, k_size_loc])  # get the number of True's

		def value_at_kmer(self, hash_index, kmer_index, k_size_loc):
			return self.vals[hash_index, kmer_index, k_size_loc]

		def process_seq(self, seq):
			for k_size_loc in range(len(k_range)):  # could do this more efficiently by putting this in the inner loop
				ksize = k_range[k_size_loc]
				for i in range(len(seq) - ksize + 1):
					kmer = seq[i:i + ksize]
					if kmer in all_kmers_set:
						tuples = tree.prefix_match(kmer, None, tst.TupleListAction())  # get all the k-mers whose prefix matches
						all_kmers_set.remove(kmer)  # Drop from the pre-filter list, since it will subsequently not be used
						hash_loc_kmer_loc = []
						# get the location of the found kmers in the counters
						for item in tuples:
							split_string = item[0].split('x')
							hash_loc_kmer_loc.append((int(split_string[1]), int(
								split_string[2])))  # first is the hash location, second is which k-mer
						if len(hash_loc_kmer_loc) > 0:
							for hash_index, kmer_index in hash_loc_kmer_loc:
								self.increment(hash_index, kmer_index, k_size_loc)

	# helper function
	def q_func(queue, counter):
		# Worker function to process the reads in the queue
		while True:
			record = queue.get()
			if record is False:  # In case I need to pass a poison pill
				return
			else:
				counter.process_seq(record)

	# Initialize the counters
	counter = Counters()
	# Start the q
	queue = multiprocessing.Queue()
	ps = list()
	for i in range(num_threads):
		p = multiprocessing.Process(target=q_func, args=(queue, counter))
		p.daemon = True
		p.start()
		ps.append(p)

	# populate the queue
	fid = khmer.ReadParser(query_file)  # This is faster than screed
	i = 0
	for record in fid:
		seq = record.sequence
		queue.put(seq)
		i += 1
		if i % 100000 == 0:
			print("Read in %d sequences" % i)

	# Wait for everything to finish
	while True:
		if queue.empty():
			break
		else:
			print("Sequences left to process: %d" % queue.qsize())
			time.sleep(1)
	queue.close()
	queue.join_thread()

	# collect the intersection counts
	containment_indices = np.zeros((len(sketches), len(k_range)))
	for sketch_index in range(len(sketches)):
		for k_size_loc in range(len(k_range)):
			uniqe_kmer_indicies = set()
			uniqe_kmers = set()
			# get the indices of the unique k-mers
			for index in range(num_hashes):
				kmer = sketches[sketch_index]._kmers[index][0:k_range[k_size_loc]]  # get the first few characters of the kmer
				if kmer not in uniqe_kmers:
					uniqe_kmers.add(kmer)
					uniqe_kmer_indicies.add(index)
			# get the total counts for these unique k-mer indices
			count = 0
			for index in uniqe_kmer_indicies:
				count += counter.value_at_kmer(sketch_index, index, k_size_loc)
			# update the containment index
			containment_indices[sketch_index, k_size_loc] = count / float(len(uniqe_kmers))  # this is the containment index estimate
			#intersection_counts[sketch_index, k_size_loc] = counter.value(sketch_index, k_size_loc)  # this old way contained duplicate counts

	results = dict()
	for k_size_loc in range(len(k_range)):
		ksize = k_range[k_size_loc]
		key = 'k=%d' % ksize
		results[key] = containment_indices[:, k_size_loc]
	df = pd.DataFrame(results, map(os.path.basename, training_file_names))
	df = df.reindex_axis(['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
	max_key = 'k=%d' % k_range[-1]
	filtered_results = df[df[max_key] > coverage_threshold].sort_values(max_key, ascending=False)
	#filtered_results.reindex_axis(['Unnamed: 0'] + ['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
	filtered_results.to_csv(results_file, index=True, encoding='utf-8')

	# If requested, plot the results
	if args.plot_file:
		df = pd.read_csv(results_file)  # annoyingly, I have to read it back in to get the format into something I can work with
		dft = df.transpose()
		fig = plt.figure()
		for key in dft.keys():
			plt.plot(k_range, dft[key].values[1:])  # [1:] since the first entry is the name of the file
		plt.legend(dft.values[0])
		plt.xlabel('K-mer size')
		plt.ylabel('Containment Index')
		fig.savefig(plot_file)




# To check results using bash
#for sketch_index in range(len(sketches)):
#	if os.path.basename(sketches[sketch_index].input_file_name) == 'G000312105.fna':
#		print(sketch_index)
#		to_print_index = sketch_index

#with open('/home/dkoslicki/Dropbox/Repositories/CMash/CMash/data/kmers.txt','w') as fid:
#	for kmer in sketches[to_print_index]._kmers:
#		fid.write('%s\n' % kmer)



