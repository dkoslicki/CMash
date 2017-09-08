import khmer
import numpy as np
from khmer.khmer_args import optimal_size
import os
import MinHash as MH
import pandas as pd
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
from itertools import *

num_threads = multiprocessing.cpu_count()
query_file = os.path.abspath('../data/1Mil.fastq')
node_graph_out = query_file + ".NodeGraph"
training_data = os.path.abspath('../data/AllSketches.h5')
results_file = os.path.abspath('../data/results.csv')
force = False
ksize = 21
fprate = 0.0001
coverage_threshold = 0.02  # desired coverage cutoff
confidence = 0.95  # desired confidence that you got all the organisms with coverage >= desired coverage

if not os.path.exists(node_graph_out) or force is True:
	hll = khmer.HLLCounter(0.01)
	hll.consume_seqfile(query_file)
	res = optimal_size(hll.estimate_cardinality(), fp_rate=fprate)
	sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
	sample_kmers.consume_seqfile(query_file)
	# Save the sample_kmers
	sample_kmers.save(node_graph_out)
else:
	sample_kmers = khmer.load_nodegraph(node_graph_out)
	node_ksize = sample_kmers.ksize()
	if node_ksize != ksize:
		raise Exception("Node graph %s has wrong k-mer size of %d (input was %d). Try --force or change -k." % (node_graph_out, node_ksize, ksize))

num_sample_kmers = sample_kmers.n_unique_kmers()

# Import all the training data
sketches = MH.import_multiple_from_single_hdf5(training_data)
# Check for issues with the sketches (can also check if all the kmers make sense (i.e. no '' or non-ACTG characters))
num_hashes = len(sketches[0]._kmers)
for i in range(len(sketches)):
	sketch = sketches[i]
	if len(sketch._kmers) != num_hashes:
		raise Exception("Unequal number of hashes for sketch of %s" % sketch.input_file_name)

training_file_names = []
for i in range(len(sketches)):
	training_file_names.append(sketches[i].input_file_name)

# Helper function that uses equations (2.1) and (2.7) that tells you where you need
# to set the threshold to ensure (with confidence 1-t) that you got all organisms
# with coverage >= c
def threshold_calc(k, c, p, confidence):
	delta = c*(1-np.sqrt(-2*np.log(1 - confidence)*(c+p) / float(c**2 * k)))
	if delta < 0:
		delta = 0
	return delta

def compute_indicies(sketch, num_sample_kmers):
	num_hashes = len(sketch._kmers)
	num_training_kmers = sketch._true_num_kmers
	count = 0
	adjust = 0
	for kmer in sketch._kmers:
		if kmer != '':
			count += sample_kmers.get(kmer)
		else:
			adjust += 1
	intersection_cardinality = count  # Should adjust this by the FP rate of the NodGraph, but I don't think I can get that...
	containment_index = count / float(num_hashes - adjust)
	jaccard_index = num_training_kmers * containment_index / float(num_training_kmers + num_sample_kmers - num_training_kmers * containment_index)
	return intersection_cardinality, containment_index, jaccard_index


def unwrap_compute_indicies(arg):
	return compute_indicies(*arg)


pool = Pool(processes=num_threads)
res = pool.map(unwrap_compute_indicies, zip(sketches, repeat(num_sample_kmers)))

# Gather up the results in a nice form
intersection_cardinalities = np.zeros(len(sketches))
containment_indexes = np.zeros(len(sketches))
jaccard_indexes = np.zeros(len(sketches))
for i in range(len(res)):
	(intersection_cardinality, containment_index, jaccard_index) = res[i]
	intersection_cardinalities[i] = intersection_cardinality
	containment_indexes[i] = containment_index
	jaccard_indexes[i] = jaccard_index

d = {'intersection': intersection_cardinalities,
	 'containment index': containment_indexes,
	 'jaccard index': jaccard_indexes}

df = pd.DataFrame(d, index=map(os.path.basename, training_file_names))

# Only get the rows above a certain threshold
if coverage_threshold <= 0:
	est_threshold = 0
else:
	est_threshold = threshold_calc(num_hashes, coverage_threshold, fprate, confidence)
filtered_results = df[df['containment index'] > est_threshold].sort_values('containment index', ascending=False)
# Export the results
filtered_results.to_csv(results_file, index=True, encoding='utf-8')
