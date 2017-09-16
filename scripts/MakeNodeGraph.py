#! /usr/bin/env python
# This script will create node graph for a given k-mer size and query file (can be used as input to QueryDNADatabase.py)
import khmer
from khmer.khmer_args import optimal_size
import os
import argparse
import screed
import threading
import multiprocessing

def restricted_float(x):
	x = float(x)
	if x < 0.0 or x >= 1.0:
		raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0)" % (x,))
	return x

def main():
	parser = argparse.ArgumentParser(
		description="This script will create node graph for a given k-mer size and query file (can be used as input to QueryDNADatabase.py)",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-fp', '--fp_rate', type=restricted_float, help="False positive rate.", default=0.0001)
	parser.add_argument('-i', '--intersect_nodegraph',
						help="Location of Node Graph. Will only insert query k-mers in bloom filter if they appear anywhere in the training"
							 " database. Note that the Jaccard estimates will now be "
							 "J(query intersect union_i training_i, training_i) instead of J(query, training_i), "
							 "but will use significantly less space (unfortunately will also disable threading).")
	parser.add_argument('-k', '--k_size', type=int, help="K-mer size", default=21)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
	parser.add_argument('in_file', help="Input file: FASTQ/A file (can be gzipped).")
	parser.add_argument('out_dir', help='Output directory')

	# Parse and check args
	args = parser.parse_args()
	query_file = os.path.abspath(args.in_file)
	ksize = args.k_size
	num_threads = args.threads
	node_graph_out = os.path.join(os.path.abspath(args.out_dir), os.path.basename(query_file) + ".NodeGraph.K" + str(ksize))
	if args.intersect_nodegraph is not None:
		intersect_nodegraph_file = args.intersect_nodegraph
	else:
		intersect_nodegraph_file = None
	intersect_nodegraph = None
	if intersect_nodegraph_file is not None:
		if not os.path.exists(intersect_nodegraph_file):
			raise Exception("Intersection nodegraph does not exist. Please re-run MakeDNADatabase.py with the -i flag.")
		try:
			intersect_nodegraph = khmer.load_nodegraph(intersect_nodegraph_file)
			if intersect_nodegraph.ksize() != ksize:
				raise Exception("Given intersect nodegraph %s has K-mer size %d while the database K-mer size is %d"
								% (intersect_nodegraph_file, intersect_nodegraph.ksize(), ksize))
		except:
			raise Exception("Could not load given intersect nodegraph %s" % intersect_nodegraph_file)
	fprate = args.fp_rate
	hll = khmer.HLLCounter(0.01, ksize)
	hll.consume_seqfile(query_file)
	full_kmer_count_estimate = hll.estimate_cardinality()
	res = optimal_size(full_kmer_count_estimate, fp_rate=fprate)
	if intersect_nodegraph is None:  # If no intersect list was given, just populate the bloom filter
		sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
		#sample_kmers.consume_seqfile(query_file)
		rparser = khmer.ReadParser(query_file)
		threads = []
		for _ in range(num_threads):
			cur_thrd = threading.Thread(target=sample_kmers.consume_seqfile_with_reads_parser, args=(rparser,))
			threads.append(cur_thrd)
			cur_thrd.start()
		for thread in threads:
			thread.join()
	else:  # Otherwise, only put a k-mer in the bloom filter if it's in the intersect list
		# (WARNING: this will cause the Jaccard index to be calculated in terms of J(query\intersect hash_list, training)
		#  instead of J(query, training)
		# (TODO: fix this after khmer is updated)
		#intersect_nodegraph_kmer_count = intersect_nodegraph.n_unique_kmers()  # Doesnt work due to khmer bug
		intersect_nodegraph_kmer_count = intersect_nodegraph.n_occupied()  # Doesnt work due to khmer bug
		if intersect_nodegraph_kmer_count < full_kmer_count_estimate:  # At max, we have as many k-mers as in the union of the training database (But makes this always return 0)
			res = optimal_size(intersect_nodegraph_kmer_count, fp_rate=fprate)
			sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
		else:
			sample_kmers = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
		for record in screed.open(query_file):
			seq = record.sequence
			for i in range(len(seq) - ksize + 1):
				kmer = seq[i:i + ksize]
				if intersect_nodegraph.get(kmer) > 0:
					sample_kmers.add(kmer)
	# Save the sample_kmers
	sample_kmers.save(node_graph_out)

if __name__ == "__main__":
	main()
