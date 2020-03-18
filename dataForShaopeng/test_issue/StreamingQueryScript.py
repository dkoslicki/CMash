#! /usr/bin/env python
import khmer
import marisa_trie as mt
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
import multiprocessing
import pandas as pd
import argparse
from argparse import ArgumentTypeError
import re
import matplotlib.pyplot as plt
from hydra import WritingBloomFilter, ReadingBloomFilter
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse import save_npz
from scipy.io import savemat
import timeit
from itertools import islice

t00 = timeit.default_timer()

# TODO: export hit matrices


def parseNumList(input):
    """Thank you stack overflow"""
    m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', input)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise ArgumentTypeError(
            "'" + input + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
    start = int(m.group(1))
    end = int(m.group(2))
    if m.group(3):
        increment = int(m.group(3))
    else:
        increment = 1
    return list(range(start, end + 1, increment))


# read in the arguments
my_range = "50-61-1"
k_range = parseNumList(my_range)
if k_range is None:
    raise Exception("The --range argument is required, no matter what the help menu says.")
training_data = "/home/dkoslicki/Desktop/CMash/dataForShaopeng/test_issue/TrainingDatabase_k_61.h5"
query_file = "/home/dkoslicki/Desktop/CMash/dataForShaopeng/small_data/taxid_1909294_104_genomic.fna.gz"
results_file = f"/home/dkoslicki/Desktop/CMash/dataForShaopeng/test_issue/out_{my_range.replace('-','_')}.csv"
npz_file = os.path.splitext(results_file)[0] + "_hit_matrix.npz"
num_threads = 12
location_of_thresh = -1
coverage_threshold = 0
streaming_database_file = os.path.splitext(training_data)[0] + ".tst"  # name of the tst training file
streaming_database_file = os.path.abspath(streaming_database_file)
hydra_file = ""
verbose = True
num_reads_per_core = 100000
sensitive = True
if not os.path.exists(streaming_database_file):
    streaming_database_file = None

# Import data and error checking
# Query file
if not os.path.exists(query_file):
    raise Exception("Query file %s does not exist." % query_file)
if not os.path.exists(training_data):
    raise Exception("Training/reference file %s does not exist." % training_data)

# Training data
if verbose:
    print("Reading in sketches")
    t0 = timeit.default_timer()
sketches = MH.import_multiple_from_single_hdf5(training_data)
if sketches[0]._kmers is None:
    raise Exception(
        "For some reason, the k-mers were not saved when the database was created. Try running MakeStreamingDNADatabase.py again.")
num_hashes = len(
    sketches[0]._kmers)  # note: this is relying on the fact that the sketches were properly constructed
max_ksize = sketches[0].ksize


def keyfunction(item):
    return os.path.basename(item.input_file_name)


sketches = sorted(sketches, key=keyfunction)  # sort the sketches by the basename of input file

# adjust the k-range if necessary
k_range = [val for val in k_range if val <= max_ksize]

# adjust location of thresh if necessary
if location_of_thresh:
    if location_of_thresh >= len(k_range):
        print("Warning, k_range is of length %d, reducing location of threshold from %d to %d" % (
        len(k_range), location_of_thresh, len(k_range)))
        location_of_thresh = len(k_range) - 1

# Get names of training files for use as rows in returned tabular data
training_file_names = []
for i in range(len(sketches)):
    training_file_names.append(str(sketches[i].input_file_name.decode('utf-8')))

training_file_names = sorted(training_file_names, key=os.path.basename)  # sort based on base name

if verbose:
    print("Finished reading in sketches")
    t1 = timeit.default_timer()
    print("Time: %f" % (t1 - t0))

if verbose:
    print("Reading in/creating ternary search tree")
    t0 = timeit.default_timer()
# Make the Marissa tree
if streaming_database_file is None:
    streaming_database_file = os.path.splitext(training_data)[0] + ".tst"
    streaming_database_file = os.path.abspath(streaming_database_file)
    print(
        "It appears a tst training file has not been created (did you remember to use MakeStreamingDNADatabase.py?).")
    print("I'm creating one anyway at: %s" % streaming_database_file)
    print("This may take a while...")
    to_insert = set()
    for i in range(len(sketches)):
        for kmer_index in range(len(sketches[i]._kmers)):
            kmer = sketches[i]._kmers[kmer_index]
            to_insert.add(
                kmer + 'x' + str(i) + 'x' + str(kmer_index))  # format here is kmer+x+hash_index+kmer_index
    tree = mt.Trie(to_insert)
    tree.save(streaming_database_file)
else:
    tree = mt.Trie()
    tree.load(streaming_database_file)

# all the k-mers of interest in a set (as a pre-filter)
if not hydra_file:  # create one
    try:
        all_kmers_bf = WritingBloomFilter(len(sketches) * len(k_range) * num_hashes * 2, 0.01)
        for sketch in sketches:
            for kmer in sketch._kmers:
                for ksize in k_range:
                    all_kmers_bf.add(kmer[0:ksize])  # put all the k-mers and the appropriate suffixes in
                    # all_kmers_bf.add(khmer.reverse_complement(kmer[0:ksize]))  # also add the reverse complement
                    all_kmers_bf.add(khmer.reverse_complement(kmer[0:ksize]))  # also add the reverse complement
    except IOError:
        print("No such file or directory/error opening file: %s" % hydra_file)
        sys.exit(1)
else:  # otherwise read it in
    try:
        all_kmers_bf = ReadingBloomFilter(hydra_file)
    except IOError:
        print("No such file or directory/error opening file: %s" % hydra_file)
        sys.exit(1)
if verbose:
    print("Finished reading in/creating ternary search tree")
    t1 = timeit.default_timer()
    print("Time: %f" % (t1 - t0))
# Seen k-mers (set of k-mers that already hit the trie, so don't need to check again)
seen_kmers = set()


# shared object that will update the intersection counts
class Counters(object):
    # This class is basically an array of counters (on the same basis as the sketches)
    # it's used to keep track (in a parallel friendly way) of which streamed k-mers went into the training file sketches
    def __init__(self):
        pass

    def return_matches(self, input_kmer, k_size_loc):
        """ Get all the matches in the trie with the kmer prefix"""
        match_info = set()
        to_return = []
        for kmer in [input_kmer, khmer.reverse_complement(input_kmer)]:  # FIXME: might need to break if one of them matches
            # for kmer in [input_kmer]:
            #if kmer not in all_kmers_bf:  # TODO: but for some reason it works here
            #    return [], False
            prefix_matches = tree.keys(kmer)  # get all the k-mers whose prefix matches
            #if not kmer in all_kmers_bf:
            #    return [], False
            #if prefix_matches and not kmer in all_kmers_bf:
            #    print(f"prefix: {prefix_matches}, is in BF: {kmer in all_kmers_bf}")
            # match_info = set()
            # get the location of the found kmers in the counters
            for item in prefix_matches:
                split_string = item.split('x')  # first is the hash location, second is which k-mer
                hash_loc = int(split_string[1])
                kmer_loc = int(split_string[2])
                match_info.add((hash_loc, k_size_loc, kmer_loc))
            # to_return = []
            saw_match = False
            if match_info:
                saw_match = True
                for tup in match_info:
                    to_return.append(tup)
            if saw_match:  # Only need to see a match to the original kmer or the reverse complement, don't return both otherwise you over-count
                break
        return to_return, saw_match

    def process_seq(self, seq):
        #  start with small kmer size, if see match, then continue looking for longer k-mer sizes, otherwise move on
        small_k_size = k_range[0]  # start with the small k-size
        to_return = []
        for i in range(len(seq) - small_k_size + 1):  # look at all k-mers
            kmer = seq[i:i + small_k_size]
            possible_match = False
            #if kmer not in seen_kmers:  # if we should process it
            if True:
                #if kmer in all_kmers_bf:  # if we should process it  # FIXME: problem might be here since if I remove the bloom filter, everything appears to work just fine....
                if True:
                    match_list, saw_match = self.return_matches(kmer, 0)
                    if saw_match:  # TODO: note, I *could* add all the trie matches and their sub-kmers to the seen_kmers
                        seen_kmers.add(kmer)  # FIXME: might also be able to add the reverse complements in here, instead of adjusting the division down near line 332
                        seen_kmers.add(khmer.reverse_complement(kmer))
                        to_return.extend(match_list)
                    possible_match = True
            # TODO: note: I could (since it'd only be for a single kmer size, keep a set of *all* small_kmers I've tried and use this as another pre-filter
            else:
                possible_match = True
            # start looking at the other k_sizes, don't overhang len(seq)
            if possible_match:
                #prev_match = True
                for other_k_size in [x for x in k_range[1:] if i + x <= len(seq)]:
                    kmer = seq[i:i + other_k_size]
                    #if kmer in all_kmers_bf:
                    #print("True")
                    #if kmer not in all_kmers_bf and tree.keys(kmer):
                    #    print(f"{kmer}\n")
                    if True:
                        k_size_loc = k_range.index(other_k_size)
                        match_list, saw_match = self.return_matches(kmer, k_size_loc)
                        if saw_match:
                            to_return.extend(match_list)
                        #if saw_match and kmer not in all_kmers_bf:  # TODO: this doesn't return a single "bad thing happened", yet I can't use `if kmer in all_kmers_bf` at the top
                        #    print("bad things happened")
                        #else:
                        #    break  # TODO: Interesting! If you break after not seeing a match, it goes back to the bad results. So something is wrong with the logic of the for loop
                        # TODO: What I would like to do is see what the kmers look like when you don't see a match for a smaller k-mer size, but somehow do see a match for the larger kmer size
                        #if saw_match and not prev_match:
                            #print(f"{kmer}\n")
                        #else:
                        #    prev_match = saw_match
                        # FIXME: ok, so here's the problem: longer kmer matches while shorter does not: for kmer='CATGTCTTTCAGGCTGGAACCGGAGGCGACCAATGCCGACCGGCTGCTTTAT', tree.keys(kmer)==[] and tree.keys(khmer.reverse_complement(kmer))==match, yet tree.keys(khmer.reverse_complement(kmer[0:-1]))==[] and tree.keys(kmer[0:-1])==[]
                        # FIXME: but curiously, tree.keys(khmer.reverse_complement(kmer)[0:-1])==match, so it depends on when the reverse complement is being applied!!!! <<<<<-------------!!!!!!!!!!!!!
                        # FIXME: Additionally, kmer in all_kmers_bf == TRUE yet, kmer[0:-1] in all_kmers_bf == FALSE
                        # FIXME: so the problem is, that reverse complement is REVERSE complement, you need to FLIP THE KMER AROUND if you are doing prefix lookups in the tree. MAYBE?!
                        # MAYBE: instead of just passing the kmer, and then returning_matches using the kmer and it's reverse complement, I should pass the kmer,

                        # FIXME: I *might* be able to circumvent all of this just by using canonical k-mers....
                    else:
                        break

        return to_return


# Initialize the counters
# TODO: note, I could be doing a partial dedup here, just to reduce the memory usage...
counter = Counters()


def map_func(sequence):
    return counter.process_seq(sequence)


pool = multiprocessing.Pool(processes=num_threads)

if verbose:
    print("Start streaming")
    t0 = timeit.default_timer()
# populate the queue
fid = khmer.ReadParser(query_file)  # This is faster than screed
match_tuples = []
# num_reads_per_core = 100000
num_reads_per_chunk = num_reads_per_core * num_threads
to_proc = [record.sequence for record in islice(fid, num_reads_per_chunk)]
i = 0
while to_proc:
    i += len(to_proc)
    if verbose:
        print("Read in %d sequences" % i)
    res = pool.map(map_func, to_proc, chunksize=int(max(1, min(num_reads_per_core, len(to_proc) / num_threads))))
    #res = map(map_func, to_proc)
    flattened_res = [item for sublist in res if sublist for item in sublist]
    flattened_res = list(set(flattened_res))  # dedup it
    match_tuples.extend(flattened_res)
    to_proc = [record.sequence for record in islice(fid, num_reads_per_chunk)]
fid.close()
pool.close()
# print(match_tuples)
if verbose:
    print("Finished streaming")
    t1 = timeit.default_timer()
    print("Time: %f" % (t1 - t0))

if verbose:
    print("Forming hit matrix")
    t0 = timeit.default_timer()
# print("Len matches: %d" % len(match_tuples))
# create k_range spare matrices. Rows index by genomes (sketch/hash index), columns index by k_mer_loc
row_ind_dict = dict()
col_ind_dict = dict()
value_dict = dict()
unique_kmers = dict()  # this will keep track of the unique k-mers seen in each genome (sketch/hash loc)
for k_size in k_range:
    row_ind_dict[k_size] = []
    col_ind_dict[k_size] = []
    value_dict[k_size] = []

match_tuples = set(match_tuples)  # uniquify, so we don't make the row/col ind dicts too large

for hash_loc, k_size_loc, kmer_loc in match_tuples:
    if hash_loc not in unique_kmers:
        unique_kmers[hash_loc] = set()
    k_size = k_range[k_size_loc]
    kmer = sketches[hash_loc]._kmers[kmer_loc][:k_size]
    if kmer not in unique_kmers[hash_loc]:  # if you've seen this k-mer before, don't add it. NOTE: this makes sure we don't over count
        row_ind_dict[k_size].append(hash_loc)
        col_ind_dict[k_size].append(kmer_loc)
        value_dict[k_size].append(1)
        unique_kmers[hash_loc].add(kmer)

hit_matrices = []

for k_size in k_range:
    mat = csc_matrix((value_dict[k_size], (row_ind_dict[k_size], col_ind_dict[k_size])),
                     shape=(len(sketches), num_hashes))
    hit_matrices.append(mat)
if verbose:
    print("Finished forming hit matrix")
    t1 = timeit.default_timer()
    print("Time: %f" % (t1 - t0))

if verbose:
    print("Computing containment indicies")
    t0 = timeit.default_timer()
containment_indices = np.zeros((len(sketches), len(k_range)))  # TODO: could make this thing sparse, or do the filtering for above threshold here
for k_size_loc in range(len(k_range)):
    containment_indices[:, k_size_loc] = (hit_matrices[k_size_loc].sum(axis=1).ravel())  # /float(num_hashes))

for k_size_loc in range(len(k_range)):
    k_size = k_range[k_size_loc]
    for hash_loc in np.where(containment_indices[:, k_size_loc])[0]:  # find the genomes with non-zero containment
        unique_kmers = set()
        for kmer in sketches[hash_loc]._kmers:  # FIXME: problem is right here since I am not counting the revcomps, though the Counters() *do* use revcomps
            unique_kmers.add(kmer[:k_size])  # find the unique k-mers
        # unique_kmers.add(khmer.reverse_complement(kmer[:k_size]))  # also add the rev-comps since I streaming queried them
        containment_indices[hash_loc, k_size_loc] /= float(len(unique_kmers))  # divide by the unique num of k-mers
if verbose:
    print("Finished computing containment indicies")
    t1 = timeit.default_timer()
    print("Time: %f" % (t1 - t0))

results = dict()
for k_size_loc in range(len(k_range)):
    ksize = k_range[k_size_loc]
    key = 'k=%d' % ksize
    results[key] = containment_indices[:, k_size_loc]
df = pd.DataFrame(results, map(os.path.basename, training_file_names))
df = df.reindex(labels=['k=' + str(k_size) for k_size in k_range], axis=1)  # sort columns in ascending order
sort_key = 'k=%d' % k_range[location_of_thresh]
max_key = 'k=%d' % k_range[-1]
filtered_results = df[df[sort_key] > coverage_threshold].sort_values(max_key,
                                                                     ascending=False)  # only select those where the highest k-mer size's count is above the threshold

if sensitive:
    if verbose:
        print("Exporting results")
        t0 = timeit.default_timer()
    filtered_results.to_csv(results_file, index=True, encoding='utf-8')

t11 = timeit.default_timer()
print(f"Total time: {t11 - t00}")