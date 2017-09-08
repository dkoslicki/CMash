"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import khmer
import screed
import h5py
import numpy as np
import os
import tempfile
import multiprocessing
from multiprocessing import Pool
import re
from itertools import *
import collections
from blist import *  # note, the import functions import the _mins etc. as lists, and the CE class imports them as blists.
# This shouldn't cause an issue, but will lead to slow performance if a CE is imported, then additional things are added.
# I.e. If you import a CE, don't add new elements, or you might have a bad day (or at least a long one).
import bisect
import scipy.optimize
import ctypes
import warnings
import subprocess
import filecmp
import shutil
import traceback
import random
warnings.simplefilter("ignore", RuntimeWarning)

# To Do:
# Implement hash_murmur3 to leave the khmer package. Need to implement reverse complement myself, etc.
# After that point, can use amino acids
# SNAP paired or single reads
# Get DIAMOND implemented
# Make the snap streaming automatically chunk the index_dirs if there are too many (can get max command len with xargs --show-limits)
# Make the count vector over a shared array (just like the kmer matricies)
# Use sam tools to partition the reads into aligned and unaligned (be careful with mate pairs)
# Implement paired reads for snap-aligner



notACTG = re.compile('[^ACTG]')


def unwrap_count_vector(arg):
    """
    Helper function for parallelizing the count_vector
    :param arg:
    :param kwarg:
    :return:
    """
    return arg[0].jaccard_count(arg[1])


def unwrap_jaccard_vector(arg):
    """
    Helper function for parallelizing the jaccard_vector
    :param arg:
    :param kwarg:
    :return:
    """
    return arg[0].jaccard(arg[1])
    #return arg[1].jaccard(arg[0]) #Would like to test effect of using the other denominators


class CountEstimator(object):
    """
    A simple bottom n-sketch MinHash implementation.
    n is the number of sketches to keep
    Still don't know what max_prime is...
    """

    def __init__(self, n=None, max_prime=9999999999971., ksize=None, input_file_name=None, save_kmers='n', hash_list=None,
                 rev_comp=False):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception

        #if ksize % 2 == 0:
        #    raise Exception("Due to an issue with khmer, only odd ksizes are allowed")

        self.ksize = ksize
        self.hash_list = hash_list

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize sketch to size n
        #self._mins = [float("inf")]*n
        self._mins = blist([p]*n)

        # initialize the corresponding counts
        self._counts = blist([0]*n)

        # initialize the list of kmers used, if appropriate
        if save_kmers == 'y':
            self._kmers = blist(['']*n)
        else:
            self._kmers = None

        # Initialize file name (if appropriate)
        self.input_file_name = input_file_name
        if self.input_file_name:
            self.parse_file(rev_comp=rev_comp)

        # Optional container for the true number of k-mers in the genome used to populate the sketch
        self._true_num_kmers = 0

    def parse_file(self, rev_comp=False):
        """
        opens a file and populates the CountEstimator with it
        """
        for record in screed.open(self.input_file_name):
            self.add_sequence(record.sequence, rev_comp)

    def down_sample(self, h):
        """
        This will down-sample a sketch to have exactly h elements
        :param h: number of elements you wish to save
        :return: None
        """
        self._mins = self._mins[0:h]
        self._counts = self._counts[0:h]
        self._kmers = self._kmers[0:h]

    def add(self, kmer, rev_comp=False):
        """
        Add kmer into sketch, keeping sketch sorted, update counts accordingly
        """
        _mins = self._mins
        _counts = self._counts
        _kmers = self._kmers

        if rev_comp:
            h1 = khmer.hash_murmur3(kmer)
            h2 = khmer.hash_murmur3(khmer.reverse_complement(kmer))
            #h1 = hash(kmer)
            #h2 = hash(khmer.reverse_complement(kmer))
            h = min(h1, h2)
            if h == h2:
                kmer = khmer.reverse_complement(kmer)
        else:
            h = khmer.hash_murmur3(kmer)
            #h = hash(kmer)

        h = h % self.p
        if self.hash_list:  # If I only want to include hashes that occur in hash_list
            if h not in self.hash_list:  # If the kmer isn't in the hash_list, then break
                return

        if h >= _mins[-1]:
            return

        i = bisect.bisect_left(_mins, h)  # find index to insert h
        if _mins[i] == h:  # if h in mins, increment counts
            _counts[i] += 1
            return
        else:  # otherwise insert h, initialize counts to 1, and insert kmer if necessary
            _mins.insert(i, h)
            _mins.pop()
            _counts.insert(i, 1)
            _counts.pop()
            if _kmers:
                _kmers.insert(i, np.string_(kmer))
                _kmers.pop()
            return

        assert 0, "should never reach this"

    def add_sequence(self, seq, rev_comp=False):
        """
         Sanitize and add a sequence to the sketch.
        """
        # seq = seq.upper().replace('N', 'G')
        seq = notACTG.sub('G', seq.upper())  # more intelligent sanatization?
        # seq = seq.upper()
        for kmer in kmers(seq, self.ksize):
            #if notACTG.search(kmer) is None:  # If it's free of non-ACTG characters
            self.add(kmer, rev_comp)

    def jaccard_count(self, other):
        """
        Jaccard index weighted by counts
        """
        truelen = len(self._mins)
        while truelen and self._mins[truelen - 1] == self.p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
            truelen -= 1
        if truelen == 0:
            raise ValueError

        (total1, total2) = self.common_count(other)
        return (total2 / float(sum(other._counts)), total1 / float(sum(self._counts)))
        # The entries here are returned as (A_{CE1,CE2}, A_{CE2,CE1})

    def jaccard(self, other):
        """
        Jaccard index
        """
        truelen = len(self._mins)
        while truelen and self._mins[truelen - 1] == self.p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
            truelen -= 1
        if truelen == 0:
            raise ValueError

        return self.common(other) / float(truelen)
    #similarity = jaccard_count

    def common_count(self, other):
        """
        Calculate number of common k-mers between two sketches, weighted by their counts
        """
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")
        if self.p != other.p:
            raise Exception("different primes - cannot compare")

        common1 = 0
        common2 = 0
        for (count1, count2) in _yield_count_overlaps(self._mins, other._mins, self._counts, other._counts):
            common1 += count1  # The unpopulated hashes have count 0, so we don't have to worry about that here
            common2 += count2
        return (common1, common2)

    def common(self, other):
        """
        Calculate number of common k-mers between two sketches.
        """
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")
        if self.p != other.p:
            raise Exception("different primes - cannot compare")

        common = 0
        for val in _yield_overlaps(self._mins, other._mins):
            if val != self.p:  # Make sure not to include the un-populated hashes p
                common += 1
        return common

    def _truncate(self, n):
        self._mins = self._mins[:n]

    def export(self, export_file_name):
        """
        This function will export the CountEstimator using hdf5
        """
        fid = h5py.File(export_file_name, 'w')
        grp = fid.create_group("CountEstimator")
        mins_data = grp.create_dataset("mins", data=self._mins)
        counts_data = grp.create_dataset("counts", data=self._counts)
        if self._kmers:
            kmer_data = grp.create_dataset("kmers", data=self._kmers)

        grp.attrs['class'] = np.string_("CountEstimator")
        grp.attrs['filename'] = np.string_(self.input_file_name)
        grp.attrs['ksize'] = self.ksize
        grp.attrs['prime'] = self.p
        grp.attrs['true_num_kmers'] = self._true_num_kmers
        fid.close()

    def count_vector(self, other_list):
        """
        Function that returns the Y vector of MetaPalette. That is, the vector where Y[i] = Jaccard_count(self, other_CEs[i]
        :param other_list: a list of count estimator classes
        :return: a numpy vector with the same basis as other_list giving the jaccard_count of self with other[i]
        """
        Y = np.zeros(len(other_list))

        pool = Pool(processes=multiprocessing.cpu_count())
        Y_tuple = pool.map(unwrap_count_vector, zip([self] * len(other_list), other_list))
        pool.terminate()
        for i in range(len(other_list)):
            Y[i] = Y_tuple[i][1]  # Gotta make sure it's not [1] (one's the CKM vector, the other is the "coverage")

        return Y

    def jaccard_vector(self, other_list):
        """
        Function that returns the Y vector of Jaccard values. That is, the vector where Y[i] = Jaccard(self, other_CEs[i]
        :param other_list: a list of count estimator classes
        :return: a numpy vector with the same basis as other_list giving the jaccard of self with other[i]
        """
        pool = Pool(processes=multiprocessing.cpu_count())
        Y = np.array(pool.map(unwrap_jaccard_vector, zip([self]*len(other_list), other_list)))
        pool.terminate()

        return Y


def import_single_hdf5(file_name):
    """
    This function will read an HDF5 file and populate the CountEstimator class accordingly
    :param file_name: input file name of HDF5 file created by CountEstimator.export(file_name)
    :return: CountEstimator
    """
    fid = h5py.File(file_name, 'r')  # This automatically handles non-existent files for me
    grp = fid["CountEstimator"]
    file_name = grp.attrs['filename']
    ksize = grp.attrs['ksize']
    prime = grp.attrs['prime']
    true_num_kmers = grp.attrs['true_num_kmers']
    mins = grp["mins"][...]  # For some reason, slicing is WAY slower than using ... in this case.
    counts = grp["counts"][...]
    CE = CountEstimator(n=len(mins), max_prime=3, ksize=ksize)
    CE.p = prime
    CE._mins = mins
    CE._counts = counts
    CE._true_num_kmers = true_num_kmers
    CE.input_file_name = file_name
    if "kmers" in grp:
        CE._kmers = grp["kmers"][...]
    else:
        CE._kmers = None

    fid.close()
    return CE


def import_multiple_hdf5(input_files_list):
    """
    Import a bunch of HDF5 Count Estimators from a given list of HDF5 files
    :param file_names: List of HDF5 file names of Count Estimators
    :return: list of Count Estimators
    """
    CEs = list()
    pool = Pool(processes=multiprocessing.cpu_count())
    CEs = pool.map(import_single_hdf5, input_files_list, chunksize=144)
    pool.terminate()

    return CEs


def export_multiple_hdf5(CEs, out_folder):
    """
    Exports a list of Count Estimators to a bunch of HDF5 files in a certain folder
    :param CEs: a list of Count Estimators
    :return: None
    """
    for CE in CEs:
        if CE.input_file_name == None:
            raise Exception("This function only works when count estimator were formed from files (i.e. CE.input_filename != None")

    for CE in CEs:
        CE.export(os.path.join(out_folder, os.path.basename(CE.input_file_name) + ".CE.h5"))

    return


def export_multiple_to_single_hdf5(CEs, export_file_name):
    """
    This will take a list of count estimators and export them to a single, large HDF5 file
    :param CEs: list of Count Estimators
    :param file_name: output file name
    :return: None
    """
    fid = h5py.File(export_file_name, 'w')
    grp = fid.create_group("CountEstimators")
    for CE in CEs:
        try:
            subgrp = grp.create_group(os.path.basename(CE.input_file_name))  # the key of a subgroup is the basename (not the whole file)
            mins_data = subgrp.create_dataset("mins", data=CE._mins)
            counts_data = subgrp.create_dataset("counts", data=CE._counts)
            if CE._kmers:
                kmer_data = subgrp.create_dataset("kmers", data=CE._kmers)

            subgrp.attrs['class'] = np.string_("CountEstimator")
            subgrp.attrs['filename'] = np.string_(CE.input_file_name)  # But keep the full file name on hand
            subgrp.attrs['ksize'] = CE.ksize
            subgrp.attrs['prime'] = CE.p
            subgrp.attrs['true_num_kmers'] = CE._true_num_kmers
        except ValueError:
            raise Exception("It appears that the training file name %s exists twice in the input data. Please make sure all names are unique (i.e. remove duplicates) and tray again." % CE.input_file_name)

    fid.close()


def import_multiple_from_single_hdf5(file_name, import_list=None):
    """
    This function will import multiple count estimators stored in a single HDF5 file.
    :param file_name: file name for the single HDF5 file
    :param import_list: List of names of files to import
    :return: a list of Count Estimators
    """
    CEs = list()
    fid = h5py.File(file_name, 'r')
    if "CountEstimators" not in fid:
        raise Exception("This function imports a single HDF5 file containing multiple sketches."
                        " It appears you've used it on a file containing a single sketch."
                        "Try using import_single_hdf5 instead")

    grp = fid["CountEstimators"]
    if import_list:
        iterator = [os.path.basename(item) for item in import_list]
    else:
        iterator = grp.keys()

    for key in iterator:
        if key not in grp:
            raise Exception("The key " + key + " is not in " + file_name)

        subgrp = grp[key]
        file_name = subgrp.attrs['filename']
        ksize = subgrp.attrs['ksize']
        prime = subgrp.attrs['prime']
        mins = subgrp["mins"][...]
        counts = subgrp["counts"][...]
        true_num_kmers = subgrp.attrs["true_num_kmers"]
        CE = CountEstimator(n=len(mins), max_prime=3, ksize=ksize)
        CE.p = prime
        CE._mins = mins
        CE._counts = counts
        CE._true_num_kmers = true_num_kmers
        CE.input_file_name = file_name
        if "kmers" in subgrp:
            CE._kmers = subgrp["kmers"][...]
        else:
            CE._kmers = None

        CEs.append(CE)

    fid.close()
    return(CEs)


class CE_map(object):
    """
    Helper function for mapping CountEstimator class over multiple input_file arguments
    """
    def __init__(self, n, max_prime, ksize, save_kmers):
        self.n = n
        self.max_prime = max_prime
        self.ksize = ksize
        self.save_kmers = save_kmers

    def __call__(self, input_file):
        return CountEstimator(n=self.n, max_prime=self.max_prime, ksize=self.ksize, input_file_name=input_file, save_kmers=self.save_kmers)


def compute_multiple(n=None, max_prime=9999999999971., ksize=None, input_files_list=None, save_kmers='n', num_threads=multiprocessing.cpu_count()):
    """
    Batch compute Count Estimators from a given list of file names.
    :param n: number of hashes to keep
    :param max_prime:
    :param ksize: kmer size to use
    :param input_files_list: list of input genomes (fasta/fastq)
    :param save_kmers: flag if you want to save kmers or not ('y' or 'n')
    :return: a list of Count Estimators
    """
    if n is None:
        raise Exception
    if ksize is None:
        raise Exception
    if input_files_list is None:
        raise Exception

    pool = Pool(processes=num_threads)
    CEs = pool.map(CE_map(n, max_prime, ksize, save_kmers), input_files_list)
    pool.close()
    return CEs


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


def jaccard_count(ij):
    """
    Clone of jaccard_count from the count_estimator class, just so I can use shared memory arrays
    :param ij: a tuple of indicies to use in the global shared_mins and shared_counts
    :return: entries of the CKM matrix
    """
    mins1 = shared_mins[ij[0]]
    mins2 = shared_mins[ij[1]]
    counts1 = shared_counts[ij[0]]
    counts2 = shared_counts[ij[1]]
    truelen = len(mins1)
    while truelen and mins1[truelen - 1] == p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
        truelen -= 1
    if truelen == 0:
        raise ValueError
    common1 = 0
    common2 = 0
    for (count1, count2) in _yield_count_overlaps(mins1, mins2, counts1, counts2):
        common1 += count1  # The unpopulated hashes have count 0, so we don't have to worry about that here
        common2 += count2
    return (common2 / float(sum(counts2)), common1 / float(sum(counts1)))


def form_jaccard_count_matrix(all_CEs):
    """
    Forms the jaccard count kmer matrix when given a list of count estimators
    :param all_CEs: a list of count estimators
    :return: a numpy array of the jaccard count matrix
    """
    A = np.zeros((len(all_CEs), len(all_CEs)), dtype=np.float64)
    # Can precompute all the indicies
    indicies = []
    for i in range(len(all_CEs)):
        for j in range(len(all_CEs)):
            indicies.append((i, j))

    shared_mins_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._mins))
    global shared_mins
    shared_mins = np.ctypeslib.as_array(shared_mins_base.get_obj())
    shared_mins = shared_mins.reshape(len(all_CEs), len(all_CEs[0]._mins))
    shared_counts_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._counts))
    global shared_counts
    shared_counts = np.ctypeslib.as_array(shared_counts_base.get_obj())
    shared_counts = shared_counts.reshape(len(all_CEs), len(all_CEs[0]._counts))
    global p
    p = all_CEs[0].p
    for i in range(len(all_CEs)):
        shared_mins[i] = all_CEs[i]._mins
        shared_counts[i] = all_CEs[i]._counts

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    chunk_size = np.floor(len(indicies)/float(multiprocessing.cpu_count()))
    if chunk_size < 1:
        chunk_size = 1
    res = pool.imap(jaccard_count, indicies, chunksize=chunk_size)
    for (i, j), val in zip(indicies, res):
        A[i, j] = val[0]
        A[j, i] = val[1]

    pool.terminate()
    return A


def jaccard(ij):
    """
    Clone of jaccard_count from the count_estimator class, just so I can use shared memory arrays
    :param ij: a tuple of indicies to use in the global shared_mins and shared_counts
    :return: entries of the CKM matrix
    """
    mins1 = shared_mins[ij[0]]
    mins2 = shared_mins[ij[1]]
    truelen = len(mins1)
    while truelen and mins1[truelen - 1] == p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
        truelen -= 1
    if truelen == 0:
        raise ValueError

    common = 0
    for val in _yield_overlaps(mins1, mins2):
        if val != p:  # Make sure not to include the un-populated hashes p
            common += 1
    return common/float(truelen)


def form_jaccard_matrix(all_CEs):
    """
    Forms the jaccard count kmer matrix when given a list of count estimators
    :param all_CEs: a list of count estimators
    :return: a numpy array of the jaccard count matrix
    """
    A = np.zeros((len(all_CEs), len(all_CEs)), dtype=np.float64)
    # Can precompute all the indicies
    indicies = []
    for i in range(len(all_CEs)):
        for j in range(len(all_CEs)):
            indicies.append((i, j))

    shared_mins_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._mins))
    global shared_mins
    shared_mins = np.ctypeslib.as_array(shared_mins_base.get_obj())
    shared_mins = shared_mins.reshape(len(all_CEs), len(all_CEs[0]._mins))
    global p
    p = all_CEs[0].p
    for i in range(len(all_CEs)):
        shared_mins[i] = all_CEs[i]._mins

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    chunk_size = np.floor(len(indicies)/float(multiprocessing.cpu_count()))
    if chunk_size < 1:
        chunk_size = 1
    res = pool.imap(jaccard, indicies, chunksize=chunk_size)
    for (i, j), val in zip(indicies, res):
        A[i, j] = val
        A[j, i] = val

    pool.terminate()
    return A


def _yield_count_overlaps(mins1, mins2, counts1, counts2):
    """
    Return (\sum_{i \in indicies(mins1\cap min2)} counts1[i], \sum_{i \in indicies(mins1\cap min2)} counts2[i])
    """
    i = 0
    j = 0
    processed = 0
    try:
        while processed <= min(len(mins1), len(mins2)):
            while mins1[i] < mins2[j]:
                i += 1
                processed += 1
            while mins1[i] > mins2[j]:
                j += 1
                processed += 1
            if mins1[i] == mins2[j]:
                yield (counts1[i], counts2[j])
                i += 1
                j += 1
                processed += 1
    except IndexError:
        return


def _yield_overlaps(x1, x2):
    """yield common hash values while iterating over two sorted lists of hashes
    To properly compute the estimate, I need this to only process min(len(x1), len(x2)) elements
    Returns an iterable object
    """
    i = 0
    j = 0
    processed = 0
    try:
        while processed <= min(len(x1), len(x2)):
            while x1[i] < x2[j]:
                i += 1
                processed += 1
            while x1[i] > x2[j]:
                j += 1
                processed += 1
            if x1[i] == x2[j]:
                yield x1[i]
                i += 1
                j += 1
                processed += 1
    except IndexError:
        return


def kmers(seq, ksize):
    """yield all k-mers of len ksize from seq.
    Returns an iterable object
    """
    for i in range(len(seq) - ksize + 1):
        yield seq[i:i+ksize]


# taken from khmer 2.0; original author Jason Pell.
def is_prime(number):
    """Check if a number is prime."""
    if number < 2:
        return False
    if number == 2:
        return True
    if number % 2 == 0:
        return False
    for _ in range(3, int(number ** 0.5) + 1, 2):
        if number % _ == 0:
            return False
    return True


def get_prime_lt_x(target):
    """Backward-find a prime smaller than (or equal to) target.

    Step backwards until a prime number (other than 2) has been
    found.

    Arguments: target -- the number to step backwards from
    """
    if target == 1:
        return 1

    i = int(target)
    if i % 2 == 0:
        i -= 1
    while i > 0:
        if is_prime(i):
            return i
        i -= 2

    if i <= 0:
        raise RuntimeError("unable to find a prime number < %d" % (target))


def cluster_matrix(A_eps, A_indicies, taxonomy, cluster_eps=.01):
    """
    This function clusters the indicies of A_eps such that for a given cluster, there is another element in that cluster
    with similarity (based on A_eps) >= cluster_eps for another element in that same cluster. For two elements of
    distinct clusters, their similarity (based on A_eps) < cluster_eps.
    :param A_eps: The jaccard or jaccard_count matrix containing the similarities
    :param A_indicies: The basis of the matrix A_eps (in terms of all the CEs)
    :param cluster_eps: The similarity threshold to cluster on
    :return: (a list of sets of indicies defining the clusters, LCAs of the clusters)
    """
    #A_indicies_numerical = np.where(A_indicies == True)[0]
    A_indicies_numerical = A_indicies
    # initialize the clusters
    clusters = []
    for A_index in range(len(A_indicies_numerical)):
        # Find nearby elements
        nearby = set(np.where(A_eps[A_index, :] >= cluster_eps)[0]) | set(np.where(A_eps[:, A_index] >= cluster_eps)[0])
        in_flag = False
        in_counter = 0
        in_indicies = []
        for i in range(len(clusters)):
            if nearby & clusters[i]:
                clusters[i].update(nearby)  # add the nearby indicies to the cluster
                in_counter += 1  # keep track if the nearby elements belong to more than one of the previously formed clusters
                in_indicies.append(i)  # which clusters nearby shares elements with
                in_flag = True  # tells if it forms a new cluster
        if not in_flag:  # if new cluster, then append to clusters
            clusters.append(set(nearby))
        if in_counter > 1:  # If it belongs to more than one cluster, merge them together
            merged_cluster = set()
            for in_index in in_indicies[::-1]:
                merged_cluster.update(clusters[in_index])
                del clusters[in_index]  # delete the old clusters (now merged)
            clusters.append(merged_cluster)  # append the newly merged clusters
    clusters_full_indicies = []
    for cluster in clusters:
        cluster_full_indicies = set()
        for item in cluster:
            cluster_full_indicies.add(A_indicies_numerical[item])
        clusters_full_indicies.append(cluster_full_indicies)
    # Check to make sure the clustering didn't go wrong
    if sum([len(item) for item in clusters_full_indicies]) != len(A_indicies_numerical):  # Check the correct number of indicies
        raise Exception("For some reason, the total number of indicies in the clusters doesn't equal the number of indicies you started with")
    if set([item for subset in clusters_full_indicies for item in subset]) != set(A_indicies_numerical):  # Make sure no indicies were missed or added
        raise Exception("For some reason, the indicies in all the clusters doesn't match the indicies you started with")
    return clusters_full_indicies, cluster_LCAs(clusters_full_indicies, taxonomy)


def cluster_LCAs(clusters, taxonomy):
    """
    This function returns the lowest common ancestor in each one of the input clusters
    :param clusters: input clusters
    :param taxonomy: input taxonomy
    :return: a list with the ith element being the lowest common ancestor of the ith cluster
    """
    LCAs = []
    for cluster in clusters:
        found_LCA = False
        if len(cluster) == 1:
            LCAs.append(taxonomy[list(cluster)[0]].split()[2].split('|')[-1])
            found_LCA = True
            continue
        cluster_taxonomy = []
        for index in cluster:
            cluster_taxonomy.append(taxonomy[index])
        for rank in range(7, -1, -1):
            rank_names = []
            dummy_name = 0
            for tax_path in cluster_taxonomy:
                split_taxonomy = tax_path.split()[2].split('|')
                if len(split_taxonomy) < rank + 1:
                    rank_names.append(str(dummy_name))
                else:
                    rank_names.append(split_taxonomy[rank])
            if len(set(rank_names)) == 1 and "0" not in rank_names:
                LCAs.append(rank_names[0])
                found_LCA = True
                break
        if not found_LCA:
            LCAs.append('sk__-1_microorganism')  # In case they don't even have the kingdom in common
    return LCAs


def _write_single_cluster(tup):
    """
    Helper function. Writes a single fast file consisting of all the sequences in input_file_names
    :param tup: input tuple (out_dir, LCA, cluster_index, input_file_names)
    :return: the name of the created file
    """
    out_dir = tup[0]
    LCA = tup[1]
    cluster_index = tup[2]
    input_file_names = tup[3]
    out_file_name = os.path.join(out_dir, LCA + "_" + str(cluster_index) + "_" + ".fa")  # put the cluster index in the name in case there are shared LCAs
    out_file = open(out_file_name, 'w')
    i = 0
    for file_name in input_file_names:
        for record in screed.open(file_name):
            out_file.write(">" + LCA + "_" + str(i))
            out_file.write("\n")
            out_file.write(record.sequence)
            out_file.write("\n")
            i += 1
    out_file.close()
    return out_file_name

def make_cluster_fastas(out_dir, LCAs, clusters, CEs, threads=multiprocessing.cpu_count()):
    """
    This function will write a single fasta file for each of the clusters
    :param out_dir: the output directory in which to write the fasta files
    :param LCAs: the least common ancestors (from cluster_LCAs())
    :param clusters: the clusters (from cluster_matrix())
    :param CEs: The list of count estimators
    :param threads: number of threads to use
    :return: a list of files created (to be used in build_reference())
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    pool = multiprocessing.Pool(processes=threads)
    file_names = pool.map(_write_single_cluster, zip(repeat(out_dir), LCAs, range(len(LCAs)), [[CEs[i].input_file_name for i in cluster] for cluster in clusters]), chunksize=1)
    pool.close()
    #pool.terminate()
    #pool.join()
    return file_names


##########################################################################
# Tests

def test_jaccard_1():
    E1 = CountEstimator(n=0, ksize=21)
    E2 = CountEstimator(n=0, ksize=21)

    E1._mins = [1, 2, 3, 4, 5]
    E2._mins = [1, 2, 3, 4, 6]
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/5.0


def test_jaccard_2_difflen():
    E1 = CountEstimator(n=0, ksize=21)
    E2 = CountEstimator(n=0, ksize=21)

    E1._mins = [1, 2, 3, 4, 5]
    E2._mins = [1, 2, 3, 4]
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/4.0


def test_yield_overlaps():
    x1 = [1, 3, 5]
    x2 = [2, 4, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 0


def test_yield_overlaps_2():
    x1 = [1, 3, 5]
    x2 = [1, 2, 4, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 1
    assert len(list(_yield_overlaps(x2, x1))) == 1


def test_yield_overlaps_3():
    x1 = [1, 3, 6]
    x2 = [1, 2, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 2
    assert len(list(_yield_overlaps(x2, x1))) == 2


def test_CountEstimator():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=1)
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=1)
    sequence1 = "AAAAAAAA"  # 100% of the 1-mers of this sequence show up in the other
    sequence2 = "AAAACCCCCCCC"  # 4/12ths of the 1-mers in this sequence show up in the other
    CE1.add_sequence(sequence1)
    CE2.add_sequence(sequence2)
    assert CE1.jaccard_count(CE2) == (4/12., 1.0)
    assert CE2.jaccard_count(CE1) == (1.0, 4/12.)
    assert CE1.jaccard(CE2) == 1.0  # all of the unique kmers in seq1 show up in seq2
    assert CE2.jaccard(CE1) == 0.5  # half of the unique kmers in seq2 show up in seq1


def test_import_export():
    CE1 = CountEstimator(n=5, max_prime=9999999999971., ksize=1)
    CE2 = CountEstimator(n=5, max_prime=9999999999971., ksize=1)
    sequence1 = "AAAA"
    sequence2 = "AAAACCCC"
    CE1.add_sequence(sequence1)
    CE2.add_sequence(sequence2)
    temp_file = tempfile.mktemp()  # Make temporary file
    CE1.export(temp_file)  # Export the CountEstimator to temp file
    CE_Import = import_single_hdf5(temp_file)  # Read in the data
    os.remove(temp_file)  # Remove the temporary file
    assert CE_Import.jaccard_count(CE2) == CE1.jaccard_count(CE2)  # Make sure the results of the import are the same as the original CountEstimator


def test_hash_list():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1='acgtagtctagtctacgtagtcgttgtattataaaatcgtcgtagctagtgctat'
    CE1.add_sequence(seq1)
    hash_list = {424517919, 660397082}
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, hash_list=hash_list, save_kmers='y')
    CE2.add_sequence(seq1)
    assert CE1.jaccard(CE2) == 0.4
    assert CE1.jaccard_count(CE2) == (1.0, 2/7.)


def test_vector_formation():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq2 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq3 = 'ttgatactcaatccgcatgcatgcatgacgatgcatgatgtacgactgatgcatgatcgaactgatgcactcgtgatgczxerqwewdfhg'
    CE1.add_sequence(seq1)
    CE2.add_sequence(seq2)
    CE3.add_sequence(seq3)
    Y = CE1.count_vector([CE1, CE2, CE3])
    assert (Y == np.array([1.,1.,0.5625])).all()
    Y2 = CE1.jaccard_vector([CE1, CE2, CE3])
    assert (Y2 == np.array([1.,1.,0.4])).all()


def form_matrices_test():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq2 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq3 = 'ttgatactcaatccgcatgcatgcatgacgatgcatgatgtacgactgatgcatgatcgaactgatgcactcgtgatgczxerqwewdfhg'
    CE1.add_sequence(seq1)
    CE2.add_sequence(seq2)
    CE3.add_sequence(seq3)
    A = form_jaccard_count_matrix([CE1, CE2, CE3])
    assert (A == np.array([[1., 1., 0.80952380952380953], [1., 1., 0.80952380952380953], [0.5625, 0.5625, 1.]])).all()
    B = form_jaccard_matrix([CE1, CE2, CE3])
    assert (B == np.array([[1., 1., 0.4], [1., 1., 0.4], [0.4, 0.4, 1.]])).all()


def test_suite():
    """
    Runs all the test functions
    :return:
    """
    from sys import platform as _platform
    test_jaccard_1()
    test_jaccard_2_difflen()
    test_yield_overlaps()
    test_yield_overlaps_2()
    test_yield_overlaps_3()
    test_CountEstimator()
    test_import_export()
    test_hash_list()
    test_vector_formation()
    form_matrices_test()
