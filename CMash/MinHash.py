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
from blist import *  # note, the import functions import the _mins etc. as lists, and the CE class imports them as blists.
# This shouldn't cause an issue, but will lead to slow performance if a CE is imported, then additional things are added.
# I.e. If you import a CE, don't add new elements, or you might have a bad day (or at least a long one).
import bisect
import ctypes
import warnings
from six import string_types
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
            if CE._kmers is not None:
                kmer_data = subgrp.create_dataset("kmers", data=CE._kmers)

            subgrp.attrs['class'] = np.string_("CountEstimator")
            subgrp.attrs['filename'] = np.string_(CE.input_file_name)  # But keep the full file name on hand
            subgrp.attrs['ksize'] = CE.ksize
            subgrp.attrs['prime'] = CE.p
            subgrp.attrs['true_num_kmers'] = CE._true_num_kmers
        except ValueError:
            fid.close()
            raise Exception("It appears that the training file name %s exists twice in the input data. Please make sure all names are unique (i.e. remove duplicates) and try again." % CE.input_file_name)

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
        fid.close()
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
            fid.close()
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


def delete_from_database(database_location, delete_list):
    """
    This function will delete specified entries from count estimators stored in a single HDF5 file.
    :param database_location: file name for the single HDF5 database file
    :param delete_list: List of names of files to delete (or a single name)
    :return: None
    """
    file_name = os.path.abspath(database_location)
    fid = h5py.File(file_name, 'r')  # Let's take a peak at what's already in there
    grp = fid["CountEstimators"]
    #all_keys = grp.keys()  # All the current data sets in there
    all_keys = [item for item in grp.keys()]  # All the current data sets in there (python 3 loads it as a view)
    fid.close()
    fid = h5py.File(file_name, "a")  # Open in append mode to modify in place
    if isinstance(delete_list, string_types):  # If it's just a single name, delete it
        keys_to_delete = delete_list  # Yeah it's bad naming (should be key_to_delete) but I need the print command below to show up
        if keys_to_delete in all_keys:  # If it's actually one of the keys in there
            del fid["CountEstimators"][keys_to_delete]  # Delete it
    else:
        keys_to_delete = list(set(map(os.path.basename, delete_list)))  # Uniqueify and get base names
        for key_to_delete in keys_to_delete:  # If it's a whole list
            if key_to_delete in all_keys:  # If it's actually one of the keys in there
                del fid["CountEstimators"][key_to_delete]  # Delete it
    fid.close()
    print("The following entries have been deleted from %s:" % file_name)
    print(keys_to_delete)


def insert_to_database(database_location, insert_list):
    """
    This function will insert specified FASTA/Q files into the HDF5 database at database_location
    :param database_location: location of HDF5 database
    :param insert_list: list of (full paths) to FASTA/Q files to insert
    :return: None
    """
    if isinstance(insert_list, string_types):
        insert_list = [insert_list]  # Make into a list if it's only one guy
    insert_list = list(set(insert_list))  # Uniqueify it
    # Get the info from one of the existing guys
    fid = h5py.File(database_location, 'r')
    grp = fid["CountEstimators"]
    #orig_keys = grp.keys()
    orig_keys = [item for item in grp.keys()]  # Python3 h5py imports it as a view, not the data itself
    fid.close()
    temp_CE = import_multiple_from_single_hdf5(database_location, import_list=[orig_keys[0]])[0]
    k_size = temp_CE.ksize
    p = temp_CE.p
    n = len(temp_CE._mins)
    if temp_CE._kmers is not None:
        save_kmers = 'y'
    else:
        save_kmers = 'n'
    rev_comp = False
    # Now create the new sketches and stick them in the database
    fid = h5py.File(database_location, 'a')
    grp = fid["CountEstimators"]
    for insert_file_name in insert_list:
        if os.path.basename(insert_file_name) not in orig_keys:
            # Populate the hash and find the true number of k-mers
            kmers = set()
            to_insert_CE = CountEstimator(n=n, max_prime=p, ksize=k_size, save_kmers=save_kmers)
            for record in screed.open(insert_file_name):
                seq = record.sequence
                for i in range(len(seq) - k_size + 1):
                    kmer = seq[i:i + k_size]
                    # No need to care about revcomp since nodegraph ID's them
                    kmers.add(kmer)
                    to_insert_CE.add(kmer)
            to_insert_CE._true_num_kmers = len(kmers)
            to_insert_CE.input_file_name = insert_file_name
            # Now stick it in the database
            subgrp = grp.create_group(
                os.path.basename(insert_file_name))  # the key of a subgroup is the basename (not the whole file)
            mins_data = subgrp.create_dataset("mins", data=to_insert_CE._mins)
            counts_data = subgrp.create_dataset("counts", data=to_insert_CE._counts)
            if to_insert_CE._kmers:
                kmer_data = subgrp.create_dataset("kmers", data=to_insert_CE._kmers)
            subgrp.attrs['class'] = np.string_("CountEstimator")
            subgrp.attrs['filename'] = np.string_(to_insert_CE.input_file_name)  # But keep the full file name on hand
            subgrp.attrs['ksize'] = to_insert_CE.ksize
            subgrp.attrs['prime'] = to_insert_CE.p
            subgrp.attrs['true_num_kmers'] = to_insert_CE._true_num_kmers
    fid.close()


def union_databases(database1_file_name, database2_file_name, out_file):
    """
    This funtion will union the two training databases and export to the specified HDF5 file.
    :param database1_file_name: Input HDF5 file of one of the input databases.
    :param database2_file_name: Input HDF5 file of the other input database.
    :param out_file: File name for the exported HDF5 file
    :return: None
    """
    out_file = os.path.abspath(out_file)
    database1_file_name = os.path.abspath(database1_file_name)
    database2_file_name = os.path.abspath(database2_file_name)
    # Read in the existing database
    CEs1 = import_multiple_from_single_hdf5(database1_file_name)
    CEs2 = import_multiple_from_single_hdf5(database2_file_name)
    # Sanity check the existing databases
    if len(set([item.ksize for item in CEs1]).union(set([item.ksize for item in CEs2]))) > 1:
        raise Exception("Incompatible k-mer lengths.")
    if len(set([item.p for item in CEs1]).union(set([item.p for item in CEs2]))) > 1:
        raise Exception("Incompatible primes. Re-run with same -p value for both databases")
    # Join up the two existing databases
    all_CEs = list(set(CEs1).union(set(CEs2)))
    all_input_names = set([os.path.basename(item.input_file_name) for item in CEs1]).union(
        set([os.path.basename(item.input_file_name) for item in CEs2]))
    # For some odd reason, even if the data is identical, it's identifying the same CEs (loaded from different sources) as the same
    # so instead, let's just pick one of them to include
    included_names = set()
    to_include_CEs = list()
    for CE in all_CEs:
        if CE.input_file_name not in included_names:
            included_names.add(CE.input_file_name)
            to_include_CEs.append(CE)
    #if len(all_input_names) != len(all_CEs):
    #    raise Exception("Number of items in union does not match total number of unique input FASTQ/A file names")
    # export the union of the databases
    export_multiple_to_single_hdf5(to_include_CEs, out_file)


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


def test_form_matrices():
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

def test_delete_from_database():
    seq1 = "ATCGTATGAGTATCGTCGATGCATGCATCGATGCATGCTACGTATCGCATGCATG"
    seq2 = "ATCTACTCAACATTAACTACTCATATTAACTCACATTCATATCCATACTACTCGT"
    seq3 = "ACTCATGTTAGATCGATATTGACTGATGACTCGTTGCACTGCATGCTGCATGATGC"
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE1.add_sequence(seq1)
    CE2.add_sequence(seq2)
    CE3.add_sequence(seq3)
    # CE's must have names to delete them
    CE1.input_file_name = "seq1"
    CE2.input_file_name = "seq2"
    CE3.input_file_name = "seq3"
    temp_file = tempfile.mktemp()
    export_multiple_to_single_hdf5([CE1, CE2, CE3], temp_file)
    fid = h5py.File(temp_file, 'r')
    assert len(fid["CountEstimators"].keys()) == 3
    fid.close()
    # delete one
    delete_from_database(temp_file, "seq1")
    # Check if it was deleted
    fid = h5py.File(temp_file, 'r')
    assert len(fid["CountEstimators"].keys()) == 2
    fid.close()
    # Delete two
    delete_from_database(temp_file, ["seq2", "seq3"])
    fid = h5py.File(temp_file, 'r')
    assert len(fid["CountEstimators"].keys()) == 0
    fid.close()
    os.remove(temp_file)

def test_insert_to_database():
    try:
        import CMash
        file1 = CMash.get_data("PRJNA67111.fna")
        file2 = CMash.get_data("PRJNA32727.fna")
        file3 = CMash.get_data("PRJNA298068.fna")
    except ImportError:
        file1 = os.path.join(os.path.dirname(__file__), "data", "PRJNA67111.fna")
        file2 = os.path.join(os.path.dirname(__file__), "data", "PRJNA32727.fna")
        file3 = os.path.join(os.path.dirname(__file__), "data", "PRJNA298068.fna")
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y', input_file_name=file1)
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y', input_file_name=file2)
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y', input_file_name=file3)
    temp_file = tempfile.mktemp()
    export_multiple_to_single_hdf5([CE1], temp_file)
    insert_to_database(temp_file, file2)  # Tests if only a single file is inserted
    CEs = import_multiple_from_single_hdf5(temp_file)
    assert len(CEs) == 2
    assert len(CEs[0]._kmers) == len(CEs[1]._kmers)
    insert_to_database(temp_file, [file2, file3])  # Tests if a list of files is inserted (with some already in the database)
    CEs = import_multiple_from_single_hdf5(temp_file)
    assert len(CEs) == 3
    assert len(CEs[0]._kmers) == len(CEs[2]._kmers)
    os.remove(temp_file)


def test_union_databases():
    try:
        import CMash
        file1 = CMash.get_data("PRJNA67111.fna")
        file2 = CMash.get_data("PRJNA32727.fna")
        file3 = CMash.get_data("PRJNA298068.fna")
    except ImportError:
        file1 = os.path.join(os.path.dirname(__file__), "data", "PRJNA67111.fna")
        file2 = os.path.join(os.path.dirname(__file__), "data", "PRJNA32727.fna")
        file3 = os.path.join(os.path.dirname(__file__), "data", "PRJNA298068.fna")
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y', input_file_name=file1)
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y', input_file_name=file2)
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y', input_file_name=file3)
    temp_file1 = tempfile.mktemp()
    temp_file2 = tempfile.mktemp()
    temp_file3 = tempfile.mktemp()
    export_multiple_to_single_hdf5([CE1], temp_file1)
    export_multiple_to_single_hdf5([CE1, CE2, CE3], temp_file2)
    try:
        union_databases(temp_file1, temp_file2, temp_file3)
    except:
        raise Exception("Unioning databases test did not succeed")
    all_3 = import_multiple_from_single_hdf5(temp_file3)
    assert len(all_3) == 3
    assert len(set([item.input_file_name for item in all_3])) == 3
    os.remove(temp_file1)
    os.remove(temp_file2)
    os.remove(temp_file3)

def test_suite():
    """
    Runs all the test functions
    :return:
    """
    #from sys import platform as _platform
    test_jaccard_1()
    test_jaccard_2_difflen()
    test_yield_overlaps()
    test_yield_overlaps_2()
    test_yield_overlaps_3()
    test_CountEstimator()
    test_import_export()
    test_hash_list()
    test_vector_formation()
    test_form_matrices()
    test_delete_from_database()
    test_insert_to_database()
    test_union_databases()
    print("All tests successful!")

