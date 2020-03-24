from CMash import MinHash as MH
import marisa_trie as mt
import khmer
import timeit
import os
import h5py


class MakeTSTNew:
    def __init__(self, training_database_file_name: str, TST_export_file_name: str):
        self.training_database_file_name = training_database_file_name
        self.TST_export_file_name = TST_export_file_name

    @staticmethod
    def yield_trie_items_to_insert_no_import(file_name):
        fid = h5py.File(file_name, 'r')
        if "CountEstimators" not in fid:
            fid.close()
            raise Exception("This function imports a single HDF5 file containing multiple sketches."
                            " It appears you've used it on a file containing a single sketch."
                            "Try using import_single_hdf5 instead")

        grp = fid["CountEstimators"]
        iterator = grp.keys()

        iterator = sorted(iterator, key=os.path.basename)  # sort so that we know the order of the input

        for (i, key) in enumerate(iterator):
            if key not in grp:
                fid.close()
                raise Exception("The key " + key + " is not in " + file_name)

            subgrp = grp[key]
            if "kmers" not in subgrp:
                raise Exception("Kmers were not saved when creating the count estimators. Please make sure save_kmers='y' "
                                "when creating the count estimators.")
            else:
                temp_kmers = subgrp["kmers"][...]
                kmers = [kmer.decode('utf-8') for kmer in temp_kmers]
                for (kmer_index, kmer) in enumerate(kmers):
                    # add both the original k-mer and the reverse complement, as the MinHashes were created without reverse complement
                    if kmer:
                        yield kmer + 'x' + str(i) + 'x' + str(kmer_index)  # format here is kmer+x+hash_index+kmer_index
                        # rev-comp kmer
                        kmer_rc = khmer.reverse_complement(kmer)
                        yield kmer_rc + 'x' + str(i) + 'x' + str(kmer_index)  # format here is kmer+x+hash_index+kmer_index

    def make_TST(self):
        tree = mt.Trie(self.yield_trie_items_to_insert_no_import(self.training_database_file_name))
        tree.save(self.TST_export_file_name)


class MakeTSTOld:
    def __init__(self, genome_sketches: list, TST_export_file_name: str):
        self.genome_sketches = genome_sketches
        self.TST_export_file_name = TST_export_file_name

    @classmethod
    def make_TST(self):
        genome_sketches = self.genome_sketches
        to_insert = set()
        # add both the original k-mer and the reverse complement, as the MinHashes were created without reverse complement
        for i in range(len(genome_sketches)):
            for kmer_index in range(len(genome_sketches[i]._kmers)):
                # normal kmer
                kmer = genome_sketches[i]._kmers[kmer_index]
                # only insert the kmer if it's actually non-empty
                if kmer:
                    to_insert.add(
                        kmer + 'x' + str(i) + 'x' + str(kmer_index))  # format here is kmer+x+hash_index+kmer_index
                    # rev-comp kmer
                    kmer = khmer.reverse_complement(genome_sketches[i]._kmers[kmer_index])
                    to_insert.add(
                        kmer + 'x' + str(i) + 'x' + str(kmer_index))  # format here is kmer+x+hash_index+kmer_index

        # export the TST
        tree = mt.Trie(to_insert)
        tree.save(self.TST_export_file_name)
