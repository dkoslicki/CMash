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

def main():
    import timeit
    from CMash import MinHash as MH
    small_database_file = "/home/dkoslicki/Desktop/CMash/tests/TempData/cmash_db_n5000_k60_1000.h5"
    TST_export_file_new = "/home/dkoslicki/Desktop/CMash/tests/TempData/cmash_db_n5000_k60_new.tst"
    TST_export_file_old = "/home/dkoslicki/Desktop/CMash/tests/TempData/cmash_db_n5000_k60_old.tst"

    # new way
    t0 = timeit.default_timer()
    M = MakeTSTNew(small_database_file, TST_export_file_new)
    M.make_TST()
    t1 = timeit.default_timer()
    print(f"New timing: {t1 - t0}")

    # old way
    CEs = MH.import_multiple_from_single_hdf5(small_database_file)
    t0 = timeit.default_timer()
    M = MakeTSTOld(CEs, TST_export_file_old)
    M.make_TST()
    t1 = timeit.default_timer()
    print(f"Old timing: {t1 - t0}")


# simple way to do testing
if __name__ == "__main__":
    main()
