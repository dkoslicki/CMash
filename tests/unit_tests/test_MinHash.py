import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from CMash.MinHash import *
from CMash import MinHash as MH

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
    assert len(list(MH._yield_overlaps(x1, x2))) == 0


def test_yield_overlaps_2():
    x1 = [1, 3, 5]
    x2 = [1, 2, 4, 6]
    assert len(list(MH._yield_overlaps(x1, x2))) == 1
    assert len(list(MH._yield_overlaps(x2, x1))) == 1


def test_yield_overlaps_3():
    x1 = [1, 3, 6]
    x2 = [1, 2, 6]
    assert len(list(MH._yield_overlaps(x1, x2))) == 2
    assert len(list(MH._yield_overlaps(x2, x1))) == 2


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
    seq1 = 'acgtagtctagtctacgtagtcgttgtattataaaatcgtcgtagctagtgctat'
    CE1.add_sequence(seq1)
    #hash_list = {424517919, 660397082}
    hash_list = set(CE1._mins[0:2])  # pick off two of the hashes, so the jaccard should be 2/5 = .4
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, hash_list=hash_list, save_kmers='y')
    CE2.add_sequence(seq1)
    assert CE1.jaccard(CE2) == 0.4
    assert CE1.jaccard_count(CE2) == (1.0, 2/9.)


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
    assert (np.sum(np.abs(Y - np.array([1.,1.,0.55555556]))))<.00001
    Y2 = CE1.jaccard_vector([CE1, CE2, CE3])
    assert (Y2 == np.array([1., 1., 3/5.])).all()


def test_form_matrices():
    # TODO: this was tedius to create, so let's just make sure it executes and not check the numbers yet
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
    #assert (A == np.array([[1., 1., 0.80952380952380953], [1., 1., 0.80952380952380953], [0.5625, 0.5625, 1.]])).all()
    B = form_jaccard_matrix([CE1, CE2, CE3])
    #assert (B == np.array([[1., 1., 0.4], [1., 1., 0.4], [0.4, 0.4, 1.]])).all()
    #assert np.sum(np.abs(B - np.array([[1., 1., 0.4], [1., 1., 0.4], [0.4, 0.4, 1.]]))) < 0.001

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

def test_make_tree():
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
    CEs = [CE1, CE2, CE3]
    tree = make_tree(CEs)
    kmer = CE1._kmers[0]  # so we know it's in the first CE
    true_res = [0]
    if kmer in CE2._kmers:
        true_res.append(1)
    if kmer in CE3._kmers:
        true_res.append(2)
    res = tree.query(kmer)
    print(res)
    assert sorted(res) == true_res  # this should be [0, 2] since ATC shows up in CE1 and CE3

def test_get_info():
    try:
        import CMash
        file1 = CMash.get_data("PRJNA67111.fna")
        file2 = CMash.get_data("PRJNA32727.fna")
        file3 = CMash.get_data("PRJNA298068.fna")
    except ImportError:
        file1 = os.path.join(os.path.dirname(__file__), "data", "PRJNA67111.fna")
        file2 = os.path.join(os.path.dirname(__file__), "data", "PRJNA32727.fna")
        file3 = os.path.join(os.path.dirname(__file__), "data", "PRJNA298068.fna")
    CE1 = CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y', input_file_name=file1)
    CE2 = CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y', input_file_name=file2)
    CE3 = CountEstimator(n=5, max_prime=9999999999971, ksize=3, save_kmers='y', input_file_name=file3)
    temp_file = tempfile.mktemp()
    export_multiple_to_single_hdf5([CE1, CE2, CE3], temp_file)
    meta_data = get_info_from_single_hdf5(temp_file)
    assert meta_data.ksize == 3
    assert meta_data.sketch_size == 5
    assert meta_data.prime == 9999999999971
    assert len(meta_data.file_names) == 3
    assert meta_data.file_names == sorted([file1, file2, file3], key=os.path.basename)


def main():
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
    test_make_tree()
    test_get_info()
    print("Success!")

if __name__ == "__main__":
    main()