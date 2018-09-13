# find the unique k-mers
# TP's should include at least some unique k-mer hits
# possibly include some with overwhelming amount of non-unique hits?
# or sets of ones where there is no unique hits among them
# scipy.sparse.save_npz

from CMash import MinHash as MH
import numpy as np
import scipy as sp
from scipy.io import loadmat
import pandas as pd


# First, read in the sketches of just the appropriate critters
#cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.csv'
cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/test.csv'
training_base_name = '/nfs1/Koslicki_Lab/koslickd/RepoPhlAn-7-24-18/out/microbes_24072018/fna/'
#training_hdf_file = '/home/dkoslicki/Data/MiCOPMinHash/AllBacteria.hd5'
training_hdf_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/Test.h5'
hit_matrices_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/test_hit_matrix.npz'
coverage_threshold = 0
sort_key = 'k=5'

# read in the file and sort as needed
df = pd.read_csv(cmash_out_file, index_col=0)
max_key = df.keys()[-1]
if not sort_key:
	sort_key = max_key
cmash_thresh = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)
names_passed_thresh = list(cmash_thresh.index)
names_passed_thresh_with_path = []
for name in names_passed_thresh:
	names_passed_thresh_with_path.append(training_base_name + name)
CEs = MH.import_multiple_from_single_hdf5(training_hdf_file, import_list=names_passed_thresh_with_path)

# impor the hit matrices
#hit_matrices_dict = np.load(hit_matrices_file)
hit_matrices_dict = loadmat(hit_matrices_file)

# now, for each one of the sketches, look for unique k-mer in it
k_range = sorted([int(i.split('=')[1]) for i in cmash_thresh.keys()])

