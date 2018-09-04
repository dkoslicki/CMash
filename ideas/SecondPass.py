# ok, so to increase the specificity, I could
# 1. Use MiCOP/alignment approach (a pain)
# 2. Take a look at the hit_matrices to find informative k-mers (a-la Kraken, but in a much more intelligent way)
# 3. Take a second pass: use the non-zero found guys to make new sketches (with many more k-mers, perhaps a larger k-mer size), and do it again

# TODO: Note: this, oddly enough, INCREASES the number of FP's (note it needs a much higher cutoff
import pandas as pd
import os

# let's try #3 here

# first, read in an example, and find the ones for which to make a reduced database
cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.csv'
file_names_file = '/home/dkoslicki/Data/MiCOPMinHash/GoodFileNamesAbsolute.txt'
out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/reduced_file_names.txt'
coverage_threshold = 0
sort_key = None

# read in the file and sort as needed
df = pd.read_csv(cmash_out_file, index_col=0)
max_key = df.keys()[-1]
if not sort_key:
	sort_key = max_key
cmash_thresh = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)
names_passed_thresh = list(cmash_thresh.index)

# read in all the reference files
ref_file_names_to_loc = dict()
fid = open(file_names_file, 'r')
for line in fid.readlines():
	line = line.strip()
	base_name = os.path.basename(line)
	ref_file_names_to_loc[base_name] = line

fid.close()

# get the locations of the reference sequences to build
out_fid = open(out_file, 'w')
ref_loc_to_build = []
for name in names_passed_thresh:
	ref_loc_to_build.append(ref_file_names_to_loc[name])  # if using the bz2 extension

with open(out_file, 'w') as fid:
	for line in ref_loc_to_build:
		fid.write(line)
		fid.write("\n")

