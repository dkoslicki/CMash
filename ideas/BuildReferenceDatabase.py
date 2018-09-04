# This script will build a reference database, given a CMash classified file
import pandas as pd
import os
import sys
import screed
import khmer
from screed import Record
from khmer.utils import write_record

# small test data
#cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/Test.csv'
#file_names_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/FileNames.txt'
#file_names_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/Fused/FileNames.txt'  # using the bbmap fuse thingy

# on a real sample
cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.csv'
file_names_file = '/home/dkoslicki/Data/MiCOPMinHash/GoodFileNamesAbsolute.txt'

out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/ref.fa'
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
	#ref_loc_to_build.append(ref_file_names_to_loc[os.path.splitext(name)[0]])  # since bbmap stuff changed the file extension


# This uses khmer to merge the contigs and put them in one fasta file
for loc in ref_loc_to_build:
	print(os.path.basename(loc))
	fid = khmer.ReadParser(loc)
	seq = ""
	i = 0
	for record in fid:
		if i == 0:
			header = record.name
		seq += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
		seq += record.sequence
		i += 1
	print("there")
	record = Record(name=header, sequence=seq)
	write_record(record, out_fid)
	fid.close()
out_fid.close()

# This relies on using bbmap to do the contig merging, and then will use cat to concatenate them
# ls *.bz2 | xargs -P 4 -I{} ~/Documents/bbmap/./fuse.sh in={} out=Fused/{}
# cd Fused
# ls *.bz2 | xargs -P 4 -I{} bzip2 -d {}


for loc in ref_file_names_to_loc.values():#ref_loc_to_build:
	with open(loc) as infile:
		for line in infile:
			out_fid.write(line)
out_fid.close()


