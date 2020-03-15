"""
This module computes the intersection of the kmers in the input sample
and the kmers in the training sample prior to testing each training sample
for membership in the input sample. This eliminates the need to pre-screen kmers
for presence in the TST via the prefilter since every kmer in the input must be
in both the input and the training samples.
"""
import argparse, math, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/CMash")
import MinHash as MH
import itertools

__location__ = "your KMC location here"
numThreads=8
#cmashBaseName="cmash_db_n1000_k60"
cmashBaseName="TrainingDatabase"
cmashDatabase=cmashBaseName + ".h5"
cmashDump=cmashBaseName + "_dump.fa"
db_21mers_loc = 'TrainingDatabase_dump'


kmc_loc = __location__ + 'kmc'
kmc_dump_loc = __location__ + 'kmc_dump'
kmc_tools_loc = __location__ + 'kmc_tools'

"""Parse command line arguments."""
def parse_args():    # handle user arguments
    parser = argparse.ArgumentParser(description="Intersect kmers from training samples and input sample")
    parser.add_argument('reads', help='Path to reads file.')
    parser.add_argument('--input_type', default='AUTO',
            choices=['fastq', 'fasta', 'AUTO'],
            help='Type of input file (fastq/fasta). Default: try to auto-determine')
    parser.add_argument('--temp_dir', default = 'AUTO/',
            help='Directory to write temporary files to.')
    parser.add_argument('--threads', type=int, default=4,
            help='How many compute threads for KMC to use. Default: 4')
    args = parser.parse_args()
    return args

"""
This code is based on dump_kmers.py in the Metalign repository.
Dumps the kmers from the input's CountEstimator to a fasta file.
"""
def dump_training_kmers():
    print("dumping training k-mers")

    with open("/dev/null", 'w') as f:
        subprocess.Popen(['rm', cmashDump], stderr=f).wait() 

    training_database = cmashDatabase  # first input is the training file name
    dump_file = cmashDump  # second input is the desired output dump file
    CEs = MH.import_multiple_from_single_hdf5(training_database)
    fid = open(dump_file, 'w')
    i = 0
    for CE in CEs:
        for kmer in CE._kmers:
            fid.write('>seq%d\n' % i)
            fid.write('%s\n' % kmer)
            i += 1

    fid.close()

"""
Opens a KMC process to count the dumped training kmers.
-ci1 excludes all kmers which appear less than one time (excludes no kmers).
"""
def count_training_kmers():
    with open("/dev/null", 'w') as f:
        subprocess.Popen(['rm', cmashBaseName+".kmc_pre"], stderr=f).wait() 
        subprocess.Popen(['rm', cmashBaseName+".kmc_suf"], stderr=f).wait() 

    subprocess.Popen([kmc_loc, '-v', '-k21', '-fa', '-ci1', \
            '-t'+str(numThreads), '-jlogsample', cmashDump,\
            cmashBaseName+'_dump', '.']).wait()


"""
Opens a KMC process to count the dumped kmers from the input sample.
-ci1 excludes all kmers which appear less than one time (excludes no kmers).
"""
def count_input_kmers(args):
    #db_21mers_loc = args.data + 'cmash_db_n1000_k60_dump'
    if args.input_type == 'fastq':
            type_arg = '-fq'
    else:
            type_arg = '-fa'

    #count kmers in the input sample (not the training db)
    subprocess.Popen([kmc_loc, '-v', '-k21', type_arg, '-ci1',
            '-t' + str(args.threads), '-jlog_sample', args.reads,
            args.temp_dir + 'reads_21mers', args.temp_dir]).wait()

"""
This source is from the Metalign repository.
Opens a new KMC process to intersect the kmers
from the input sample and the training database.
The kmers in the intersection are written to a file as-is. Then,
the contents of the file are read and rewritten in FASTA format. 
"""
def intersect(args):
    #intersect kmers
    print("Intersecting input sample & training sample...")
    subprocess.Popen([kmc_tools_loc, 'simple', db_21mers_loc,
            args.temp_dir + 'reads_21mers', 'intersect',
            args.temp_dir + '21mers_intersection']).wait()

    #dump intersection
    print ("dumping intersection to FASTA file")
    subprocess.Popen([kmc_dump_loc, args.temp_dir + '21mers_intersection',
            args.temp_dir + '21mers_intersection_dump']).wait()

    #read intersection & rewrite in fasta format
    with(open(args.temp_dir + '21mers_intersection_dump', 'r')) as infile:
            with(open(args.temp_dir + '21mers_intersection_dump.fa', 'w')) as fasta:
                    for line in infile:
                            seq = line.split()[0]
                            fasta.write('>seq' + '\n' + seq + '\n')



if __name__ == "__main__":
    args = parse_args()
    dump_training_kmers()
    count_training_kmers()
    count_input_kmers(args)
    intersect(args)
