import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+"/CMash/CMash")
import MinHash as MH
import itertools
training_database = sys.argv[1]  # first input is the training file name
dump_file = sys.argv[2]  # second input is the desired output dump file
CEs = MH.import_multiple_from_single_hdf5(training_database)
fid = open(dump_file, 'w')
i = 0
for CE in CEs:
    for kmer in CE._kmers:
        fid.write('>seq%d\n' % i)
        fid.write('%s\n' % kmer)
        i += 1

fid.close()
