#!/bin/bash
#set -e
# this script will demonstrate how to re-train CMash on a custom database

trainingFiles="training_files.txt"  # location of the full path names to the training genomes (fasta/q, compressed or not)
numThreads=8  # number of threads to utilize for the re-training
cmashBaseName="cmash_db_n1000_k60"  # base name you wish for the output training database
outputDir="test"  # name of the output directory you would like the results to be placed into. THIS DIRECTORY MUST ALREADY EXIST # FIXME
kmerRange="30-60-10"  # range of k-mer sizes you want to be able to query with. This translates to [30, 40, 50, 60]


# don't change these unless you know what you're doing
cmashDatabase="${cmashBaseName}.h5"
cmashDump="${cmashBaseName}_dump.fa"
prefilterName="${cmashBaseName}_${kmerRange}.bf"
maxKsize=$(echo ${kmerRange} | cut -d'-' -f2)
numHashes=1000


# re-train CMash
echo "re-training CMash"
rm ${cmashDatabase} 2> /dev/null
python ../../scripts/MakeStreamingDNADatabase.py ${trainingFiles} ${outputDir}/${cmashDatabase} -n ${numHashes} -k ${maxKsize} -v

# make streaming pre-filter
echo "making streaming pre-filter"
rm ${prefilterName} 2> /dev/null
python ../../scripts/MakeStreamingPrefilter.py ${outputDir}/${cmashDatabase} ${outputDir}/${prefilterName} ${kmerRange}
