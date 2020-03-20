#!/bin/bash
#set -e
# this script will demonstrate how to re-train CMash on a custom database

trainingFiles="training_files.txt"  # location of the full path names to the training genomes (fasta/q, compressed or not)
numThreads=8  # number of threads to utilize for the re-training
cmashBaseName="cmash_db_n1000_k60"  # base name you wish for the output training database
outputDir="test"  # name of the output directory you would like the results to be placed into
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
python ../../scripts/MakeStreamingPrefilter.py ${cmashDatabase} ${outputDir}/${prefilterName} ${kmerRange}

: << 'END'
# dump all the k-mers in the new training database
echo "dumping training k-mers"
rm ${cmashDump} 2> /dev/null
python dump_kmers.py ${cmashDatabase} ${cmashDump}

# do the kmc counting of the k-mers
echo "running kmc"
rm "${cmashBaseName}.kmc_pre" 2> /dev/null
rm "${cmashBaseName}.kmc_suf" 2> /dev/null
# count all the k-mers in the training database, -ci0 means include them all, -cs3 use small counters
../KMC/bin/kmc -v -k60 -fa -ci0 -cs3 -t"${numThreads}" -jlogsample ${cmashDump} "${cmashBaseName}_dump" .


# Next, make a silly mock community
echo "mock community organisms:"
head -n 2 ${trainingFiles}
rm $bbmapRef 2> /dev/null
head -n 2 ${trainingFiles} | xargs -I{} zcat {} >> "${bbmapRef}_temp"
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "${bbmapRef}_temp" | tail -n +2 > "${bbmapRef}"
rm "${bbmapRef}_temp"
eval ${bbmapRandomReads} ref=${bbmapRef} out="${mockCommunity}_temp.fa" adderrors=false simplenames=t minlength=250 maxlength=250 maxsnps=0 coverage=20 maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 maxns=0 maxinslen=0 maxdellen=0 maxsublen=0 maxnlen=0 mininslen=0 mindellen=0 minsublen=0 minnlen=0
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "${mockCommunity}_temp.fa" | tail -n +2 > "${mockCommunity}"
rm "${mockCommunity}_temp.fa"

# Run metalign on the new data
#cd ..
#python metalign.py local_tests/mock_community.fa local_tests/ --output local_tests/test_out.csv --verbose --keep_temp_files --temp_dir local_tests/metalign_temp

# then take a look at local_tests/test_out.csv and make sure everything worked. Should show two strains corresponding to the ones used to form the mock community

END