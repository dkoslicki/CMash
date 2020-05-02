#!/bin/bash
#set -e
# A basic end-to-end workflow to make sure things are working correctly (esp wrt integration with CMash)
# basic steps:
# make small training database
# re-train Cmash
# make a mock community
# use metalign to profile the mock community (results should be two strains corresponding to the two selected organisms)

# requires bbmap: https://sourceforge.net/projects/bbmap/

bbmapRandomReads="~/Documents/bbmap/randomreads.sh"
trainingFiles="training_files.txt"
numThreads=8
cmashBaseName="cmash_db_n1000_k60"
cmashDatabase="${cmashBaseName}.h5"
cmashDump="${cmashBaseName}_dump.fa"
prefilterName="cmash_filter_n1000_k60_30-60-10.bf"
bbmapRef='mock_reference.fa'
mockCommunity="mock_community.fa"



# copy some of the training data over
echo "getting organism files"
rm -rf organism_files 2> /dev/null
mkdir organism_files
find ../data/organism_files/ -type f -name "*.fna.gz" | head -n 25 | xargs -I{} cp {} organism_files/

# Pull off the taxonomy information
echo "getting taxonomy"
rm db_info.txt 2> /dev/null
rm headers.txt 2> /dev/null
for file in `find organism_files/ -type f -name "*.fna.gz"`
do
  zcat ${file} | grep '>' | cut -d' ' -f1 | sed 's/>//g' >> headers.txt  # pull off accessions
done
LC_ALL=C grep -F -w -f headers.txt ../data/db_info.txt > db_info.txt  # look for the taxonomy info
rm headers.txt

# store the names of the training files
rm $trainingFiles 2> /dev/null
fullPath="$(pwd)/organism_files"
find $fullPath -type f -name "*.fna.gz" > ${trainingFiles}

# re-train CMash
echo "re-training CMash"
rm ${cmashDatabase} 2> /dev/null
python ../CMash/scripts/MakeStreamingDNADatabase.py ${trainingFiles} ${cmashDatabase} -n 1000 -k 60 #-v

# make streaming pre-filter
echo "making streaming pre-filter"
rm ${prefilterName} 2> /dev/null
python ../CMash/scripts/MakeStreamingPrefilter.py ${cmashDatabase} ${prefilterName} 30-60-10

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

