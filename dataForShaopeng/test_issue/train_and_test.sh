#!/bin/bash
# make the training database
#/home/dkoslicki/anaconda3/envs/CMash/bin/python ../../scripts/MakeStreamingDNADatabase.py file_names.txt TrainingDatabase_k_61.h5 -n 1000 -k 61 -v

# remove all output files
rm out*

# test a couple of ranges
inputFile="/home/dkoslicki/Desktop/CMash/dataForShaopeng/small_data/taxid_1909294_104_genomic.fna.gz"
max=61

# problem one
min=61
step=1
/home/dkoslicki/anaconda3/envs/CMash/bin/python ../../scripts/StreamingQueryDNADatabase.py ${inputFile} TrainingDatabase_k_61.h5 out_${min}_${max}_${step}.csv ${min}-${max}-${step} -c 0 -l 0 --sensitive -v


min=60
max=61
/home/dkoslicki/anaconda3/envs/CMash/bin/python ../../scripts/StreamingQueryDNADatabase.py ${inputFile} TrainingDatabase_k_61.h5 out_${min}_${max}_${step}.csv ${min}-${max}-${step} -c 0 -l 0 --sensitive -v

min=50
max=61
/home/dkoslicki/anaconda3/envs/CMash/bin/python ../../scripts/StreamingQueryDNADatabase.py ${inputFile} TrainingDatabase_k_61.h5 out_${min}_${max}_${step}.csv ${min}-${max}-${step} -c 0 -l 0 --sensitive -v

min=40
max=61
/home/dkoslicki/anaconda3/envs/CMash/bin/python ../../scripts/StreamingQueryDNADatabase.py ${inputFile} TrainingDatabase_k_61.h5 out_${min}_${max}_${step}.csv ${min}-${max}-${step} -c 0 -l 0 --sensitive -v


# sort them all for easy comparison
ls out* | xargs -I{} sh -c 'sort -o {} {}'