#!/bin/bash

### read parameters
while getopts k:r:s:n:o:c:h opts
do case "$opts" in
	k) maxk="$OPTARG";; # one large k value to check CI
	r) ref="$OPTARG";; # candidate database path
	s) shuf_seed="$OPTARG";; # seed used in shuffle cmd
	n) num="$OPTARG";; # number of genomes picked
	o) out_dir="$OPTARG";; # output dir
	c) cmash_dir="$OPTARG";; # the CMash locol copy to use
	h) echo "
Description:
Randomly pick n samples from the NCBI bacteria database and build the CI matrix.

Usage: bash <pipe> 

Optional parameters:
-k <maxk>: the k-mer length, default 61
-s <seed>: the seed in shuffle cmd, default 100
-n <num>: number of ref genomes to pick, default 100
-r <ref_path>: location where the refs are stored, default in Shaopeng's folder
-o <out_dir>: output dir, default current location
-c <cmash_dir>: locol copy of CMash to use, default the master branch
"
exit;;
[?]) echo "Use -h for help information"
exit;;
esac
done



### check parameters
# local variables in CSE server
conda_path="/data/sml6467/software/miniconda3"
python_exec="python3.8"
. ${conda_path}/etc/profile.d/conda.sh
conda activate CMash-env
threads=48
num_hashes=2000

# input parameters
if [ -z "$maxk" ]; then
	maxk=61
fi

if [ -z "$ref" ]; then
	ref="/data/sml6467/resource/NCBI_Genbank_bacteria_genome_20200525/NCBI_Genbank_bacteria_genome/genome_files"
fi

if [ -z "$shuf_seed" ]; then
	shuf_seed=100
fi

if [ -z "$num" ]; then
	num=100   # use 10 for test, then change to 100
fi

if [ -z "$out_dir" ]; then
	out_dir=$PWD
fi

if [ -z "$cmash_dir" ]; then
       cmash_dir="/data/sml6467/github/CMash_master"
fi

time_tag=$(date +"%H-%M_%m_%d_%y")
mkdir output_${time_tag}
cd output_${time_tag}
out_dir=$PWD

date > running_log_${time_tag}.log
echo -e "\n
start running
maxk:	$maxk
seed:	$shuf_seed
ref:	$ref
num:	$num
Cmash:	$cmash_dir
out_dir:\n$out_dir\n"  >> running_log_${time_tag}.log
cat running_log_${time_tag}.log



### random sample
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# store absolute path of selected files
find ${ref} -name "*.fna.gz" | shuf --random-source=<(get_seeded_random $shuf_seed) | head -${num} > sampled_files.txt



### ran CMash to get estimated CI
# setup 
ref_file="sampled_files.txt"
query_file="sampled_files.txt"
export PYTHONPATH=${cmash_dir}:$PYTHONPATH
testFile=$(python -c "from CMash import MinHash as MH; print(MH.__file__)")
correctFile=${cmash_dir}/CMash/MinHash.py
[ "$testFile" == "$correctFile" ] && echo "Files are correct" || exit 1
scriptDir=${cmash_dir}/scripts
moduleDir=${cmash_dir}/CMash
ltime="/usr/bin/time -av -o temp_runLog"
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 

# build training DB
${ltime} ${python_exec} ${scriptDir}/MakeStreamingDNADatabase.py -k ${maxk} -v -t ${threads} -n ${num_hashes} ${ref_file} TrainingDB_k${maxk}.h5 
mv temp_runLog TrainingDB_k${maxk}.log 2> /dev/null

# get estimated CI
mkdir estimated_CI
cd estimated_CI
for file in `cat ../${query_file}`
do
	echo $file
	name=`echo ${file##*/}`
	${ltime} ${python_exec} ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../TrainingDB_k${maxk}.h5 estimated_CI_k${maxk}_${name}_results.csv ${maxk}-${maxk}-1 -v -c 0 -l 0 -t ${threads} --sensitive
	mv temp_runLog estimated_CI_k${maxk}_${name}_results.log
done

# merge results
sed '1d' $(ls *_results.csv | head -1) | sort -k1,1V | cut -d"," -f 1 | cat <(echo "file") - > _temp_0row.txt
for file in `ls *_results.csv`
do
	name=`echo ${file##estimated_CI_k${maxk}_}`
	name=`echo ${name%_results.csv}`
	sed '1d' $file | sort -k1,1V | cut -d"," -f 2 | cat <(echo $name) - > _temp_${name}.txt 
done

paste _temp_* > merged_CI_table.txt
mv merged_CI_table.txt ../
rm _temp_*
cd ..
conda deactivate

echo "pipe done"
date
