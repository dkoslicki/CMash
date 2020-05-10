#!/bin/bash
date
echo "pipe start"

### input parameters to be set
ref_file="/data/sml6467/projects/202002_CMash_test/src/filepath.txt"
query_file="/data/sml6467/projects/202002_CMash_test/src/querypath.txt"
input_range="5-60-5"
CMash="/data/sml6467/github/CMash_master"
conda_path="/data/sml6467/software/miniconda3"
python_exec="python3.8"
outdir="/data/sml6467/projects/202002_CMash_test/results/r1_metagenome_against_ref_165_species"
threads=48
num_hashes=2000

### parameter check
# check input range: 1) can't start from 1; 2) maxk should be covered
[ -z "$input_range" ] && echo "Please specify your testing range in the format start-end-gap" && exit 1
temp_range=`echo $input_range | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
temp_range=${r_adj_start}-${r_end}-${r_gap}
maxk=${r_end}
# check outdir exist
outdir=`realpath ${outdir}`
[ ! -d ${outdir} ] && mkdir -p ${outdir} 


### dependent paramebers (no need to touch)
# activate conda env inside bash
. ${conda_path}/etc/profile.d/conda.sh  
conda activate CMash-env
# simplify time cmd
ltime="/usr/bin/time -av -o temp_runLog"
# dir of pipe folder
pipe_path="$(dirname "`pwd`")"
summary_code="${pipe_path}/rep1_summary.py"
# test CMash python env match
export PYTHONPATH=${CMash}:$PYTHONPATH
testFile=$(python -c "from CMash import MinHash as MH; print(MH.__file__)")
correctFile=${CMash}/CMash/MinHash.py
[ "$testFile" == "$correctFile" ] && echo "Files are correct" || exit 1
scriptDir=${CMash}/scripts
moduleDir=${CMash}/CMash
time_tag=`date +"%m_%d_%H-%M"`


### Running the py code
cd ${outdir} 
mkdir output_${time_tag}
cd output_${time_tag}

# Step1, build the ref database for all ks in range
mkdir ref_db_${temp_range}
cd ref_db_${temp_range}
for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "generating reference for k=${i}"
	${ltime} ${python_exec} ${scriptDir}/MakeStreamingDNADatabase.py -k ${i} -v -t ${threads} -n ${num_hashes}  ${ref_file} TrainingDB_k${i}.h5
        mv temp_runLog TrainingDB_k${i}.log	2> /dev/null
done
cd ..


# Step2, truncated CI for all k-point
### temp adjust to test code
mkdir truncated_CI
cd truncated_CI
for file in `cat ${query_file}`
do
	echo $file
	name=`echo ${file##*/}`
	${ltime} ${python_exec} ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../ref_db_${temp_range}/TrainingDB_k${maxk}.h5 truncation_${name}_results.csv ${temp_range} -v -c 0 -l 0 -t ${threads}  --sensitive
	mv temp_runLog truncation_${name}_results.log
done
cd ..

# Step3, estimated CI and ground truth CI
for i in $(seq ${r_adj_start} ${r_gap} ${maxk})
do
	echo "running on k ${i}"
	for file in `cat ${query_file}`
	do
		echo $file
		name=`echo ${file##*/}`
		# estimated CI
		${ltime} ${python_exec} ${scriptDir}/StreamingQueryDNADatabase.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5 estimated_CI_k${i}_${name}_results.csv ${i}-${i}-1 -v -c 0 -l 0 -t ${threads} --sensitive 
		mv temp_runLog estimated_CI_k${i}_${name}_results.log
		# ground truth CI
		${ltime} ${python_exec} ${moduleDir}/GroundTruth.py ${file} ./ref_db_${temp_range}/TrainingDB_k${i}.h5  ground_truth_CI_k${i}_${name}_results.csv ${i}-${i}-1 -c 0 -l 0 
		mv temp_runLog ground_truth_CI_k${i}_${name}_results.log
	done
done

mkdir ground_truth_CI
mv ground_truth_CI_k* ground_truth_CI
mkdir estimated_CI
mv estimated_CI_k* estimated_CI


# Step3, summaize results by py
mkdir summary
${python_exec} ${summary_code} ${temp_range} ${query_file}

conda deactivate
date

