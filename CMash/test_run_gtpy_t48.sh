# Manual pipe: check time usage for GT.py with all threads limited to 48 (only works for 48)
# The CMash repo is switched to "Update_GT_thread" branch!!!!!! Otherwise the GT.py will fail due to -t parameter

date
cd /data/sml6467/projects/202002_CMash_test/results/20200719_conda_test/CMASH-reproducibles/task1_k-mer_truncation/example/high_depth_test/Output_CMash_BBMap_6_file
CMash="/data/sml6467/github/CMash_master"
conda activate CMash-env
maxk=60
temp_range="10-60-5"
threads=48
num_hashes=2000
ltime="/usr/bin/time -av -o temp_runLog"
export PYTHONPATH=${CMash}:$PYTHONPATH
scriptDir=${CMash}/scripts
moduleDir=${CMash}/CMash


query="/data/sml6467/projects/202002_CMash_test/results/20200719_conda_test/CMASH-reproducibles/task1_k-mer_truncation/example/high_depth_test/all_BB_path.txt"
file="/data/sml6467/projects/202002_CMash_test/results/20200719_conda_test/CMASH-reproducibles/task1_k-mer_truncation/example/output_07_20_02-47/BBMap_simu/BBMap_simulated_meta_18x.fq"


### single test 
# mkdir single_test
# cd single_test/
# for i in $(seq 10 5 60); do
# 	${ltime} python ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../ref_db_10-60-5/TrainingDB_k${i}.h5 est_CI_k${i}_BB18x.csv ${i}-${i}-1 -v -c 0 -l 0 -t ${threads} --sensitive
# 	mv temp_runLog  est_CI_k${i}_BB18x.log
# 	${ltime} python ${moduleDir}/GroundTruth.py ${file} ../ref_db_10-60-5/TrainingDB_k${i}.h5 GT_CI_k${i}_BB18x.csv ${i}-${i}-1 -c 0 -l 0 -t ${threads}
# 	mv  temp_runLog GT_CI_k${i}_BB18x.log
# done


### re run all analysis 
mkdir rerun_all_w_GT_t48
cd rerun_all_w_GT_t48
# truncated CI
for file in `cat ${query}`; do
	echo $file
	name=`echo ${file##*/}`
	${ltime} python ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../ref_db_${temp_range}/TrainingDB_k${maxk}.h5 truncation_${name}_results.csv ${temp_range} -v -c 0 -l 0 -t ${threads}  --sensitive
	mv temp_runLog truncation_${name}_results.log
done
mkdir truncated_CI
mv truncation_* ./truncated_CI

# est and GT CI
for i in $(seq 10 5 60); do
	for file in `cat ${query}`; do
		echo $file
		name=`echo ${file##*/}`
		${ltime} python ${scriptDir}/StreamingQueryDNADatabase.py ${file} ../ref_db_${temp_range}/TrainingDB_k${i}.h5 estimated_CI_k${i}_${name}_results.csv ${i}-${i}-1 -v -c 0 -l 0 -t ${threads} --sensitive
		mv temp_runLog estimated_CI_k${i}_${name}_results.log
		${ltime} python ${moduleDir}/GroundTruth.py ${file} ../ref_db_${temp_range}/TrainingDB_k${i}.h5 ground_truth_CI_k${i}_${name}_results.csv ${i}-${i}-1 -c 0 -l 0 -t 48
		mv temp_runLog ground_truth_CI_k${i}_${name}_results.log
	done
done

mkdir ground_truth_CI
mv ground_truth_CI_k* ground_truth_CI
mkdir estimated_CI
mv estimated_CI_k* estimated_CI



### copy the time grep code from the previous pipeline 
function cal_sec() {
        grep_word=$1
        file=$2
        temp_time=$(grep ${grep_word} $file | cut -d")" -f 3 | cut -d" " -f 2)
        out_sec=$(echo ${temp_time} | awk -F: '{ if (NF == 1) {print $NF} else if (NF == 2) {print $1 * 60 + $2} else if (NF==3) {print $1 * 3600 + $2 * 60 + $3} }')
        echo ${out_sec}
}


cd truncated_CI
grep_word="System"
[ -f time_truncated.txt ] && rm time_truncated.txt
for file in $(ls truncation_*_results.log); do
        name=$(echo ${file#truncation_})
        name=$(echo ${name%_results.log})
        kk_time=$(grep ${grep_word} $file | cut -d":" -f 2 | cut -d" " -f 2)
        echo -e "$name\t$kk_time" >> time_truncated.txt
done
mv time_truncated.txt ../Sys_truncated.txt

grep_word="User"
for file in $(ls truncation_*_results.log); do
        name=$(echo ${file#truncation_})
        name=$(echo ${name%_results.log})
        kk_time=$(grep ${grep_word} $file | cut -d":" -f 2 | cut -d" " -f 2)
        echo -e "$name\t$kk_time" >> time_truncated.txt
done
mv time_truncated.txt ../User_truncated.txt

grep_word="Elapsed"
for file in $(ls truncation_*_results.log); do
        name=$(echo ${file#truncation_})
        name=$(echo ${name%_results.log})
        kk_time=$(cal_sec $grep_word $file)
        echo -e "$name\t$kk_time" >> time_truncated.txt
done
mv time_truncated.txt ../Elapsed_truncated.txt




cd ../estimated_CI
grep_word="System"
[ -f time_truncated.txt ] && rm time_truncated.txt
for name in $(cut -f 1 ../Sys_truncated.txt); do
	[ -f time_Est_${name}.txt ] && rm time_Est_${name}.txt
	for file in $(ls estimated_CI_k*_${name}_results.log); do
                kk=$(echo ${file#estimated_CI_k})
                kk=$(echo ${kk%_${name}_results.log})
                kk_time=$(grep ${grep_word} $file | cut -d":" -f 2 | cut -d" " -f 2)
                echo -e "$kk\t$kk_time" >> time_Est_${name}.txt
        done
        cat <(echo -e "k_size\ttime_sec")  <(sort -k1,1n time_Est_${name}.txt) > temp_time && mv temp_time time_Est_${name}.txt && mv time_Est_${name}.txt ../Sys_Est_${name}.txt
done

grep_word="User"
for name in $(cut -f 1 ../Sys_truncated.txt); do
        [ -f time_Est_${name}.txt ] && rm time_Est_${name}.txt
        for file in $(ls estimated_CI_k*_${name}_results.log); do
                kk=$(echo ${file#estimated_CI_k})
                kk=$(echo ${kk%_${name}_results.log})
                kk_time=$(grep ${grep_word} $file | cut -d":" -f 2 | cut -d" " -f 2)
                echo -e "$kk\t$kk_time" >> time_Est_${name}.txt
        done
        cat <(echo -e "k_size\ttime_sec")  <(sort -k1,1n time_Est_${name}.txt) > temp_time && mv temp_time time_Est_${name}.txt && mv time_Est_${name}.txt ../User_Est_${name}.txt
done

grep_word="Elapsed"
for name in $(cut -f 1 ../Sys_truncated.txt); do
        [ -f time_Est_${name}.txt ] && rm time_Est_${name}.txt
        for file in $(ls estimated_CI_k*_${name}_results.log); do
                kk=$(echo ${file#estimated_CI_k})
                kk=$(echo ${kk%_${name}_results.log})
                kk_time=$(cal_sec $grep_word $file)
                echo -e "$kk\t$kk_time" >> time_Est_${name}.txt
        done
        cat <(echo -e "k_size\ttime_sec")  <(sort -k1,1n time_Est_${name}.txt) > temp_time && mv temp_time time_Est_${name}.txt && mv time_Est_${name}.txt ../Elapsed_Est_${name}.txt
done




cd ../ground_truth_CI
grep_word="System"
for name in $(cut -f 1 ../Sys_truncated.txt); do
	[ -f time_GT_${name}.txt ] && rm time_GT_${name}.txt
	for file in $(ls ground_truth_CI_k*_${name}_results.log); do
                kk=$(echo ${file#ground_truth_CI_k})
                kk=$(echo ${kk%_${name}_results.log})
                kk_time=$(grep ${grep_word} $file | cut -d":" -f 2 | cut -d" " -f 2)
                echo -e "$kk\t$kk_time" >> time_GT_${name}.txt
        done
        cat <(echo -e "k_size\ttime_sec")  <(sort -k1,1n time_GT_${name}.txt) > temp_time && mv temp_time time_GT_${name}.txt && mv time_GT_${name}.txt ../Sys_GT_${name}.txt
done

grep_word="User"
for name in $(cut -f 1 ../Sys_truncated.txt); do
        [ -f time_GT_${name}.txt ] && rm time_GT_${name}.txt
	for file in $(ls ground_truth_CI_k*_${name}_results.log); do
                kk=$(echo ${file#ground_truth_CI_k})
                kk=$(echo ${kk%_${name}_results.log})
                kk_time=$(grep ${grep_word} $file | cut -d":" -f 2 | cut -d" " -f 2)
                echo -e "$kk\t$kk_time" >> time_GT_${name}.txt
        done
        cat <(echo -e "k_size\ttime_sec")  <(sort -k1,1n time_GT_${name}.txt) > temp_time && mv temp_time time_GT_${name}.txt && mv time_GT_${name}.txt ../User_GT_${name}.txt
done



grep_word="Elapsed"
for name in $(cut -f 1 ../Sys_truncated.txt); do
	[ -f time_GT_${name}.txt ] && rm time_GT_${name}.txt
        for file in $(ls ground_truth_CI_k*_${name}_results.log); do
                kk=$(echo ${file#ground_truth_CI_k})
                kk=$(echo ${kk%_${name}_results.log})
                kk_time=$(cal_sec $grep_word $file)
                echo -e "$kk\t$kk_time" >> time_GT_${name}.txt
        done
        cat <(echo -e "k_size\ttime_sec")  <(sort -k1,1n time_GT_${name}.txt) > temp_time && mv temp_time time_GT_${name}.txt && mv time_GT_${name}.txt ../Elapsed_GT_${name}.txt
done
conda deactivate

### merge all results
cd ..
mkdir time_summary
mv *txt time_summary
cd time_summary

function merge_tb() {
	prefix=$1
	echo "the prefix used is $prefix"
	cut -f 1 $(ls ${prefix}*.txt | head -1) > _start.txt
	for file in $(ls ${prefix}*.txt); do
		echo $file
		name=$(echo ${file#${prefix}})
		name=$(echo ${name%.txt})
		cut -f 2 $file | awk -v out_name=${name} '{if (NR==1){print out_name;} else {print $0} }' > _temp_${name}.txt
	done
	paste _start.txt _temp_*.txt > merged_${prefix}time.txt
	rm _start.txt _temp_*.txt
}

for prefix in Elapsed_Est_ Elapsed_GT_ Sys_Est_ Sys_GT_ User_Est_ User_GT_; do
	merge_tb $prefix
done


echo "pipe done"
date















