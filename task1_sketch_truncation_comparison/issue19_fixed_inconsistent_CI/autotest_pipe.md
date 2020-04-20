# Auto-test pipe for issue1

In case we may need to modify both the python script "MakeStreamingDNADatabase.py" and "StreamingQueryDNADatabase.py". It might be better to run the auto-test pipe from generating trainingDB step, given that it's pretty fast.



### Pseudo code:

```bash
### Usage
bash ${pipe} -r ${range1},${range2},${range3}... -o ${out_dir}

### Parameters
range= start-end-gap #end is the max_k
filenames.txt #file that record full path of the 165 files
k_size=${end}

### Running
### for each input range: run the following step:
###### 1. creating ref DB
MakeStreamingDNADatabase.py -k ${k_size} -v filenames.txt TrainingDB_k${k_size}.h5
###### 2. get true CI
for file in filenames.txt
  StreamingQueryDNADatabase.py ${file} TrainingDB_k${k_size}.h5 true_CI_out_${file}.csv ${k_size}-${k_size}-1 -v -c 0 -l 0 --sensitive
  StreamingQueryDNADatabase.py ${file} TrainingDB_k${k_size}.h5 truncation_CI_out_${file}.csv ${range} -v -c 0 -l 0 --sensitive
END for
###### 3. merge results and compare
merge all true_CI_out_${file}.csv => merged_true_CI_${k_size}.csv
merge the column with k=${k_size} of truncation_CI_out_${file}.csv => merged_trunc_CI_${k_size}.csv
compare the 2 merged files (maybe a plot or matrix substraction)
```



### Real script

```bash
#script location: /gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/src/a4.0_CMash_auto_test_pipe.sh

#!/bin/bash

### read parameters
while getopts r:o:f:h opts
do case "$opts" in
r) input_range="$OPTARG";; #the input range for CMash, separated by comma
o) our_dir="$OPTARG";; #output dir for auto-test pipe
f) genome_files="$OPTARG";; #input genome files to compare by CMash, here we already have them
h) echo "
Testing CMash codes for true CI with different input ranges.

Usage: bash <pipe> -r <range1,range2,...> -o <out_dir>"
exit;;
[?]) echo "Use -h for help information"
exit;;
esac
done

# edit this var if chaning new files
genome_files="/gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/src/filenames.txt"

### Check ENVs and parameters
source /gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/src/a0_env.sh
echo "check ENVs"

###### check out_dir
if [ -z "$out_dir" ]; then
  echo "NO output dir specified, the results will be stored at current dir:"
  echo "$PWD"
  out_dir=${PWD}
else
  echo "The output files would be stored at:"
  echo " ${out_dir}"
fi
time_tag=`date +"%F_%H-%M"`
cd ${out_dir}
mkdir test_CMash_true_CI_${time_tag}
cd test_CMash_true_CI_${time_tag}
echo -e "bingo!\n"

###### check input range to make sure the max_k are same for all range
[ -z "$input_range" ] && echo "Please specify your testing range in the format start1-end1-gap1,start2-end2-gap2,..." && exit
range_num=`echo $input_range | awk -F"," '{print NF}'`
prev_maxk=`echo $input_range | cut -d"," -f 1 | cut -d"-" -f 2` #check if the max_k is same for all range in the running loop
date > running_record.log
echo "The genome files are stored at: " >> running_record.log
echo ${genome_files} >> running_record.log
echo "The input range to test are: " >> running_record.log
echo ${input_range} >> running_record.log
echo " " >> running_record.log
for ((i=1; i<=$range_num; i+=1));do
  cur_maxk=`echo $input_range | cut -d"," -f $i | cut -d"-" -f 2`
  if (($prev_maxk != $cur_maxk)); then
    echo "The max_k in input_range are NOT same, please check INPUT!!!!!!"
    exit 1
  fi
done


### Start to run data
echo -e " "
echo "Start running CMash"
k_size=${prev_maxk} # the max_k for true CI

###### build trainingDB
${ltime} ${python_exec} ${CMash_scripts}/MakeStreamingDNADatabase.py  -k ${k_size} -v ${genome_files} TrainingDB_k${k_size}.h5
mv temp_runLog  TrainingDB_k${k_size}.log

###### get CI from k-k-1 input
mkdir CI_from_${k_size}-${k_size}-1
cd CI_from_${k_size}-${k_size}-1
for file in `cat ${genome_files}`; do
  echo "processing $file"
  name=`echo ${file##*/}`
  ${ltime} ${python_exec} ${CMash_scripts}/StreamingQueryDNADatabase.py ${file} ../TrainingDB_k${k_size}.h5 true_CI_${k_size}-${k_size}-1_${name}_results.csv ${k_size}-${k_size}-1 -v -c 0 -l 0 --sensitive
  mv temp_runLog true_CI_${k_size}-${k_size}-1_${name}_results.log
done
# Collect results
Rscript ${pipe_path}/a4.1_merge_true_CI_kk1_results.R  && \
  mv temp_out.csv ../merged_${k_size}-${k_size}-1_results.csv
cd ..

###### run all ranges staring from END (so max_k can be covered by whatever gap used)
### to add new range manually:
### before this loop, add input_range, range_num, k_size, genome_files
### and source the a0 env

for ((i=1; i<=$range_num; i+=1));do
  ### increase start point by 1 if START=1 (some running error while truncation to 1
  temp_range=`echo $input_range | cut -d"," -f $i | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  ### make sure all ranges will cover the max_k
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
  temp_range=${r_adj_start}-${r_end}-${r_gap}
  ### run the loop
  echo "running range for ${i}th range is $temp_range" >> running_record.log
  mkdir CI_from_${temp_range}
  cd CI_from_${temp_range}
  for file in `cat ${genome_files}`; do
    echo "processing $file"
    name=`echo ${file##*/}`
    ${ltime} ${python_exec} ${CMash_scripts}/StreamingQueryDNADatabase.py ${file} ../TrainingDB_k${k_size}.h5 truncation_CI_${temp_range}_${name}_results.csv ${temp_range} -v -c 0 -l 0 --sensitive
    mv temp_runLog truncation_CI_${temp_range}_${name}_results.log
  done
  Rscript ${pipe_path}/a4.2_merge_true_CI_from_truncation_results.R && \
    mv temp_out.csv ../merged_${temp_range}_true_CI_results.csv
  cd ..
done
```



### The merging code in R script a4.1 and a4.2

```R
library(tidyverse)

### a4.1
temp_f2 <- list.files(pattern=paste0("^true_CI_",".*_genomic.fna.gz_results.csv$"))
temp_list <- vector(mode='list')
  for (input_f in temp_f2) {
    temp_list[[input_f]] <- read.csv(input_f, as.is=T)
    rm(input_f)
} #put all files together

for (i in seq(temp_list)) {
    t_name=names(temp_list)[i]
    t_name=sub("true_CI_.*-1_", "", t_name)
    t_name=sub("_genomic.fna.gz_results.csv$", "", t_name)
    colnames(temp_list[[i]]) <- c("X", t_name)
    rm(i, t_name)
} #rename each df with the file name

merged_temp <- temp_list %>% reduce(full_join, by="X")
merged_temp$X <- sub("_genomic.fna.gz$", "", merged_temp$X)
merged_temp <- merged_temp[order(merged_temp$X), ]
row.names(merged_temp)<-merged_temp$X
merged_temp <- merged_temp[, row.names(merged_temp)]
write.csv(merged_temp, file="temp_out.csv")


### a4.2
all_files1=list.files(pattern="*_results.csv")
raw_tb <- vector(mode='list')
for ( input_f in all_files1 ){
  raw_tb[[input_f]]=read.csv(input_f, as.is=T)
  rm(input_f)
}

temp_list <- lapply(raw_tb, "[", c("X", tail(colnames(raw_tb[[1]]), n=1)))

for (i in seq(temp_list)) {
        t_name=names(temp_list)[i]
        t_name=sub("^truncation_CI_[0-9]+-[0-9]+-[0-9]+_", "", t_name)
        t_name=sub("_genomic.fna.gz_results.csv$", "", t_name)
        colnames(temp_list[[i]]) <- c("X", t_name)
        rm(i, t_name)
}

merged_temp <- temp_list %>% reduce(full_join, by="X")
merged_temp$X <- sub("_genomic.fna.gz$", "", merged_temp$X)
merged_temp <- merged_temp[order(merged_temp$X), ]
row.names(merged_temp)<-merged_temp$X
merged_temp <- merged_temp[, row.names(merged_temp)]
write.csv(merged_temp, file="temp_out.csv")
```

