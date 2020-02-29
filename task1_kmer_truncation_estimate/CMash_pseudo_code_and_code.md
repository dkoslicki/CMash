# CMash running task1

### Introduction: 

[CMash record](https://github.com/ShaopengLiu1/PSU_Bioinformatics/blob/master/r3_Koslicki_group/CMash_record.md)



### Revision note:

- Update:
  - Step1 seems in trouble (overflow or due to lack of MEM or data is too big?)
  - Step2 should be ready soon (k61 results for demo)
- Q:
  - [x] MakeStreamingDNADatabase.py will generate different database when repeat (md5 check dif)? The TST is same.
  - [x] Errors: step1 record in folder
    - Code was manually tested for smaller range on single filer (bingo)
    - For full data, it already run overnight (>5h on ICS) but no results generated

### Pseudo-code

```bash
### CI: CMash Index

### i) find true CI with k=61 + all variation data on each k-size point
max_k=61
MakeStreamingDNADatabase.py <all_genome> <max_k>  -> TB_61.hs
for gi in <all_genome>:
  StreamingqueryDNADatabase.py	gi	TB_61.hs	1-61-3	-c 0 -l 0 --sensitive
# In the output matrix: col[k=61] is the true value of each gi; col[k<61] is the variation due to truncation
  
### ii) find the true CI with all k<61
for max_k in range (1,58,3):
  MakeStreamingDNADatabase.py <all_genome> <max_k>   -> TB_${max_k}.hs
  for gi in <all_genome>:
    StreamingqueryDNADatabase.py	gi	TB_${max_k}.hs	${max_k}-${max_k}-1	-c 0 -l 0 --sensitive
# there will be 19 output, each of them is a true CI at k=${max_k}

### iii) clustering
run cluster_matrix(A_eps, A_indicies, cluster_eps=.01) on the TrainingDatabase_{n}_k_60.h5
```



### Real script:

- Env check

```bash
### Adjust env check: PYTHONPATH specified to absolute path
export PYTHONPATH="/storage/home/sml6467/scratch/tools/CMash_github_master/CMash":$PYTHONPATH  #similar to $PATH, guide the dir where to look up executives
testFile=$(python -c "from CMash import MinHash as MH; print(MH.__file__)")
correctFile="/storage/home/sml6467/scratch/tools/CMash_github_master/CMash/CMash/MinHash.py"
[ "$testFile" == "$correctFile" ] && echo "Files are correct" || exit 1
```



- Step1: k61

```bash
### Step1: k61 true CI and all truncation records
# filenames.txt records all the absolute path of the 165 genome files
python_exec="/opt/aci/sw/python/3.6.3_anaconda-5.0.1/bin/python"
CMash_scripts="/storage/home/sml6467/scratch/tools/CMash_github_master/CMash/scripts"
ltime="/usr/bin/time -av -o temp_runLog"
ml python/3.6.3-anaconda5.0.1
cd /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step1_true_CI_k61_and_truncation_errors
# generate the TrainingDB for k=61
${ltime} ${python_exec} ${CMash_scripts}/MakeStreamingDNADatabase.py  -k ${k_size} -v filenames.txt TrainingDB_k${k_size}.h5
mv temp_runLog  TrainingDB_k${k_size}.log
# generate the trunction output with the TDB_k61.h5
for file in `cat filenames.txt`
do
  echo "processing $file"
  name=`echo ${file##*/}`
  ${ltime} ${python_exec} ${CMash_scripts}/StreamingQueryDNADatabase.py ${file} TrainingDB_k${k_size}.h5 truncation_${name}_results.csv 1-61-3 -v -c 0 -l 0 --sensitive
  mv temp_runLog truncation_${name}_results.log  2> /dev/null
  unset name file
done
```

- Step2: k= seq 1 3 61

```bash
### Step2: true CI for all k<61
all_files="/gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step1_true_CI_k61_and_truncation_errors/filenames.txt"
# loop through
cd /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step2_true_CI_for_all_short_k_mers
for ((k_size=61; k_size>=1; k_size-=3))
do
  ${ltime} ${python_exec} ${CMash_scripts}/MakeStreamingDNADatabase.py  -k ${k_size} -v ${all_files} TrainingDB_k${k_size}.h5
  mv temp_runLog  TrainingDB_k${k_size}.log
  for file in `cat ${all_files}`
  do
    echo "processing $file"
    name=`echo ${file##*/}`
    ${ltime} ${python_exec} ${CMash_scripts}/StreamingQueryDNADatabase.py ${file} TrainingDB_k${k_size}.h5 true_CI_k${k_size}_${name}_results.csv ${k_size}-${k_size}-1 -v -c 0 -l 0 --sensitive
    mv temp_runLog true_CI_k${k_size}_results.log 
    unset file name
  done
done 
```

- Step3: cluster for ref

```bash
### Step3: (waiting for step1 output)
```





