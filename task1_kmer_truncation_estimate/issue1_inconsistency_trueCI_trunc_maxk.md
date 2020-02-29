# Issue1: inconsistent CI output (github issue #19)

### Issue statement: 

The calculated true CI are different when input parameter is 61-61-1 vs 4-61-3. Namely,

```bash
code1: StreamingQueryDNADatabase.py 61-61-1 -v -c 0 -l 0 --sensitive
code2: StreamingQueryDNADatabase.py 4-61-3 -v -c 0 -l 0 --sensitive
2 codes above generate different CI for k=61.
```

Possible reason:

```bash
Bloom filter size depends on the number of ks in the parameter, this may cause the difference.
```



### Data location

Current output for CMash test of k= 4-61-3:

```bash
/storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step3_results_summary
```

1. col and row of all tables are pre-sorted in same order, so each table is in symmatric design.
2. each column is one output file from "StreamingQueryDNADatabase.py", so $f_{*j}$ means the CI calculated for file j vs all other files.



### Merging code for results in step1,2 

```R
library(tidyverse)
out1_dir="/storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step1_true_CI_k61_and_truncation_errors/finished"
out2_dir="/gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step2_true_CI_for_all_short_k_mers/finished"
out3_dir="/storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step3_results_summary"

# merge step1: truncation output
setwd(out1_dir)
### store all files in a list
all_files1=list.files(pattern="*_results.csv")
print("The total files of in put is: ")
print(length(all_files1))
raw_tb <- vector(mode='list')
for ( input_f in all_files1 ){
  print(paste0("processing ", input_f))
  raw_tb[[input_f]]=read.csv(input_f, as.is=T)
  rm(input_f)
}
### merge results for each k
for (k in seq(4,61,3)) {
  temp_list <- lapply(raw_tb, "[", c("X", paste0("k.", k)))
  for (i in seq(temp_list)) {
    t_name=names(temp_list)[i]
    t_name=sub("^truncation_", "", t_name)
    t_name=sub("_genomic.fna.gz_results.csv$", "", t_name)
    colnames(temp_list[[i]]) <- c("X", t_name)
    rm(i, t_name)
  }
  merged_temp <- temp_list %>% reduce(full_join, by="X")
  merged_temp$X <- sub("_genomic.fna.gz$", "", merged_temp$X)
  merged_temp <- merged_temp[order(merged_temp$X), ]
  row.names(merged_temp)<-merged_temp$X
  merged_temp <- merged_temp[, row.names(merged_temp)]
  write.csv(merged_temp, file=paste0(out3_dir,"/merged_truncated_CI_k",k,".csv"))
  print(paste0("Finishing merging files for k=",k))
  rm(k, temp_list, merged_temp)
}

# merge step2: true CI output
setwd(out2_dir)
for (k in seq(4,61,3)) {
  temp_f2 <- list.files(pattern=paste0("^true_CI_k",k,".*_genomic.fna.gz_results.csv$"))
  temp_list <- vector(mode='list')
  for (input_f in temp_f2) {
    temp_list[[input_f]] <- read.csv(input_f, as.is=T)
    rm(input_f)
  } #put all files together
  for (i in seq(temp_list)) {
    t_name=names(temp_list)[i]
    t_name=sub(paste0("^true_CI_k",k,"_"), "", t_name)
    t_name=sub("_genomic.fna.gz_results.csv$", "", t_name)
    colnames(temp_list[[i]]) <- c("X", t_name)
    rm(i, t_name)
  } #rename each df with the file name
  merged_temp <- temp_list %>% reduce(full_join, by="X")
  merged_temp$X <- sub("_genomic.fna.gz$", "", merged_temp$X)
  merged_temp <- merged_temp[order(merged_temp$X), ]
  row.names(merged_temp)<-merged_temp$X
  merged_temp <- merged_temp[, row.names(merged_temp)]
  write.csv(merged_temp, file=paste0(out3_dir,"/true_CI_k",k,".csv"))
  print(paste0("Finishing merging files for k=",k))
  rm(k, temp_list, merged_temp, temp_f2)
}
```



### Follow up:

1. setup auto-test pipe for CMash results
2. try to track the code "StreamingQueryDNADatabase.py" and check the results consistency by auto-test pipe

