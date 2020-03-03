# Issue1: inconsistent CI output (github issue #19)

### Issue statement: 

The calculated true CI are different when input ranges are different (i.e. 61-61-1 vs 4-61-3). Namely,

```bash
code1: StreamingQueryDNADatabase.py 61-61-1 -v -c 0 -l 0 --sensitive
code2: StreamingQueryDNADatabase.py 4-61-3 -v -c 0 -l 0 --sensitive
2 codes above generate different CI for k=61.
```

Possible reason:

```bash
Bloom filter size depends on the number of ks in the parameter, this may cause the difference.
```



---

### Progress:

1. auto-test pipe is done ready, [click here](https://github.com/dkoslicki/CMash/blob/shaopeng/task1_kmer_truncation_estimate/issue19_inconsistent_CI/autotest_pipe.md) for more details

2. test multiple ranges with max_k=20, results:

   - [corrplot of 20-20-1 vs 12-20-4 / 16-20-2 / 14-20-3 / 10-20-5 vs 8-20-2](https://drive.google.com/open?id=1cM3X06e9OEHySM3SEBRD7QoQq1quQb0W)

   - with md5sum check, all results except **20-20-1** are same

   ```bash
   61cc414ab45330d610eaeb9c6702ff07  merged_10-20-5_true_CI_results.csv
   61cc414ab45330d610eaeb9c6702ff07  merged_12-20-4_true_CI_results.csv
   61cc414ab45330d610eaeb9c6702ff07  merged_14-20-3_true_CI_results.csv
   61cc414ab45330d610eaeb9c6702ff07  merged_16-20-2_true_CI_results.csv
   0af72845361fa00a6786c7ace9a51f7b  merged_20-20-1_results.csv
   61cc414ab45330d610eaeb9c6702ff07  merged_8-20-2_true_CI_results.csv
   ```

   - conclusion:

   ```bash
   #1: k-k-1 results are different with other start-end-gap results
   #2: ALL other start-end-gap results are identical
   #3: possible reason does make sense, but when the size get high it won't influence the results, we may need more details for the affect of size.
   ```

   

3. repeat multiple range test with max_k=40 (running):





---

### Follow up:

1. Use a higher k (like 40), repeat the finding:
   - [ ] use different range size, 40-40-1(1), 38-40-2(2), 36-40-2(3), ...... to see if it performs the same
2. try to track the code "StreamingQueryDNADatabase.py" and check the results consistency by auto-test pipe



---

### Data location

###### Note: 

1. col and row of all tables are pre-sorted in same order, so each table is in symmatric design.
2. each column is one output file from "StreamingQueryDNADatabase.py", so $f_{*j}$ means the CI calculated for file j vs all other files.

###### Data:

1. range= 61-61-1 vs 4-61-1 (CMash manual processing data)

   ```bash
   /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200217_CMash_task1_trunction_kmer/step3_results_summary
   ```

2. range= 20-20-1, 12-20-4, 16-20-2, 14-20-3 (auto-test pipe)

   ```bash
   /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/20200302_auto_test_pipe/test_CMash_true_CI_2020-03-02_13-05
   ```

   



