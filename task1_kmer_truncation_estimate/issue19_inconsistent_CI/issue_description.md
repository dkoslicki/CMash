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
No idea now.
```

Current update:

```bash
# Computationally confirmed:
1. the main problem is the difference between max-max-1 and start-max-gap, e.g. with maxk=61, the difference can be 0.6 vs 0.4 for closely related species.
2. while there are tiny differences among different start-max-gap combinations (with same maxk), they look like normal flunctuations due to random sampling process and can be neglectable.
3. from the py code, the len(range_k) would only affect a) BF size (164); b) location to apply threshold (123) which manually set to 0; c) output df (320). So the BF size should NOT be an issue.
```



Test pipe:

```bash
# usage:
bash <pipe> -h #for help information
bash <pipe> -r range1,range2,...  #it will automatically run max-max-1, so don't need to include it in the ranges

Auto test pipe with original CMash:
/gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/src/a4.0_CMash_auto_test_pipe.sh

Auto test pipe for "modified" CMash:
# this's a local copy of CMash/shaopeng branch
# I edited the "StreamingQueryDNADatabase.py" file directly
/gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/src/a5.1_auto_test_pipe_for_edited_CMash_script.sh

```



---

### Follow up:

1. try to track the code "StreamingQueryDNADatabase.py" and check the results consistency by auto-test pipe

---

### Running record:

1. auto-test pipe is done ready, [click here](https://github.com/dkoslicki/CMash/blob/shaopeng/task1_kmer_truncation_estimate/issue19_inconsistent_CI/autotest_pipe.md) for more details

2. test multiple ranges with max_k=20, results:

   - location: /gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/results/20200302_auto_test_pipe/test_CMash_true_CI_2020-03-02_13-05

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

   - location:
     - 1st run: /gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/results/20200302_auto_test_pipe/test_CMash_true_CI_2020-03-04_00-58
     - 2nd run for validation: /gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/results/20200302_auto_test_pipe/test_CMash_true_CI_2020-03-09_15-53

   - md5 check output of 1st run:

   ```bash
   c24bca8164ab98ba8b54b658f5a24ee9  merged_26-40-2_true_CI_results.csv
   c24bca8164ab98ba8b54b658f5a24ee9  merged_28-40-2_true_CI_results.csv
   c24bca8164ab98ba8b54b658f5a24ee9  merged_30-40-2_true_CI_results.csv
   c24bca8164ab98ba8b54b658f5a24ee9  merged_30-40-5_true_CI_results.csv
   c24bca8164ab98ba8b54b658f5a24ee9  merged_32-40-2_true_CI_results.csv
   c24bca8164ab98ba8b54b658f5a24ee9  merged_34-40-2_true_CI_results.csv
   c24bca8164ab98ba8b54b658f5a24ee9  merged_36-40-2_true_CI_results.csv
   0a97ed393a5e4a5b2d494bdb744a93f3  merged_38-40-2_true_CI_results.csv
   96d6f0dfccbcfbca8ac30e7770bbbd70  merged_40-40-1_results.csv
   ```

   - summary:
     - All results are same except k_range size =1 or 2 
     - However, when checking the results, k_range size=2 is almost same to size > 2. There is only 1 record with <1% difference (no idea why only 1 record is different)
   - conclusion:

   ```bash
   #4: On top of k_range size, the max_k itself may also influence the results
   ```



4. try to edit BF size and rerun the process
   - location: /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_CMash_test/results/step5_20200308_try_editing_CMash_scripts
   - what I've tried,
     - a3: times another 10 in WriteBF function to make sure the BF size is large enough & larger than all previous running
     - a4: manually fix len(range_k) to 10 for BF size so that all BF are large enough & same & larger than previous running
   - output files:
     - https://drive.google.com/open?id=1KcIiWhull92j7yEPS-91p40rEIP8dYyt
   - results:
     - Increasing BF size does NOT influence the results
     - Fix BF size does NOT influence the results
     - The difference of results remains same with different BF size (len_krange)

---

### Note:

1. col and row of all tables are pre-sorted in same order, so each table is in symmatric design.
2. each column is one output file from "StreamingQueryDNADatabase.py", so $f_{*j}$ means the CI calculated for file j vs all other files.



