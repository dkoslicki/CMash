# Test the CMash after been fixed on 20200319

### Conclusion:

1. Bingo now the results are consist : )
2. There are some unexpected NAs.
3. The local autotest pipe is ready for fixed code, so does my own branch.



### Results:

1. The CMash index (CI) matrix from different k_range are **EXACT same** except there are some NAs in different location. (So I don't put the plot here as it's not necessary)
2. However, the missing values are unexpected. Issue description:
   - Based on several manual check, missing values happen on low CI pairs. For example: taxid_1909294_134 vs taxid_1909294_219 is (NA and 0.0006, with one of them act as query correspondingly).
   - As I added "-c 0" parameter, it's supposed to report all 165 records. This can be confirmed that there are many records with <0.1 (the default) CI
   - based on the number of NAs for each k_range, it looks like this is due to output strategy (seems CMash will omit the record if the whole row are NAs)



### Next step:

- [ ] Come back to the original design to get the results of true CI vs truncated CI
- [ ] maybe explore why there are missing values



### Running code:

1. Run the new CMash code by auto-test pipe

```
# the CMash code on ICS has been merged with your new update
bash /gpfs/group/dmk333/default/shaopeng/projects/202002_CMash_test/src/a4.0_CMash_auto_test_pipe.sh -r 51-61-5,53-61-2,53-61-4,55-61-2,55-61-3,56-61-5,57-61-2,57-61-4,58-61-3,59-61-2
```

Output files md5sum:

(The differences comes from missing values, but not the CI)

```bash
459e45d13084dda1c2be04b7d3b3b29e  merged_51-61-5_true_CI_results.csv
459e45d13084dda1c2be04b7d3b3b29e  merged_53-61-2_true_CI_results.csv
459e45d13084dda1c2be04b7d3b3b29e  merged_53-61-4_true_CI_results.csv
459e45d13084dda1c2be04b7d3b3b29e  merged_55-61-2_true_CI_results.csv
459e45d13084dda1c2be04b7d3b3b29e  merged_55-61-3_true_CI_results.csv
9d3a5a831e8d1edb783d5408b172931c  merged_56-61-5_true_CI_results.csv
9d3a5a831e8d1edb783d5408b172931c  merged_57-61-2_true_CI_results.csv
9d3a5a831e8d1edb783d5408b172931c  merged_57-61-4_true_CI_results.csv
eb22f189bc8fce79196f900dd37db642  merged_58-61-3_true_CI_results.csv
7966f880c85b181f926b0bed81d2cf66  merged_59-61-2_true_CI_results.csv
7966f880c85b181f926b0bed81d2cf66  merged_61-61-1_results.csv
```



2. Checking results by R

```R
full_file=read.csv("merged_61-61-1_results.csv", row.names=1)

cutoff=1e-07

for (temp_file in list.files(pattern=".*results.csv")) {
  print(paste0("Processing ", temp_file))
  temp_csv <- read.csv(temp_file, row.names=1)
  temp_dif <- sum(full_file-temp_csv > cutoff, na.rm=T)
  temp_exact <- sum(full_file-temp_csv == 0, na.rm=T)
  print(paste0("The number of cells with changes more than ",cutoff," is ", temp_dif))
  print(paste0("The number of exact matched cells is ",temp_exact))
}
```

Output:

```
[1] "Processing merged_51-61-5_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_53-61-2_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_53-61-4_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_55-61-2_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_55-61-3_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_56-61-5_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_57-61-2_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_57-61-4_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_58-61-3_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_59-61-2_true_CI_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
[1] "Processing merged_61-61-1_results.csv"
[1] "The number of cells with changes more than 1e-07 is 0"
[1] "The number of exact matched cells is 27192"
```

Briefly, 

1. all cells (if they are not missing) in the CI matrix are same.
2. but the place where missing values located are slightly different



3. List number of NAs

```R
for (temp_file in list.files(pattern=".*results.csv")) {
  print(paste0("Processing ", temp_file))
  temp_csv <- read.csv(temp_file, row.names=1)
  print(sum(is.na(temp_csv)))
}
```

Output

```R
[1] "Processing merged_51-61-5_true_CI_results.csv"
[1] 0
[1] "Processing merged_53-61-2_true_CI_results.csv"
[1] 0
[1] "Processing merged_53-61-4_true_CI_results.csv"
[1] 0
[1] "Processing merged_55-61-2_true_CI_results.csv"
[1] 0
[1] "Processing merged_55-61-3_true_CI_results.csv"
[1] 0
[1] "Processing merged_56-61-5_true_CI_results.csv"
[1] 22
[1] "Processing merged_57-61-2_true_CI_results.csv"
[1] 22
[1] "Processing merged_57-61-4_true_CI_results.csv"
[1] 22
[1] "Processing merged_58-61-3_true_CI_results.csv"
[1] 30
[1] "Processing merged_59-61-2_true_CI_results.csv"
[1] 33
[1] "Processing merged_61-61-1_results.csv"
[1] 33
```

