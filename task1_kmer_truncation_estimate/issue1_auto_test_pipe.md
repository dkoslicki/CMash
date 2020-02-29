# Auto-test pipe for issue1

In case we may need to modify both the python script "MakeStreamingDNADatabase.py" and "StreamingQueryDNADatabase.py". It might be better to run the auto-test pipe from generating trainingDB step, given that it's pretty fast.



### Pseudo code:

```bash
### Parameters
range= start-end-gap #end is the max_k
filenames.txt #file that record full path of the 165 files
k_size=${end}

### Running
###### creating ref DB
MakeStreamingDNADatabase.py -k ${k_size} -v filenames.txt TrainingDB_k${k_size}.h5
###### get true CI
for file in filenames.txt
  StreamingQueryDNADatabase.py ${file} TrainingDB_k${k_size}.h5 true_CI_out_${file}.csv ${k_size}-${k_size}-1 -v -c 0 -l 0 --sensitive
  StreamingQueryDNADatabase.py ${file} TrainingDB_k${k_size}.h5 truncation_CI_out_${file}.csv ${range} -v -c 0 -l 0 --sensitive
END for
###### merge results and compare
merge all true_CI_out_${file}.csv => merged_true_CI_${k_size}.csv
merge the column with k=${k_size} of truncation_CI_out_${file}.csv => merged_trunc_CI_${k_size}.csv
compare the 2 merged files (maybe a plot or matrix substraction)
```



### Real script

