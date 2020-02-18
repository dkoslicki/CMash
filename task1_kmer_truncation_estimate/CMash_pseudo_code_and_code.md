# CMash running task1

### Introduction: 

[CMash record](https://github.com/ShaopengLiu1/PSU_Bioinformatics/blob/master/r3_Koslicki_group/CMash_record.md)



### Pseudo-code

```bash
### CI: CMash Index

### i) find true CI with k=61 + all variation data on each k-size point
max_k=61
TrainStreamingDNADatabase.py <all_genome> <max_k>  -> TB_61.hs
for gi in <all_genome>:
  StreamingqueryDNADatabase.py	gi	TB_61.hs	1-61-3	-c 0 -l 0 
# In the output matrix: col[k=61] is the true value of each gi; col[k<61] is the variation due to truncation
  
### ii) find the true CI with all k<61
for max_k in range (1,58,3):
  TrainStreamingDNADatabase.py <all_genome> <max_k>  -> TB_${max_k}.hs
  for gi in <all_genome>:
    StreamingqueryDNADatabase.py	gi	TB_${max_k}.hs	${max_k}-${max_k}-1	-c 0 -l 0 
# there will be 19 output, each of them is a true CI at k=${max_k}
```



### Real script:



