#! /usr/bin/env python
import argparse
import os
import glob
import pandas as pd
import numpy as np

### import parameters
parser = argparse.ArgumentParser(description="Summarize the results of Groundtruth_CI, estimated_CI and truncated_CI.")
parser.add_argument('range', type=str, help="Range of k-mer sizes in the formate <start>-<end>-<increment>.")
parser.add_argument('query', type=str, help="Query file of analysis shown as absolute path.")

args = parser.parse_args()
input_range= args.range
start=input_range.split("-")[0]
end=input_range.split("-")[1]
gap=input_range.split("-")[2]
query_file= args.query

### merge all results by each metagenome input
f=open(query_file, 'r')
query_list=[ x.strip() for x in list(f)]
f.close()

### summarize groundtruth_CI
os.chdir("ground_truth_CI")
for line in query_list:
    name=line.split("/")[-1]
    print("Processing ground truth of file: " + name)
    merged_gt_CI=pd.concat([ pd.read_csv(x, header=0, index_col=0) for x in glob.glob("ground_truth_CI_k*"+name+"_results.csv") ] , axis=1, join="inner", sort=True)
    merged_gt_CI.to_csv("merged_Ground_truth_CI_"+name+".csv", index=True)
    
os.popen("mv merged_Ground_truth_CI_* ../summary")
os.chdir("../estimated_CI/")

### summary estimated_CI
for line in query_list:
    name=line.split("/")[-1]
    print("Processing estimated CI of file: " + name)
    merged_est_CI=pd.concat([ pd.read_csv(x, header=0, index_col=0) for x in glob.glob("estimated_CI_k*"+name+"_results.csv")], axis=1, join="inner", sort=True)
    merged_est_CI.to_csv("merged_Estimated_CI_"+name+".csv", index=True)

os.popen("mv merged_Estimated_CI_* ../summary")

### go to plot
os.popen("cp ../truncated_CI/truncation_*_results.csv ../summary")
os.chdir("../summary/")
for line in query_list: #single plot for each metagenome
    name=line.split("/")[-1]




