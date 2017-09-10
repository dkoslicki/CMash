# These are a collection of functions that may be useful to include in MinHash.py in future releases
import numpy as np
import os
import screed
import multiprocessing
from itertools import *

def cluster_matrix(A_eps, A_indicies, taxonomy, cluster_eps=.01):
    """
    This function clusters the indicies of A_eps such that for a given cluster, there is another element in that cluster
    with similarity (based on A_eps) >= cluster_eps for another element in that same cluster. For two elements of
    distinct clusters, their similarity (based on A_eps) < cluster_eps.
    :param A_eps: The jaccard or jaccard_count matrix containing the similarities
    :param A_indicies: The basis of the matrix A_eps (in terms of all the CEs)
    :param cluster_eps: The similarity threshold to cluster on
    :return: (a list of sets of indicies defining the clusters, LCAs of the clusters)
    """
    #A_indicies_numerical = np.where(A_indicies == True)[0]
    A_indicies_numerical = A_indicies
    # initialize the clusters
    clusters = []
    for A_index in range(len(A_indicies_numerical)):
        # Find nearby elements
        nearby = set(np.where(A_eps[A_index, :] >= cluster_eps)[0]) | set(np.where(A_eps[:, A_index] >= cluster_eps)[0])
        in_flag = False
        in_counter = 0
        in_indicies = []
        for i in range(len(clusters)):
            if nearby & clusters[i]:
                clusters[i].update(nearby)  # add the nearby indicies to the cluster
                in_counter += 1  # keep track if the nearby elements belong to more than one of the previously formed clusters
                in_indicies.append(i)  # which clusters nearby shares elements with
                in_flag = True  # tells if it forms a new cluster
        if not in_flag:  # if new cluster, then append to clusters
            clusters.append(set(nearby))
        if in_counter > 1:  # If it belongs to more than one cluster, merge them together
            merged_cluster = set()
            for in_index in in_indicies[::-1]:
                merged_cluster.update(clusters[in_index])
                del clusters[in_index]  # delete the old clusters (now merged)
            clusters.append(merged_cluster)  # append the newly merged clusters
    clusters_full_indicies = []
    for cluster in clusters:
        cluster_full_indicies = set()
        for item in cluster:
            cluster_full_indicies.add(A_indicies_numerical[item])
        clusters_full_indicies.append(cluster_full_indicies)
    # Check to make sure the clustering didn't go wrong
    if sum([len(item) for item in clusters_full_indicies]) != len(A_indicies_numerical):  # Check the correct number of indicies
        raise Exception("For some reason, the total number of indicies in the clusters doesn't equal the number of indicies you started with")
    if set([item for subset in clusters_full_indicies for item in subset]) != set(A_indicies_numerical):  # Make sure no indicies were missed or added
        raise Exception("For some reason, the indicies in all the clusters doesn't match the indicies you started with")
    return clusters_full_indicies, cluster_LCAs(clusters_full_indicies, taxonomy)


def cluster_LCAs(clusters, taxonomy):
    """
    This function returns the lowest common ancestor in each one of the input clusters
    :param clusters: input clusters
    :param taxonomy: input taxonomy
    :return: a list with the ith element being the lowest common ancestor of the ith cluster
    """
    LCAs = []
    for cluster in clusters:
        found_LCA = False
        if len(cluster) == 1:
            LCAs.append(taxonomy[list(cluster)[0]].split()[2].split('|')[-1])
            found_LCA = True
            continue
        cluster_taxonomy = []
        for index in cluster:
            cluster_taxonomy.append(taxonomy[index])
        for rank in range(7, -1, -1):
            rank_names = []
            dummy_name = 0
            for tax_path in cluster_taxonomy:
                split_taxonomy = tax_path.split()[2].split('|')
                if len(split_taxonomy) < rank + 1:
                    rank_names.append(str(dummy_name))
                else:
                    rank_names.append(split_taxonomy[rank])
            if len(set(rank_names)) == 1 and "0" not in rank_names:
                LCAs.append(rank_names[0])
                found_LCA = True
                break
        if not found_LCA:
            LCAs.append('sk__-1_microorganism')  # In case they don't even have the kingdom in common
    return LCAs


def _write_single_cluster(tup):
    """
    Helper function. Writes a single fast file consisting of all the sequences in input_file_names
    :param tup: input tuple (out_dir, LCA, cluster_index, input_file_names)
    :return: the name of the created file
    """
    out_dir = tup[0]
    LCA = tup[1]
    cluster_index = tup[2]
    input_file_names = tup[3]
    out_file_name = os.path.join(out_dir, LCA + "_" + str(cluster_index) + "_" + ".fa")  # put the cluster index in the name in case there are shared LCAs
    out_file = open(out_file_name, 'w')
    i = 0
    for file_name in input_file_names:
        for record in screed.open(file_name):
            out_file.write(">" + LCA + "_" + str(i))
            out_file.write("\n")
            out_file.write(record.sequence)
            out_file.write("\n")
            i += 1
    out_file.close()
    return out_file_name

def make_cluster_fastas(out_dir, LCAs, clusters, CEs, threads=multiprocessing.cpu_count()):
    """
    This function will write a single fasta file for each of the clusters
    :param out_dir: the output directory in which to write the fasta files
    :param LCAs: the least common ancestors (from cluster_LCAs())
    :param clusters: the clusters (from cluster_matrix())
    :param CEs: The list of count estimators
    :param threads: number of threads to use
    :return: a list of files created (to be used in build_reference())
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    pool = multiprocessing.Pool(processes=threads)
    file_names = pool.map(_write_single_cluster, zip(repeat(out_dir), LCAs, range(len(LCAs)), [[CEs[i].input_file_name for i in cluster] for cluster in clusters]), chunksize=1)
    pool.close()
    #pool.terminate()
    #pool.join()
    return file_names

