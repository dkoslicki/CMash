import numpy as np
from CMash import MinHash as MH
import seaborn
import matplotlib.pyplot as plt
from collections import Counter
# Data prep

# In bash:
#mkdir dataForShaopeng
#cd dataForShaopeng/
#mkdir data
#cd data
#wget https://ucla.box.com/shared/static/c1g8xjc9glh68oje9e549fjqj0y8nc17.gz && tar -zxvf c1g8xjc9glh68oje9e549fjqj0y8nc17.gz && rm c1g8xjc9glh68oje9e549fjqj0y8nc17.gz  # grabbed this from the Metalign setup_data.sh
#ls | xargs -I{} sh -c 'readlink -f {} >> ../file_names.txt'
#cd ..
#head -n 10 file_names.txt > file_names_10.txt
#head -n 100 file_names.txt > file_names_100.txt
#head -n 1000 file_names.txt > file_names_1000.txt
#cd ../scripts/
#python MakeStreamingDNADatabase.py ../dataForShaopeng/file_names_10.txt ../dataForShaopeng/TrainingDatabase_10_k_60.h5 -n 1000 -k 60 -v
#python MakeStreamingDNADatabase.py ../dataForShaopeng/file_names_100.txt ../dataForShaopeng/TrainingDatabase_100_k_60.h5 -n 1000 -k 60 -v
#python MakeStreamingDNADatabase.py ../dataForShaopeng/file_names_1000.txt ../dataForShaopeng/TrainingDatabase_1000_k_60.h5 -n 1000 -k 60 -v

################################################
# Simon stuff
# cd /home/dkoslicki/Data/Repos/CMash/dataForShaopeng
# rm file_names_100.txt
# cd /home/dkoslicki/Data/Repos/CMash/dataForShaopeng/organism_files/organism_files
# ls | head -n 100 | xargs -I{} sh -c 'readlink -f {} >> ../../file_names_100.txt'
# cd ../../../scripts/
# python MakeStreamingDNADatabase.py ../dataForShaopeng/file_names_100.txt ../dataForShaopeng/TrainingDatabase_100_k_20.h5 -n 1000 -k 20 -v
#################################################

def cluster_matrix(A_eps, A_indicies, cluster_eps=.01):
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
    return clusters_full_indicies#, cluster_LCAs(clusters_full_indicies, taxonomy)


n = 100
cluster_eps = .01
k = 20
dir_base_name = "/home/dkoslicki/Data/Repos/CMash/dataForShaopeng/"
CEs = MH.import_multiple_from_single_hdf5(f"{dir_base_name}TrainingDatabase_{n}_k_{k}.h5")
mat = MH.form_jaccard_matrix(CEs)
clusters_full_indicies = cluster_matrix(mat, range(len(CEs)), cluster_eps=cluster_eps)
cluster_sizes = [len(x) for x in clusters_full_indicies]
max_cluster_loc = np.argmax(cluster_sizes)
max_cluster_indicies = list(clusters_full_indicies[max_cluster_loc])
print(len(max_cluster_indicies))
sub_mat = mat[max_cluster_indicies,:][:,max_cluster_indicies]
sub_CEs = [CEs[x] for x in max_cluster_indicies]
out_file_names = [x.input_file_name.decode('utf-8') for x in sub_CEs]
fid = open(f"{dir_base_name}to_select.txt", 'w')
for name in out_file_names:
    fid.write(f"{name}\n")
fid.close()
#seaborn.heatmap(sub_mat)
#plt.show()
clustergrid = seaborn.clustermap(sub_mat)
plt.show()
# to check the kinds of organisms
#cat to_select.txt  | xargs -I{} sh -c 'zcat {} | head -n 1'

####################################################################
# Simon stuff
# Get all the k_mers in all these sketches
all_kmers = set()
for CE in sub_CEs:
    all_kmers.update(CE._kmers)
print(len(all_kmers))

# also get a list of all their counts
all_counts = list()
for CE in sub_CEs:
    all_counts.extend(CE._counts)

print(Counter(all_counts))

# reordered indicies
re_ordered_indicies = clustergrid.dendrogram_row.reordered_ind
sub_mat_reordered = sub_mat[re_ordered_indicies, :][:, re_ordered_indicies]
seaborn.heatmap(sub_mat_reordered)
plt.show()

# then print out these file names
sub_CEs_reordered = [sub_CEs[i] for i in re_ordered_indicies]
out_file_names = [x.input_file_name.decode('utf-8') for x in sub_CEs_reordered]
fid = open(f"{dir_base_name}to_select.txt", 'w')
for name in out_file_names:
    fid.write(f"{name}\n")
fid.close()

# Then print these out into something readable by matlab
# Basically, with the 100, find a good cluster, write the to to_select.txt
# python MakeStreamingDNADatabase.py ../dataForShaopeng/to_select.txt ../dataForShaopeng/TrainingDatabase_100_k_20.h5 -n 1000 -k 20 -v
# Then read them in and make the A matrices
n = 100
k = 50
dir_base_name = "/home/dkoslicki/Data/Repos/CMash/dataForShaopeng/"
import_list = []
with open(f"{dir_base_name}to_select.txt", 'r') as fid:
    for line in fid.readlines():
        import_list.append(line.strip())

CEs = MH.import_multiple_from_single_hdf5(f"{dir_base_name}TrainingDatabase_{n}_k_{k}.h5")
# reorder them since they are auto-imported alphabetically, not in the order of the to_select.txt input for training
CEs_reordered = []
for file_name in import_list:
    for CE in CEs:
        if CE.input_file_name.decode('utf-8') == file_name:
            CEs_reordered.append(CE)
CEs = CEs_reordered
mat = MH.form_jaccard_matrix(CEs)
seaborn.heatmap(mat)  # sanity check

# get the basis (all kmers in all CEs)
all_kmers = set()
for CE in CEs:
    all_kmers.update(CE._kmers)

basis_kmers = list(all_kmers)
A_mat = np.zeros((len(basis_kmers), len(CEs)))
for j in range(len(CEs)):
    CE = CEs[j]
    CE_kmers = set(CE._kmers)
    CE_kmers_to_counts = dict()
    for (kmer, count) in zip(CE._kmers, CE._counts):
        CE_kmers_to_counts[kmer] = count
    for i in range(len(basis_kmers)):
        kmer = basis_kmers[i]
        if kmer in CE_kmers:
            A_mat[i, j] = CE_kmers_to_counts[kmer]

seaborn.heatmap(A_mat)
# dump it
write_file = f"/home/dkoslicki/Dropbox/SimonFoucart/March 3 2020 Visit/TestDataForUnknownSpecies/A_{k}.csv"
np.savetxt(write_file, A_mat, delimiter=",", fmt="%d")

