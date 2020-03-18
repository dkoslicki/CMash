import sys
import khmer
import screed
sys.path.append("/home/dkoslicki/Desktop/CMash/")
from CMash import MinHash as MH
database = "/home/dkoslicki/Desktop/CMash/dataForShaopeng/test_issue/TrainingDatabase_k_61.h5"
import_list = ['/home/dkoslicki/Desktop/CMash/dataForShaopeng/small_data/taxid_1909294_104_genomic.fna.gz',
               '/home/dkoslicki/Desktop/CMash/dataForShaopeng/small_data/taxid_1909294_109_genomic.fna.gz']
CEs = MH.import_multiple_from_single_hdf5(database, import_list=import_list)
print(f"Jaccard matrix @ 61: {MH.form_jaccard_matrix(CEs)}")
print(f"Raw intersection size @ 61: {len(set(CEs[0]._kmers).intersection(set(CEs[1]._kmers)))}")
print(f"Jaccard @ 61: {len(set(CEs[0]._kmers).intersection(set(CEs[1]._kmers))) / float(len(set(CEs[0]._kmers).union(set(CEs[1]._kmers))))}")
print(f"Jaccard adjust sketch size @ 61: {len(set(CEs[0]._kmers).intersection(set(CEs[1]._kmers))) / float(len(CEs[0]._kmers))}")
print(f"Contain 104 @ 61: {len(set(CEs[0]._kmers).intersection(set(CEs[1]._kmers))) / float(len(set(CEs[0]._kmers)))}")
print(f"Contain 109 @ 61: {len(set(CEs[0]._kmers).intersection(set(CEs[1]._kmers))) / float(len(set(CEs[1]._kmers)))}")
# FIXME: problem is definitely 60-60-1 as it reports a containment of 0.493 between these, yet the real containment both ways is 0.145

def reduce_to(kmers, k_size):
    """
    small helper function that takes in the kmers, truncates them to the k_size
    """
    return_kmers = set()
    for kmer in kmers:
        return_kmers.add(kmer[0:k_size])
    return return_kmers

def reduce_to_with_revcomp(kmers, k_size):
    """
    small helper function that takes in the kmers, truncates them to the k_size
    """
    return_kmers = set()
    for kmer in kmers:
        return_kmers.add(kmer[0:k_size])
        return_kmers.add(khmer.reverse_complement(kmer[0:k_size]))
    return return_kmers

k_size = 40
print(f"Contain 104 @ {k_size}: {len(reduce_to(CEs[0]._kmers,k_size).intersection(reduce_to(CEs[1]._kmers,k_size))) / float(len(reduce_to(CEs[0]._kmers,k_size)))}")
print(f"Contain 109 @ {k_size}: {len(reduce_to(CEs[0]._kmers,k_size).intersection(reduce_to(CEs[1]._kmers,k_size))) / float(len(reduce_to(CEs[1]._kmers,k_size)))}")

print(f"With revcomp contain 104 @ {k_size}: {len(set(reduce_to_with_revcomp(CEs[0]._kmers,k_size)).intersection(set(reduce_to_with_revcomp(CEs[1]._kmers,k_size)))) / float(len(set(reduce_to_with_revcomp(CEs[0]._kmers,k_size))))}")
print(f"With revcomp contain 109 @ {k_size}: {len(set(reduce_to_with_revcomp(CEs[0]._kmers,k_size)).intersection(set(reduce_to_with_revcomp(CEs[1]._kmers,k_size)))) / float(len(set(reduce_to_with_revcomp(CEs[1]._kmers,k_size))))}")

print(f"With revcomp, normalize without contain 104 @ {k_size}: {len(set(reduce_to_with_revcomp(CEs[0]._kmers,k_size)).intersection(set(reduce_to_with_revcomp(CEs[1]._kmers,k_size)))) / float(len(set(reduce_to(CEs[0]._kmers,k_size))))}")
print(f"With revcomp, normalize without contain 109 @ {k_size}: {len(set(reduce_to_with_revcomp(CEs[0]._kmers,k_size)).intersection(set(reduce_to_with_revcomp(CEs[1]._kmers,k_size)))) / float(len(set(reduce_to(CEs[1]._kmers,k_size))))}")

# FIXME: need to check what the containment looks like when I allow either revcomp or non-revcomp matches

# TODO: when setting MakeStreamingDNADatabase to no revcomp and removing all revcomps from StreamingQueryDNADatabase, the results are MUCH more consistent


# FIXME: might also be able to fix this by making alread_seen_kmers include the reverse complement of hits

# FIXME: problem remains in the very first column of the output

# Let's check what the actual containments are so we have a ground truth to check against
max_k = 60
kmer_sets = []
for i, input_file in enumerate(import_list):
    kmer_sets.append(set())
    for record in screed.open(input_file):
        seq = record.sequence
        for j in range(len(seq) - max_k + 1):
            kmer_sets[i].add(seq[j:j + max_k])

print(f"True containment 104 containment at, plain kmers {max_k}: {len(kmer_sets[0].intersection(kmer_sets[1])) / float(len(kmer_sets[0]))}")
print(f"True containment 109 containment at, plain kmers {max_k}: {len(kmer_sets[0].intersection(kmer_sets[1])) / float(len(kmer_sets[1]))}")

# and check the containments with rev-comps
max_k = 40
kmer_sets_rev = []
for i, input_file in enumerate(import_list):
    kmer_sets_rev.append(set())
    for record in screed.open(input_file):
        seq = record.sequence
        for j in range(len(seq) - max_k + 1):
            kmer_sets_rev[i].add(seq[j:j + max_k])
            kmer_sets_rev[i].add(khmer.reverse_complement(seq[j:j + max_k]))


print(f"True containment 104 at {max_k} with rev-comp: {len(kmer_sets_rev[0].intersection(kmer_sets_rev[1])) / float(len(kmer_sets_rev[0]))}")
print(f"True containment 109 at {max_k} with rev-comp: {len(kmer_sets_rev[0].intersection(kmer_sets_rev[1])) / float(len(kmer_sets_rev[1]))}")



# check containment with canonical k-mers
max_k = 40
kmer_sets_can = []
for i, input_file in enumerate(import_list):
    kmer_sets_can.append(set())
    for record in screed.open(input_file):
        seq = record.sequence
        for j in range(len(seq) - max_k + 1):
            kmer = seq[j:j + max_k]
            kmer_rc = khmer.reverse_complement(kmer)
            if kmer_rc < kmer:
                kmer = kmer_rc
            kmer_sets_can[i].add(kmer)


print(f"True containment 104 at {max_k}, canonical: {len(kmer_sets_can[0].intersection(kmer_sets_can[1])) / float(len(kmer_sets_can[0]))}")
print(f"True containment 109 at {max_k}, canonical: {len(kmer_sets_can[0].intersection(kmer_sets_can[1])) / float(len(kmer_sets_can[1]))}")


######
# also check the largest size
# and check the containments with rev-comps
max_k = 61
kmer_sets_rev = []
for i, input_file in enumerate(import_list):
    kmer_sets_rev.append(set())
    for record in screed.open(input_file):
        seq = record.sequence
        for j in range(len(seq) - max_k + 1):
            kmer_sets_rev[i].add(seq[j:j + max_k])
            kmer_sets_rev[i].add(khmer.reverse_complement(seq[j:j + max_k]))


print(f"True containment 104 at {max_k} with rev-comp: {len(kmer_sets_rev[0].intersection(kmer_sets_rev[1])) / float(len(kmer_sets_rev[0]))}")
print(f"True containment 109 at {max_k} with rev-comp: {len(kmer_sets_rev[0].intersection(kmer_sets_rev[1])) / float(len(kmer_sets_rev[1]))}")

# check containment with canonical k-mers
max_k = 61
kmer_sets_can = []
for i, input_file in enumerate(import_list):
    kmer_sets_can.append(set())
    for record in screed.open(input_file):
        seq = record.sequence
        for j in range(len(seq) - max_k + 1):
            kmer = seq[j:j + max_k]
            kmer_rc = khmer.reverse_complement(kmer)
            if kmer_rc < kmer:
                kmer = kmer_rc
            kmer_sets_can[i].add(kmer)


print(f"True containment 104 at {max_k}, canonical: {len(kmer_sets_can[0].intersection(kmer_sets_can[1])) / float(len(kmer_sets_can[0]))}")
print(f"True containment 109 at {max_k}, canonical: {len(kmer_sets_can[0].intersection(kmer_sets_can[1])) / float(len(kmer_sets_can[1]))}")


# FIXME: problem is with the bloom filter. See StreamingQueryScript.py lines 207 ish, don't use the bloom filter and things go fine

