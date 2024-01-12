# CMash

# Please Note: CMash has largely been supplanted by [Sourmash](https://github.com/sourmash-bio/sourmash). While Sourmash technically does not possess the ability to change k-mer sizes, we have decided to adopt Sourmash and incorporate our other results (eg. estimating ANI and AAI) into the Sourmash code base. Consider this repo depreciated.

CMash is a fast and accurate way to estimate the similarity of two sets. This is a probabilisitic data analysis approach, and uses containment min hashing. Please see the [associated paper]( https://doi.org/10.1101/2021.12.06.471436) for further details (and please cite if you use it):
```
Liu, S., & Koslicki, D. (2021). CMash: fast, multi-resolution estimation of k-mer-based Jaccard and containment indices. bioRxiv.  https://doi.org/10.1101/2021.12.06.471436
```
# Be aware, this is a work in progress and isn't guaranteed to be functional

## Installation
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/cmash/README.html)

The best way to install is via bioconda:
```bash
conda install cmash
```

Alternatively: use [virtualenv](https://virtualenv.pypa.io/en/stable/):
```bash
virtualenv CMashVE
source CMashVE/bin/activate
pip install -U pip
git clone https://github.com/dkoslicki/CMash.git
cd CMash
pip install -r requirements.txt
```

Note: if installing on Ubuntu and you get a ``no package 'cairo' found`` error, this can be fixed by running ``sudo apt install libcairo2-dev`` and installing again.


# The below is very outdated and will be changed upon v1.0.0 release


## Usage
The basic paradigm is to create a reference/training database, form a sample bloom filter, and then query the database.

#### Forming a reference/training database
Say you have three reference fasta/q file: ``ref1.fa``, ``ref2.fa`` and ``ref3.fa``. In a file (here called ``FileNames.txt``), place the absolute paths pointing to the fasta/q files:
```bash
cat FileNames.txt
# /abs/path/to/ref1.fa
# /abs/path/to/ref2.fa
# /abs/path/to/ref3.fa
```
Then you can create the training database via:
```bash
MakeDNADatabase.py FileNames.txt TrainingDatabase.h5
```
See ``MakeDNADatabase.py -h`` for more options when forming a database.

#### Creating a sample bloom filter
Given a (large) query fasta/q file ``Metagenome.fa``, you can *optionally* create a bloom filter via ``MakeNodeGraph.py Metagenome.fa .``. 
See ``MakeNodeGraph.py -h`` for more details about this function.

This step is not strictly necessary (as the next step automatically forms a nodegraph/bloom filter if you didn't already create one). 
However, I've provided this script in case you want to pre-process a bunch of metagenomes.

#### Query the database
To get containment and Jaccard index estimates of the references files in your query file ``Metagenome.fa``, use something like the following ``QueryDNADatabase.py Metagenome.fa TrainingDatabase.h5 Output.csv``.

There are a bunch of options available: ``QueryDNADatabase.py -h``. The output file is a CSV file with rows corresponding (in this case) to ``ref1.fa``, ``ref2.fa``, and ``ref3.fa`` and columns corresponding to the containment index estimate, intersection cardinality, and Jaccard index estimate.

## Streaming

There are streaming versions of database formation and querying. These are contained in the scripts:
``MakeStreamingDNADatabase.py`` and ``StreamingQueryDNADatabase.py``. Both of these are similar to the bloom filter counterparts (given above) except for the fact that:
1. They stream the data (so no bloom filters are required)
2. They allow multiple k-mer sizes to be used simultaneously. This allows for visualizing the containment index as a function of ``k`` (ala Figure 4a in [the MetaPalette publication](http://msystems.asm.org/content/1/3/e00020-16)).

#### Other functionality
The module ``MinHash`` (imported in python via ``from CMash import MinHash as MH``) has a bunch more functionality, including (but not limited to!):
1. Fast updates to the training databases (via ``help(MH.delete_from_database)``, ``help(MH.insert_to_database)``, ``help(MH.union_databases)``)
2. Ability to form a matrix of Jaccard indexes (for comparison of all pairwise Jaccard indexes of organisms in the training database). This is useful for identifying redundances/patterns/structure in your training database: ``help(MH.form_jaccard_count_matrix)`` and ``help(MH.form_jaccard_matrix)``.
3. Access to the k-mers that MinHash randomly selected (see the class ``CountEstimator`` and the associated ``_kmers`` data structure.)

I'd encourage you to poke through the source code of ``MinHash.py`` and take a look at the scripts as well.

Protein databases (and for that matter, arbitrary K-length strings) coming soon...
