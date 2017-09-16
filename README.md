# CMash
CMash is a fast and accurate way to estimate the similarity of two sets. This is a probabilisitic data analysis approach, and uses containment min hashing. Please see the [associated paper](http://www.biorxiv.org/content/early/2017/09/04/184150) for further details (and please cite if you use it):
>Improving Min Hash via the Containment Index with applications to Metagenomic Analysis
>David Koslicki, Hooman Zabeti
>bioRxiv 184150; doi: https://doi.org/10.1101/184150

## Installation
The easiest way to install this is to use [virtualenv](https://virtualenv.pypa.io/en/stable/):
```bash
virtualenv -p python3 CMashVE  # or python3 -m venv CMashVE
source CMashVE/bin/activate
pip install -U pip
pip3 install CMash
```
You can also just use ``pip3 install CMash`` if you don't want to create a virtual environment.

To get the absolute latest edition of CMash, then you can build from the Github repository via:
```bash
virtualenv -p python3 CMashVE
source CMashVE/bin/activate
pip install -U pip
git clone https://github.com/dkoslicki/CMash.git
cd CMash
pip3 install -r requirements.txt
```

Note that this repository itself is python2 and python3 compatible, but the dependency ``khmer`` requires python3 (though curiously enough, it appears ``khmer`` version ``2.1.1`` runs just fine in python2.)
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

#### Other functionality
The module ``MinHash`` (imported in python via ``from CMash import MinHash as MH``) has a bunch more functionality, including (but not limited to!):
1. Fast updates to the training databases (via ``help(MH.delete_from_database)``, ``help(MH.insert_to_database)``, ``help(MH.union_databases)``)
2. Ability to form a matrix of Jaccard indexes (for comparison of all pairwise Jaccard indexes of organisms in the training database). This is useful for identifying redundances/patterns/structure in your training database: ``help(MH.form_jaccard_count_matrix)`` and ``help(MH.form_jaccard_matrix)``.
3. Access to the k-mers that MinHash randomly selected (see the class ``CountEstimator`` and the associated ``_kmers`` data structure.)

I'd encourage you to poke through the source code of ``MinHash.py`` and take a look at the scripts as well.

Protein databases (and for that matter, arbitrary K-length strings) coming soon...
