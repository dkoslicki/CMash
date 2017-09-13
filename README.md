# CMash
Fast and accurate set similarity estimation via containment min hash

In progress... Full documentation soon.

Until further notice:
1. ``MakeDNADatabase.py --help`` to create a reference database.
2. ``MakeNodeGraph.py --help`` to prep a fasta/q file for comparison against a reference database.
3. ``QueryDNADatbase.py --help`` to query a fasta/q file against a reference database. Forgot step #2? No problem! The file will be prepped automatically!

Install via traditional ``git clone`` (with ``pip install khmer`` and ``argparse``, ``blist``, ``h5py``, ``numpy``, ``pandas``, ``matplotlib``, and ``screed`` installed similarly).

Protein databases (and for that matter, arbitrary K-length strings) coming soon...
