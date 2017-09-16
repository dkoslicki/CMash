import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

SCRIPTS = []
SCRIPTS.extend([os.path.join("scripts", script)
                for script in os.listdir(os.path.join(os.path.dirname(__file__), "scripts"))
                if script.endswith(".py")])

HERE = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(HERE, 'README.md'), 'r') as fid:
    LONG_DESCRIPTION = fid.read()

setup(
	name="CMash",
	version="0.2.0",
	author="David Koslicki",
	author_email="dmkoslicki@gmail.com",
	description=("Fast and accurate set similarity estimation via containment min hash (for genomic datasets)."),
	long_description=LONG_DESCRIPTION,
	#license="BSD-3-Clause",  # see classifiers
	keywords="jaccard min hash containment genomics metagenomics",
	url="https://github.com/dkoslicki/CMash",
	packages=['CMash'],
	install_requires=[
		'khmer>=2.1.1',
		'screed',
		'h5py',
		'numpy',
		'blist',
		'argparse',
		'pandas',
		'setuptools>=24.2.0',
		'six'
	],
	python_requires='>=3',
	zip_safe=False,
	package_data={'CMash': ['data/*.fna']},
	scripts=SCRIPTS,
	classifiers=[
		"Development Status :: 3 - Alpha",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
		"Topic :: Scientific/Engineering :: Mathematics",
		"Programming Language :: Python :: 3.5",
		"License :: OSI Approved :: BSD License",
		"Intended Audience :: Science/Research",
		"Programming Language :: Python :: 3.5",
		"Natural Language :: English",
		"Operating System :: MacOS :: MacOS X",
		"Operating System :: POSIX :: Linux",
	],
)
