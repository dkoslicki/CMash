#!/bin/bash

# This is to compare the results of CMash against a brute force calculated containment indicies.
# WARNING: this uses a LOT of memory and is CPU intensive as well, so use with caution (and use on a server)

# In case you have multiple versions installed (eg. Metalign as well as CMash), make sure python is looking in the right place:
export PYTHONPATH="$(dirname $(dirname "`pwd`"))":$PYTHONPATH

#Make sure the correct CMash is being pulled from
testFile=$(python -c "from CMash import MinHash as MH; print(MH.__file__)")
parentDir=`dirname $PWD`
parentDir=`dirname ${parentDir}`
correctFile="${parentDir}/CMash/MinHash.py"
if [ "$testFile" == "$correctFile" ];
then
  echo "Files are correct"
else
  echo "Files are not correct"
  exit 1
fi

testOrganism="../Organisms/taxid_1192839_4_genomic.fna.gz"
maxK=22
kSizes="10-${maxK}-2"
numHashes=2000
containmentThresh=0
locationOfThresh=-1

scriptsDir="${parentDir}/scripts"
modulesDir="${parentDir}/CMash"

# make the training database
echo "Training on data"
rm TrainingDatabase.h5 2> /dev/null
rm TrainingDatabase.tst 2> /dev/null
/usr/bin/time python ${scriptsDir}/MakeStreamingDNADatabase.py filenames.txt TrainingDatabase.h5 -n ${numHashes} -k ${maxK}
if test -f TrainingDatabase.h5; then
  if test -f TrainingDatabase.tst; then
    echo "Training file successfully created"
  else
    echo "SOMETHING WENT WRONG!!!!"
    exit 1
  fi
  else
    echo "SOMETHING WENT WRONG!!!!"
    exit 1
fi

echo "Classifying sample with CMash"
outName="est_results.csv"
rm ${outName} 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ${scriptsDir}/StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 ${outName} $kSizes --sensitive -l $locationOfThresh -c $containmentThresh
if test -f ${outName}; then
  echo "Estimate successful:"
  #cat ${outName}
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "Computing ground truth containment indicies"
outName="true_results.csv"
rm ${outName} 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ${modulesDir}/GroundTruth.py ${testOrganism} TrainingDatabase.h5 ${outName} $kSizes -l $locationOfThresh -c $containmentThresh
if test -f ${outName}; then
  echo "Ground truth successful:"
  #cat ${outName}
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

# then compare the results
python -c "import pandas as pd
import numpy as np
df = pd.read_csv('true_results.csv', index_col=0)
df.sort_index(inplace=True)
df2 = pd.read_csv('est_results.csv', index_col=0)
df2.sort_index(inplace=True)
diff = np.abs(df - df2)
diff.sort_index(inplace=True)
print('True:');
print(df)
print('')
print('CMash:');
print(df2)
print('')
print('|true-CMash|:');
print(diff)
print('')
print('Total error per k-mer size:')
print(diff.sum())
print('')
print('Median error per k-mer size:')
print(diff.median())
print('')
print('Mean error per k-mer size:')
print(diff.mean())
"