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
maxK=5
kSizes="4-${maxK}-1"
numHashes=10
containmentThresh=.01
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
  echo "sensitive classify successful"
  cat ${outName}
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
  echo "sensitive classify successful"
  cat ${outName}
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi
