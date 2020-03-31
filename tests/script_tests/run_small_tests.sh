#!/bin/bash

# For manual tests. Run this in the tests/script_tests folder

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

scriptsDir="${parentDir}/scripts"
testOrganism="../Organisms/taxid_1192839_4_genomic.fna.gz"

# make the training database
echo "Training on data"
rm TrainingDatabase.h5 2> /dev/null
rm TrainingDatabase.tst 2> /dev/null
/usr/bin/time python ${scriptsDir}/MakeStreamingDNADatabase.py filenames.txt TrainingDatabase.h5 -k 10
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

echo "Classifying sample, sensitive settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ${scriptsDir}/StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-10-1 --sensitive
if test -f results.csv; then
  echo "sensitive classify successful"
  cat results.csv
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

: << 'END1'
echo "Classifying sample, specific settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ${scriptsDir}/StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-21-2
if test -f results.csv; then
  echo "specific classify successful"
  cat results.csv
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi
END1

# intersection tests

echo "Classifying sample, sensitive settings, with KMC intersect"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ${scriptsDir}/StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-10-1 --sensitive --intersect
if test -f results.csv; then
  echo "sensitive classify successful"
  cat results.csv
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

: << 'END'
echo "Classifying sample, specific settings, with KMC intersect"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ${scriptsDir}/StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-21-2 --intersect
if test -f results.csv; then
  echo "specific classify successful"
  cat results.csv
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi
END