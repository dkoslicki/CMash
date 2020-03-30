#!/bin/bash

# For manual tests. Run this in the tests/script_tests folder

testOrganism="../Organisms/taxid_1192839_4_genomic.fna.gz"

# make the training database
echo "Training on data"
rm TrainingDatabase.h5 2> /dev/null
rm TrainingDatabase.tst 2> /dev/null
/usr/bin/time MakeStreamingDNADatabase.py filenames.txt TrainingDatabase.h5
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

echo "Training on data, verbose"
rm TrainingDatabase.h5 2> /dev/null
rm TrainingDatabase.tst 2> /dev/null
/usr/bin/time MakeStreamingDNADatabase.py filenames.txt TrainingDatabase.h5 --verbose
if test -f TrainingDatabase.h5; then
  if test -f TrainingDatabase.tst; then
    echo "Training file verbose successfully created"
  else
    echo "SOMETHING WENT WRONG!!!!"
    exit 1
  fi
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "Classifying sample, default settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-21-2
if test -f results.csv; then
  echo "default classify successful"
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "Classifying sample, verbose settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-21-2 -v
if test -f results.csv; then
  echo "verbose classify successful"
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "Classifying sample, sensitive settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-21-2 --sensitive
if test -f results.csv; then
  echo "sensitive classify successful"
  cat results.csv
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "Create streaming prefilter"
rm prefilter.bf 2> /dev/null
# make a streaming pre-filter
/usr/bin/time MakeStreamingPrefilter.py TrainingDatabase.h5 prefilter.bf 10-21-2
if test -f prefilter.bf; then
  echo "streaming prefilter creation successful"
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "Classifying sample, default settings with prefilter"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time StreamingQueryDNADatabase.py ${testOrganism} TrainingDatabase.h5 results.csv 10-21-2 -f prefilter.bf --sensitive
if test -f results.csv; then
  echo "default classify with prefilter successful"
  echo "results with prefilter:"
  cat results.csv
else
  echo "SOMETHING WENT WRONG!!!!"
  exit 1
fi

echo "ALL TESTS SUCCESSFULLY COMPLETED!"
