#!/bin/bash

# For manual tests. Run this in the tests folder

# make the training database
echo "Training on data"
rm TrainingDatabase.h5 2> /dev/null
rm TrainingDatabase.tst 2> /dev/null
/usr/bin/time python ../scripts/MakeStreamingDNADatabase.py filenames.txt TrainingDatabase.h5
if test -f TrainingDatabase.h5; then
  if test -f TrainingDatabase.tst; then
    echo "Training file successfully created"
  else
    echo "SOMETHING WENT WRONG!!!!"
  fi
  else
    echo "SOMETHING WENT WRONG!!!!"
fi

echo "Training on data, verbose"
rm TrainingDatabase.h5 2> /dev/null
rm TrainingDatabase.tst 2> /dev/null
/usr/bin/time python ../scripts/MakeStreamingDNADatabase.py filenames.txt TrainingDatabase.h5 --verbose
if test -f TrainingDatabase.h5; then
  if test -f TrainingDatabase.tst; then
    echo "Training file verbose successfully created"
  else
    echo "SOMETHING WENT WRONG!!!!"
  fi
else
  echo "SOMETHING WENT WRONG!!!!"
fi

echo "Classifying sample, default settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ../scripts/StreamingQueryDNADatabase.py Organisms/taxid_1192839_4_genomic.fna.gz TrainingDatabase.h5 results.csv 5-21-2
if test -f results.csv; then
  echo "default classify successful"
else
  echo "SOMETHING WENT WRONG!!!!"
fi

echo "Classifying sample, verbose settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ../scripts/StreamingQueryDNADatabase.py Organisms/taxid_1192839_4_genomic.fna.gz TrainingDatabase.h5 results.csv 5-21-2 -v
if test -f results.csv; then
  echo "verbose classify successful"
else
  echo "SOMETHING WENT WRONG!!!!"
fi

echo "Classifying sample, sensitive settings"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ../scripts/StreamingQueryDNADatabase.py Organisms/taxid_1192839_4_genomic.fna.gz TrainingDatabase.h5 results.csv 5-21-2 --sensitive
if test -f results.csv; then
  echo "sensitive classify successful"
else
  echo "SOMETHING WENT WRONG!!!!"
fi

echo "Create streaming prefilter"
rm prefilter.bf 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ../scripts/MakeStreamingPrefilter.py TrainingDatabase.h5 prefilter.bf 5-21-2
if test -f prefilter.bf; then
  echo "streaming prefilter creation successful"
else
  echo "SOMETHING WENT WRONG!!!!"
fi

echo "Classifying sample, default settings with prefilter"
rm results.csv 2> /dev/null
# make a streaming pre-filter
/usr/bin/time python ../scripts/StreamingQueryDNADatabase.py Organisms/taxid_1192839_4_genomic.fna.gz TrainingDatabase.h5 results.csv 5-21-2 -f prefilter.bf
if test -f results.csv; then
  echo "default classify with prefilter successful"
else
  echo "SOMETHING WENT WRONG!!!!"
fi
