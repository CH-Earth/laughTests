#!/bin/bash

# create the paths for the output files
mkdir -p test_cases/output/actors/celia1990
mkdir -p test_cases/output/non-actors/celia1990

mkdir -p test_cases/output/actors/colbeck1976
mkdir -p test_cases/output/non-actors/colbeck1976

mkdir -p test_cases/output/actors/miller1998
mkdir -p test_cases/output/non-actors/miller1998

mkdir -p test_cases/output/actors/mizoguchi1990
mkdir -p test_cases/output/non-actors/mizoguchi1990

mkdir -p test_cases/output/actors/vanderborght2005
mkdir -p test_cases/output/non-actors/vanderborght2005

mkdir -p test_cases/output/actors/wigmosta1999
mkdir -p test_cases/output/non-actors/wigmosta1999


# # modify the paths in the model input file
# # we create a new directories to preserve copies of the original files in case
# # something goes wrong
cd ..
SUMMADIR=`pwd`
export SUMMADIR
cd Laugh-Tests
BASEDIR=`pwd`
export BASEDIR
TEST_CASE_SETTINGS_DIR=test_cases/settings/syntheticTestCases
for DIR in ${TEST_CASE_SETTINGS_DIR}/*; do
    for FILE in `grep -l '<BASEDIR>' -R ${DIR}`; do
        RUN_FILE="$FILE"
        # Create a file to use for testing while preserving the orignal fileManager file
        sed "s|<BASEDIR>|${BASEDIR}/test_cases|" $FILE > "${RUN_FILE/./_test.}"
    done
done
