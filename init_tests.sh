#!/bin/bash

# create the paths for the output files
mkdir -p test_cases/output/actors/celia1990
mkdir -p test_cases/output/non_actors/celia1990

mkdir -p test_cases/output/actors/colbeck1976
mkdir -p test_cases/output/non_actors/colbeck1976

mkdir -p test_cases/output/actors/miller1998
mkdir -p test_cases/output/non_actors/miller1998

mkdir -p test_cases/output/actors/mizoguchi1990
mkdir -p test_cases/output/non_actors/mizoguchi1990

mkdir -p test_cases/output/actors/vanderborght2005
mkdir -p test_cases/output/non_actors/vanderborght2005

mkdir -p test_cases/output/actors/wigmosta1999
mkdir -p test_cases/output/non_actors/wigmosta1999

mkdir -p test_cases/output/actors/umpqua_snowIncept
mkdir -p test_cases/output/non_actors/umpqua_snowIncept

mkdir -p test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist
mkdir -p test_cases/output/non_actors/reynolds_canopySrad_windPrfile_stomResist

mkdir -p test_cases/output/actors/reynolds_groundwatr
mkdir -p test_cases/output/non_actors/reynolds_groundwatr

mkdir -p test_cases/output/actors/senator_alb_method
mkdir -p test_cases/output/non_actors/senator_alb_method


# # modify the paths in the model input file
# # we create a new directories to preserve copies of the original files in case
# # something goes wrong
cd ..
SUMMADIR=`pwd`
export SUMMADIR
cd laugh_tests
BASEDIR=`pwd`
export BASEDIR
TEST_CASE_SETTINGS_DIR=test_cases/settings/syntheticTestCases
for DIR in ${TEST_CASE_SETTINGS_DIR}/*; do
    for FILE in `grep -l '<BASEDIR>' -R ${DIR}`; do
        RUN_FILE="$FILE"
        # Create a file to use for testing while preserving the original fileManager file
        sed "s|<BASEDIR>|${BASEDIR}/test_cases|" $FILE > "${RUN_FILE/./_test.}"
    done
done
WRR_CASE_SETTINGS_DIR=test_cases/settings/wrrPaperTestCases
for DIR in ${WRR_CASE_SETTINGS_DIR}/*; do
    for FILE in `grep -l '<BASEDIR>' -R ${DIR}`; do
        RUN_FILE="$FILE"
        # Create a file to use for testing while preserving the original fileManager file
        sed "s|<BASEDIR>|${BASEDIR}/test_cases|" $FILE > "${RUN_FILE/./_test.}"
    done
done
