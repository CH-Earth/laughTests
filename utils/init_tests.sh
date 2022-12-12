#!/bin/bash

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
