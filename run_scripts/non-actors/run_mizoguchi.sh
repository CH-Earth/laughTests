#! /bin/bash

$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/mizoguchi1990/non_actors/summa_fileManager_mizoguchi_test.txt > ../../test_cases/output/non_actors/mizoguchi1990/mizoguchi.txt

echo "Mizoguchi Test Finished. Output located in $BASEDIR/test_cases/output/non_actors/mizoguchi1990/"