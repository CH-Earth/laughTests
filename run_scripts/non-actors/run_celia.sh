#! /bin/bash
$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/celia1990/non_actors/summa_fileManager_celia1990_test.txt > ../../test_cases/output/non_actors/celia1990/celia.txt
echo "Celia Test Finished, output is located in $BASEDIR/test_cases/output/non_actors/celia1990/"