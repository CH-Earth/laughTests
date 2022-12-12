#! /bin/bash

$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/miller1998/non_actors/summa_fileManager_millerClay_test.txt > ../../test_cases/output/non_actors/miller1998/miller_clay.txt
echo "Miller 1 Finished"
$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/miller1998/non_actors/summa_fileManager_millerLoam_test.txt > ../../test_cases/output/non_actors/miller1998/miller_loam.txt
echo "Miller 2 Finished"
$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/miller1998/non_actors/summa_fileManager_millerSand_test.txt > ../../test_cases/output/non_actors/miller1998/miller_sand.txt
echo "Miller 3 Finished"

echo "All Miller Tests Finished. Output located in $BASEDIR/test_cases/output/non_actors/miller1998/"