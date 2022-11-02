#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases//settings/syntheticTestCases/miller1998/actors/clay/Summa_Actors_Settings.json > ../../test_cases/output/actors/miller1998/miller_clay.txt
# echo "Miller 1 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/miller1998/actors/loam/Summa_Actors_Settings.json > ../../test_cases/output/actors/miller1998/miller_loam.txt
echo "Miller 2 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/miller1998/actors/sand/Summa_Actors_Settings.json > ../../test_cases/output/actors/miller1998/miller_sand.txt
# echo "Miller 3 Finished"

echo "All Miller Tests Finished. Output located in $BASEDIR/test_cases/output/actors/miller1998/"