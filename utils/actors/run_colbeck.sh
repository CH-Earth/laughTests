#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/colbeck1976/actors/exp1/Summa_Actors_Settings.json > ../../test_cases/output/actors/colbeck1976/colbeck_1.txt
echo "Colbeck 1 has finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/colbeck1976/actors/exp2/Summa_Actors_Settings.json > ../../test_cases/output/actors/colbeck1976/colbeck_2.txt
echo "Colbeck 2 has finsihed"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/colbeck1976/actors/exp3/Summa_Actors_Settings.json > ../../test_cases/output/actors/colbeck1976/colbeck_3.txt
echo "Colbeck 3 has finished"

echo "All Colbeck Test Have Finished. Output files are located in $BASEDIR/test_cases/output/actors/colbeck1976/"