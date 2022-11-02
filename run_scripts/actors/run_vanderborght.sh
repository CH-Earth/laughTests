#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/vanderborght2005/actors/exp1/Summa_Actors_Settings.json > ../../test_cases/output/actors/vanderborght2005/vanderborght_1.txt
echo "Vanderborght 1 Finished"

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/vanderborght2005/actors/exp2/Summa_Actors_Settings.json > ../../test_cases/output/actors/vanderborght2005/vanderborght_2.txt
echo "Vanderborght 2 Finished"

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/syntheticTestCases/vanderborght2005/actors/exp3/Summa_Actors_Settings.json > ../../test_cases/output/actors/vanderborght2005/vanderborght_3.txt
echo "Vanderborght 3 Finished"

echo "All Vanderborght Tests Finished. Ouptut located in $BASEDIR/test_cases/output/actors/vanderborght2005"
