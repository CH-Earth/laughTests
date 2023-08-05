#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/umpqua_snowIncept/actors/summa_fileManager_umpquaStickySnow_test.txt > ../../test_cases/output/actors/umpqua_snowIncept/umpquaStickySnow.txt
echo "Umpqua 1 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/umpqua_snowIncept/actors/summa_fileManager_umpquaLightSnow_test.txt > ../../test_cases/output/actors/umpqua_snowIncept/umpquaLightSnow.txt
echo "Umpqua 2 Finished"

echo "All Umpqua Tests Finished, output is located in $BASEDIR/test_cases/output/actors/umpqua_snowIncept/"