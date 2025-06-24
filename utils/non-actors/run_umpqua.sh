#! /bin/bash

$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/wrrPaperTestCases/umpqua_snowIncept/non_actors/summa_fileManager_umpquaStickySnow_test.txt > ../../test_cases/output/non_actors/umpqua_snowIncept/umpquaStickySnow.txt
echo "Umpqua 1 Finished"
$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/wrrPaperTestCases/umpqua_snowIncept/non_actors/summa_fileManager_umpquaLightSnow_test.txt > ../../test_cases/output/non_actors/umpqua_snowIncept/umpquaLightSnow.txt
echo "Umpqua 2 Finished"

echo "All Umpqua Tests Finished, output is located in $BASEDIR/test_cases/output/non_actors/umpqua_snowIncept/"
