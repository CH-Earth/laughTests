#! /bin/bash

$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/wrrPaperTestCases/umpqua1996/non_actors/summa_fileManager_umpquaAndrea_test.txt > ../../test_cases/output/non_actors/umpqua1996/umpquaAndrea.txt
echo "Umpqua 1 Finished"
$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/wrrPaperTestCases/umpqua1996/non_actors/summa_fileManager_umpquaHedPom_test.txt > ../../test_cases/output/non_actors/umpqua1996/umpquaHedpom.txt
echo "Umpqua 2 Finished"

echo "All Umpqua Tests Finished, output is located in $BASEDIR/test_cases/output/non_actors/umpqua1996/"