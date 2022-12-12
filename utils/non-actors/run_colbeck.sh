#! /bin/bash

$SUMMADIR/bin/summa.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/colbeck1976/fileManagers/summa_fileManager_colbeck1976-exp1_test.txt > ../../test_cases/output/colbeck1976/colbeck_1.txt
echo "Colbeck 1 has finished"
$SUMMADIR/bin/summa.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/colbeck1976/fileManagers/summa_fileManager_colbeck1976-exp2_test.txt > ../../test_cases/output/colbeck1976/colbeck_2.txt
echo "Colbeck 2 has finsihed"
$SUMMADIR/bin/summa.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/colbeck1976/fileManagers/summa_fileManager_colbeck1976-exp3_test.txt > ../../test_cases/output/colbeck1976/colbeck_3.txt
echo "Colbeck 3 has finished"

echo "All Colbeck Test Have Finished. Output files are located in $BASEDIR/test_cases/output/colbeck1976/"