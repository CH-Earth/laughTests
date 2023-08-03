#! /bin/bash

$SUMMADIR/bin/summa.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/vanderborght2005/fileManagers/summa_fileManager_vanderborght2005_exp1_test.txt > ../../test_cases/output/vanderborght2005/vanderborght_1.txt
echo "Vanderborght 1 Finished"

$SUMMADIR/bin/summa.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/vanderborght2005/fileManagers/summa_fileManager_vanderborght2005_exp2_test.txt > ../../test_cases/output/vanderborght2005/vanderborght_2.txt
echo "Vanderborght 2 Finished"

$SUMMADIR/bin/summa.exe -g 1 1 -m $BASEDIR/test_cases/settings/syntheticTestCases/vanderborght2005/fileManagers/summa_fileManager_vanderborght2005_exp3_test.txt > ../../test_cases/output/vanderborght2005/vanderborght_3.txt
echo "Vanderborght 3 Finished"

echo "All Vanderborght Tests Finished. Ouptut located in $BASEDIR/test_cases/output/non_actors/vanderborght2005"