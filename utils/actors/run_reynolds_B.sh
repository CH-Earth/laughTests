#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_groundwatr/actors/summa_fileManager_reynoldsNoXplict_test.txt > ../../test_cases/output/actors/reynolds_groundwatr/reynoldsNoXplict.txt
echo "Reynolds B 1 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_groundwatr/actors/summa_fileManager_reynoldsLumpedQTopmodel_test.txt > ../../test_cases/output/actors/reynolds_groundwatr/reynoldsLumpedQTopmodel.txt
echo "Reynolds B 2 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_groundwatr/actors/summa_fileManager_reynoldsDistributedQTopmodel_test.txt > ../../test_cases/output/actors/reynolds_groundwatr/reynoldsDistributedQTopmodel.txt
echo "Reynolds B 3 Finished"

echo "All Reynolds B Tests Finished, output is located in $BASEDIR/test_cases/output/actors/reynolds_groundwatr/"