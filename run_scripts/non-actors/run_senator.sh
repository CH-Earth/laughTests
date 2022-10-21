#! /bin/bash

$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/wrrPaperTestCases/senator_alb_method/non_actors/summa_fileManager_senatorConDecay_test.txt > ../../test_cases/output/non_actors/senator_alb_method/senatorConDecay.txt
echo "Senator 1 Finished"
$SUMMADIR/bin/summa_sundials.exe -g 1 1 -m $BASEDIR/test_cases/settings/wrrPaperTestCases/senator_alb_method/non_actors/summa_fileManager_senatorVarDecay_test.txt > ../../test_cases/output/non_actors/senator_alb_method/senatorVarDecay.txt
echo "Senator 2 Finished"

echo "All Senator Tests Finished, output is located in $BASEDIR/test_cases/output/non_actors/senator_alb_method/"