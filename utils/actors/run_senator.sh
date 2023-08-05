#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/senator_alb_method/actors/summa_fileManager_senatorConDecay_test.txt > ../../test_cases/output/actors/senator_alb_method/senatorConDecay.txt
echo "Senator 1 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/senator_alb_method/actors/summa_fileManager_senatorVarDecay_test.txt > ../../test_cases/output/actors/senator_alb_method/senatorVarDecay.txt
echo "Senator 2 Finished"

echo "All Senator Tests Finished, output is located in $BASEDIR/test_cases/output/actors/senator_alb_method/"