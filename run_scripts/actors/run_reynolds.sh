#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad/actors/summa_fileManager_reynoldsUEB2stream_test.txt > ../../test_cases/output/actors/reynolds_canopySrad/reynoldsUEB2stream.txt
echo "Reynolds 1 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad/actors/summa_fileManager_reynoldsNLscatter_test.txt > ../../test_cases/output/actors/reynolds_canopySrad/reynoldsNLscatter.txt
echo "Reynolds 2 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c  $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad/actors/summa_fileManager_reynoldsCLM2stream_test.txt > ../../test_cases/output/actors/reynolds_canopySrad/reynoldsCLM2stream.txt
echo "Reynolds 3 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad/actors/summa_fileManager_reynoldsBeersLaw_test.txt > ../../test_cases/output/actors/reynolds_canopySrad/reynoldsBeersLaw.txt
echo "Reynolds 4 Finished"

echo "All Reynolds Tests Finished, output is located in $BASEDIR/test_cases/output/actors/reynolds_canopySrad/"