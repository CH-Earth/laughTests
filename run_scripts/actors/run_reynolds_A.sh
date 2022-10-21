#! /bin/bash

$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad_windPrfile_stomResist/actors/summa_fileManager_reynoldsDefault_test.txt > ../../test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/reynoldsDefault.txt
echo "Reynolds A 1 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad_windPrfile_stomResist/actors/summa_fileManager_reynoldsUEB2stream_test.txt > ../../test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/reynoldsUEB2stream.txt
echo "Reynolds A 2 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad_windPrfile_stomResist/actors/summa_fileManager_reynoldsNLscatter_test.txt > ../../test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/reynoldsNLscatter.txt
echo "Reynolds A 3 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad_windPrfile_stomResist/actors/summa_fileManager_reynoldsCLM2stream_test.txt > ../../test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/reynoldsCLM2stream.txt
echo "Reynolds A 4 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad_windPrfile_stomResist/actors/summa_fileManager_reynoldsExponential_test.txt > ../../test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/reynoldsExponential.txt
echo "Reynolds A 5 Finished"
$SUMMADIR/bin/summaMain -g 1 -n 1 -c $BASEDIR/test_cases/settings/wrrPaperTestCases/reynolds_canopySrad_windPrfile_stomResist/actors/summa_fileManager_reynoldsJarvis_test.txt > ../../test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/reynoldsJarvis.txt
echo "Reynolds A 6 Finished"

echo "All Reynolds A Tests Finished, output is located in $BASEDIR/test_cases/output/actors/reynolds_canopySrad_windPrfile_stomResist/"