#!/bin/bash

# This script runs the non-actor test cases. To use this script you have to set the SUMMA_EXE variable to where the summa.exe is located.
# This script runs one test as specified by the user from the input. To run multiple tests use the python script or create your own script that calls the cases within this script.

case $TEST in
    celia1990)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_celia1990_test.txt > $OUTPUT_PATH/celia.txt
        ;;
    colbeck1976)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_colbeck1976-exp1_test.txt > $OUTPUT_PATH/summa_colbeck_1.txt
        echo "Colbeck 1 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_colbeck1976-exp2_test.txt > $OUTPUT_PATH/summa_colbeck_2.txt
        echo "Colbeck 2 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_colbeck1976-exp3_test.txt > $OUTPUT_PATH/summa_colbeck_3.txt
        echo "Colbeck 3 Finished"
        ;;
    miller1998)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_millerClay_test.txt > $OUTPUT_PATH/summa_miller_clay.txt
        echo "Miller Clay Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_millerLoam_test.txt > $OUTPUT_PATH/summa_miller_loam.txt
        echo "Miller Loam Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_millerSand_test.txt > $OUTPUT_PATH/summa_miller_sand.txt
        echo "Miller Sand Finished"
        ;;

    mizoguchi1990)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_mizoguchi1990_test.txt > $OUTPUT_PATH/summa_mizoguchi.txt
        ;;

    vanderborght2005)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_vanderborght2005-exp1_test.txt > $OUTPUT_PATH/summa_vanderborght_1.txt
        echo "Vanderborght 1 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_vanderborght2005-exp2_test.txt > $OUTPUT_PATH/summa_vanderborght_2.txt
        echo "Vanderborght 2 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_vanderborght2005-exp3_test.txt > $OUTPUT_PATH/summa_vanderborght_3.txt
        echo "Vanderborght 3 Finished"
        ;;
esac