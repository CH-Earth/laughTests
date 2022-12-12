#!/bin/bash

# This script runs the non-actor test cases. To use this script you have to set the SUMMA_EXE variable to where the summa.exe is located.
# This script runs one test as specified by the user from the input. To run multiple tests use the python script or create your own script that calls the cases within this script.

case $TEST in
    celia1990)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_celia1990_test.txt > $OUTPUT_PATH/celia.txt
        ;;
    colbeck1976)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_colbeck1976-exp1_test.txt > $OUTPUT_PATH/colbeck_1.txt
        echo "Colbeck 1 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_colbeck1976-exp2_test.txt > $OUTPUT_PATH/colbeck_2.txt
        echo "Colbeck 2 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_colbeck1976-exp3_test.txt > $OUTPUT_PATH/colbeck_3.txt
        echo "Colbeck 3 Finished"
        ;;
    miller1998)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_millerClay_test.txt > $OUTPUT_PATH/miller_clay.txt
        echo "Miller Clay Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_millerLoam_test.txt > $OUTPUT_PATH/miller_loam.txt
        echo "Miller Loam Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_millerSand_test.txt > $OUTPUT_PATH/miller_sand.txt
        echo "Miller Sand Finished"
        ;;

    mizoguchi1990)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/summa_fileManager_mizoguchi1990_test.txt > $OUTPUT_PATH/mizoguchi.txt
        ;;

    vanderborght2005)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH//summa_fileManager_vanderborght2005_exp1_test.txt > $OUTPUT_PATH//vanderborght_1.txt
        echo "Vanderborght 1 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH//summa_fileManager_vanderborght2005_exp2_test.txt > $OUTPUT_PATH//vanderborght_2.txt
        echo "Vanderborght 2 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH//summa_fileManager_vanderborght2005_exp3_test.txt > $OUTPUT_PATH//vanderborght_3.txt
        echo "Vanderborght 3 Finished"
        ;;
esac