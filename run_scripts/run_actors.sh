#!/bin/bash

# This script runs the actor test cases. To use this script you have to set the SUMMA_EXE variable to where the summa_actors executable is located.
# This script runs one test as specified by the user from the input. To run multiple tests use the python script or create your own script that calls the cases within this script.

case $TEST in
    celia1990)
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_celia1990.json > $OUTPUT_PATH/celia_actors.txt
        ;;
    colbeck1976)
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_colbeck1976-exp1.json > $OUTPUT_PATH/colbeck1_actors.txt
        echo "Colbeck 1 Finished"
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_colbeck1976-exp2.json > $OUTPUT_PATH/colbeck2_actors.txt
        echo "Colbeck 2 Finished"
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_colbeck1976-exp3.json > $OUTPUT_PATH/colbeck3_actors.txt
        echo "Colbeck 3 Finished"
        ;;
    miller1998)
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_millerClay.json > $OUTPUT_PATH/millerClay_actors.txt
        echo "Miller Clay Finished"
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_millerLoam.json > $OUTPUT_PATH/millerLoam_actors.txt
        echo "Miller Loam Finished"
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_millerSand.json > $OUTPUT_PATH/millerSand_actors.txt
        echo "Miller Sand Finished"
        ;;

    mizoguchi1990)
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_mizoguchi1990.json > $OUTPUT_PATH/mizoguchi_actors.txt
        ;;

    vanderborght2005)
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_vanderborght2005-exp1.json > $OUTPUT_PATH/vanderborght1_actors.txt
        echo "Vanderborght 1 Finished"
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_vanderborght2005-exp2.json > $OUTPUT_PATH/vanderborght2_actors.txt
        echo "Vanderborght 2 Finished"
        $SUMMA_EXE -g 1 -n 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_vanderborght2005-exp3.json > $OUTPUT_PATH/vanderborght3_actors.txt
        echo "Vanderborght 3 Finished"
        ;;
esac