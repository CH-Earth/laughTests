#!/bin/bash

# This script runs the actor test cases. To use this script you have to set the SUMMA_EXE variable to where the summa_actors executable is located.
# This script runs one test as specified by the user from the input. To run multiple tests use the python script or create your own script that calls the cases within this script.

case $TEST in
    celia1990)
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_celia1990.json > $OUTPUT_PATH/celia_actors.txt
        ;;

    colbeck1976)
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_colbeck1976-exp1.json > $OUTPUT_PATH/colbeck1_actors.txt
        echo "Colbeck 1 Finished"
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_colbeck1976-exp2.json > $OUTPUT_PATH/colbeck2_actors.txt
        echo "Colbeck 2 Finished"
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_colbeck1976-exp3.json > $OUTPUT_PATH/colbeck3_actors.txt
        echo "Colbeck 3 Finished"
        ;;

    miller1998)
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_millerClay.json > $OUTPUT_PATH/millerClay_actors.txt
        echo "Miller Clay Finished"
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_millerLoam.json > $OUTPUT_PATH/millerLoam_actors.txt
        echo "Miller Loam Finished"
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_millerSand.json > $OUTPUT_PATH/millerSand_actors.txt
        echo "Miller Sand Finished"
        ;;

    mizoguchi1990)
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_mizoguchi1990.json > $OUTPUT_PATH/mizoguchi_actors.txt
        ;;

    vanderborght2005)
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_vanderborght2005-exp1.json > $OUTPUT_PATH/vanderborght1_actors.txt
        echo "Vanderborght 1 Finished"
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_vanderborght2005-exp2.json > $OUTPUT_PATH/vanderborght2_actors.txt
        echo "Vanderborght 2 Finished"
        $SUMMA_EXE -g 1 1 -c $FILE_MANAGER_PATH/Summa_Actors_Settings_vanderborght2005-exp3.json > $OUTPUT_PATH/vanderborght3_actors.txt
        echo "Vanderborght 3 Finished"
        ;;

    wigmosta1999)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManagers/fileManager_wigmosta1999-exp1_test.txt > $OUTPUT_PATH/wigmosta_1_actors.txt
        echo "Wigmosta 1 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManagers/fileManager_wigmosta1999-exp2_test.txt > $OUTPUT_PATH/wigmosta_2_actors.txt
        echo "Wigmosta 2 Finished"
        ;;

    reynolds_canopySrad_windPrfile_stomResist)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsDefault.json > $OUTPUT_PATH/summa_reynolds_Default.txt
        echo "Reynolds Default (canopySrad BeersLaw, windPrfile logBelowCanopy, stomResist BallBerry) Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsCLM2stream.txt > $OUTPUT_PATH/summa_reynolds_CLM_2stream.txt
        echo "Reynolds canopySrad CLM_2stream Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsUEB2stream.json > $OUTPUT_PATH/summa_reynolds_UEB_2stream.txt
        echo "Reynolds canopySrad UEB_2stream Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsNLscatter.json > $OUTPUT_PATH/summa_reynolds_NL_scatter.txt
        echo "Reynolds canopySrad NL_scatter Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsExponential.json > $OUTPUT_PATH/summa_reynolds_exponential.txt
        echo "Reynolds windPrfile exponential Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsJarvis.json > $OUTPUT_PATH/summa_reynolds_Jarvis.txt
        echo "Reynolds stomResist Jarvis Finished"
        ;;

    reynolds_groundwatr)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsLumpedQTopmodel.json > $OUTPUT_PATH/summa_reynolds_lumped_qTopmodl.txt
        echo "Reynolds groundwatr qTopmodl lumped Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsDistributedQTopmodel.json > $OUTPUT_PATH/summa_reynolds_distributed_qTopmodl.txt
        echo "Reynolds groundwatr qTopmodl distributed Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsNoXplict.json > $OUTPUT_PATH/summa_reynolds_noXplict.txt
        echo "Reynolds groundwatr noXplict Finished"
        ;;

    senator_alb_method)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_senatorConDecay.json > $OUTPUT_PATH/summa_senator_conDecay.txt
        echo "Senator alb_method conDecay Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_senatorVarDecay.json > $OUTPUT_PATH/summa_senator_varDecay.txt
        echo "Senator alb_method varDecay Finished"
        ;;

    umpqua_snowIncept)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_fileManager_umpquaStickySnow.json > $OUTPUT_PATH/summa_umpqua_stickySnow.txt
        echo "Umpqua snowIncept stickySnow Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_fileManager_umpquaLightSnow.json > $OUTPUT_PATH/summa_umpqua_lightSnow.txt
        echo "Umpqua snowIncept lightSnow Finished"
        ;;

esac

# MISSING WIGMOSTA AND WRR
