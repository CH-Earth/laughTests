#!/bin/bash

# This script runs the non-actor test cases. To use this script you have to set the SUMMA_EXE variable to where the summa executable is located.
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

    wigmosta1999)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_wigmosta1999-exp1_test.txt > $OUTPUT_PATH/summa_wigmosta_1.txt
        echo "Wigmosta 1 Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_wigmosta1999-exp2_test.txt > $OUTPUT_PATH/summa_wigmosta_2.txt
        echo "Wigmosta 2 Finished"
        ;;

    reynolds_canopySrad_windPrfile_stomResist)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsDefault.txt > $OUTPUT_PATH/summa_reynolds_Default.txt
        echo "Reynolds Default (canopySrad BeersLaw, windPrfile logBelowCanopy, stomResist BallBerry) Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsCLM2stream.txt > $OUTPUT_PATH/summa_reynolds_CLM_2stream.txt
        echo "Reynolds canopySrad CLM_2stream Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsUEB2stream.txt > $OUTPUT_PATH/summa_reynolds_UEB_2stream.txt
        echo "Reynolds canopySrad UEB_2stream Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsNLscatter.txt > $OUTPUT_PATH/summa_reynolds_NL_scatter.txt
        echo "Reynolds canopySrad NL_scatter Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsExponential.txt > $OUTPUT_PATH/summa_reynolds_exponential.txt
        echo "Reynolds windPrfile exponential Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsJarvis.txt > $OUTPUT_PATH/summa_reynolds_Jarvis.txt
        echo "Reynolds stomResist Jarvis Finished"
        ;;

    reynolds_groundwatr)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsLumpedQTopmodel.txt > $OUTPUT_PATH/summa_reynolds_lumped_qTopmodl.txt
        echo "Reynolds groundwatr qTopmodl lumped Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsDistributedQTopmodel.txt > $OUTPUT_PATH/summa_reynolds_distributed_qTopmodl.txt
        echo "Reynolds groundwatr qTopmodl distributed Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_reynoldsNoXplict.txt > $OUTPUT_PATH/summa_reynolds_noXplict.txt
        echo "Reynolds groundwatr noXplict Finished"
        ;;

    senator_alb_method)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_senatorConDecay.txt > $OUTPUT_PATH/summa_senator_conDecay.txt
        echo "Senator alb_method conDecay Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_senatorVarDecay.txt > $OUTPUT_PATH/summa_senator_varDecay.txt
        echo "Senator alb_method varDecay Finished"
        ;;

    umpqua_snowIncept)
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_fileManager_umpquaStickySnow.txt > $OUTPUT_PATH/summa_umpqua_stickySnow.txt
        echo "Umpqua snowIncept stickySnow Finished"
        $SUMMA_EXE -g 1 1 -m $FILE_MANAGER_PATH/fileManager_fileManager_umpquaLightSnow.txt > $OUTPUT_PATH/summa_umpqua_lightSnow.txt
        echo "Umpqua snowIncept lightSnow Finished"
        ;;

esac

# MISSING WIGMOSTA AND WRR

