#! /bin/bash

/SUMMA/bin/summa_sundials.exe -g 1 1 -m /SUMMA/Laugh-Tests/settings/syntheticTestCases/wigmosta1999/summa_fileManager_wigmosta1999_exp1_test.txt > wigmosta_1.txt &
/SUMMA/bin/summa_sundials.exe -g 1 1 -m /SUMMA/Laugh-Tests/settings/syntheticTestCases/wigmosta1999/summa_fileManager_wigmosta1999_exp2_test.txt > wigmosta_2.txt &