# Laugh-Tests
Testing Framework for SUMMA

We include some verification data to check the results of SUMMA simulations after changes to the code are made. This testing frameworks can also be used to compare two versions SUMMA together. Some python scripts for comparing the output of SUMMA are provided in the `post_scripts/` directory.

## Installing the Laugh-Tests
To install the laugh tests, clone this git repo into the top level of your summa project. This is the directory that contains `bin/` and `build/`.

Make sure to add the `laugh_tests/` to the summa gitignore file to prevent the `laugh_tests/` from being commited to your summa repo.

Directory paths are set up by running the `init_test.sh` script. To get the directories to take effect run the script with `source init_tests.sh`.

The tests should now be set up and can be run with the various scripts provided in `run_scripts/non-actors/`.

## Running the Laugh-Tests
Make sure to run the laugh tests with the same terminal you ran `source init_tests.sh` with to ensure the bash variables take effect.

Inside `run_scripts/non-actors/` exists a script for running each of the laugh tets. There is also a `run_all.sh` file for running all the tests. This will take a couple hours to complete. There is also a `run_fast.sh` script for running a series of tests that should take around 5-10 seconds to complete.

All output from the tests should go to `/test_cases/output/non-actors/`



