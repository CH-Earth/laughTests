# Laugh-Tests
Testing Framework for SUMMA

Included is a `laugh_test.py` python program that controls the testing framework. This includes setting up the directories and running the tests. There is also tests included to verify the output of the tests.

## Installing the Laugh-Tests


## Installing the Laugh-Tests
To install the laugh tests, clone this git repo into the top level of your summa project. This is the directory that contains `bin/` and `build/`.

The types of tests that can be run. There are automated tests for normal `backward euler` that does not include sundials. There are automated tests for `sundials`. There are also automated tests setup for using `single_precision` forcing data and `double_precision` forcing data. These tests can be set up with the python script.

Make sure to add the `laugh_tests/` to the summa gitignore file to prevent the `laugh_tests/` from being commited to your summa repo.

Directory paths are set up by running the `init_test.sh` script. To get the directories to take effect run the script with `source init_tests.sh`.

The tests should now be set up and can be run with the various scripts provided in `run_scripts/non_actors/`.

## Running the Laugh-Tests
Make sure to run the laugh tests with the same terminal you ran `source init_tests.sh` with to ensure the bash variables take effect.

Inside `run_scripts/non_actors/` exists a script for running each of the laugh tets. There is also a `run_all.sh` file for running all the tests. This will take a couple hours to complete. There is also a `run_fast.sh` script for running a series of tests that should take around 5-10 seconds to complete.

All output from the tests should go to `/test_cases/output/non_actors/`



