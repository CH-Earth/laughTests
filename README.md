# Laugh-Tests
Testing Framework for SUMMA

Included is a `laugh_test.py` python program that controls the testing framework. This includes setting up the directories and running the tests. There is also tests included to verify the output of the tests.

## Installing the Laugh-Tests
To get a copy of the laugh tests, simply clone the repo to your system you wish to run the laugh_tests on. In previous configurations of this testing framework you were required to put this directory at the top level of your summa directory. THIS IS NO LONGER THE CASE. This repo can be cloned anywhere on your system. There is a configuration file `test_settings.json` that is used for controlling the testing framework.

## Running the Laugh-Tests
To run the laugh tests first you must install and compile SUMMA. This can be done with the documentation found here https://summa.readthedocs.io/en/latest/SUMMA_documentation/.

Once you have installed summa you can now set up the laugh tests. The laugh tests are configured through the `test_settings.json` file included in the laugh_test repository. This file will need to be modified to work in your environment and to run the specific tests you desire. The settings are configured by the user by modifying the value of the key-value pairs in the `test_settings.json` file. Detailed explanations of the key-value pairs are provided below:

 - Test_List: This is a the list of test you wish to run. The directories and file manager files will automatically be set up by running the `laugh_test.py` file which instructions can be found in the `laugh_test.py` section below
    - Example: "Test_List": ["celia1990", "colbeck1976", "miller1998", "mizoguchi1990", "vanderborght2005", "wigmosta1999"]

 - Precision: Either use "single" or "double". This specifies what precision of forcing data to use for the laugh tests.
    - Example: "Precision": "single",

 - Solver: Either use "numrec", "ida" for V4 and "v3" for older implementaitons of SUMMA. This specifies what solver to use for the solution to summa. The file manager file will automatically be set up from your selection here.
    - Example: "Solver": "numrec",

 - Summa_exe: Provide a full path to the summa executable that you compiled.
    - "Summa_Exe": "/absolute/path/to/summa.exe"

 - Version: actors or non-actors versions of the programs, they have some differences in setting up
    - "Version": "non-actors"

### Using `laugh_tests.py`
Once you have set up the `test_settings.json` file. You then have to configure your setup with the `laugh_tests.py` file. You run it with the following `python3 laugh_tests.py [options]`. Below are a list of options:
- init
- run
- reset
- clean

The first command to use is `python3 laugh_tests.py init`. This will set up all of the output directories and the file managers. After running init the next command to run is `python3 laugh_tests.py run`. This will run all of the tests that you specified in the Test_List key in the `test_settings.json` file.

If you plan to run more tests with multiple configurations use the `reset` option followed by `init` again. If you wish to delete everything including your generated output use `clean`

### EXAMPLE TEST_SETTINGS.json
```json
{
    "Test_List": ["celia1990", "colbeck1976", "miller1998", "mizoguchi1990", "vanderborght2005", "wigmosta1999"],
    
    "Precision": "single",
    
    "Solver": "numrec",

    "Summa_Exe": "/absolute/path/to/summa.exe"

    "Version": "non-actors"
}
```





