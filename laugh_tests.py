import os
import sys
import subprocess
import xarray as xr
import json

"""
Create all of the output directories for the test cases
"""
def initOutputDirs(test_cases):
    outputDir = "./test_cases/output/"
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    for dir in test_cases:
        path = os.path.join(outputDir, dir)
        if not os.path.exists(path):
            os.mkdir(path)

"""
Clean up the output directories for the test cases
"""
def cleanOutputDirs(test_cases):
    outputDir = "./test_cases/output/"
    for dir in test_cases:
        path = os.path.join(outputDir, dir)
        if os.path.exists(path):
            os.rmdir(path)

"""
Clean up the file managers for the test cases
This cleans up the generated file managers for the test cases
"""
def cleanFileManagers(test_cases):
    settings_dir = "./test_cases/settings/syntheticTestCases/"
    for dir in test_cases:
        path = os.path.join(settings_dir, dir)
        if os.path.exists(path):
            file_managers = os.path.join(path, "fileManagers")
            for file_manager in os.listdir(file_managers):
                if file_manager.endswith("_test.txt"):
                    file_manager = os.path.join(file_managers, file_manager)
                    os.remove(file_manager)

"""
Setup the file managers for the test cases
"""
def setupBE(test_cases, precision, solver):
    basedir = os.getcwd()
    settings_dir = "./test_cases/settings/syntheticTestCases/"
    if precision == "single":
        precision = "single_precision.txt"
    else:
        precision = "double_precision.txt"
    for dir in test_cases:
        path = os.path.join(settings_dir, dir)
        print(path)
        if os.path.exists(path):
            print("Setting up " + dir)
            file_managers = os.path.join(path, "fileManagers")

            for file_manager in os.listdir(file_managers):
                
                file_manager = os.path.join(file_managers, file_manager)
                with open(file_manager, 'r') as f:
                    data = f.read()
                    data = data.replace('<SOLVER>', solver)
                    data = data.replace('<BASEDIR>', basedir)
                    data = data.replace('<PRECISION>', precision)
                
                file_manager = os.path.splitext(file_manager)[0]+'_test.txt'
                with open(file_manager, 'w') as f:
                    f.write(data)
"""
Run the test cases
"""
def runTest(test_cases, summa_exe):
    base_dir = os.getcwd()
    run_script = os.path.join(base_dir, "run_scripts/non-actors/run_non_actors.sh")
    output_path = os.path.join(base_dir, "test_cases/output")
    settings_dir = os.path.join(base_dir, "test_cases/settings/syntheticTestCases")

    for dir in test_cases:
        print("Starting Test: " + dir)
        test_output = os.path.join(output_path, dir)
        path = os.path.join(settings_dir, dir)
        if os.path.exists(path):
            file_manager_path = os.path.join(path, "fileManagers")
            os.environ["SUMMA_EXE"] = summa_exe
            os.environ["FILE_MANAGER_PATH"] = file_manager_path
            os.environ["OUTPUT_PATH"] = test_output
            os.environ["TEST"] = dir
            subprocess.call(run_script)
            print("Finished Test: " + dir)


# Start of main
with open("test_settings.json", "r") as test_settings_file:
    test_settings = json.load(test_settings_file)
    test_cases = test_settings["Test_List"]
    precision = test_settings["Precision"]
    solver = test_settings["Solver"]
    summa_exe = test_settings["Summa_Exe"]

if len(args) < 2:
    if args[1] == "init":
        initOutputDirs(test_cases)
        setupBE(test_cases, precision, solver)
    elif args[1] == "clean":
        cleanOutputDirs(test_cases)
        cleanFileManagers(test_cases)
    elif args[1] == "reset":
        cleanFileManagers(test_cases)
    elif args[1] == "run":
        runTest(test_cases)
    else:
        print("Invalid argument")
        print("Valid arguments are: init, run, reset, clean")