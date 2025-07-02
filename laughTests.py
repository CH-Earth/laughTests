import os
import sys
import subprocess
import json
import re
import shutil

"""
Create all of the output directories for the test cases
"""
def initOutputDirs(test_cases):
    outputDir = "./test_cases/output/"
    cleanOutputDirs(test_cases)
    for test in test_cases:
        path = os.path.join(outputDir, test["name"])
        os.makedirs(path, exist_ok=True)

"""
Clean up the output directories for the test cases
"""
def cleanOutputDirs(test_cases):
    print("Cleaning output directories")
    outputDir = "./test_cases/output/"
    for test in test_cases:
        path = os.path.join(outputDir, test["name"])
        if os.path.exists(path):
            shutil.rmtree(path)  # Recursively delete directory and contents
            os.makedirs(path)    # Recreate the empty directory

"""
Clean up the file managers for the test cases
This cleans up the generated file managers for the test cases
"""
def cleanFileManagers(test_cases):
    print("Cleaning fileManagers")
    settings_dir = "./test_cases/settings/"
    for test in test_cases:
        path = os.path.join(settings_dir, test["type"])
        path = os.path.join(path, test["name"])
        if os.path.exists(path):
            file_managers = os.path.join(path, "fileManagers")
            for file_manager in os.listdir(file_managers):
                if file_manager.endswith("_test.txt"):
                    file_manager = os.path.join(file_managers, file_manager)
                    os.remove(file_manager)
"""
Setup the file managers for the test cases
"""
def setup(test_cases, precision, solver):
    # Form full path to the test_cases dir
    basedir = os.getcwd() + "/test_cases"
    settings_dir = "./test_cases/settings/"
    cleanFileManagers(test_cases)

    for test in test_cases:
        test_settings_dir = os.path.join(settings_dir, test["type"])
        test_settings_dir = os.path.join(test_settings_dir, test["name"])
        if not os.path.exists(test_settings_dir):
            print("Test case does not exist: " + test["name"])
            continue
        print("Setting up " + test["name"])
        file_managers_path = os.path.join(test_settings_dir, "fileManagers")
        for file_manager in os.listdir(file_managers_path):
            if file_manager.startswith('.'):
                continue
            sub_test = re.findall("fileManager_([^._]+).txt", file_manager)
            file_manager = os.path.join(file_managers_path, file_manager)
            print("\tSetting sub-test " + sub_test[0])
            with open(file_manager, 'r') as f:
                data = f.read()
                data = data.replace('<SOLVER>', solver)
                data = data.replace('<BASEDIR>', basedir)
                data = data.replace('<PRECISION>', precision + '_precision.txt')
            file_manager = os.path.splitext(file_manager)[0]+'_test.txt'
            with open(file_manager, 'w') as f:
                f.write(data)
                            
"""
Run the test cases
"""
def runTest(test_cases, summa_exe):
    print("Running test cases")
    base_dir = os.getcwd()
    output_path = os.path.join(base_dir, "test_cases/output")
    print(base_dir)
    print(output_path)
    for exe_counter, exe in enumerate(summa_exe):
        print("\n\nRunning executable", exe)
        suffix = "exe_" + str(exe_counter)
        for test in test_cases:
            num_gru = "25" if test["type"] == "multiGruTestCases" else "1"
            settings_dir = os.path.join(base_dir, "test_cases/settings")
            settings_dir = os.path.join(settings_dir, test["type"])
            settings_dir = os.path.join(settings_dir, test["name"])
            print("Running", test["name"])
            file_managers_path = os.path.join(settings_dir, "fileManagers")
            for file_manager in os.listdir(file_managers_path):
                if file_manager.endswith("_test.txt"):
                    file_manager = os.path.join(file_managers_path, file_manager)
                    print("\tRunning file manager " + file_manager)
                    sub_test = re.findall("fileManager_([^_]+)_test.txt", file_manager)
                    log_file_path = os.path.join(output_path, test["name"], "log_" + sub_test[0] + "_" + suffix + ".txt")
                    with open(log_file_path, 'w') as log_file:
                        result = subprocess.run([exe, "-s", suffix, "-g", "1", num_gru, "-m", file_manager], stdout=log_file, check=False)
                    print(f"Test {test['name']} exit code was: {result.returncode}")
                        
def main():
    with open("settings.json", "r") as settings_file:
        settings = json.load(settings_file)
        if 'Precision' not in settings:
            settings['Precision'] = 'double'
        if 'all' in settings['Test_List']:
            with open("test_inventory.json", "r") as inventory_file:
                inventory = json.load(inventory_file)
                settings['Test_List'] = inventory['tests']
        else:
            with open("test_inventory.json", "r") as inventory_file:
                inventory = json.load(inventory_file)
                new_test_list = []
                for test in settings['Test_List']:
                    if test not in ['syntheticTestCases', 'multiGruTestCases', 'wrrPaperTestCases']:
                        for inventory_test in inventory['tests']:
                            if test == inventory_test['name']:
                                new_test_list.append(inventory_test)
                    else:
                        for inventory_test in inventory['tests']:
                            if inventory_test['type'] == test:
                                if inventory_test not in new_test_list:
                                    new_test_list.append(inventory_test)
                settings['Test_List'] = new_test_list
        
    if len(sys.argv) == 2:
        if sys.argv[1] == "init":
            print("Initializing test cases")
            initOutputDirs(settings['Test_List'])
            setup(settings['Test_List'], settings['Precision'], settings['Solver'])
        elif sys.argv[1] == "clean":
            cleanOutputDirs(settings['Test_List'])
            cleanFileManagers(settings['Test_List'])
        elif sys.argv[1] == "run":
            runTest(settings['Test_List'], settings['Executables'])
        else:
            print("Invalid argument")
            print("Valid arguments are: init, run, clean")
    
if __name__ == "__main__":
    main()