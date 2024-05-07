import os
import sys
import subprocess
import json
import re

"""
Create all of the output directories for the test cases
"""
def initOutputDirs(test_cases):
    outputDir = "./test_cases/output/"
    for test in test_cases:
        path = os.path.join(outputDir, test["name"])
        os.makedirs(path, exist_ok=True)

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

def determinePrecisionParameters(precision):
    if precision == "single":
        precision = "single_precision.txt"
        output = "single_precision"
    elif precision == "double":
        precision = "double_precision.txt"
        output = "double_precision"
    else:
        throw("Invalid precision")
    return precision, output

def getVersionParameters(version):
    if version == "non-actors":
        version = "summa"
    elif version == "actors":
        version = "summaActors"
    else:
        throw("Invalid version")
    return version

def setActorConfig(file_manager, version, path, specific_expirement):
    print(specific_expirement)
    if version == "non-actors":
        return
    elif version == "actors":
        for settings_file in os.listdir(path):
            if settings_file.endswith("{}.json".format(specific_expirement[0])):
                settings_file = os.path.join(path, settings_file)
                with open(settings_file, 'r') as f:
                    data = json.load(f)
                    data["Job_Actor"]["file_manager_path"] = os.path.abspath(file_manager)
                with open(settings_file, 'w') as f:
                    json.dump(data, f, indent=4)
    else:
        throw("Invalid version")



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
            sub_test = re.findall("fileManager_([^._]+).txt", file_manager)
            file_manager = os.path.join(file_managers_path, file_manager)
            print("\tSetting sub-test " + sub_test[0])
            with open(file_manager, 'r') as f:
                data = f.read()
                data = data.replace('<SOLVER>', solver)
                data = data.replace('<BASEDIR>', test_settings_dir)
                data = data.replace('<PRECISION>', precision)
            file_manager = os.path.splitext(file_manager)[0]+'_test.txt'
            with open(file_manager, 'w') as f:
                f.write(data)
                            

def getRunScript(version, base_dir):
    if version == "non-actors":
        run_script = os.path.join(base_dir, "run_scripts/run_non_actors.sh")
    elif version == "actors":
        run_script = os.path.join(base_dir, "run_scripts/run_actors.sh")
    else:
        print("Invalid version")
        sys.exit()
    return run_script

"""
Run the test cases
"""
def runTest(test_cases, summa_exe, version):
    base_dir = os.getcwd()
    run_script = getRunScript(version, base_dir)
    output_path = os.path.join(base_dir, "test_cases/output")
    settings_dir = os.path.join(base_dir, "test_cases/settings/syntheticTestCases")

    for dir in test_cases:
        print("Starting Test: " + dir)
        test_output = os.path.join(output_path, dir)
        path = os.path.join(settings_dir, dir)
        if os.path.exists(path):
            file_manager_path = os.path.join(path, "fileManagers")
            os.environ["SUMMA_EXE"] = summa_exe
            if version == "actors":
                os.environ["FILE_MANAGER_PATH"] = path
            else:
                os.environ["FILE_MANAGER_PATH"] = file_manager_path
            os.environ["OUTPUT_PATH"] = test_output
            os.environ["TEST"] = dir
            subprocess.call(run_script)
            print("Finished Test: " + dir)


def main():
    with open("test_settings.json", "r") as settings_file:
        settings = json.load(settings_file)
        if 'Precision' not in settings:
            settings['Precision'] = 'double'
        if 'all' in settings['Test_List']:
            with open("test_inventory.json", "r") as inventory_file:
                inventory = json.load(inventory_file)
                settings['Test_List'] = inventory['tests']
        
    # print(settings)

    if len(sys.argv) == 2:
        if sys.argv[1] == "init":
            print("Initializing test cases")
            initOutputDirs(settings['Test_List'])
            setup(settings['Test_List'], settings['Precision'], settings['Solver'])
        elif sys.argv[1] == "clean":
            print("Cleaning output directories")
            cleanFileManagers(settings['Test_List'])
        else:
            print("Invalid argument")
            print("Valid arguments are: init, run, reset, clean")
    
             

        # test_cases = test_settings["Test_List"]
        # precision = test_settings["Precision"]
        # solver = test_settings["Solver"]
        # summa_exe = test_settings["Summa_Exe"]
        # version = test_settings["Version"]

#     if len(sys.argv) == 2:
#         if sys.argv[1] == "init":
#             initOutputDirs(test_cases)
#             setup(test_cases, precision, solver, version)
#         elif sys.argv[1] == "clean":
#             cleanFileManagers(test_cases)
#         elif sys.argv[1] == "reset":
#             cleanFileManagers(test_cases)
#         elif sys.argv[1] == "run":
#             runTest(test_cases, summa_exe, version)
#         elif sys.argv[1] == "check":
#             checkOutput()
#     else:
#         print("Invalid argument")
#         print("Valid arguments are: init, run, reset, clean")

if __name__ == "__main__":
    main()