import os
import sys
import subprocess
import xarray as xr


def initOutputDirs(test_cases):
    outputDir = "./test_cases/output/"
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    for dir in test_cases:
        path = os.path.join(outputDir, dir)
        if not os.path.exists(path):
            os.mkdir(path)


def cleanOutputDirs(test_cases):
    outputDir = "./test_cases/output/"
    for dir in test_cases:
        path = os.path.join(outputDir, dir)
        if os.path.exists(path):
            os.rmdir(path)


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


def setupBE(test_cases, precision, summa_dir):
    settings_dir = "./test_cases/settings/syntheticTestCases/"

    for dir in test_cases:
        path = os.path.join(settings_dir, dir)
        print(path)
        if os.path.exists(path):
            print("Setting up " + dir)
            file_managers = os.path.join(path, "fileManagers")
            if precision == "single":
                precision = "single_precision.txt"
            else:
                precision = "double_precision.txt"
            for file_manager in os.listdir(file_managers):
                
                file_manager = os.path.join(file_managers, file_manager)
                with open(file_manager, 'r') as f:
                    data = f.read()
                    data = data.replace('<SOLVER>', 'be')
                    data = data.replace('<BASEDIR>', summa_dir)
                    data = data.replace('<PRECISION>', precision)
                
                file_manager = os.path.splitext(file_manager)[0]+'_test.txt'
                with open(file_manager, 'w') as f:
                    f.write(data)

                

syntheticTestCases = ["celia1990", "colbeck1976", "miller1998", "mizoguchi1990", "vanderborght2005"]
wrrPaperTestCases = ["figure01", "figure02", "figure03", "figure05", "figure06", "figure07", "figure08", "figure09"]

os.environ["SUMMA_DIR"] = "/home/summa/summa/build/summa.exe"
# args = sys.argv
# if len(args) < 2:
#     print("Usage: python init_tests.py [init,clean]")
#     sys.exit(1)
# if args[1] == "clean":
#     # cleanOutputDirs()
#     cleanFileManagers(syntheticTestCases)
# elif args[1] == "init":
#     initOutputDirs(syntheticTestCases)

# if args[1] == "be" and args[2] == "single":
#     setupBE(syntheticTestCases, "single", args[3])

# if args[1] == "be" and args[2] == "double":
#     setupBE(syntheticTestCases, "double", args[3])




# initOutputDirs(syntheticTestCases)
# initOutputDirs(wrrPaperTestCases)

# cleanOutputDirs(syntheticTestCases)
# cleanOutputDirs(wrrPaperTestCases)