import os
import sys
import subprocess
import xarray as xr


def initOutputDirs(test_cases):
    # List of tests to create output directories for
    outputDir = "./test_cases/output/"
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    for dir in test_cases:
        path = os.path.join(outputDir, dir)
        if not os.path.exists(path):
            os.mkdir(path)

def cleanOutputDirs(test_cases):
    # List of tests to create output directories for
    outputDir = "./test_cases/output/"
    for dir in test_cases:
        path = os.path.join(outputDir, dir)
        if os.path.exists(path):
            os.rmdir(path)



syntheticTestCases = ["celia1990", "colbeck1976", "miller1998", "mizoguchi1990", "vanderborght2005", "wigmosta1999"]
wrrPaperTestCases = ["figure01", "figure02", "figure03", "figure05", "figure06", "figure07", "figure08", "figure09"]

args = sys.argv
if len(args) < 2:
    print("Usage: python init_tests.py [init,clean]")
    sys.exit(1)
if args[1] == "clean":
    cleanOutputDirs()
elif args[1] == "init":
    initOutputDirs(syntheticTestCases)
else:
    print("Usage: python init_tests.py [init,clean]")


convertToFloat(syntheticTestCases)
# initOutputDirs(syntheticTestCases)
# initOutputDirs(wrrPaperTestCases)

# cleanOutputDirs(syntheticTestCases)
# cleanOutputDirs(wrrPaperTestCases)