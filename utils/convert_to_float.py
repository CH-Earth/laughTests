import os
import sys
import subprocess
import xarray as xr


"""
Converts double precision netCDF files to single precision
"""
def convertToFloat(test_cases):
    inputDir = "../test_cases/input_data/"
    for dir in test_cases:
        path = os.path.join(inputDir, dir)
        if os.path.exists(path):
            for file in os.listdir(path):
                new_file = file.replace("_double_precision", "_single_precision")
                temp_file = file.replace("_double_precision", "_temp")
                if file.endswith("_double_precision.nc"):                    
                    subprocess.call(["cdo", "-b", "f32", "copy", os.path.join(path, file), os.path.join(path, temp_file)])
                    dataset = xr.open_dataset(os.path.join(path, temp_file))
                    if dataset["hruId"].values.dtype == "float32":
                        subprocess.call(["ncap2", "-s", "hruId=int(hruId)", os.path.join(path, temp_file), os.path.join(path, new_file)])
                        os.remove(os.path.join(path, temp_file))

syntheticTestCases = ["celia1990", "colbeck1976", "miller1998", "mizoguchi1990", "vanderborght2005", "wigmosta1999"]

convertToFloat(syntheticTestCases)