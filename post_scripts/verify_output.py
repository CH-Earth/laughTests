from os import listdir
from os.path import isfile, join
from pathlib import Path
import xarray as xr
import numpy as np


def verify_data(verified_data, compare_data, num_hru, output_variables):
    verified_dataset = xr.open_dataset(verified_data)
    to_compare_dataset = xr.open_dataset(compare_data)
    total_errors = 0
    for i_hru in range(0, num_hru):
        print("Checking HRU ", i_hru, " of ", num_hru, " for errors")
        verified_hru = verified_dataset.isel(hru=i_hru).copy()
        hru_to_compare = to_compare_dataset.isel(hru=i_hru).copy()
        for var in output_variables:
            error_counter = 0
            
            if len(verified_hru[var].values) != len(hru_to_compare[var].values):
                print("ERROR: output variable", var, "does not contain the same amount of data")
                print("     verified_hru = ", len(verified_hru[var].values))
                print("     hru_to_compare = ", len(hru_to_compare[var].values))
    
        
            verified_data_list = []
            to_verify_data = []
            if (verified_hru[var].values.ndim > 1):
                # 2D output case
                for list in verified_hru[var].values:
                    for data in list:
                        verified_data_list.append(data)
                
                for list in hru_to_compare[var].values:
                    for data in list:
                        to_verify_data.append(data)
            else:
                # 1D output case
                for data in verified_hru[var].values:
                    verified_data_list.append(data)
                
                for data in hru_to_compare[var].values:
                    to_verify_data.append(data)

                                
            # check length
            if len(to_verify_data) != len(to_verify_data):
                print("ERROR: output variable", var, "does not contain the same amount of data")
                print("     verified_hru = ", len(verified_data_list))
                print("     hru_to_compare = ", len(to_verify_data))

            # check values
            for elem in range(0, len(verified_data_list)):
                if abs(verified_data_list[elem] - to_verify_data[elem]) != 0:
                    error_counter += 1
            
            print("    Var", var, "has", error_counter, "errors")
            total_errors += error_counter
                
        
    return total_errors


"""
    This script is used to verify the output variables of the model.
"""
def get_output_vars(model_output_file):
    model_output_vars = []
    open_file = open(model_output_file, "r")
    lines = open_file.readlines()

    for line in lines:
        var = line.split(" ")[0]
        if var != "!" and var != "\n":
            model_output_vars.append(var)
    
    return model_output_vars




num_hru = 1
model_output_file = "/Summa-Actors/laugh_tests/test_cases/settings/meta/Model_Output.txt"
celia_non_actors = "/Summa-Actors/laugh_tests/test_cases/output/celia1990/summa_celia1990_double_precision_G1-1_timestep.nc"
celia_actors = "/Summa-Actors/laugh_tests/test_cases/output/celia1990/summaActors_celia1990_double_precisionGRU1-1_timestep.nc"
output_vars = get_output_vars(model_output_file)
total_errors_celia = verify_data(celia_non_actors, celia_actors, num_hru, output_vars)
print()

model_output_file = "/Summa-Actors/laugh_tests/test_cases/settings/syntheticTestCases/colbeck1976/summa_defineModelOutput.txt"
output_vars = get_output_vars(model_output_file)

colbeck1_non_actors = "/Summa-Actors/laugh_tests/test_cases/output/colbeck1976/summa_colbeck1976-exp1_double_precision_G1-1_timestep.nc"
colbeck1_actors = "/Summa-Actors/laugh_tests/test_cases/output/colbeck1976/summaActors_colbeck1976-exp1_double_precisionGRU1-1_timestep.nc"
total_errors_colbeck1 = verify_data(colbeck1_non_actors, colbeck1_actors, num_hru, output_vars)

colbeck2_non_actors = "/Summa-Actors/laugh_tests/test_cases/output/colbeck1976/summa_colbeck1976-exp2_double_precision_G1-1_timestep.nc"
colbeck2_actors = "/Summa-Actors/laugh_tests/test_cases/output/colbeck1976/summaActors_colbeck1976-exp2_double_precisionGRU1-1_timestep.nc"
total_errors_colbeck2 = verify_data(colbeck2_non_actors, colbeck2_actors, num_hru, output_vars)

colbeck3_non_actors = "/Summa-Actors/laugh_tests/test_cases/output/colbeck1976/summa_colbeck1976-exp3_double_precision_G1-1_timestep.nc"
colbeck3_actors = "/Summa-Actors/laugh_tests/test_cases/output/colbeck1976/summaActors_colbeck1976-exp3_double_precisionGRU1-1_timestep.nc"
total_errors_colbeck3 = verify_data(colbeck3_non_actors, colbeck3_actors, num_hru, output_vars)



print()
print("Finished Checking Output")
print("     Total Errors Celia = ", total_errors_celia)
print("     Total Errors Colbeck1 = ", total_errors_colbeck1)
print("     Total Errors Colbeck2 = ", total_errors_colbeck2)
print("     Total Errors Colbeck3 = ", total_errors_colbeck3)


