from os import listdir
from os.path import isfile, join
from pathlib import Path
import xarray as xr
import numpy as np


def verify_data(verified_data, compare_data, num_hru, output_variables):
    try:
        verified_dataset = xr.open_dataset(verified_data)
        to_compare_dataset = xr.open_dataset(compare_data)
    except FileNotFoundError:
        print("Check the variables \'verified_data\' and \'compare_data\'. They may not point to the correct output files or the output filenames may have changed.")
        print("Verified_path =", verified_data)
        print("Compare_path =", compare_data)
        exit()
    total_errors = 0
    for i_hru in range(0, num_hru):
        verified_hru = verified_dataset.isel(hru=i_hru).copy()
        hru_to_compare = to_compare_dataset.isel(hru=i_hru).copy()
        for var in output_variables:
            error_counter = 0
            try:
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
                if len(verified_data_list) != len(to_verify_data):
                    print("ERROR: output variable", var, "does not contain the same amount of data")
                    print("     verified_hru = ", len(verified_data_list))
                    print("     hru_to_compare = ", len(to_verify_data))

                # check values
                for elem in range(0, len(verified_data_list)):
                    if abs(verified_data_list[elem] - to_verify_data[elem]) != 0:
                        # print("variable -",var, "has different values at", elem)
                        # print("     verified_hru = ", verified_data[elem])
                        # print("     hru_to_compare = ", to_verify_data[elem])
                        error_counter += 1
                
                print("    Var", var, "has", error_counter, "errors")
                total_errors += error_counter
                        
            except (TypeError, KeyError) as error:
                print(error)
                # print("variable - ", var, "Cannot be compared with len")
                # print("     verified_hru = ",verified_hru[var].values)
                # print("     hru_to_compare = ", hru_to_compare[var].values)
        
        return total_errors


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
actors_output_path = "/Summa-Actors/laugh_tests/test_cases/output/actors"
non_actors_output_path = "/SUMMA/laugh_tests/test_cases/output/non-actors"
verification_dir_path = "../validation_data/m2_mac_Sept23"


print("Checking Celia")
celia_non_actors = "/SUMMA/laugh_tests/test_cases/output/non-actors/celia1990/celia1990_G1-1_timestep.nc"
celia_actors = "/Summa-Actors/laugh_tests/test_cases/output/actors/celia1990/celia1990GRU1-1_timestep.nc"

model_output_file = "../test_cases/settings/meta/Model_Output.txt"

output_vars = get_output_vars(model_output_file)
celia_errors = verify_data(celia_non_actors, celia_actors, num_hru, output_vars)
print()


print("Checking Colbeck_1")
colbeck_1_non_actors = "/SUMMA/laugh_tests/test_cases/output/non-actors/colbeck1976/colbeck1976-exp1_G1-1_timestep.nc"
colbeck_1_actors = "/Summa-Actors/laugh_tests/test_cases/output/actors/colbeck1976/colbeck1976-exp1GRU1-1_timestep.nc"
model_output_file = "../test_cases/settings/syntheticTestCases/colbeck1976/common/summa_defineModelOutput.txt"
output_vars = get_output_vars(model_output_file)
colbeck_1_errors = verify_data(colbeck_1_non_actors, colbeck_1_actors, num_hru, output_vars)
print()


print("Checking Colbeck_2")
colbeck_2_non_actors = "/SUMMA/laugh_tests/test_cases/output/non-actors/colbeck1976/colbeck1976-exp2_G1-1_timestep.nc"
colbeck_2_actors = "/Summa-Actors/laugh_tests/test_cases/output/actors/colbeck1976/colbeck1976-exp2GRU1-1_timestep.nc"
model_output_file = "../test_cases/settings/syntheticTestCases/colbeck1976/common/summa_defineModelOutput.txt"
output_vars = get_output_vars(model_output_file)
colbeck_2_errors = verify_data(colbeck_2_non_actors, colbeck_2_actors, num_hru, output_vars)
print()


print("Checking Colbeck_3")
colbeck_3_non_actors = "/SUMMA/laugh_tests/test_cases/output/non-actors/colbeck1976/colbeck1976-exp3_G1-1_timestep.nc"
colbeck_3_actors = "/Summa-Actors/laugh_tests/test_cases/output/actors/colbeck1976/colbeck1976-exp3GRU1-1_timestep.nc"
model_output_file = "../test_cases/settings/syntheticTestCases/colbeck1976/common/summa_defineModelOutput.txt"
output_vars = get_output_vars(model_output_file)
colbeck_3_errors = verify_data(colbeck_3_non_actors, colbeck_3_actors, num_hru, output_vars)
print()


print("Checking Miller_loam")
miller_loam_non_actors = "/SUMMA/laugh_tests/test_cases/output/non-actors/miller1998/millerLoam_G1-1_timestep.nc"
miller_loam_actors = "/Summa-Actors/laugh_tests/test_cases/output/actors/miller1998/millerLoamGRU1-1_timestep.nc"
model_output_file = "../test_cases/settings/meta/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
miller_errors = verify_data(verified_path, compare_path, num_hru, output_vars)
print()


print("Checking Mizoguchi")
mizoguchi = "/mizoguchi1990/mizoguchi1990_G1-1_timestep.nc"
verified_path = verification_dir_path + mizoguchi
compare_path = actors_output_path + mizoguchi
model_output_file = "./test_cases/settings/meta/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
mizoguchi_errors = verify_data(verified_path, compare_path, num_hru, output_vars)
print()


print("Checking Vanderborght_1")
vanderborght_1 = "/vanderborght2005/vanderborght2005_exp1_G1-1_timestep.nc"
verified_path = verification_dir_path + vanderborght_1
compare_path = actors_output_path + vanderborght_1
model_output_file = "./test_cases/settings/syntheticTestCases/vanderborght2005/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
vanderborght_1_erros = verify_data(verified_path, compare_path, num_hru, output_vars)
print()


print("Checking Vanderborght_2")
vanderborght_2 = "/vanderborght2005/vanderborght2005_exp2_G1-1_timestep.nc"
verified_path = verification_dir_path + vanderborght_2
compare_path = actors_output_path + vanderborght_2
model_output_file = "./test_cases/settings/syntheticTestCases/vanderborght2005/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
vanderborght_2_erros = verify_data(verified_path, compare_path, num_hru, output_vars)
print()


print("Checking Vanderborght_ 2")
vanderborght_3 = "/vanderborght2005/vanderborght2005_exp3_G1-1_timestep.nc"
verified_path = verification_dir_path + vanderborght_3
compare_path = actors_output_path + vanderborght_3
model_output_file = "./test_cases/settings/syntheticTestCases/vanderborght2005/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
vanderborght_3_erros = verify_data(verified_path, compare_path, num_hru, output_vars)
print()
print("Finished Checking Output")
print("     Total Errors For Celia", celia_errors)
print("     Total Errors For Colbeck_1", colbeck_1_errors)
print("     Total Errors For Colbeck_2", colbeck_2_errors)
print("     Total Errors For Colbeck_3", colbeck_3_errors)
print("     Total Errors For Miller_Loam", miller_errors)
print("     Total Errors For Mizoguchi", mizoguchi_errors)
print("     Total Errors For Vanderborght", vanderborght_1_erros)
print("     Total Errors For Vanderborght", vanderborght_2_erros)
print("     Total Errors For Vanderborght", vanderborght_3_erros)

