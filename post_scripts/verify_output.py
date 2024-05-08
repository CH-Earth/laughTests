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
            
            # print("    Var", var, "has", error_counter, "errors")
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
model_output_file = "../test_cases/settings/meta/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
celia_exe_0 = "../test_cases/output/celia1990/celia1990_exe_0_G1-1_timestep.nc"
celia_exe_1 = "../test_cases/output/celia1990/celia1990_exe_1_G1-1_timestep.nc"
total_errors_celia = verify_data(celia_exe_0, celia_exe_1, num_hru, output_vars)

model_output_file = "../test_cases/settings/syntheticTestCases/colbeck1976/summa_defineModelOutput.txt"
output_vars = get_output_vars(model_output_file)
colbeck_exe_0_exp1 = "../test_cases/output/colbeck1976/colbeck1976-exp1_exe_0_G1-1_timestep.nc"
colbeck_exe_1_exp1 = "../test_cases/output/colbeck1976/colbeck1976-exp1_exe_1_G1-1_timestep.nc"
total_errors_colbeck_exp1 = verify_data(colbeck_exe_0_exp1, colbeck_exe_1_exp1, num_hru, output_vars)

colbeck_exe_0_exp2 = "../test_cases/output/colbeck1976/colbeck1976-exp2_exe_0_G1-1_timestep.nc"
colbeck_exe_1_exp2 = "../test_cases/output/colbeck1976/colbeck1976-exp2_exe_1_G1-1_timestep.nc"
total_errors_colbeck_exp2 = verify_data(colbeck_exe_0_exp2, colbeck_exe_1_exp2, num_hru, output_vars)

colbeck_exe_0_exp3 = "../test_cases/output/colbeck1976/colbeck1976-exp3_exe_0_G1-1_timestep.nc"
colbeck_exe_1_exp3 = "../test_cases/output/colbeck1976/colbeck1976-exp3_exe_1_G1-1_timestep.nc"
total_errors_colbeck_exp3 = verify_data(colbeck_exe_0_exp3, colbeck_exe_1_exp3, num_hru, output_vars)

model_output_file = "../test_cases/settings/meta/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
miller_clay_exe_0 = "../test_cases/output/miller1998/millerClay_exe_0_G1-1_timestep.nc"
miller_clay_exe_1 = "../test_cases/output/miller1998/millerClay_exe_1_G1-1_timestep.nc"
total_errors_miller_clay = verify_data(miller_clay_exe_0, miller_clay_exe_1, num_hru, output_vars)

miller_loam_exe_0 = "../test_cases/output/miller1998/millerLoam_exe_0_G1-1_timestep.nc"
miller_loam_exe_1 = "../test_cases/output/miller1998/millerLoam_exe_1_G1-1_timestep.nc"
total_errors_miller_loam = verify_data(miller_loam_exe_0, miller_loam_exe_1, num_hru, output_vars)

miller_sand_exe_0 = "../test_cases/output/miller1998/millerSand_exe_0_G1-1_timestep.nc"
miller_sand_exe_1 = "../test_cases/output/miller1998/millerSand_exe_1_G1-1_timestep.nc"
total_errors_miller_sand = verify_data(miller_sand_exe_0, miller_sand_exe_1, num_hru, output_vars)

mizoguchi_exe_0 = "../test_cases/output/mizoguchi1990/mizoguchi1990_exe_0_G1-1_timestep.nc"
mizoguchi_exe_1 = "../test_cases/output/mizoguchi1990/mizoguchi1990_exe_1_G1-1_timestep.nc"
total_errors_mizoguchi = verify_data(mizoguchi_exe_0, mizoguchi_exe_1, num_hru, output_vars)

model_output_file = "../test_cases/settings/syntheticTestCases/vanderborght2005/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
vaderborght1_exe_0 = "../test_cases/output/vanderborght2005/vanderborght2005_exp1_exe_0_G1-1_timestep.nc"
vaderborght1_exe_1 = "../test_cases/output/vanderborght2005/vanderborght2005_exp1_exe_1_G1-1_timestep.nc"
total_errors_vanderborght1 = verify_data(vaderborght1_exe_0, vaderborght1_exe_1, num_hru, output_vars)

vaderborght2_exe_0 = "../test_cases/output/vanderborght2005/vanderborght2005_exp2_exe_0_G1-1_timestep.nc"
vaderborght2_exe_1 = "../test_cases/output/vanderborght2005/vanderborght2005_exp2_exe_1_G1-1_timestep.nc"
total_errors_vanderborght2 = verify_data(vaderborght2_exe_0, vaderborght2_exe_1, num_hru, output_vars)

vaderborght3_exe_0 = "../test_cases/output/vanderborght2005/vanderborght2005_exp3_exe_0_G1-1_timestep.nc"
vaderborght3_exe_1 = "../test_cases/output/vanderborght2005/vanderborght2005_exp3_exe_1_G1-1_timestep.nc"
total_errors_vanderborght3 = verify_data(vaderborght3_exe_0, vaderborght3_exe_1, num_hru, output_vars)

num_hru = 50
model_output_file = "../test_cases/settings/meta/Model_Output.txt"
output_vars = get_output_vars(model_output_file)
wigmosta_exe_0_exp_1 = "../test_cases/output/wigmosta1999/wigmosta1999_exp1_exe_0_G1-1_timestep.nc" 
wigmosta_exe_1_exp_1 = "../test_cases/output/wigmosta1999/wigmosta1999_exp1_exe_1_G1-1_timestep.nc"
total_errors_wigmosta_exp1 = verify_data(wigmosta_exe_0_exp_1, wigmosta_exe_1_exp_1, num_hru, output_vars)

wigmosta_exe_0_exp_2 = "../test_cases/output/wigmosta1999/wigmosta1999_exp2_exe_0_G1-1_timestep.nc"
wigmosta_exe_1_exp_2 = "../test_cases/output/wigmosta1999/wigmosta1999_exp2_exe_1_G1-1_timestep.nc"
total_errors_wigmosta_exp2 = verify_data(wigmosta_exe_0_exp_2, wigmosta_exe_1_exp_2, num_hru, output_vars)

num_hru = 1
reynolds_canopy_clm2stream_exe_0 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsCLM2stream_exe_0_G1-1_timestep.nc"
reynolds_canopy_clm2stream_exe_1 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsCLM2stream_exe_1_G1-1_timestep.nc"
total_errors_reynolds_canopy_clm2stream = verify_data(reynolds_canopy_clm2stream_exe_0, reynolds_canopy_clm2stream_exe_1, num_hru, output_vars)

reynolds_default_exe_0 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsDefault_exe_0_G1-1_timestep.nc"
reynolds_default_exe_1 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsDefault_exe_1_G1-1_timestep.nc"
total_errors_reynolds_default = verify_data(reynolds_default_exe_0, reynolds_default_exe_1, num_hru, output_vars)

reynolds_exponential_exe_0 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsExponential_exe_0_G1-1_timestep.nc"
reynolds_exponential_exe_1 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsExponential_exe_1_G1-1_timestep.nc"
total_errors_reynolds_exponential = verify_data(reynolds_exponential_exe_0, reynolds_exponential_exe_1, num_hru, output_vars)

reynolds_jarvis_exe_0 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsJarvis_exe_0_G1-1_timestep.nc"
reynolds_jarvis_exe_1 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsJarvis_exe_1_G1-1_timestep.nc"
total_errors_reynolds_jarvis = verify_data(reynolds_jarvis_exe_0, reynolds_jarvis_exe_1, num_hru, output_vars)

reynolds_nlscatter_exe_0 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsNLscatter_exe_0_G1-1_timestep.nc"
reynolds_nlscatter_exe_1 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsNLscatter_exe_1_G1-1_timestep.nc"
total_errors_reynolds_nlscatter = verify_data(reynolds_nlscatter_exe_0, reynolds_nlscatter_exe_1, num_hru, output_vars)

reynolds_ueb2stream_exe_0 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsUEB2stream_exe_0_G1-1_timestep.nc"
reynolds_ueb2stream_exe_1 = "../test_cases/output/reynolds_canopySrad_windPrfile_stomResist/reynoldsUEB2stream_exe_1_G1-1_timestep.nc"
total_errors_reynolds_ueb2stream = verify_data(reynolds_ueb2stream_exe_0, reynolds_ueb2stream_exe_1, num_hru, output_vars)

reynolds_groundwatr_exe_0 = "../test_cases/output/reynolds_groundwatr/reynoldsDistributedQTopmodel_exe_0_G1-1_timestep.nc"
reynolds_groundwatr_exe_1 = "../test_cases/output/reynolds_groundwatr/reynoldsDistributedQTopmodel_exe_1_G1-1_timestep.nc"
# total_errors_reynolds_groundwatr = verify_data(reynolds_groundwatr_exe_0, reynolds_groundwatr_exe_1, num_hru, output_vars)

reynolds_lumped_exe_0 = "../test_cases/output/reynolds_groundwatr/reynoldsLumpedQTopmodel_exe_0_G1-1_timestep.nc"
reynolds_lumped_exe_1 = "../test_cases/output/reynolds_groundwatr/reynoldsLumpedQTopmodel_exe_1_G1-1_timestep.nc"
total_errors_reynolds_lumped = verify_data(reynolds_lumped_exe_0, reynolds_lumped_exe_1, num_hru, output_vars)

reynolds_no_xplict_exe_0 = "../test_cases/output/reynolds_groundwatr/reynoldsNoXplict_exe_0_G1-1_timestep.nc"
reynolds_no_xplict_exe_1 = "../test_cases/output/reynolds_groundwatr/reynoldsNoXplict_exe_1_G1-1_timestep.nc"
total_errors_reynolds_no_xplict = verify_data(reynolds_no_xplict_exe_0, reynolds_no_xplict_exe_1, num_hru, output_vars)

senator_alb_con_exe_0 = "../test_cases/output/senator_alb_method/albedoTestCon_exe_0_G1-1_timestep.nc"
senator_alb_con_exe_1 = "../test_cases/output/senator_alb_method/albedoTestCon_exe_1_G1-1_timestep.nc"
total_errors_senator_alb = verify_data(senator_alb_con_exe_0, senator_alb_con_exe_1, num_hru, output_vars)

senator_alb_var_exe_0 = "../test_cases/output/senator_alb_method/albedoTestVar_exe_0_G1-1_timestep.nc"
senator_alb_var_exe_1 = "../test_cases/output/senator_alb_method/albedoTestVar_exe_1_G1-1_timestep.nc"
total_errors_senator_alb_var = verify_data(senator_alb_var_exe_0, senator_alb_var_exe_1, num_hru, output_vars)

umpaqua_light_snow_exe_0 = "../test_cases/output/umpqua_snowIncept/umpquaLightSnow_exe_0_G1-1_timestep.nc"
umpaqua_light_snow_exe_1 = "../test_cases/output/umpqua_snowIncept/umpquaLightSnow_exe_1_G1-1_timestep.nc"
# total_errors_umpaqua_light_snow = verify_data(umpaqua_light_snow_exe_0, umpaqua_light_snow_exe_1, num_hru, output_vars)

umpaqua_sticky_snow_exe_0 = "../test_cases/output/umpqua_snowIncept/umpquaStickySnow_exe_0_G1-1_timestep.nc"
umpaqua_sticky_snow_exe_1 = "../test_cases/output/umpqua_snowIncept/umpquaStickySnow_exe_1_G1-1_timestep.nc"
# total_errors_umpaqua_sticky_snow = verify_data(umpaqua_sticky_snow_exe_0, umpaqua_sticky_snow_exe_1, num_hru, output_vars)

num_gru = 25
model_output_file = "../test_cases/settings/multiGruTestCases/northamerica2005/outputControl.txt"
output_vars = get_output_vars(model_output_file)
northamerica_exe_0 = "../test_cases/output/northamerica2005/NorthAmerica_exe_0_G001-025_timestep.nc"
northamerica_exe_1 = "../test_cases/output/northamerica2005/NorthAmerica_exe_1_G001-025_timestep.nc"
total_errors_northamerica = verify_data(northamerica_exe_0, northamerica_exe_1, num_gru, output_vars)


print()
print("Finished Checking Output")
print("\tTotal Errors Celia = ", total_errors_celia)
print("\tTotal Errors Colbeck1 = ", total_errors_colbeck_exp1)
print("\tTotal Errors Colbeck2 = ", total_errors_colbeck_exp2)
print("\tTotal Errors Colbeck3 = ", total_errors_colbeck_exp3)
print("\tTotal Errors Miller Clay = ", total_errors_miller_clay)
print("\tTotal Errors Miller Loam = ", total_errors_miller_loam)
print("\tTotal Errors Miller Sand = ", total_errors_miller_sand)
print("\tTotal Errors Mizoguchi = ", total_errors_mizoguchi)
print("\tTotal Errors Vanderborght1 = ", total_errors_vanderborght1)
print("\tTotal Errors Vanderborght2 = ", total_errors_vanderborght2)
print("\tTotal Errors Vanderborght3 = ", total_errors_vanderborght3)
print("\tTotal Errors Wigmosta Exp 1 = ", total_errors_wigmosta_exp1)
print("\tTotal Errors Wigmosta Exp 2 = ", total_errors_wigmosta_exp2)
print("\tTotal Errors Reynolds Canopy CLM2Stream = ", total_errors_reynolds_canopy_clm2stream)
print("\tTotal Errors Reynolds Default = ", total_errors_reynolds_default)
print("\tTotal Errors Reynolds Exponential = ", total_errors_reynolds_exponential)
print("\tTotal Errors Reynolds Jarvis = ", total_errors_reynolds_jarvis)
print("\tTotal Errors Reynolds NLscatter = ", total_errors_reynolds_nlscatter)
print("\tTotal Errors Reynolds UEB2Stream = ", total_errors_reynolds_ueb2stream)
# print("\tTotal Errors Reynolds Groundwatr = ", total_errors_reynolds_groundwatr)
print("\tTotal Errors Reynolds Lumped = ", total_errors_reynolds_lumped)
print("\tTotal Errors Reynolds No Xplict = ", total_errors_reynolds_no_xplict)
print("\tTotal Errors Senator Alb Con = ", total_errors_senator_alb)
print("\tTotal Errors Senator Alb Var = ", total_errors_senator_alb_var)
# print("\tTotal Errors Umpaqua Light Snow = ", total_errors_umpaqua_light_snow)
# print("\tTotal Errors Umpaqua Sticky Snow = ", total_errors_umpaqua_sticky_snow)
print("\tTotal Errors NorthAmerica = ", total_errors_northamerica)



