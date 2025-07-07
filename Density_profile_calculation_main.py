import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

import Density_profile_function as func

np.set_printoptions(threshold=np.inf)


# print the current path
func.current_path()

Work_path = input("Please input the working directory: \n")

System_info = (input("Please input the co-solvent present or not (y/n): \n")).lower()

Cosolvent = input("Please input co-solvent formula matching the file: \n")

# here require the input for the atom number of cosolvent
Atom_cosolvent = float(input("Please input the atom number of cosolvent: \n"))

# here require input the time range
Time_range = float(input("Please input the time range analysis begin from: \n"))


Save_path = input("Please input the saving directory: \n")

# Lipid_perleaflet = 402


Lipid_perleaflet = float(input("Please input the lipid amount per leaflet"))

File_name, File_dict = func.file_read(Work_path)


File_name.sort()

File_dict = func.data_interp(File_name, File_dict, Cosolvent, Atom_cosolvent)

Output_result = func.len_integra_main(
    File_dict, Time_range, System_info, Cosolvent, Lipid_perleaflet
)

file_folder_name = Work_path[(Work_path.rfind("/") + 1) :]
file_prefix_ndx = File_name[0].index("_r")
Output_result.to_csv(
    f"{Save_path}/{File_name[0][:file_prefix_ndx]}_ScriptAnalysis_{file_folder_name}.csv"
)
print(Output_result)
