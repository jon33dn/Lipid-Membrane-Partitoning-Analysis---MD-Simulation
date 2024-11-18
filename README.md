# Comprehensive Analysis of Lipid Structure and Small Molecule Partitioning from MD Simulation Trajectories

Please cite the following work. Thank you.

> Tan, L., Scott, H. L., Smith, M. D., Pingali, S. V., O’Neill, H. M., Morrell-Falvey, J. L., ... & Nickels, J. D. (2023).
> Amphiphilic Co-solvents modulate the structure of membrane domains. ACS Sustainable Chemistry & Engineering, 11(4), 1598-1609

#### Details of analysis and descriptions can be found in this paper [Link](https://pubs.acs.org/doi/10.1021/acssuschemeng.2c06876).

The goal of this script is used to analyze the membrane simulation w/o co-sovlent.

## Set up

This script use miniconda or anaecoda as package manager.
Please check out the installation page for [miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)
or [anaconda](https://docs.anaconda.com/anaconda/install/)

Run this command to install the python packages.

```shell
conda config --append channels conda-forge
bash package_install.sh
```

After finishing install all the required packages, you are good to run the analysis script.

```shell
bash density_run_script.sh
```

## Required Input

This script is required to input the following information.

- Density profile path

- system information (if co-solvent present) (y/n)

- co-solvent formula (if the system without co-sovlent, just type H2O)

  - To be noticed, this formula will be used in the script to look for the co-solvent density file

  <details>
  <summary>Example</summary>
      Butanol: C4H9OH <br>
      Ethanol: C2H5OH <br>
  </details>

- co-solvent atom numbers (if the system without co-solvent, just type 0)

- time range for analysis (used in the xy.xvg)

- Output csv path

- Lipids per leaflet

## File formatting

The input path only need to reach the EtOH folder, not the concentraion folder.
All the .xvg files are the output of **gmx density** (membrane centered by C21C31)
and **gmx energy**

Input density profile format with co-solvent EtOH example.

After running density_run_script.sh, it will go through every concentration in the EtOH folder
and output a .csv file for each concentration in the Output path.

```bash
EtOH
├── 6p0EtOH
│   ├── ${system}_6p0EtOH_r1_C21C31_along_z.xvg
│   ├── ${system}_6p0EtOH_r1_C2H5OH_along_z.xvg
│   ├── ${system}_6p0EtOH_r1_P_along_z.xvg
│   ├── ${system}_6p0EtOH_r1_water_along_z.xvg
│   ├── ${system}_6p0EtOH_r1_xy.xvg
│   ├── ${system}_6p0EtOH_r2_C21C31_along_z.xvg
│   ├── ${system}_6p0EtOH_r2_C2H5OH_along_z.xvg
│   ├── ${system}_6p0EtOH_r2_P_along_z.xvg
│   ├── ${system}_6p0EtOH_r2_water_along_z.xvg
│   ├── ${system}_6p0EtOH_r0_C21C31_along_z.xvg
│   ├── ${system}_6p0EtOH_r0_C2H5OH_along_z.xvg
│   ├── ${system}_6p0EtOH_r0_P_along_z.xvg
│   ├── ${system}_6p0EtOH_r0_water_along_z.xvg
│   ├── ${system}_6p0EtOH_r0_xy.xvg
│   └── ${system}_6p0EtOH_r2_xy.xvg
```

Input density profile format without co-solvent.

The folder need to be Wateronly

```bash
Wateronly
├── ${system}_Wateronly_r0_C21C31_along_z.xvg
├── ${system}_Wateronly_r0_P_along_z.xvg
├── ${system}_Wateronly_r0_water_along_z.xvg
├── ${system}_Wateronly_r0_xy.xvg
├── ${system}_Wateronly_r1_C21C31_along_z.xvg
├── ${system}_Wateronly_r1_P_along_z.xvg
├── ${system}_Wateronly_r1_water_along_z.xvg
├── ${system}_Wateronly_r1_xy.xvg
├── ${system}_Wateronly_r2_C21C31_along_z.xvg
├── ${system}_Wateronly_r2_P_along_z.xvg
├── ${system}_Wateronly_r2_water_along_z.xvg
└── ${system}_Wateronly_r2_xy.xvg
```
