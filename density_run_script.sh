source /opt/homebrew/Caskroom/miniforge/base/bin/activate
# activate env by source activate
#
# source /opt/homebrew/Caskroom/miniforge/base/bin/activate
conda activate denscal_env
echo please input density file path
read dens_path
echo please input the co-solvent present or not "(y/n)"
read system_info
echo please input co-solvent formula
read formula
echo please input co-solvent atom number
read atom
echo please input time range at the beginning point for analysis
read time
echo please input file saving path
read save_path
echo please input the lipid amount per leaflet
read lipid_perleaflet
num=($(eval ls "$dens_path" |wc -l))

# wc command only work for macos machine
# if using in linux, should be | wc --lines
file_pattern=($(eval ls "$dens_path"))
for((i=0;i<$num;i++));do
echo -e "${dens_path}/${file_pattern["$i"]}
${system_info}
${formula}
${atom}
${time}
${save_path}
${lipid_perleaflet}" | python Density_profile_calculation_main.py
#python Density_profile_calculation_main.py;
done

conda deactivate
