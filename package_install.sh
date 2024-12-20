current_path=$(eval pwd)
echo $current_path


conda_path=$(eval conda env list | grep -o '/.*'| head -1)

echo $conda_path

conda env create -f ./Density_profile_calculation.yml

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' '1s/^/#/' density_run_script.sh
else
    sed -i '1s/^/#/' density_run_script.sh
fi
# echo "source $conda_path" | cat - filename.txt > temp && mv temp filename.txt
echo "source $conda_path/bin/activate" | cat - density_run_script.sh > temp && mv temp density_run_script.sh
# cd $current_path

