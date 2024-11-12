current_path=$(eval pwd)
echo $current_path


conda_path=$(eval conda env list | grep -o '/.*'| head -1)

echo $conda_path

# conda create -f ./Density_profile_calculation.yml


echo "Your content here" | cat - filename.txt > temp && mv temp filename.txt
echo 
# cd $current_path

