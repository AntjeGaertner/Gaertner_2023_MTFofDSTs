#!/bin/sh


#--- This bash scripts installs python environments necessary to run analysis and figure scripts published in Gärtner et al. (submitted to Forests)


#--- Install Geopandas2 environment
conda env create --file python_env_Geopandas2.yml # Create python environment
eval "$(conda shell.bash hook)"
conda activate Geopandas2 # Activate environment
python -m ipykernel install --user --name Geopandas2 --display-name "Geopandas2" # Install Jupyter kernel


#--- Install stats environment
conda env create --file python_env_stats.yml # Create python environment
eval "$(conda shell.bash hook)"
conda activate stats # Activate environment
python -m ipykernel install --user --name stats --display-name "stats" # Install Jupyter kernel


#--- Install P3_21 environment
conda env create --file python_env_P3_21.yml # Create python environment
eval "$(conda shell.bash hook)"
conda activate P3_21 # Activate environment
python -m ipykernel install --user --name P3_21 --display-name "P3_21" # Install Jupyter kernel
