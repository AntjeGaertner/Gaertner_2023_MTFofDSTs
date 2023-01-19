#------------#
#   README   #
#------------#

This repository contains data and analysis script to recreate Figures and Tables presented in the manuscript "Mean time to fall of standing dead trees explained by temperature, tree size, plant functional type, humidity and mortality cause" submitted to Forests (January 2023).

#--- Authors (Affiliations):
Antje Gärtner (1); Anna Maria Jönsson (1); Daniel B. Metcalfe (1,2); Thomas A. M. Pugh (1,3,4); Torbern Tagesson (1,5); Anders Ahlström (1)

1 	Department of Physical Geography and Ecosystem Science, Lund University, Lund, Sweden
2 	Department of Ecology and Environmental Sciences, Umeå University, Umeå, Sweden
3 	School of Geography, Earth & Environmental Sciences, University of Birmingham, Birmingham, UK
4 	Birmingham Institute of Forest Research, University of Birmingham, Birmingham, UK;
5 	Department of Geosciences and Natural Resource Management, University of Copenhagen, Copenhagen, Denmark
*	Correspondence: antje.gartner@nateko.lu.se




#--- Data availability
All data used and presented in Gärtner et al. (submitted to Forests) is stored in the /data folder. The mean time to fall (MTF) database (MTF_database.xlsx) contains MTFs at site and species level, as well as site information, correspondingg climate and location information and full reference list of publications included in the analysis.




#--- Initialisation of python environments
Scripts use 3 different python environments (Geopandas2, stats and P3_21) to manage different package requirements and incompatibilities between them. We provide the package environments as .yml files that can be installed using anaconda or mini conda distributions (https://www.anaconda.com/products/distribution). The environments can be installed using the provided bash script "install_python_environments.sh" which also installs Jupiter kernels with the correct names to run analysis and figure scripts. Run as follows from the command line:

	bash install_python_environments.sh



#--- Running analysis presented in Gärtner et al. (submitted to Forests)
The analysis, figures and tables presented in Gärtner et al. (submitted to Forests) can be run comprehensively using the bash script "run_DB_climate_analysis.sh". Run from the command line as follows:

	bash run_DB_climate_analysis.sh

This script will print paths and statements to the command line detailing the analysis flow.