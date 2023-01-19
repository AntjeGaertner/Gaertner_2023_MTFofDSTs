#!/bin/sh
#-----------------------------#
#   Database control script   #
#-----------------------------#
exp_name=$1 #"$(date +"%y%m%d")" #$1 #$4     #"$(date +"%y%m%d")" Experimental identifier (yymmdd)+ x  

# Subsets for which we only run Species level regressions
SUB_1='Moisture'
SUB_2='WoodQ'

#--- SETTINGS

#--- Database
run_db=false                     # Run database? (true | false)

#--- Climate
run_climate=false                # Run climate? (true | false)


#--- Analysis
run_analysis_spe=true            # Run SPECIES alysis? (true | false)
run_analysis_site=true           # Run SITE analysis? (true | false)






#--- ARGUMENTS passed to Python
#- Data selection
remove_data=$2                   # Remove Ritchie2013? (true | "") - empty quotes evaluate to False using bool() in python


#- Script input from command line
alpha=$3 #0.05                   # Significance threshold (float)





echo ' '
echo "Today: ${exp_name}: " # Experimental log first line
echo "Start time: $(date +"%H:%M")"
echo ' '
# Selecting subset based on experiment string
if [[ "$exp_name" == *"$SUB_1"* ]] || [["$exp_name" == *"$SUB_2"* ]]; then
  run_analysis_spe=true           # Run SPECIES analysis? (true | false)
  run_analysis_site=false
  echo "Running only species level regressions."
fi


#------------------#
#   Run database   #
#------------------#


if [ "$run_db" = true ]; then
    echo ' '
    echo ' '
    echo '#--- Running the Database'
    bash run_nb.sh
    
    echo ' '
    echo 'Extracting locations'
    jupyter nbconvert --to notebook --inplace --execute 00_211025_snag-MTF_Site_map.ipynb
    
    echo ' '
    echo 'Calculating MTF for digitised figure data'
    jupyter nbconvert --execute --to html --template full 00_211025_snag-MTF_Master.ipynb
    mv 00_211025_snag-MTF_Master.html "${exp_name}_"snag-MTF_Master.html
    
    echo ' '
    echo 'Calculating average CRU climate based on 6-hourly data'
    jupyter nbconvert --to notebook --inplace --execute climate/220426_climate_extraction_CRUNCEPv7_obs_period.ipynb  # Extract climate for CRUNCEP obs
    jupyter nbconvert --to notebook --inplace --execute climate/220510_climate_extraction_CRUNCEPv7_climatology.ipynb # Extract climate CRUNCEPv7 climatology
    
    echo ' '
    echo 'Calcluating CDI'
    jupyter nbconvert --to notebook --inplace --execute climate/220826_climate_extraction_CDI.ipynb
fi


#------------------#
#   Climate only   #
#------------------#

if [ "$run_climate" = true ]; then
    
    echo ' '
    echo ' '
    echo '#--- Extracting climate'
    
    echo ' '
    echo 'converting count to cmass'
    jupyter nbconvert --to notebook --inplace --execute 00_220829_snag-MTF_count2cmass_regression.ipynb
    
    #--- CRUNCEP 1981-2010
    echo ' '
    echo 'CRUNCEP 1981-2010'
    jupyter nbconvert --to notebook --inplace --execute 00_220510_snag-MTF_climate_regression_CRUNCEPv7_climatology.ipynb
    NB_ARGS='True','CRUclim' jupyter nbconvert --execute --to html --template full 220201_Whittaker_biomes.ipynb
    mv 220201_Whittaker_biomes.html "spe_reg_results/${exp_name}_"CRUclim_Whittaker_biomes.html
    
    
    echo ' '
    echo 'CRUNCEPv7 obs time period'
    jupyter nbconvert --to notebook --inplace --execute 00_220413_snag-MTF_climate_regression_CRUNCEPv7_obs_period.ipynb
    NB_ARGS='True','CRUNCEPv7' jupyter nbconvert --execute --to html --template full 220201_Whittaker_biomes.ipynb
    mv 220201_Whittaker_biomes.html "spe_reg_results/${exp_name}_"CRUNCEPv7_Whittaker_biomes.html
    
    echo ' '
    echo 'CHELSA 30s'
    jupyter nbconvert --to notebook --inplace --execute 00_220412_snag-MTF_climate_regression_CHELSA30s.ipynb
    NB_ARGS='True','CHELSA30s' jupyter nbconvert --execute --to html --template full 220201_Whittaker_biomes.ipynb
    mv 220201_Whittaker_biomes.html "spe_reg_results/${exp_name}_"CHELSA30s_Whittaker_biomes.html
    
    echo ' '
    echo 'WorldClim 30s'
    jupyter nbconvert --to notebook --inplace --execute 00_220413_snag-MTF_climate_regression_WorldClim2_30s.ipynb
    NB_ARGS='True','WorldClim30s' jupyter nbconvert --execute --to html --template full 220201_Whittaker_biomes.ipynb
    mv 220201_Whittaker_biomes.html "spe_reg_results/${exp_name}_"WorldClim30s_Whittaker_biomes.html
    
    echo ' '
    echo 'WorldClim 10m'
    jupyter nbconvert --to notebook --inplace --execute 00_220413_snag-MTF_climate_regression_WorldClim2_10m.ipynb
    NB_ARGS='True','WorldClim10m' jupyter nbconvert --execute --to html --template full 220201_Whittaker_biomes.ipynb
    mv 220201_Whittaker_biomes.html "spe_reg_results/${exp_name}_"WorldClim10m_Whittaker_biomes.html
fi



#-------------------------#
#   Run the regressions   #
#-------------------------#

#--- Species
if [ "$run_analysis_spe" = true ]; then
    echo ' '
    echo ' '
    echo '#--- Multiple linear regressions'
    echo "SPECIES"
    echo "Runnning with a significance level of  ${alpha}"
    echo " "
    echo " "
    
    # Make current experiment list
    echo ' '
    echo 'Make experiment list'
    NB_ARGS="${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_generate_all_possible_covariate_models.ipynb

    #--- CRUNCEP 1981-2010
    echo ' '
    echo 'CRUNCEP 1981-2010'
    NB_ARGS='CRUclim',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates.ipynb
    mv SPECIES_multiple_regressions_climates.html  "spe_reg_results/${exp_name}_CRUclim_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"
    echo "spe_reg_results/${exp_name}_CRUclim_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"


    #--- CRUNCEPv7 obs time period
    echo ' '
    echo 'CRUNCEPv7 obs time period'
    NB_ARGS='CRUNCEPv7',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates.ipynb
    mv SPECIES_multiple_regressions_climates.html "spe_reg_results/${exp_name}_CRUNCEPv7_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"
    echo "spe_reg_results/${exp_name}_CRUNCEPv7_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"


    #--- Chelsa 30s
    echo ' '
    echo 'CHELSA 30s'
    NB_ARGS='CHELSA30s',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates.ipynb
    mv SPECIES_multiple_regressions_climates.html "spe_reg_results/${exp_name}_CHELSA30s_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"
    echo "spe_reg_results/${exp_name}_CHELSA30s_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"


    #--- WorldClim 30s
    echo ' '
    echo 'WorldClim 30s'
    NB_ARGS='WorldClim30s',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates.ipynb
    mv SPECIES_multiple_regressions_climates.html "spe_reg_results/${exp_name}_WorldClim30s_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"
    echo "spe_reg_results/${exp_name}_WorldClim30s_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"


    #--- WorldClim 10m
    echo ' '
    echo 'WorldClim 10m'
    NB_ARGS='WorldClim10m',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates.ipynb
    mv SPECIES_multiple_regressions_climates.html "spe_reg_results/${exp_name}_WorldClim10m_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"
    echo "spe_reg_results/${exp_name}_WorldClim10m_${alpha}_${remove_data}_SPECIES_multiple_regressions_climates.html"



    #--- Analyse the results
    echo ' '
    echo 'Analysing output'
    NB_ARGS="${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates_model_analysis.ipynb
    mv SPECIES_multiple_regressions_climates_model_analysis.html "spe_reg_results/${exp_name}_SPECIES_multiple_regressions_climates_model_analysis_${alpha}_${remove_data}.html"
    echo "spe_reg_results/${exp_name}_SPECIES_multiple_regressions_climates_model_analysis_${alpha}_${remove_data}.html"
    
fi


#--- Site
if [ "$run_analysis_site" = true ]; then
    echo ' '
    echo ' '
    echo '#--- Multiple linear regressions'
    echo "SITE"
    echo "Runnning with a significance level of  ${alpha}"
    echo " "
    echo " "
    
    # Make current experiment list
    echo ' '
    echo 'Make experiment list'
    NB_ARGS="${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_generate_all_possible_covariate_models.ipynb

    #--- CRUNCEP 1981-2010
    echo ' '
    echo 'CRUNCEP 1981-2010'
    NB_ARGS='CRUclim',"${alpha}","${exp_name}","${remove_data}"  jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates.ipynb
    mv SITE_multiple_regressions_climates.html "site_reg_results/${exp_name}_"CRUclim_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html
    echo "site_reg_results/${exp_name}_"CRUclim_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html


    #--- CRUNCEPv7 obs time period
    echo ' '
    echo 'CRUNCEPv7 obs time period'
    NB_ARGS='CRUNCEPv7',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates.ipynb
    mv SITE_multiple_regressions_climates.html "site_reg_results/${exp_name}_"CRUNCEPv7_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html
    echo "site_reg_results/${exp_name}_"CRUNCEPv7_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html


    #--- Chelsa 30s
    echo ' '
    echo 'CHELSA 30s'
    NB_ARGS='CHELSA30s',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates.ipynb
    mv SITE_multiple_regressions_climates.html "site_reg_results/${exp_name}_"CHELSA30s_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html
    echo "site_reg_results/${exp_name}_"CHELSA30s_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html


    #--- WorldClim 30s
    echo ' '
    echo 'WorldClim 30s'
    NB_ARGS='WorldClim30s',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates.ipynb
    mv SITE_multiple_regressions_climates.html "site_reg_results/${exp_name}_"WorldClim30s_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html
    echo "site_reg_results/${exp_name}_"WorldClim30s_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html


    #--- WorldClim 10m
    echo ' '
    echo 'WorldClim 10m'
    NB_ARGS='WorldClim10m',"${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates.ipynb
    mv SITE_multiple_regressions_climates.html "site_reg_results/${exp_name}_"WorldClim10m_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html
    echo "site_reg_results/${exp_name}_"WorldClim10m_"${alpha}_${remove_data}_"SITE_multiple_regressions_climates.html



    #--- Analyse the results
    echo ' '
    echo 'Analysing Site output'
    NB_ARGS="${alpha}","${exp_name}","${remove_data}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates_model_analysis.ipynb
    mv SITE_multiple_regressions_climates_model_analysis.html "site_reg_results/${exp_name}_"SITE_multiple_regressions_climates_model_analysis_"${alpha}_${remove_data}".html
    echo "site_reg_results/${exp_name}_"SITE_multiple_regressions_climates_model_analysis_"${alpha}_${remove_data}".html
fi


echo ''
echo "End time: $(date +"%H:%M")"
