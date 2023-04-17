#!/bin/sh

alpha=0.01

# Create folders to store intermediary results
mkdir spe_reg_results/
mkdir site_reg_results/

#--- Run Mortality selection for both site & species
echo ""
echo "#---------------------------------------#"
echo "#   Model selection by mortality group  #"
echo "#---------------------------------------#"
echo ""
echo "Running model selection for the largest data subset with values on MAT, DBH, Management and PFT."
echo "Model selection for mortality groups including all causes of mortality (M_All), dead standing trees (DSTs) that died by fire (M_Fire), DSTs that died by other causes than fire (including insects; M_NoFire) and DSTs that died by other causes than fire or insects (M_Other)."
echo ""
echo ""

exp_name="$(date +"%y%m%d")_Management"
exp_name_management="${exp_name}"

echo ""
# Remove Ricthie & all data
echo "#--- Log files include links to .html files detailing results by climate data set at Site and Species level."
echo "M_All:"
nohup bash run_DB_climate_analysis.sh "${exp_name}" true ${alpha} > "${exp_name}_true_${alpha}_run_log.txt" &
BACK_PID=$!
echo "${exp_name}_true_${alpha}_run_log.txt"
echo ""
wait ${BACK_PID}

# Fire
echo "M_Fire:"
nohup bash run_DB_climate_analysis.sh "${exp_name}_FireOnly" true ${alpha} > "${exp_name}_FireOnly_true_${alpha}_run_log.txt" &
BACK_PID=$!
echo "${exp_name}_FireOnly_true_${alpha}_run_log.txt"
echo ""
wait ${BACK_PID}

# No Fire
echo "M_NoFire:"
nohup bash run_DB_climate_analysis.sh "${exp_name}_NoFire" true ${alpha} > "${exp_name}_NoFire_true_${alpha}_run_log.txt" &
BACK_PID=$!
echo "${exp_name}_NoFire_true_${alpha}_run_log.txt"
echo ""
wait ${BACK_PID}

# No Fire & Insects
echo "M_Other:"
nohup bash run_DB_climate_analysis.sh "${exp_name}_NoInsectsFire" true ${alpha} > "${exp_name}_NoInsectsFire_true_${alpha}_run_log.txt" &
BACK_PID=$!
echo "${exp_name}_NoInsectsFire_true_${alpha}_run_log.txt"
echo ""
wait ${BACK_PID}


# Species
echo ' '
echo 'Analysing Species output:'
NB_ARGS="${alpha}","${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates_model_analysis_final_table.ipynb
mv SPECIES_multiple_regressions_climates_model_analysis_final_table.html "spe_reg_results/${exp_name}_snag-MTF_regression_multiple_species_model_analysis_final_table_${alpha}.html"
echo "spe_reg_results/${exp_name}_snag-MTF_regression_multiple_species_model_analysis_final_table_${alpha}.html"

# Site
echo ' '
echo 'Analysing Site output:'
NB_ARGS="${alpha}","${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates_model_analysis_final_table.ipynb
mv SITE_multiple_regressions_climates_model_analysis_final_table.html "site_reg_results/${exp_name}_multiple_regressions_site_climates_model_analysis_final_table_${alpha}.html"
echo "site_reg_results/${exp_name}_multiple_regressions_site_climates_model_analysis_final_table_${alpha}.html"

# Combine to best models table for Species and Site
echo ' '
echo 'Compiling Table 4 and B2:'
NB_ARGS="${alpha}","${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN Table4_Appendix_TableB2_multiple_regressions_final_table_SITE_SPECIES.ipynb
mv Table4_Appendix_TableB2_multiple_regressions_final_table_SITE_SPECIES.html "${exp_name}_Table4_B2_multiple_regressions_site_climates_model_analysis_${alpha}.html"
echo "${exp_name}_Table4_B2_multiple_regressions_site_climates_model_analysis_${alpha}.html"

echo ' '
echo 'Compiling Table 5:'
standardise_covars=""
NB_ARGS="${exp_name}","${alpha}","${standardise_covars}" jupyter nbconvert --execute --to html --template full --log-level WARN Table5_Figure6_Management_subset_covariate_values.ipynb
mv Table5_Figure6_Management_subset_covariate_values.html "${exp_name}_Table5_regression_coeff_${alpha}_False.html"
echo "${exp_name}_Table5_regression_coeff_${alpha}_False.html"

echo ' '
echo 'Plotting Figure 6:'
standardise_covars="True"
NB_ARGS="${exp_name}","${alpha}","${standardise_covars}" jupyter nbconvert --execute --to html --template full --log-level WARN Table5_Figure6_Management_subset_covariate_values.ipynb
mv Table5_Figure6_Management_subset_covariate_values.html "${exp_name}_Figure6_regression_coeff_${alpha}_True.html"
echo "${exp_name}_Figure6_regression_coeff_${alpha}_True.html"

echo ""
echo "Climate comparison, Table B1:"
NB_ARGS="${alpha}","${exp_name_full}","${exp_name_management}",true jupyter nbconvert --execute --to html --template full --log-level WARN TableB1_climate_comparison_site_species.ipynb
mv TableB1_climate_comparison_site_species.html "${exp_name}_TableB1_climate_comparison_site_species_true_${alpha}.html"
echo "${exp_name}_TableB1_climate_comparison_site_species_true_${alpha}.html"
echo ""
echo ""


#-----------------------#
#   Additional tests    #
#-----------------------#

echo ""
echo "#------------------------------------#"
echo "#   Running the temperature subset   #"
echo "#------------------------------------#"
echo ""
echo "Running the mortality groups again for the narrower temperature range of M_Fire to test if MAT remains a significant predictor variable."

exp_name="$(date +"%y%m%d")_TempSS_FMort"

# M_No Fire
echo "Temp SS M_NoFire:"
nohup bash run_DB_climate_analysis.sh "${exp_name}_NoFire" true ${alpha} > "${exp_name}_NoFire_true_${alpha}_run_log.txt" &
BACK_PID=$!
echo "${exp_name}_NoFire_true_${alpha}_run_log.txt"
echo ""
wait ${BACK_PID}

# M_No Fire & Insects
echo "Temp SS  M_Other:"
nohup bash run_DB_climate_analysis.sh "${exp_name}_NoInsectsFire" true ${alpha} > "${exp_name}_NoInsectsFire_true_${alpha}_run_log.txt" &
BACK_PID=$!
echo "${exp_name}_NoInsectsFire_true_${alpha}_run_log.txt"
echo ""
wait ${BACK_PID}


#--- Analyse the regression results
# Species
echo ' '
echo 'Temp SS, Analysing Species output:'
NB_ARGS="${alpha}","${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN SPECIES_multiple_regressions_climates_model_analysis_final_table.ipynb
mv SPECIES_multiple_regressions_climates_model_analysis_final_table.html "spe_reg_results/${exp_name}_"SPECIES_multiple_regressions_climates_model_analysis_final_table_"${alpha}.html"
echo "spe_reg_results/${exp_name}_"SPECIES_multiple_regressions_climates_model_analysis_final_table_"${alpha}.html"

# Site
echo ' '
echo 'Temp SS, Analysing Site output:'
NB_ARGS="${alpha}","${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN SITE_multiple_regressions_climates_model_analysis_final_table.ipynb
mv SITE_multiple_regressions_climates_model_analysis_final_table.html "site_reg_results/${exp_name}_"SITE_multiple_regressions_climates_model_analysis_final_table_"${alpha}".html
echo "site_reg_results/${exp_name}_"SITE_multiple_regressions_climates_model_analysis_final_table_"${alpha}".html

# Combine to best models table for Species and Site
echo ' '
echo 'Temp SS, Compiling Table B3:'
NB_ARGS="${alpha}","${exp_name}" jupyter nbconvert --execute --to html --template full --log-level WARN Table4_Appendix_TableB2_multiple_regressions_final_table_SITE_SPECIES.ipynb
mv Table4_Appendix_TableB2_multiple_regressions_final_table_SITE_SPECIES.html "${exp_name}_"Table4_Appendix_TableB2_multiple_regressions_final_table_SITE_SPECIES_"${alpha}".html
echo "${exp_name}_"Table4_Appendix_TableB2_multiple_regressions_final_table_SITE_SPECIES_"${alpha}".html


echo ""
echo "#-------------------------------------#"
echo "#   Running the Wood quality subset   #"
echo "#-------------------------------------#"
echo ""

exp_name="$(date +"%y%m%d")_WoodQ_NoInsectsFire"
#--- Wood Quality
nohup bash run_DB_climate_analysis.sh "${exp_name}" true ${alpha} > "${exp_name}_true_${alpha}_run_log.txt" &
echo "${exp_name}_true_${alpha}_run_log.txt"
echo "Table B4:"
echo "spe_reg_results/${exp_name}_SPECIES_multiple_regressions_climates_model_analysis_${alpha}_true.html"
echo ""
BACK_PID=$!
wait ${BACK_PID1}


echo ""
echo "#---------------------------------#"
echo "#   Running the Moisture subset   #"
echo "#---------------------------------#"
echo ""
exp_name="$(date +"%y%m%d")_Moisture"

#--- Moisture
nohup bash run_DB_climate_analysis.sh "${exp_name}" true ${alpha} > "${exp_name}_true_${alpha}_run_log.txt" &
echo "${exp_name}_true_${alpha}_run_log.txt"
echo "spe_reg_results/${exp_name}_SPECIES_multiple_regressions_climates_model_analysis_${alpha}_true.html"
echo ""
BACK_PID=$!
wait ${BACK_PID}

echo ""
echo ""
echo "#--- Running remaining Figures and Tables"
echo ""
echo ""
echo "Table 3  - Database overview"
jupyter nbconvert --execute --to html --template full --log-level WARN Table3_DB_overview.ipynb
echo ""
echo "Figure 3 - Study locations"
jupyter nbconvert --execute --to html --template full --log-level WARN Figure3_site_locations_Whittaker_biomes.ipynb
echo ""
echo "Figure 4 - Data distributions by category"
jupyter nbconvert --execute --to html --template full --log-level WARN Figure4_violin_plots_categories.ipynb
echo ""
echo "Figure 5 - Scatter DBH, mortality cause and MAT"
jupyter nbconvert --execute --to html --template full --log-level WARN Figure5_Scatter_DBH_Mort_cause_MAT.ipynb
echo ""
echo "Figure 7 - Species model visualisation"
jupyter nbconvert --execute --to html --template full --log-level WARN Figure7_species_MFire_MOther.ipynb
echo ""
echo "Appendix Figure A1            - Survey duration by reference"
jupyter nbconvert --execute --to html --template full --log-level WARN Appendix_FigureA1_survey_duration_references_plot.ipynb
echo ""
echo "Appendix Figure A2 & A3       - Simulated errors"
jupyter nbconvert --execute --to html --template full --log-level WARN Appendix_FigureA2_A3_simulated_errors.ipynb
echo ""
echo "Appendix Figure A4 & A5       - MTFcount to MTFsize regression"
jupyter nbconvert --execute --to html --template full --log-level WARN Appendix_FigureA4_A5_count2cmass_regression.ipynb
echo ""
echo "Appendix Figure B1            - Mortality cause differences"
jupyter nbconvert --execute --to html --template full --log-level WARN Appendix_FigureB1_mortality_cause_differences.ipynb
echo ""
echo "Appendix Table A2 & Figure B3 - Visualization of covariate ranges by mortality cause"
jupyter nbconvert --execute --to html --template full --log-level WARN Appendix_FigureB3_regression_covariates_visualisation.ipynb

echo ""
