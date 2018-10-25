
# This is to run some of the key scripts sequentially and automatically.


# scripts that extract the environmental predictors (beware: these require
# large environmental rasters that don't come as a part of this repository)
source("3_Extract_environment_for_PLOTS_and_COUNTRIES.r")
source("3_Extract_environment_for_prediction_grids.r")

# cleaning of the data
source("4.1_Data_cleaning_and_putting_PLOTS_and_COUNTRIES_together.r")
# source("5.0_GAM_testing_post-review_addition_of_insularity_predictors_and_elongation.r")

# the actual analyses
source("5.1_GAM_main_analysis_1_model_fitting.r")
source("5.1_GAM_main_analysis_2_model_inference_figures.r")

source("5.2_GAM_crossvalidation.r")
source("5.3_GAM_external_validation.r")
source("5.4_GAM_spatial_correlograms.r")
source("6.0_GAM_sensitivity_analysis_on_subset_of_data.r")
# source("7.0_GAM_deviance_partitioning_and_PCA.r")
source("8.0_GAM_make_predictions_to_regular_grids_at_two_grains.r")
source("8.0_GAM_maps_of_prediction_uncertainty.r")
source("8.0_GAM_maps_of_raw_data_and_predictions_in_data_units.r")
# source("9.0_Extract_plot_references.r")

# ------------------------------------------------------------------------------

# session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# how many lines of R code do we have altoghether?
system("find . -name '*.r' | xargs wc -l")
