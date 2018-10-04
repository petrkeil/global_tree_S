
# This is a sequence of scripts that need to be executed to re-run the actual data
# extraction and analysis. Steps 1-2 are missing since they require some
# manual adjustments.

source("3_Extract_environment_for_PLOTS_and_COUNTRIES.r")
source("3_Extract_environment_for_prediction_grids.r")
source("4_Data_cleaning_and_putting_PLOTS_and_COUNTRIES_together.r")
