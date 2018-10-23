################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Dec 7 2017
################################################################################

# Description: This is just a simple little script that sorts the list of references
# and removes duplicates.

################################################################################


PLOTS <- read.csv("../Data/PLOTS/PLOTS_all_data_raw_tab_delimited.txt",
                  sep="\t", header=TRUE)

a <- sort(unique(PLOTS$Ref_long))


write.table(a, file="../Bibliography/Sources_for_forest_plots.txt", row.names = FALSE, 
            col.names = FALSE, quote=FALSE)
