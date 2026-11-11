# Recreate every table in the paper, with its footnote, in output/tables/.
# Tables 1, 2, and S1 are recalculated from the planning formulas. The other
# tables use the paper's saved results in data/saved_results/.
#
# Run it with: Rscript make_tables.R
source("R/load.R")

export_paper_tables()
