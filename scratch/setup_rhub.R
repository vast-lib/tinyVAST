
# Directions from: https://r-hub.github.io/rhub/articles/rhubv2.html
setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )

library(rhub)

# Options
rhub_platforms()

# Setup
rhub_setup(overwrite = FALSE)

# Check
rhub::rhub_doctor()

# Run
platforms = c( "m1-san", "clang-ubsan", "valgrind" )
# 2 -- UBSAN
# 9 -- clang-ubsan
# 30 -- valgrind
rhub::rhub_check( platforms = platforms )

# Check valgrind on dev
if( FALSE ){
  rhub::rhub_check( platforms = "valgrind", branch = "dev" )
  # Or main
  rhub::rhub_check( platforms = "valgrind", branch = "main" )
}

# HOW TO CHECK `valgrind` FOR ISSUES
# Click Actions tab -> click R-hub tab on left-hand side ->
# click action -> job -> "Run r-hub/actions/run-check@v1"
# SEARCH: "Artifact download URL" to find URL to download errors
# e.g., here: https://github.com/vast-lib/tinyVAST/actions/runs/13868485634/job/38811880343
# -> Check tinyVAST-Ex.Rout and search "tinyVAST.cpp" for errors
