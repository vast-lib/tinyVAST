
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
}
