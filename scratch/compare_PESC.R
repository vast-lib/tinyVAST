
# Import and save data
if( FALSE ){
  library(VAST)
  library(usethis)
  data( PESC_example_red_grouper )

  red_grouper_diet <- example

  setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
  use_data( red_grouper_diet )
}

devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )

