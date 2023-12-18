
setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
library(TMB)
compile("tinyVAST.cpp", framework="TMBad")

devtools::document( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)')
devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)',
                         force=TRUE, dep=FALSE )
?tinyVAST::fit

setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\vignettes)' )
vignettes = list.files( pattern="Rmd")
for( vign in vignettes ){
  message( "Running " , vign)
  devtools::build_rmd( file.path(getwd(), vign) )
  rmarkdown::render( file.path(getwd(), vign), rmarkdown::pdf_document())
}
