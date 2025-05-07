
setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)' )
#setwd( R'(C:\Users\james\OneDrive\Desktop\Git\tinyVAST\src)' )

file.remove( "tinyVAST.dll" )
file.remove( "tinyVAST.o" )
file.remove( "compile_output.txt" )

sink( file = "compile_output.txt" )
sink( type = "message" )
  TMB::compile("tinyVAST.cpp", framework = "TMBad",  flags = "-pedantic -O2 -Wall" ) # eigen.disable.warning = FALSE )
sink(type = "message")
sink()
#compile_cmd <- TMB::compile("tinyVAST.cpp", framework = "TMBad",  eigen.disable.warning = FALSE, flag = TRUE)

# Instead of using TMB::compile directly, we'll run system2
# First, extract compiler command and flags
#compiler <- Sys.getenv("CXX")          # e.g., "g++"
#flags <- paste(compile_cmd$compile, collapse = " ")  # Flags and file

system2(
  command = "TMB::compile",
  args = c( "tinyVAST.cpp", "framework = 'TMBad'" )
)
