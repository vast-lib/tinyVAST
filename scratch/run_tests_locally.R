

test_dir = R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\tests\testthat)'
save_dir = R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\scratch)'

library(tinyVAST)
library(testthat)

tests = list.files(test_dir)

sink( file.path(save_dir,"test_output.txt") )
#zzfil <- tempfile(fileext=".txt", tmpdir = save_dir )
#zz <- file(zzfil, "w")  # open an output file connection
for(i in seq_along(tests) ){
  source( file.path(test_dir, tests[i]) )
}
#close(zz)
sink()

test_package( "tinyVAST" )
