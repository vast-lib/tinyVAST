
setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )

# Download Excel:
# https://afspubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fmcf2.10023&file=mcf210023-sup-0001-TableS1-S24.xlsx
# and then export tab "ST 13-16 Tot ret (bioma) 52-15.csv" as CSV

# Read CSV
CSV = read.csv( file.path( getwd(), "data-raw", "mcf210023-sup-0001-tables1-s24 -- ST 13-16 Tot ret (bioma) 52-15.csv"),
                skip=4, sep="\t")

# Exclude unknown Management.Area
Pink = CSV[,2:15]
dimnames(Pink) = list( "Year"=CSV[,1], "Region"=colnames(Pink) )
Chum = CSV[,19:32]
dimnames(Chum) = list( "Year"=CSV[,18], "Region"=colnames(Pink) )
Sockeye = CSV[,36:49]
dimnames(Sockeye) = list( "Year"=CSV[,35], "Region"=colnames(Pink) )
f = \(x,sp) data.frame( "Species"=sp, expand.grid(dimnames(x)), "Biomass"=unlist(x) )
salmon_returns = rbind( f(Pink,"pink"), f(Chum,"chum"), f(Sockeye,"sockeye") )

# Fix formatting
colnames(salmon_returns)[2:3] = c("Year","Region")
salmon_returns$Year = as.numeric(as.character(salmon_returns$Year))

#
usethis::use_data(salmon_returns)
