
library(sdmTMB)
library(tinyVAST)
set.seed(123)

year_set = c(2021,2022)
species_set = c("anchovy", "sardine", "pollock")
n_j = 5

df = expand.grid(
  species = species_set,
  sample = seq_len(n_j),
  year = year_set
)
df$year_sample = interaction( df$year, df$sample )
df$species_year = interaction( df$species, df$year )

# Fake coordinates for sdmTMB
xy_jt = matrix( rnorm(nlevels(df$year_sample)*2), ncol = 2 )
df = cbind( df, setNames(as.data.frame(xy_jt),c("X","Y")) )

b_s = rnorm(length(species_set))
b_st = rnorm( length(species_set) * length(year_set) )
b_jt = rnorm( n_j * length(year_set) )

p_i = b_s[df$species] + b_st[df$species_year] + b_jt[df$year_sample]

# Normal errors
Y_i = rnorm( 
  n = nrow(df),
  mean = p_i,
  sd = 1
)

##################
# FAILS
# Does not propertly detect corner constraint
##################

mytv = tinyVAST(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + species_year + year_sample
)

# Also doesn't detect using nested notation
mytv = tinyVAST(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + year / species + factor(year) / factor(sample)
)

##################
# Duplicate in sdmTMB ...
# same issue
##################

mesh = make_mesh(
  data = df,
  xy_cols = c("X","Y"),
  n_knots = nlevels(df$year_sample)
)

mysdmTMB = sdmTMB(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + species_year + year_sample,
  mesh = mesh,
  spatial = "off"
)


##################
# Duplicate in lm ...
# same issue
##################

mylm = lm(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + species_year + year_sample
)

##################
# WORKS
# Manually code corner constraint
##################

# Fix level
df$level = as.character(df$year_sample)
for( year in unique(df$year) ){
  which_rows = which( df$year == year )
  which_base = unique(df[which_rows,'level'])[1]
  df$level = ifelse( df$level == which_base, "base", df$level )
}
df$level = factor(df$level)
df$level = relevel( df$level, "base" )

# Run
mytv2 = tinyVAST(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + species_year + level
)

##################
# repeat in sdmTMB
##################

mysdmTMB2 = sdmTMB(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + species_year + level,
  mesh = mesh,
  spatial = "off"
)

#################
# repeat in lm
#################

mylm2 = lm(
  data = cbind( df, y = Y_i ),
  formula = y ~ 0 + species_year + level
)
 
##################
# Try new package
##################

library(fullRankMatrix)
X <- model.matrix( ~ 0 + species_year + year_sample, df )
X_full <- make_full_rank_matrix(X)$matrix

mylm3 = lm(
  formula = Y_i ~ 0 + X_full
)

