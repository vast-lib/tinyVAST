
##################
# INSTRUCTIONS
# Modified from C:\Users\James.Thorson\Desktop\Work files\Collaborations\2024 -- tinyVAST\Dryad\Reprex_code_R1.R
###############

#
if( FALSE ){
  setwd( R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2024 -- tinyVAST\Dryad\2025-04-04)' )
  combined_samples = read.csv( file="Formatted_data.csv" )
  goaai_outline = st_read( "goaai_outline.shp" )

  #
  alaska_sponge_coral_fish = list( combined_samples = combined_samples, goaai_outline = goaai_outline )
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)' )
  usethis::use_data(alaska_sponge_coral_fish)
}

library(tinyVAST)
library(fmesher)
library(sf)

#


# Define formula
Formula = Count ~ 0 + interaction(Group4,Year) + offset(log(AreaSwept))

############
# Loop
############

for( config in config_set ){
  run_dir = file.path( date_dir, config )
    dir.create(run_dir, recursive=TRUE )

  if( !("myfit.RDS" %in% list.files(run_dir)) ){
    # Define SEM
    if( config == "same-slope" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4

        Coral -> Flat_trawl, b1
        Sponge -> Flat_trawl, b2
        Coral -> Rock_trawl, b3
        Sponge -> Rock_trawl, b4
      "
    }
    if( config == "diff-slope" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4

        Coral -> Flat_trawl, d1
        Sponge -> Flat_trawl, d2
        Coral -> Rock_trawl, d3
        Sponge -> Rock_trawl, d4
      "
    }
    if( config == "no-link" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4
      "
    }
    if( config == "spatial-Qratio" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4
        Flat -> Flat_trawl, NA, 1
        Rock -> Rock_trawl, NA, 1
      "
    }
    if( config == "spatial-Qratio-hyperstable" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4
        Flat -> Flat_trawl, a1
        Rock -> Rock_trawl, a2
      "
    }
    if( config == "missing-covariate" ){
      sem = "
        Coral -> Flat_trawl, d1
        Sponge -> Flat_trawl, d2
        Coral -> Rock_trawl, d3
        Sponge -> Rock_trawl, d4
      "
    }
    if( config == "spatial-Qratio-linked" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4
        Flat -> Flat_trawl, NA, 1
        Rock -> Rock_trawl, NA, 1
        Coral -> Flat_trawl, d1
        Sponge -> Flat_trawl, d2
        Coral -> Rock_trawl, d3
        Sponge -> Rock_trawl, d4
      "
    }
    if( config == "spatial-Qratio-hyperstable-linked" ){
      sem = "
        Coral -> Flat, b1
        Sponge -> Flat, b2
        Coral -> Rock, b3
        Sponge -> Rock, b4
        Flat -> Flat_trawl, a1
        Rock -> Rock_trawl, a2
        Coral -> Flat_trawl, d1
        Sponge -> Flat_trawl, d2
        Coral -> Rock_trawl, d3
        Sponge -> Rock_trawl, d4
      "
    }
    if( config == "full_covariance" ){    # c("Coral", "Sponge", "Rock", "Flat", "Rock_trawl", "Flat_trawl")
      sem = "
        Coral <-> Sponge, c1
        Coral <-> Rock, c2
        Coral <-> Flat, c3
        Coral <-> Rock_trawl, c4
        Coral <-> Flat_trawl, c5

        Sponge <-> Rock, d1
        Sponge <-> Flat, d2
        Sponge <-> Rock_trawl, d3
        Sponge <-> Flat_trawl, d4

        Rock <-> Flat, e1
        Rock <-> Rock_trawl, e2
        Rock <-> Flat_trawl, e3

        Flat <-> Rock_trawl, f1
        Flat <-> Flat_trawl, f2

        Rock_trawl <-> Flat_trawl, g1
      "
    }
    cat(sem)

    spatial_graph = fm_mesh_2d( loc = Combo[,c('X','Y')],
                                cutoff=10 )
    control = tinyVASTcontrol( getsd = TRUE,
                               nlminb_loops = 2,
                               newton_loops = 1,
                               profile = c("alpha_j"),
                               trace = 1,
                               getJointPrecision = TRUE )
    family = list(
      Coral = tweedie(),
      Sponge = tweedie(),
      Rock = tweedie(),
      Flat = tweedie(),
      Rock_trawl = tweedie(),
      Flat_trawl = tweedie()
    )
    myfit = tinyVAST(
      data = Combo,
      formula = Formula,
      space_term = sem,
      family = family,
      space_columns = c("X", "Y"),
      variable_column = "Group4",
      variables = c("Coral", "Sponge", "Rock", "Flat", "Rock_trawl", "Flat_trawl"),
      time_column = "time",
      distribution_column = "Group4",
      spatial_domain = spatial_graph,
      control = control
    )

    # Save object
    saveRDS( myfit, file=file.path(run_dir,"myfit.RDS") )

    # Save
    myprint = print(myfit)
    myprint$AIC = AIC(myfit)
    capture.output( myprint,
                    file = file.path(run_dir,"print.txt") )
    mysummary = summary(myfit, "space_term")
    capture.output( mysummary,
                    file = file.path(run_dir,"summary.txt") )
  }
}

#################
#
#################

# Custom function for adding legend in multi-panel figure
# https://stackoverflow.com/questions/52975447/reorganize-sf-multi-plot-and-add-a-legend
add_legend <-
function( legend,
          col = sf.colors(),
          legend_x = c(0.9,  1.0),
          legend_y = c(0.05, 0.45),
          text_col = "black",
          ...){

    # Get the axis limits and calculate size
    axisLimits <- par()$usr
    xLength <- axisLimits[2] - axisLimits[1]
    yLength <- axisLimits[4] - axisLimits[3]

    xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
    xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
    yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
    yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]
    if( diff(legend_y) > diff(legend_x) ){
      align = c("lt","rb")[2]
      gradient = c("x","y")[2]
    }else{
      align = c("lt","rb")[1]
      gradient = c("x","y")[1]
    }

    # Add the legend
    plotrix::color.legend( xl = xl,
                           xr = xr,
                           yb = yb,
                           yt = yt,
                           legend = legend,
                           rect.col = col,
                           gradient="y",
                           col = text_col,
                           ... )
}



# Get polygon
library(sf)
goaai_outline = st_read( "goaai_outline" )

# Make extrapolation grid
for( p in 1:2 ){
  if(p==1) grid = st_make_grid( GOAAI, cellsize=c(20e3,20e3) )
  if(p==2) grid = st_make_grid( st_box, cellsize=c(5e3,5e3) )

  # Intersect and project
  grid = st_intersection( grid, GOAAI )
  grid = st_make_valid( grid )
  grid = st_transform( grid, crs=st_crs('+proj=natearth +lon_0=-170 +units=km') )

  # Make long-form
  loc_gz = st_coordinates(st_centroid( grid ))
  colnames(loc_gz) = c("X", "Y")

  png( file=file.path(run_dir, paste0("density-",c("full","bank")[p],".png")), width=6, height=6, res=200, units="in" )
    plots = data.frame( "Group4"=myfit$internal$variables )
      plots$Year = ifelse( plots$Group4 %in% c("Flat_trawl","Rock_trawl"), 1990, "camera" )
    if(p==1) par( mfrow=c(6,2) )
    if(p==2) par( mfrow=c(3,4) )
    for( i in seq_len(nrow(plots)) ){
      pred_orig = pred1 = predict( myfit, newdata=data.frame(loc_gz, Group4=plots$Group4[i], Year=plots$Year[i], AreaSwept=1), what="p_g", se.fit=TRUE )
      pred1$se.fit = ifelse( pred1$fit > max(pred1$fit-log(Inf),na.rm=TRUE), pred1$se.fit, NA )
      pred1$fit = ifelse( pred1$fit > max(pred1$fit-log(Inf),na.rm=TRUE), pred1$fit, NA )
      sf_plot = st_sf( grid, pred1$fit, pred1$se.fit, remove_origdata=TRUE )
      # Densities
      plot( sf_plot[,1], reset=FALSE, key.pos=NULL, border=NA, main=plots$Group4[i], pal=viridisLite::viridis )
      plot( goaai_outline, add=TRUE, lwd=2 )
      if(p==1) add_legend( legend=round(range(sf_plot[[1]],na.rm=TRUE),2), legend_x=c(0.25,0.3), legend_y=c(0.05,0.45), col=viridisLite::viridis(10) )
      if(p==2) add_legend( legend=round(range(sf_plot[[1]],na.rm=TRUE),2), legend_x=c(0.25,0.3), legend_y=c(0.55,0.95), col=viridisLite::viridis(10) )
      # Standard errors
      plot( sf_plot[,2], reset=FALSE, key.pos=NULL, border=NA, main=paste("SE[",plots$Group4[i],"]"), pal=sf.colors )
      plot( goaai_outline, add=TRUE, lwd=2 )
      if( p==2 ){
        tmp = subset( Combo, Group4==plots$Group4[i])
        points( x=tmp$X, y=tmp$Y, pch=19, cex=0.5, col="green" )
      }
      if(p==1) add_legend( legend=round(range(sf_plot[[2]]),2), legend_x=c(0.25,0.3), legend_y=c(0.05,0.45) )
      if(p==2) add_legend( legend=round(range(sf_plot[[2]]),2), legend_x=c(0.25,0.3), legend_y=c(0.55,0.95) )
    }
  dev.off()
}

#################
# Single graph with estimates
#################

#
config = "spatial-Qratio-linked"
  run_dir = file.path( date_dir, config )
  myfit = readRDS( file.path(run_dir,"myfit.rds") )

library(igraph)

layout = cbind( x=c(2,2,1,3,1,3), y=c(2,3,1,1,4,4) )
  rownames(layout) = c( "Coral", "Sponge", "Rock", "Flat", "Rock_trawl", "Flat_trawl" )
Switch = \(x) switch( x, "Coral"="C", "Sponge"="S",
                       "Rock"="R", "Flat"="F", "Rock_trawl"="R2", "Flat_trawl"="F2" )
png( file=file.path(run_dir,"estimates.png"), width=3, height=3, res=200, units="in" )
  par( mar=c(0,0,0,0) )
  #out = myfit$internal$sem_ram_output
  out = summary( myfit )
  variables = as.character(myfit$internal$variables)
  var_labels = sapply( variables, FUN=Switch )
  label = paste0( formatC(out$Estimate,digits=2,format="f"),
                  "\n(",round(out$Std_Error,2),")" )                            # formatC(out$Std_Error,digits=2,format="f")
  df <- data.frame( from = sapply(out$from,FUN=Switch),
                    to = sapply(out$to,FUN=Switch),
                    label = label,
                    x = NA,
                    y = NA )
  #df$x[c(1:4,7:10)] = c(0.6,0.6,-0.6,-0.6)
  df$x = ifelse( df$from %in% c("C","S") & df$to %in% c("F","F2"), 0.5, df$x )
  df$x = ifelse( df$from %in% c("C","S") & df$to %in% c("R","R2"), -0.5, df$x )
  df$y = ifelse( df$from=="C" & df$to %in% c("R","F"), -0.66, df$y )
  df$y = ifelse( df$from=="C" & df$to %in% c("R2","F2"), 0.33, df$y )
  df$y = ifelse( df$from=="S" & df$to %in% c("R","F"), -0.33, df$y )
  df$y = ifelse( df$from=="S" & df$to %in% c("R2","F2"), 0.66, df$y )
  df = subset( df, from!=to )
  pg <- graph_from_data_frame(d = df[,1:3], directed=TRUE, vertices=data.frame(var_labels) )
  plot.igraph( pg, vertex.shape="rectangle", vertex.size=0, vertex.size2=0, vertex.label.cex=2,
        vertex.color="grey", vertex.label.color="black", edge.label.color="black",
        edge.label.cex=1, layout=layout[as.character(variables),], xlim=1.2*c(-1,1), ylim=1.2*c(-1,1),
        edge.label.dist=0.8, edge.label.x=df$x, edge.label.y=df$y )
dev.off()


#################
# Multiple graphs
#################

#
config_set = c(
  "no-link",
  "missing-covariate",
  "same-slope",
  "spatial-Qratio",
  "spatial-Qratio-hyperstable",
  "diff-slope",
  "spatial-Qratio-linked",
  "spatial-Qratio-hyperstable-linked"
)


AIC = rep(NA, length(config_set) )
for( i in seq_along(config_set) ){
  run_dir = file.path( date_dir, config_set[i] )
  myfit = readRDS( file.path(run_dir,"myfit.rds") )
  AIC[i] = AIC(myfit)
}
AIC_table = data.frame( "config"=config_set, "deltaAIC"=AIC-min(AIC) )
write.csv( AIC_table, file=file.path(date_dir,"AIC_table.csv") )

library(igraph)

layout = cbind( x=c(2,2,1,3,1,3), y=c(2,3,1,1,4,4) )
  rownames(layout) = c( "Coral", "Sponge", "Rock", "Flat", "Rock_trawl", "Flat_trawl" )
  config_labels = sapply( config_set, FUN=switch,
                         "same-slope" = "Same response",
                         "diff-slope" = "Diff. response",
                         "no-link" = "No link",
                         "spatial-Qratio" = "Spatial Q",
                         "spatial-Qratio-hyperstable" = "Hyperstable\n spatial Q",
                         "missing-covariate" = "Missing covariate",
                         "spatial-Qratio-linked" = "Spatial Q\nand diff. response",
                         "spatial-Qratio-hyperstable-linked" = "Hyperstable spatial Q\nand diff. response")
png( file=file.path(date_dir,"graphs.png"), width=4, height=8, res=200, units="in" )
  par( mfrow=c(4,2), mar=c(1,0,2,0) )
  for( i in seq_along(config_set) ){
    run_dir = file.path( date_dir, config_set[i] )
    myfit = readRDS( file.path(run_dir,"myfit.rds") )
    out = myfit$internal$space_term_ram$output
    variables = as.character(myfit$internal$variables)
    var_labels = sapply( variables, FUN=switch, "Coral"="C", "Sponge"="S",
                         "Rock"="R", "Flat"="F", "Rock_trawl"="R2", "Flat_trawl"="F2" )
    df <- data.frame(from = var_labels[out$ram$from], to = var_labels[out$ram$to], label = out$model[,2])
    df$label = ifelse( df$from==df$to, "", df$label )
    df = subset( df, from!=to )
    pg <- graph_from_data_frame(d = df, directed=TRUE, vertices=data.frame(var_labels) )
    plot( pg, vertex.shape="rectangle", vertex.size=0, vertex.size2=0, vertex.label.cex=2,
          vertex.color="grey", vertex.label.color="black", edge.label.color="black",
          edge.label.cex=1.5, layout=layout[as.character(variables),], xlim=1.2*c(-1,1), ylim=1.2*c(-1,1) )
    title( config_labels[i] )
    legend( "top", bty="n", legend=paste0("Params = ",attr(logLik(myfit),"df")) ) # length(myfit$opt$par)+length(myfit$internal$parlist$alpha_j)) )
    legend( "bottom", bty="n", legend=paste0("AIC = ",round(AIC_table[i,'deltaAIC'],1)) )
    box()
  }
dev.off()

