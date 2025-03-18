 #
    library(tinyVAST)

    # Simulate a seperable two-dimensional AR1 spatial process
     #set.seed(101) # No errors
     set.seed(123)  # HAS ERRORS
     n_x = n_y = 25
     n_w = 10
     R_xx = exp(-0.4 * abs(outer(1:n_x, 1:n_x, FUN="-")) )
     R_yy = exp(-0.4 * abs(outer(1:n_y, 1:n_y, FUN="-")) )
     z = mvtnorm::rmvnorm(1, sigma=kronecker(R_xx,R_yy) )

     # Simulate nuissance parameter z from oscillatory (day-night) process
     w = sample(1:n_w, replace=TRUE, size=length(z))
     Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), w=w, z=as.vector(z) + cos(w/n_w*2*pi))
     Data$n = Data$z + rnorm(nrow(Data), sd=1)

     # Add columns for multivariate and/or temporal dimensions
     Data$var = "n"

     # make SPDE mesh for spatial term
     mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], n=100 )

     # fit model with cyclic confounder as GAM term
     out = tinyVAST( data = Data,
                     formula = n ~ s(w),
                     control = tinyVASTcontrol(
                       calculate_deviance_explained = FALSE,
                       nlminb_loops = FALSE,
                       getsd = 0
                     ) )
     #rep = out$obj$report()

     #out = tinyVAST( data = Data,
     #                formula = n ~ s(w) )
     #out
