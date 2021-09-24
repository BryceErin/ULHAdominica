########################
## Data and libraries ##
########################
library(INLA)
landslides <- read.csv("ULHAlandslides.csv")
coords     <- read.csv("ULHAcoords.csv")

######################
## Data Preparation ##
######################
## Scale variables
landslides$SU_Area                   <- scale(landslides$SU_Area)
landslides$SU_Perimeter              <- scale(landslides$SU_Perimeter)
landslides$SU_P.A                    <- scale(landslides$SU_P.A)
landslides$SU_P.sqrt.A               <- scale(landslides$SU_P.sqrt.A.)
landslides$SU_Max.Distance           <- scale(landslides$SU_Max.Distance)
landslides$SU_D.A                    <- scale(landslides$SU_D.A)
landslides$SU_D.sqrt.A               <- scale(landslides$SU_D.sqrt.A.)
landslides$Dist2Stream..MEAN.        <- scale(landslides$Dist2Stream..MEAN.)
landslides$Dist2Stream..STDDEV       <- scale(landslides$Dist2Stream..STDDEV.)
landslides$Eastness..MEAN.           <- scale(landslides$Eastness..MEAN.)
landslides$Eastness..STDDEV.         <- scale(landslides$Eastness..STDDEV.)
landslides$Elevation..MEAN.          <- scale(landslides$Elevation..MEAN.)
landslides$Elevation..STDDEV.        <- scale(landslides$Elevation..STDDEV.)
landslides$Northness..MEAN.          <- scale(landslides$Northness..MEAN.)
landslides$Northness..STDDEV.        <- scale(landslides$Northness..STDDEV.)
landslides$PlanarCurvature..MEAN.    <- scale(landslides$PlanarCurvature..MEAN.)
landslides$PlanarCurvature..STDDEV.  <- scale(landslides$PlanarCurvature..STDDEV.)
landslides$ProfileCurvature..MEAN.   <- scale(landslides$ProfileCurvature..MEAN.)
landslides$ProfileCurvature..STDDEV. <- scale(landslides$ProfileCurvature..STDDEV.)
landslides$RSP..MEAN.                <- scale(landslides$RSP..MEAN.)
landslides$RSP..STDDEV.              <- scale(landslides$RSP..STDDEV.)
landslides$Slope..MEAN.              <- scale(landslides$Slope..MEAN.)
landslides$Slope..STDDEV.            <- scale(landslides$Slope..STDDEV.)
landslides$TWI..MEAN.                <- scale(landslides$TWI..MEAN.)
landslides$TWI..STDDEV.              <- scale(landslides$TWI..STDDEV.)

## Create grouped version of variables (needed for RW1 covariates)
landslides$SU_Area.levels                <- inla.group(landslides$SU_Area, method = "quantile", n = 20)
landslides$SU_Max.Distance.levels        <- inla.group(landslides$SU_Max.Distance, method="quantile", n = 20)
landslides$Dist2Stream..MEAN.levels      <- inla.group(landslides$Dist2Stream..MEAN., method = "quantile", n = 20)
landslides$Elevation..MEAN.levels        <- inla.group(landslides$Elevation..MEAN., method = "quantile", n = 20)
landslides$ProfileCurvature..MEAN.levels <- inla.group(landslides$ProfileCurvature..MEAN., method = "quantile", n = 20)
landslides$Slope..MEAN.levels            <- inla.group(landslides$Slope..MEAN., method = "quantile", n = 20)
landslides$TWI..MEAN.levels              <- inla.group(landslides$TWI..MEAN., method = "quantile", n = 20)

## Preparing INLA data frame
data.inla                                       <- landslides
data.inla$y.size                                <- log(data.inla$LS_Maria2017)
data.inla$y.size[is.infinite(data.inla$y.size)] <- NA
data.inla$y.pres                                <- ifelse(data.inla$LS_Maria2017>0, 1, 0)

data.inla <- cbind(coords, 
                   data.inla[match("y.size", names(data.inla))], 
                   data.inla[match("y.pres", names(data.inla))],
                   data.inla[-match(c("SU_ID", "y.size", "y.pres"), names(data.inla))])

#######################
## Model Preparation ##
#######################
## Projection matrix A, spde object and mesh
mesh <- inla.mesh.2d(loc=  data.inla[,c('Lat', 'Lon')], max.edge = c(2, 5), cutoff   = 0.001)
A    <- inla.spde.make.A(mesh, loc = as.matrix(data.inla[,c('Lat', 'Lon')]))
spde <- inla.spde2.pcmatern(mesh        = mesh, 
                            alpha       = 2, # alpha = 2 corresponds to smoothness = 1 (see Equation 1 in Lindgren at al., 2011, JRSS-B)
                            prior.range = c(25, .5), # P(range < 25) = 0.5
                            prior.sigma = c(0.5, .5)) # (P(sd > 0.5) = 0.5)

mesh.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

## Stack: observations and covariates
effects <- list(c(mesh.index, list(Intercept = 1)),
                list(area     = data.inla$SU_Area.levels,
                     maxd     = data.inla$SU_Max.Distance,
                     maxd.l   = data.inla$SU_Max.Distance.levels,
                     d.a      = data.inla$SU_D.A,
                     d.sqrt.a = data.inla$SU_D.sqrt.A,
                     dist     = data.inla$Dist2Stream..MEAN.,
                     dist.l   = data.inla$Dist2Stream..MEAN.levels,
                     distsd   = data.inla$Dist2Stream..STDDEV,
                     eastm    = data.inla$Eastness..MEAN.,
                     northm   = data.inla$Northness..MEAN.,
                     eastsd   = data.inla$Eastness..STDDEV.,
                     northsd  = data.inla$Northness..STDDEV.,
                     elevm    = data.inla$Elevation..MEAN.levels,
                     elevsd   = data.inla$Elevation..STDDEV.,
                     planm    = data.inla$PlanarCurvature..MEAN.,
                     plansd   = data.inla$PlanarCurvature..STDDEV.,
                     profm    = data.inla$ProfileCurvature..MEAN.,
                     profm.l  = data.inla$ProfileCurvature..MEAN.levels,
                     profsd   = data.inla$ProfileCurvature..STDDEV.,
                     rspm     = data.inla$RSP..MEAN.,
                     rspsd    = data.inla$RSP..STDDEV.,
                     slopem   = data.inla$Slope..MEAN.levels,
                     slopesd  = data.inla$Slope..STDDEV.,
                     twim     = data.inla$TWI..MEAN.,
                     twim.l   = data.inla$TWI..MEAN.levels,
                     twisd    = data.inla$TWI..STDDEV.,
                     YPV      = data.inla$YPV, 
                     BAY      = data.inla$BAY,
                     OPL      = data.inla$OPL,
                     RGA      = data.inla$RGA,
                     YPC      = data.inla$YPC, 
                     YPI      = data.inla$YPI, 
                     IGO      = data.inla$IGO,
                     IPL      = data.inla$IPL, 
                     IGY      = data.inla$IGY, 
                     PBA      = data.inla$PBA, 
                     PPD      = data.inla$PPD))

stack.size <- inla.stack(tag     = 'data.stack',
                         data    = list(y = data.inla$y.size),
                         A       = list(A, 1),
                         effects = effects)

stack.pres <- inla.stack(tag     = 'data.stack',
                         data    = list(y = data.inla$y.pres),
                         A       = list(A, 1),
                         effects = effects)

## PC priors and INLA formula
u.area   <- 5
u.area   <- sd(data.inla$SU_Area.levels)
u.maxd   <- 0.5
u.dist   <- sd(data.inla$Dist2Stream..MEAN.levels)
u.twim   <- 0.1
u.slopem <- sd(data.inla$Slope..MEAN.levels)
u.elevm  <- sd(data.inla$Elevation..MEAN.levels)
u.profm  <- sd(data.inla$ProfileCurvature..MEAN.levels)

alpha.area   <- 0.01
alpha.maxd   <- 0.01
alpha.dist   <- 0.01
alpha.twim   <- 0.01
alpha.slopem <- 0.01
alpha.elevm  <- 0.01
alpha.profm  <- 0.01

hyper.area   <- list(theta = list(prior = "pc.prec", param = c(u.area, alpha.area))) 
hyper.maxd   <- list(theta = list(prior = "pc.prec", param = c(u.maxd, alpha.maxd)))
hyper.dist   <- list(theta = list(prior = "pc.prec", param = c(u.dist, alpha.dist)))
hyper.twim   <- list(theta = list(prior = "pc.prec", param = c(u.twim, alpha.twim)))
hyper.slopem <- list(theta = list(prior = "pc.prec", param = c(u.slopem, alpha.slopem)))
hyper.elevm  <- list(theta = list(prior = "pc.prec", param = c(u.elevm, alpha.elevm)))
hyper.profm  <- list(theta = list(prior = "pc.prec", param = c(u.profm, alpha.profm)))

cyclic = TRUE

formula.sizes = y ~ -1 + Intercept + 
  d.a + d.sqrt.a + distsd + eastm + eastsd + northm + northsd + 
  elevsd + planm + plansd + profm + profsd + rspm + rspsd +
  slopesd + twisd + YPV + BAY + OPL +
  f(area, model = "rw1", hyper = hyper.area, cyclic = cyclic) + 
  f(maxd.l, model = "rw1", hyper = hyper.maxd, cyclic = cyclic) + 
  f(dist.l, model = "rw1", hyper = hyper.dist, cyclic = cyclic) + 
  f(elevm, model = "rw2", hyper = hyper.elevm, cyclic = cyclic) +
  f(slopem, model = "rw2",hyper = hyper.slopem, cyclic = cyclic) +
  f(twim.l, model = "rw1", hyper = hyper.twim, cyclic = cyclic) +
  f(spatial.field, model = spde)

formula.pres = y ~ -1 + Intercept + 
  maxd  + d.a + d.sqrt.a + dist + distsd + eastm + eastsd + northm + northsd +
  elevsd + planm + plansd + profsd + rspm + rspsd + slopesd + 
  twim + twisd + RGA + YPC + YPV + YPI + BAY + IGO + IPL + IGY + OPL + PBA + PPD + 
  f(area, model = "rw1", hyper = hyper.area, cyclic = cyclic) +
  f(elevm, model = "rw2", hyper = hyper.elevm, cyclic = cyclic) +
  f(slopem, model = "rw2", hyper = hyper.slopem, cyclic = cyclic) +
  f(profm.l, model = "rw1", hyper = hyper.profm, cyclic = cyclic) +
  f(spatial.field, model = spde)

## Setting inla control arguments shared by both model
control.compute   <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)


####################
## Size model fit ##
####################
sd(data.inla$y.size, na.rm = T)
fit.size <- inla(formula.sizes, 
                 data              = inla.stack.data(stack.size), 
                 family            = "gaussian",
                 control.family    = list(hyper=list(theta = list(prior="pc.prec", param=c(2, 0.01)))), # P(measurement error sd > sd data) = 0.1
                 control.predictor = list(A = inla.stack.A(stack.size), compute = TRUE),
                 control.compute   = control.compute,
                 control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb"),
                 verbose           = F)


########################
## Presence model fit ##
########################
fit.pres <- inla(formula.pres, 
                 data              = inla.stack.data(stack.pres), 
                 family            = "binomial",
                 control.predictor = list(A = inla.stack.A(stack.pres), compute = TRUE),
                 control.compute   = control.compute,
                 control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb"),
                 verbose           = TRUE)

