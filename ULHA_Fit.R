###################################################################################################################
##                                                                                                               ##
##----------------- Unified landslide hazard assessment: a case study in the Island of Dominica -----------------##
## Code file to fit the landslide sizes and presence model (log-Gaussian amd Bernoulli likelihood, respectively) ##
## Erin Bryce - School of Mathematics and Statistics, University of Glasgow, UK                                  ##
##                                                                                                               ##
###################################################################################################################

########################
## Data and libraries ##
########################
library(INLA)
data <- read.csv("ULHAdominica.csv")
y.size <- data$y.logsize
y.pres <- data$y.pres
coords = cbind(data$lon, data$lat)

#######################
## Model preparation ##
#######################
## Projection matrix A, spde object and mesh
mesh <- inla.mesh.2d(loc= coords, 
                     max.edge = c(2, 5),
                     cutoff   = 0.001)

A <- inla.spde.make.A(mesh, loc = coords)
spde <- inla.spde2.pcmatern(mesh        = mesh, 
                            alpha       = 2, # alpha = 2 corresponds to smoothness = 1 (see Equation 1 in Lindgren at al., 2011, JRSS-B)
                            prior.range = c(25, .5), # P(range < 25) = 0.5
                            prior.sigma = c(0.5, .5)) # (P(sd > 0.5) = 0.5)
mesh.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

## Stack: observations and covariates
effects <- list(c(mesh.index, list(Intercept = 1)),
                list(area     = data$scaledSU_Area.levels,
                     maxd     = data$scaledSU_Max.Distance.levels,
                     d.a      = data$scaledSU_D.A,
                     d.sqrt.a = data$scaledSU_D.sqrt.A,
                     dist     = data$scaledDist2Stream..MEAN.levels,
                     distsd   = data$scaledDist2Stream..STDDEV,
                     eastm    = data$scaledEastness..MEAN.,
                     northm   = data$scaledNorthness..MEAN.,
                     eastsd   = data$scaledEastness..STDDEV.,
                     northsd  = data$scaledNorthness..STDDEV.,
                     elevm    = data$scaledElevation..MEAN.levels,
                     elevsd   = data$scaledElevation..STDDEV.,
                     planm    = data$scaledPlanarCurvature..MEAN.,
                     plansd   = data$scaledPlanarCurvature..STDDEV.,
                     profm    = data$scaledProfileCurvature..MEAN.,
                     profsd   = data$scaledProfileCurvature..STDDEV.,
                     rspm     = data$scaledRSP..MEAN.,
                     rspsd    = data$scaledRSP..STDDEV.,
                     slopem   = data$scaledSlope..MEAN.levels,
                     slopesd  = data$scaledSlope..STDDEV.,
                     twim     = data$scaledTWI..MEAN.levels,
                     twisd    = data$scaledTWI..STDDEV.,
                     YPV      = data$YPV, 
                     BAY      = data$BAY,
                     OPL      = data$OPL,
                     RGA      = data$RGA,
                     YPC      = data$YPC, 
                     YPI      = data$YPI, 
                     IGO      = data$IGO,
                     IPL      = data$IPL, 
                     IGY      = data$IGY, 
                     PBA      = data$PBA, 
                     PPD      = data$PPD))

stack.size <- inla.stack(tag     = 'data.stack',
                         data    = list(y = data$y.logsize),
                         A       = list(A, 1),
                         effects = effects)

stack.pres <- inla.stack(tag     = 'data.stack',
                         data    = list(y = data$y.pres),
                         A       = list(A, 1),
                         effects = effects)

## PC priors and INLA formula
u.area   <- 5
u.maxd   <- 0.5
u.dist   <- sd(data$scaledDist2Stream..MEAN.levels)
u.twim   <- 0.1
u.slopem <- sd(data$scaledSlope..MEAN.levels)
u.elevm  <- sd(data$scaledElevation..MEAN.levels)
u.profm  <- sd(data$scaledProfileCurvature..MEAN.levels)

alpha.area   <- 0.01
alpha.maxd   <- 0.01
alpha.dist   <- 0.01
alpha.twim   <- 0.9
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
  d.a + d.sqrt.a + distsd + eastm + northm + eastsd + northsd + 
  elevsd + planm + plansd + profm + profsd + rspm + rspsd +
  slopesd + twisd + YPV + BAY + OPL +
  f(area, model = "rw1", hyper = hyper.area, cyclic = cyclic) + 
  f(maxd, model = "rw1", hyper = hyper.maxd, cyclic = cyclic) + 
  f(dist, model = "rw1", hyper = hyper.dist, cyclic = cyclic) + 
  f(elevm, model = "rw2", hyper = hyper.elevm, cyclic = cyclic) +
  f(slopem, model = "rw2",hyper = hyper.slopem, cyclic = cyclic) +
  f(twim, model = "rw1", hyper = hyper.twim, cyclic = cyclic) +
  f(spatial.field, model = spde)

formula.pres = y ~ -1 + Intercept + 
  maxd  + d.a + d.sqrt.a + dist + distsd + eastm + northm + eastsd + northsd +
  elevsd + planm + plansd + profsd + rspm + rspsd + slopesd + 
  twim + twisd + RGA + YPC + YPV + YPI + BAY + IGO + IPL + IGY + OPL + PBA + PPD + 
  f(area, model = "rw1", hyper = hyper.area, cyclic = cyclic) +
  f(elevm, model = "rw2", hyper = hyper.elevm, cyclic = cyclic) +
  f(profm, model = "rw1", hyper = hyper.profm, cyclic = cyclic) +
  f(slopem, model = "rw2", hyper = hyper.slopem, cyclic = cyclic) +
  f(spatial.field, model = spde)



####################
## Size model fit ##
####################
fit.size <- inla(formula.sizes, 
                 data              = inla.stack.data(stack.size), 
                 family            = "gaussian",
                 control.family    = list(hyper=list(theta = list(prior="pc.prec", param=c(2, 0.01)))), # P(measurement error sd > sd data) = 0.1
                 control.predictor = list(A = inla.stack.A(stack.size), compute = TRUE),
                 control.compute   = list(cpo = TRUE, dic = TRUE, config=TRUE),
                 control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb"),
                 verbose           = F)

save(fit.size, file = 'FitSize13092021_INLA200629_R20200229.Rdata')

########################
## Presence model fit ##
########################
fit.pres <- inla(formula.pres, 
                 data              = inla.stack.data(stack.pres), 
                 family            = "binomial",
                 control.predictor = list(A = inla.stack.A(stack.pres), compute = TRUE),
                 control.compute   = list(cpo = TRUE, dic = TRUE, config=TRUE),
                 control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb"),
                 verbose           = TRUE)

# if complains that locations are too close for profm, run the following lines and rerun the fit
m = get("inla.models", INLA:::inla.get.inlaEnv())
m$latent$rw1$min.diff = NULL
assign("inla.models", m, INLA:::inla.get.inlaEnv())

save(fit.pres, file = 'FitPres13092021_INLA200629_R20200229.Rdata') # takes 1.38hours...


