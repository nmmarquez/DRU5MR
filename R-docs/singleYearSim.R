.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(PointPolygon)
library(dplyr)
library(tibble)
library(sp)
library(ggplot2)
library(TMB)
set.seed(123)

load("~/Documents/DRU5MR/data/prepData.rda")

# Im dumb and need to fix this in the creation code
fullDF <- fullDF %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_"))

nPSU <- polyDF %>%
    filter(age_group == "NN") %>% 
    group_by(strat, psu) %>%
    summarize(n=sum(N)) %>%
    summarise(nPSU=n(), births=round(mean(n)))

stratSPDF <- maptools::unionSpatialPolygons(spDF, spDF$strat)
stratSPDF$strat <- names(stratSPDF)
stratSPDF$polyid <- 1:nrow(stratSPDF@data)
stratList <- lapply(1:length(stratSPDF), function(i) stratSPDF[i,])

drSim <- list()

drSim$mesh <- INLA::inla.mesh.2d(
    loc=bbCoords(spDF), 
    offset = 1.,
    max.edge = c(.1, .2))

drSim$AprojField <- INLA::inla.spde.make.A(
    mesh=drSim$mesh, 
    loc=as.matrix(fullDF[,c("long", "lat")]))

drSim$bound <- rgeos::gUnaryUnion(spDF, id=rep(1, length(spDF)))
drSim$bound <- sp::SpatialPolygonsDataFrame(
    drSim$bound, data.frame(Bound=1, row.names="1"))

drSim$spde <- INLA::inla.spde2.matern(drSim$mesh)

sigmaE <- .12
rangeE <- .7
kappaE <- sqrt(8) / rangeE
tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)

Qspde <- tauE**2 * 
    (kappaE**4 * drSim$spde$param.inla$M0 +
         2 * kappaE**2 *drSim$spde$param.inla$M1 + drSim$spde$param.inla$M2)

x_ <- as.vector(ar.matrix::sim.AR(n=1, Qspde))
x <- x_ - mean(x_)
z <- as.vector(drSim$AprojField %*% x)

drSim$spdf <- sp::SpatialPointsDataFrame(
    select(fullDF, long, lat), proj4string=spDF@proj4string,
    data = fullDF %>%
        rename(x=long, y=lat) %>%
        left_join(
            yearWDF %>%
                filter(year == 2015 & !is.na(strat)) %>%
                select(-year)) %>%
        mutate(id=(1:n())-1, Bound=1, V0=1, z=z) %>%
        mutate(theta=arm::invlogit(-3.5 + z)) %>%
        as.data.frame
)

drSim$latent <- x
drSim$betas <- -3.5

drSim$spdf@data %>%
    ggplot(aes(x, y, fill=theta)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    ggtitle("")

simPolyDF <- bind_rows(lapply(nPSU$strat, function(s){
    snPSU <- nPSU$nPSU[nPSU$strat == s]
    snSamp <- nPSU$births[nPSU$strat == s] * 10
    subDF <- drSim$spdf@data %>%
        filter(strat == s)
    sampleDF <- tibble(id = sample(subDF$id, snPSU, TRUE, subDF$popW)) %>%
        left_join(subDF, by="id") %>%
        mutate(trials=snSamp, obs=rbinom(n(), trials, theta)) %>%
        mutate(polyid=which(s == nPSU$strat))
    sampleDF
}))

drSim$spdf@data %>%
    ggplot(aes(x, y, fill=theta)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    ggtitle("") +
    geom_point(aes(x=x, y=y, fill=NULL), data=simPolyDF, size=.1, alpha=.3)

drSim$spdf@data %>%
    ggplot(aes(x, y, fill=Population)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    ggtitle("Population") +
    geom_point(aes(x=x, y=y, fill=NULL), data=simPolyDF, size=.1, alpha=.3)


subPointDF <- filter(simPolyDF, urban == 0)
subPolyDF <- filter(simPolyDF, urban == 1)
modelPoint <- runFieldModel(drSim, simPolyDF, verbose = T, moption = 0)
subModelPoint <- runFieldModel(drSim, subPointDF, verbose = T, moption = 0)

save.image(file = "~/Documents/Test/intermediate.Rdata")

mpCI <- list(
    model1 = simulateFieldCI(drSim, modelPoint),
    model2 = simulateFieldCI(drSim, subModelPoint))

ggFieldEst(drSim, mpCI, F)
ggFieldEst(drSim, mpCI, T)

inputs <- buildModelInputs(
  drSim, 
  subPointDF, 
  polyDF=subPolyDF)

inputs$Data$AprojPoly <- Matrix::Matrix(
  0,
  nrow = nrow(drSim$spdf@data),
  ncol = length(unique(subPolyDF$strat)),
  sparse = T)

for(i in 1:length(unique(subPolyDF$strat))){
  s <- unique(subPolyDF$strat)[i]
  subDF <- filter(drSim$spdf@data, strat == s)
  idx <- subDF$id + 1
  inputs$Data$AprojPoly[idx, i] <- subDF$popW
}

inputs$Data$idPoly <- as.numeric(as.factor(subPolyDF$polyid)) - 1

verbose <- T
model <- "PointPolygon2"
compile("~/Documents/Test/PointPolygon2.cpp")
dyn.load(dynlib("~/Documents/Test/PointPolygon2"))
symbolic <- TRUE
control = list(eval.max = 10000, iter.max = 10000)
startTime <- Sys.time()
Obj <- TMB::MakeADFun(data = inputs$Data, parameters = inputs$Params, 
                      DLL = model, random = "z", silent = !verbose)
Obj$env$tracemgc <- verbose
Obj$env$inner.control$trace <- verbose
if (symbolic) {
  nah <- utils::capture.output(TMB::runSymbolicAnalysis(Obj))
}
Opt <- stats::nlminb(start = Obj$par, objective = Obj$fn, 
                     gradient = Obj$gr, control = control)
sdrep <- TMB::sdreport(Obj, getJointPrecision = TRUE)
runtime <- Sys.time() - startTime
modelMix <- list(obj = Obj, opt = Opt, runtime = runtime, 
     moption = 0, stack = NULL, sd=sdrep)


save.image("~/Documents/DRU5MR/data/testRun.Rdata")

mpCI <- list(
  allPoints = simulateFieldCI(drSim, modelPoint),
  censoredPoints = simulateFieldCI(drSim, subModelPoint),
  mixtureModel = simulateFieldCI(drSim, modelMix))

ggFieldEst(drSim, mpCI, F)

ggFieldEst(drSim, mpCI, F) +
  geom_point(aes(x=x, y=y, fill=NULL), data=simPolyDF, size=.01, alpha=.1)

ggFieldEst(drSim, mpCI, T) +
  geom_point(aes(x=x, y=y, fill=NULL), data=simPolyDF, size=.01, alpha=.1)

ggFieldEst(drSim, mpCI, F) +
  geom_point(aes(x=x, y=y, fill=NULL), data=subPointDF, size=.15, alpha=.2)

ggFieldEst(drSim, mpCI, T) +
  geom_point(aes(x=x, y=y, fill=NULL), data=subPointDF, size=.15, alpha=.2)


sapply(mpCI, function(x){
  sqrt(mean((x$mu - x$trueValue)^2))
})

sapply(mpCI, function(x){
  mean(x$lwr <= x$trueValue  & x$upr >= x$trueValue)
})

sapply(mpCI, function(x){
  mean(x$upr - x$lwr)
})
