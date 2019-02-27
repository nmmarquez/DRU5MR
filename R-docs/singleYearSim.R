rm(list=ls())

library(DRU5MR)
library(PointPolygon)
library(tidyverse)
library(sp)

# Im dumb and need to fix this in the creation code
fullDF <- fullDF %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_"))

maxYear <- polyDF %>%
    filter(age_group == "NN") %>%
    group_by(year) %>%
    summarize(n=sum(N)) %>%
    arrange(-n) %>%
    select(year) %>%
    unlist %>%
    unname %>%
    .[1]

nPSU <- polyDF %>%
    filter(age_group == "NN") %>% 
    group_by(strat, psu) %>%
    summarize(n=sum(N)) %>%
    summarise(nPSU=n(), briths=round(mean(n)))

stratSPDF <- maptools::unionSpatialPolygons(spDF, spDF$strat)
stratSPDF$strat <- names(stratSPDF)
stratSPDF$polyid <- 1:nrow(stratSPDF@data)
stratList <- lapply(1:length(stratSPDF), function(i) stratSPDF[i,])

unitSim <- PointPolygon::simField(
    N = 100, # use 500 squares as base field just as in example
    sigmaE = .12, # unit variance for spde process
    rangeE = .7,
    shape = spDF, # null shape which by default creates unit square
    beta0 = -3.5, # intercept 
    link = arm::invlogit,
    offset = 1.,
    max.edge = c(.15,.3))

drSim <- list()

drSim$mesh <- INLA::inla.mesh.2d(
    loc=bbCoords(spDF), 
    offset = 1.,
    max.edge = c(.15, .3))

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
        mutate(id=id-1, Bound=1, V0=1, z=z) %>%
        mutate(theta=arm::invlogit(-3.5 + z))
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

# for each stratification sample population weighted 

test <- samplePolygons(drSim, 10, polygonList = stratList)