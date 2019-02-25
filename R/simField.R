rm(list=ls())

library(haven)
library(rgdal)
library(sp)
library(maptools)
library(PointPolygon)
library(rgeos)
library(raster)
library(tidyverse)
library(INLA)
library(ar.matrix)

load("./Data/prepData.Rdata")

plot(meshDR <- inla.mesh.2d(
    loc=bbCoords(spDF), 
    offset = c(1),
    max.edge = c(.3,.6)))

AprojDR <- INLA::inla.spde.make.A(
    mesh=meshDR, 
    loc=as.matrix(fullDF[,c("long", "lat")]))

M <- 9
rangeE <- .4
sigmaE <- .04

kappaE <- sqrt(8) / rangeE
tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)
spde <- INLA::inla.spde2.matern(meshDR)
Qspde <- tauE**2 * (kappaE**4 * spde$param.inla$M0 +
                    2 * kappaE**2 *spde$param.inla$M1 + spde$param.inla$M2)
Qar1 <- Q.AR1(M, 1, .9) 
Q <- kronecker(Qar1, Qspde)
x_ <- as.vector(ar.matrix::sim.AR(n=1, Q))
x <- x_ - mean(x_)

for(i in 1:M){
    fullDF[,paste0("x", i)] <- as.vector(
        AprojDR %*% x[(meshDR$n * (i-1) + 1):(meshDR$n * i)])
}

pAgeDF <- tibble(
    ageB = arm::logit(c(.007, .006, .005, .004, .003, .002, .001) + .00025),
    age  = 1:7)

pTimeDF <- tibble(
    time = 1:M,
    rwTime = c(r.AR1(1, M, .001, .9999) %>% `-`(mean(.)))
)

fullDF <- fullDF %>%
    gather(key="time", value="x", -long, -lat, -reg, -urban, -id, -strat) %>%
    mutate(time=as.integer(gsub("x", "", time))) %>%
    right_join(expand.grid(id=1:max(fullDF$id), time=1:M, age=1:7)) %>%
    left_join(pTimeDF) %>%
    left_join(pAgeDF) %>%
    mutate(deltaU=ifelse(urban == 1, -.01, 0)) %>%
    mutate(p = arm::invlogit(ageB + deltaU + rwTime + x))

fqzDF <- fullDF %>%
    group_by(long, lat, id, reg, urban, strat, time) %>%
    summarize(fqz=1-prod(1-p)) %>%
    ungroup

fqzDF %>%
    mutate(long=round(long, 5), lat=round(lat, 5)) %>%
    ggplot(aes(x=long, y=lat, fill=fqz)) +
    geom_raster() +
    facet_wrap(~time) +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

summary(fqzDF)

rasterFull <- stack(lapply(1:M, function(i){
    spFullDF <- fqzDF %>%
        filter(time == i) %>%
        select(long, lat, fqz)
    coordinates(spFullDF) <- ~ long + lat
    proj4string(spFullDF) <- proj4string(spDF)
    gridded(spFullDF) <- TRUE
    rasterDF <- raster(spFullDF)
    rasterDF
}))

plot(rasterFull, legend=FALSE)

simField <- function(N=60, sigmaE=1, rangeE=.3, shape=NULL,
                     beta0=0, betaList=list(), link=arm::invlogit, ...){
    if(is.null(shape)){
        shape <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(
            matrix(c(0,1,1,0,0,0,1,1), ncol=2))), 1)))
    }
    bb <- shape@bbox
    shape$isPresent <- TRUE
    baseRaster <- raster::raster(
        ncols=N, nrows=round(N*hRatio(shape)),
        xmn=bb[1,1], xmx=bb[1,2], ymn=bb[2,1], ymx=bb[2,2],
        crs=shape@proj4string)
    shapeRaster <- raster::rasterize(shape, baseRaster, field=1)
    shapeExtPointsDF <- sp::SpatialPointsDataFrame(
        sp::coordinates(shapeRaster),
        data=as.data.frame(sp::coordinates(shapeRaster)),
        proj4string = shape@proj4string)
    
    validIndex <- !is.na(sp::over(shapeExtPointsDF, shape)$isPresent)
    shapePointsDF <- shapeExtPointsDF[validIndex,]
    mesh <- INLA::inla.mesh.2d(loc=bbCoords(shape), ...)
    AprojField <- INLA::inla.spde.make.A(mesh=mesh, loc=shapePointsDF)
    
    kappaE <- sqrt(8) / rangeE
    tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)
    spde <- INLA::inla.spde2.matern(mesh)
    Q <- tauE**2 * (kappaE**4 * spde$param.inla$M0 +
                        2 * kappaE**2 *spde$param.inla$M1 + spde$param.inla$M2)
    x_ <- as.vector(ar.matrix::sim.AR(n=1, Q))
    x <- x_ - mean(x_)
    N <- nrow(shapePointsDF@data)
}

