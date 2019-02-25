#' Simulate A continous 2D spatial field evolving over time 
#'
#' @description Calculates a continous underlying spatial field for the
#' Dominican Republic over time and across ages for probability of mortality
#' for one of 7 age groups. The final observed field is a linear additive model
#' as described in Wakefield 2017.
#' 
#' @param M Number of years
#' @param sigmaE spatial variance
#' @param rangeE spatial range
#' @param rhoAR Yearly autoregression
#' @param sigmaRW Variation of independent temporal component
#' @param ageB age specfifc effects
#' @param deltaUrban urban effect
#' @param offset mesh creation parameter
#' @param max.edge mesh creation parameter
#'
#' @return field object, contains multiple items. Spatial points data frame with 
#' raster values for transformed linear combination, and other values. A mesh 
#' that was used to create the latent field. The latent field
#' itself. A bounding shape where all observations take place. A projection 
#' matrix from the mesh to the entire field of interest. The spde for the matern
#' approximation.
#'
#' @examples
#' \dontrun{
#' require(tidyverse)
#'
#' fieldDR <- simField(4, rangeE = 1.3)
#'  
#' fieldDR$fqzDF %>%
#'     ggplot(aes(x=long, y=lat, fill=fqz)) +
#'     geom_raster() +
#'     facet_wrap(~year) +
#'     coord_equal() +
#'     theme_void() +
#'     scale_fill_distiller(palette = "Spectral")
#'
#' @export


rm(list=ls())

load("./Data/prepData.Rdata")

simField <- function(
    M,
    rangeE = .8,
    sigmaE = .04,
    rhoAR = .9,
    sigmaRW = .001,
    ageB = arm::logit((7:1 / 1000) + .00025),
    deltaUrban = -.1,
    offset = 1.,
    max.edge = c(.15,.3)){
    
    library(rgdal)
    library(sp)
    library(maptools)
    library(PointPolygon)
    library(rgeos)
    library(tidyverse)
    library(ar.matrix)
    
    meshDR <- INLA::inla.mesh.2d(
        loc=bbCoords(spDF), 
        offset = offset,
        max.edge = max.edge)
    
    AprojDR <- INLA::inla.spde.make.A(
        mesh=meshDR, 
        loc=as.matrix(fullDF[,c("long", "lat")]))
    
    kappaE <- sqrt(8) / rangeE
    tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)
    spde <- INLA::inla.spde2.matern(meshDR)
    Qspde <- tauE**2 * 
        (kappaE**4 * spde$param.inla$M0 +
             2 * kappaE**2 *spde$param.inla$M1 + spde$param.inla$M2)
    Qar1 <- Q.AR1(M, 1, rhoAR) 
    Q <- kronecker(Qar1, Qspde)
    x_ <- as.vector(ar.matrix::sim.AR(n=1, Q))
    x <- x_ - mean(x_)
    
    for(i in 1:M){
        fullDF[,paste0("x", i)] <- as.vector(
            AprojDR %*% x[(meshDR$n * (i-1) + 1):(meshDR$n * i)])
    }
    
    pAgeDF <- tibble(
        ageB = ageB,
        age  = 1:7)
    
    pTimeDF <- tibble(
        time = 1:M,
        rwTime = c(r.AR1(1, M, sigmaRW, .9999) %>% `-`(mean(.)))
    )
    
    idDF <- as_tibble(expand.grid(id=1:max(fullDF$id), time=1:M, age=1:7))
    
    fullDF <- fullDF %>%
        gather("time", "x", -long, -lat, -reg, -urban, -id, -strat) %>%
        mutate(time=as.integer(gsub("x", "", time))) %>%
        right_join(idDF, by=c("id", "time")) %>%
        left_join(pTimeDF, by="time") %>%
        left_join(pAgeDF, by="age") %>%
        mutate(deltaU=ifelse(urban == 1, deltaUrban, 0)) %>%
        mutate(p = arm::invlogit(ageB + deltaU + rwTime + x)) %>%
        mutate(year = time + 1999) %>%
        select(-time)
    
    fqzDF <- fullDF %>%
        group_by(long, lat, id, reg, urban, strat, year) %>%
        summarize(fqz=1-prod(1-p)) %>%
        ungroup
    
    boundShape <- rgeos::gUnaryUnion(spDF, id =  spDF@data$Bound)
    boundShape <- sp::SpatialPolygonsDataFrame(
        boundShape, data.frame(Bound = 1, row.names = "1"))
    
    field <- list(spdf = fullDF, mesh = meshDR, latent = x, fqzDF=fqzDF,
                  bound = boundShape, AprojField = AprojDR, spde = spde, 
                  deltaUrban, rangeE = rangeE, sigmaE = sigmaE, rhoAR = rhoAR,
                  sigmaRW = sigmaRW, ageB = ageB)
    class(field) <- "field"
    field
}
