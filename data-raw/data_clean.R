### Create data for use with DR simulation code.
# This script takes data from DHS 2007 & 2013 and MICS 2014 and combines them
# into a single data set. In addition we create as raster file with rural and
# urban distinctions as pulled from the 2010 census. 

rm(list=ls())
set.seed(123)
library(haven)
library(rgdal)
library(sp)
library(maptools)
library(PointPolygon)
library(rgeos)
library(raster)
library(tidyverse)

years <- 2000:2015
ageGroups <- c(0, NN=1/12, PNN1=1/2, PNN2=1, `1yr`=2, `2yr`=3, `3yr`=4, `4yr`=5)

# Modified resample function to get sum of small areas

smartResample <- function(x, y, method="bilinear", filename="", ...)  {
    
    # to do: compare projections of x and y
    
    ln <- names(x)
    nl <- nlayers(x)
    if (nl == 1) {
        y <- raster(y)
    } else {
        y <- brick(y, values=FALSE, nl=nl)
    }
    
    if (!hasValues(x)) {
        return(y)
    }	
    
    if (!method %in% c('bilinear', 'ngb')) {
        stop('invalid method') 
    }
    if (method == 'ngb') method <- 'simple'
    
    skipaggregate <- isTRUE(list(...)$skipaggregate)
    if (!skipaggregate) {
        rres <- res(y) / res(x)
        resdif <- max(rres)
        if (resdif > 2) {
            ag <- pmax(1, floor(rres-1))
            if (max(ag) > 1) {
                if (method == 'bilinear') {
                    x <- aggregate(x, ag, 'sum')
                } else {  
                    x <- aggregate(x, ag, modal)
                }
            }
        }
    }
    
    e <- raster:::.intersectExtent(x, y, validate=TRUE)
    
    filename <- trim(filename)
    if (canProcessInMemory(y, 4*nl)) {
        inMemory <- TRUE
        v <- matrix(NA, nrow=ncell(y), ncol=nlayers(x))
    } else {
        inMemory <- FALSE
        y <- writeStart(y, filename=filename, ... )
    }
    
    
    if (raster:::.doCluster()) {
        
        cl <- getCluster()
        on.exit( returnCluster() )
        
        nodes <- min(ceiling(y@nrows/10), length(cl)) # at least 10 rows per node
        
        message('Using cluster with ', nodes, ' nodes')
        utils::flush.console()
        
        tr <- blockSize(y, minblocks=nodes, n=nl*4*nodes)
        pb <- pbCreate(tr$n, label='resample', ...)
        
        clFun <- function(i) {
            #r <- tr$row[i]:(tr$row[i]+tr$nrows[i]-1)
            xy <- xyFromCell(y, cellFromRowCol(y, tr$row[i], 1) : cellFromRowCol(y, tr$row[i]+tr$nrows[i]-1, ncol(y)) ) 
            .xyValues(x, xy, method=method)
        }
        
        parallel::clusterExport(cl, c('x', 'y', 'tr', 'method'), envir=environment())
        .sendCall <- eval( parse( text="parallel:::sendCall") )
        for (ni in 1:nodes) {
            .sendCall(cl[[ni]], clFun, list(ni), tag=ni)
        }
        
        if (inMemory) {
            for (i in 1:tr$n) {
                d <- .recvOneData(cl)
                if (! d$value$success) {
                    stop('cluster error')
                }
                start <- cellFromRowCol(y, tr$row[d$value$tag], 1)
                end <- cellFromRowCol(y, tr$row[d$value$tag]+tr$nrows[d$value$tag]-1, y@ncols)
                v[start:end, ] <- d$value$value
                
                ni <- ni + 1
                if (ni <= tr$n) {
                    .sendCall(cl[[d$node]], clFun, list(ni), tag=ni)
                }
                pbStep(pb)
            }
            y <- setValues(y, v)
            if (filename != '') {
                writeRaster(y, filename, ...)
            }
            
        } else {
            
            for (i in 1:tr$n) {
                d <- .recvOneData(cl)
                y <- writeValues(y, d$value$value, tr$row[d$value$tag])
                ni <- ni + 1
                if (ni <= tr$n) {
                    .sendCall(cl[[d$node]], clFun, list(ni), tag=ni)
                }
                pbStep(pb)
            }
            y <- writeStop(y)	
        }	
        
    } else {
        
        tr <- blockSize(y, n=nl*4)
        pb <- pbCreate(tr$n, label='resample', ...)
        
        if (inMemory) {
            for (i in 1:tr$n) {
                #r <- tr$row[i]:(tr$row[i]+tr$nrows[i]-1)
                xy <- xyFromCell(y, cellFromRowCol(y, tr$row[i], 1) : cellFromRowCol(y, tr$row[i]+tr$nrows[i]-1, ncol(y)) ) 
                vals <- raster:::.xyValues(x, xy, method=method)
                
                start <- cellFromRowCol(y, tr$row[i], 1)
                end <- cellFromRowCol(y, tr$row[i]+tr$nrows[i]-1, y@ncols)
                v[start:end, ] <- vals
                
                pbStep(pb, i)
            }
            v <- setValues(y, v)
            if (filename != '') {
                writeRaster(v, filename, ...)
                
            }
            pbClose(pb)
            names(v) <- ln
            return(v)
            
        } else {
            for (i in 1:tr$n) {
                xy <- xyFromCell(y, cellFromRowCol(y, tr$row[i], 1) : cellFromRowCol(y, tr$row[i]+tr$nrows[i]-1, ncol(y)) ) 
                vals <- raster:::.xyValues(x, xy, method=method)
                
                y <- writeValues(y, vals, tr$row[i])
                
                pbStep(pb, i)
            }
            y <- writeStop(y)	
        }
    }
    
    pbClose(pb)
    names(y) <- ln
    return(y)
}


# dhsDF <- "./Data/DRBR52FL.SAV" %>%
#     read_spss

DF <- "./data-extended/bh.sav"  %>%
    read_spss() %>%
    rename(Region = HH7, PSU = HH1, `Birth(CMC)` = BH4C) %>%
    mutate(Urban = HH6=="1") %>%
    mutate(strat = paste0(sprintf("%02d", Region), "_", HH6)) %>%
    mutate(Alive = BH5 == "1") %>%
    mutate(`Age at Death` = c(1/365, 1/12, 1)[BH9U] * BH9N) %>%
    mutate(`Year of Birth` = BH4Y) %>%
    mutate(`Year at Death` = floor(`Age at Death`) + `Year of Birth`) %>%
    filter(`Year of Birth` <= 2018 & `Year of Birth` >= 1996) %>%
    mutate(`Age Group at Death` = cut(
        `Age at Death`, ageGroups, right=F, labels=names(ageGroups)[-1]))

# make an individual age group level analysis datset where each individual
# contributes one row per number of age groups that they survived to. Location
# Remains constant for the individual and time is increased according to
# year of birth and age.
micsDF <- bind_rows(lapply(1:nrow(DF), function(i){
    sigDF <- DF[i,]
    allLevels <- levels(DF$`Age Group at Death`)
    addLevels <- cumsum(as.numeric((ageGroups > 1)[-1]))
    names(addLevels) <- allLevels
    if(is.na(sigDF$`Age Group at Death`)){
        expDF <- tibble(
            id = i,
            died = 0,
            age_group = allLevels,
            strat = sigDF$strat,
            year = as.numeric(sigDF$`Year of Birth`) + unname(addLevels),
            psu = sigDF$PSU
        )
    }
    else{
        deathIDX <- which(allLevels == sigDF$`Age Group at Death`)
        subLevels <- allLevels[1:deathIDX]
        expDF <- tibble(
            id = i,
            died = c(rep(0, deathIDX-1), 1),
            age_group = subLevels,
            strat = sigDF$strat,
            year = as.numeric(sigDF$`Year of Birth`) + 
                unname(addLevels[subLevels]),
            psu = sigDF$PSU
        )
    }
    expDF})) %>%
    select(-id) %>%
    group_by(age_group, year, strat, psu) %>%
    summarize(N=n(), died=sum(died)) %>%
    mutate(lat=NA, long=NA, source="MICS_2014") %>%
    ungroup %>%
    filter(year >= min(years) & year <= 2014) %>%
    mutate(psu = as.character(psu)) %>%
    mutate(source = 
               "Dominican Republic Multiple Indicator Cluster Survey 2014") %>%
    mutate(urban = str_sub(strat, -1, -1) == "1") %>%
    filter(!(year == 2014 & (age_group %in% c("PNN1", "PNN2"))))

sources <- c(
    "Dominican Republic Demographic and Health Survey 2007",
    "Dominican Republic Demographic and Health Survey 2013",
    "Dominican Republic Multiple Indicator Cluster Survey 2014"
)

ihmeDF <- "./data-extended/DRdata.csv" %>%
    read_csv(col_types="ccdcidcddccddidcdc") %>%
    mutate(strat=NA) %>%
    select(
        age_group=ab, year, strat, psu=geo_unit, 
        N, died, lat=latnum, long=longnum, source=Title) %>%
    filter(source %in% sources)


# theres a huge difference here need to ask why that is
sum(ihmeDF$N[grepl("Cluster Survey", ihmeDF$source) & ihmeDF$year < 2015]) - 
    sum(micsDF[micsDF$year < 2015,]$N)

# load in the shape and use this to get the urban rural of points
spDF <- readOGR("./data-extended/SECCenso2010.dbf") %>%
    spTransform(CRS("+proj=longlat"))
spDF$strat <- paste0(spDF$REG, "_", spDF$ZONA)
spDF$urban <- spDF$ZONA == "1"

# because of the jitter we get some NA so lets recode those to the nearest poly
assignLoc <- function(locsDF, shape, ZONA=TRUE){
    locMatSPDF <- select(locsDF, long, lat) %>%
        as.matrix %>%
        SpatialPoints(CRS("+proj=longlat"))
    resultsDF <- locMatSPDF %>%
        over(spTransform(shape, CRS("+proj=longlat")))
    if(ZONA){
        results <- resultsDF$ZONA == "1"
    }
    else{
        results <- resultsDF$id
    }
    
    for(i in 1:length(results)){
        if(is.na(results[i])){
            m_ <- which.min(gDistance(
                locMatSPDF[i,], 
                spTransform(shape, CRS("+proj=longlat")), byid=T))
            if(ZONA){
                results[i] <- shape$ZONA[m_] == "1"
            }
            else{
                results[i] <- shape$id[m_]
            }
        }
    }
    results
}

if(!file.exists("./data-extended/dhsDF.Rds")){
    dhsDF <- ihmeDF %>%
        filter(!grepl("Cluster Survey", ihmeDF$source)) %>%
        mutate(urban=assignLoc(., spDF))
    
    saveRDS(dhsDF, "./data-extended/dhsDF.Rds")
}

dhsDF <- read_rds("./data-extended/dhsDF.Rds")

# the last step is to convert to a sufficiently detailed raster
# we will then convert this raster to point data as used in the PointPolygon
# package.
rbase <- raster(ncol=750, nrow=750)
extent(rbase) <- extent(spDF)
rasterRegSP <- rasterize(spDF, rbase, 'REG')
rasterUrbanSP <- rasterize(spDF, rbase, 'urban')

fullRaster <- rasterToPolygons(rasterRegSP)
fullRaster$reg <- fullRaster$layer
fullRasterUrban <- rasterToPolygons(rasterUrbanSP)
fullRaster$urban <- fullRasterUrban$layer
fullRaster$id <- 1:nrow(fullRaster@data)
fullRaster$ZONA <- if_else(fullRaster$urban==1, "1", "2")
fullRaster$strat <- paste0(
    sprintf("%02d", fullRaster$reg), "_", fullRaster$ZONA)

fullDF <- as_tibble(sp::coordinates(fullRaster)) %>%
    rename(long=V1, lat=V2) %>%
    bind_cols(select(fullRaster@data, reg, urban, id)) %>%
    mutate(strat = paste0(sprintf("%02d", reg), "_", 2-urban)) %>%
    mutate(long=round(long, 5), lat=round(lat, 5))

# now that we have assigned all of the points lets assign each of the points in
# the DHS to their corresponding id

if(!file.exists("./data-extended/dhsClusterDF.Rds")){
    dhsClusterDF <- dhsDF %>%
        select(psu, lat, long) %>%
        unique %>%
        mutate(id = assignLoc(., fullRaster, ZONA=FALSE))
    
    saveRDS(dhsClusterDF, "./data-extended/dhsClusterDF.Rds")
}

dhsClusterDF <- read_rds("./data-extended/dhsClusterDF.Rds") %>%
    select(psu, id)

idSubs <- mapply(function(i, j){
    suDF <- subset(fullRaster@data, strat == paste0(sprintf("%02d", i), "_", j))
    suDF$id
}, rep(1:10, each=2), rep(1:2, 10))

polyDF <- tibble(
    strat = mapply(
        function(i,j) paste0(sprintf("%02d", i), "_", j), 
        rep(1:10, each=2), rep(1:2, 10)),
    id = idSubs) %>%
    right_join(micsDF, by="strat") %>%
    mutate(point=FALSE)

pointDF <- dhsDF %>%
    left_join(dhsClusterDF, by="psu") %>%
    mutate(point=TRUE)

# Next we are going to want to build the population rasters

wgetRaster <- function(year, age, sex, loc="DOM"){
    baseURL <- "ftp://ftp.worldpop.org.uk/"
    projURL <- paste0(baseURL, "GIS/AgeSex_structures/Global_2000_2020/")
    specURL <- paste0(projURL, year, "/", toupper(loc), "/", tolower(loc), "_",
                      str_sub(tolower(sex), 1, 1), "_", age, "_", year, ".tif")
    wpopRaster <- smartResample(raster(specURL), rasterUrbanSP)
    wpopRaster
}

wgetRasterYear <- function(year, loc="DOM"){
    rasterList <- mapply(
        function(i, j) wgetRaster(year, i, j, loc),
        c(0, 1, 0, 1),
        c("Male", "Male", "Female", "Female"))
    Reduce("+", rasterList)
}


if(!file.exists("./data-extended/rasterYearList.Rds")){
    rasterYearList <- lapply(years, wgetRasterYear)
    saveRDS(rasterYearList, "./data-extended/rasterYearList.Rds")
}

rasterYearList <- read_rds("./data-extended/rasterYearList.Rds")
names(rasterYearList) <- years

popYearDFList <- lapply(years, function(x) NULL)
names(popYearDFList) <- as.character(years)

# the number of filled points dont match the urban rural rasters so this isnt
# going to work, We need to change this step
for(y in years){
    rDF <- rasterToPolygons(rasterYearList[[as.character(y)]])
    rDF$Population <- rDF$layer
    rDF$layer <- NULL
    rDF$id <- 1:nrow(rDF@data)
    rDF@data <- left_join(rDF@data, select(fullRaster@data, id, strat)) %>%
        group_by(strat) %>%
        mutate(popW=Population/sum(Population)) %>%
        ungroup %>%
        mutate(year=y)
    popYearDFList[[as.character(y)]] <- rDF
    rm(list=c("rDF"))
}

popYearDFList <- lapply(years, function(y){
    rDF <- rasterToPolygons(rasterYearList[[as.character(y)]])
    rDF$Population <- rDF$layer
    rDF$layer <- NULL
    rDF$id <- 1:nrow(rDF@data)
    rDF@data <- left_join(rDF@data, select(fullRaster@data, id, strat)) %>%
        group_by(strat) %>%
        mutate(popW=Population/sum(Population)) %>%
        ungroup %>%
        mutate(year=y)
    rDF
})

yearWDF <- bind_rows(lapply(popYearDFList, function(x)x@data)) %>% 
    arrange(year, id)

yearWMat <- matrix(
    c(yearWDF$popW), 
    nrow=nrow(popYearDFList[[1]]@data), 
    ncol = length(years))

save(yearWMat, yearWDF, polyDF, pointDF, spDF, fullDF,
     file="./Data/prepData.Rdata")
