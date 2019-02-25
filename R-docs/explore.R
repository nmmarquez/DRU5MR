rm(list=ls())
set.seed(123)
library(tidyverse)
library(haven)
library(rgdal)
library(sp)
library(maptools)
library(PointPolygon)

ageGroups <- c(0, NN=1/12, PNN1=1/2, PNN2=1, `1yr`=2, `2yr`=3, `3yr`=4, `4yr`=5)

dhsDF <- "./Data/DRBR52FL.SAV" %>%
    read_spss

DF <- "~/Documents/DRExplore/Data/bh.sav"  %>%
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
aggDF <- bind_rows(lapply(1:nrow(DF), function(i){
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
    mutate(lat, lon)


# we are going to need to restructure this data so it works for modeling
# we just need to sample a number of individuals equal to whats observed in the
# data and not worry about the outcome for now.
DF %>%
    filter(`Age at Death` <= 5) %>%
    ggplot(aes(x=`Age at Death`)) +
    geom_density() +
    theme_classic()

spDF <- readOGR("~/Documents/DRExplore/Data/SECCenso2010.dbf")
plot(spDF[spDF$ZONA == "2",], main="Rural Areas")
plot(spDF[spDF$ZONA == "1",], main="Urban Areas") 

# plot(spDF[spDF$REG == "01",], main="01") # Norte
# plot(spDF[spDF$REG == "02",], main="02") # Sur
# plot(spDF[spDF$REG == "03",], main="03") # Nordeste
# plot(spDF[spDF$REG == "04",], main="04") # Noroeste
# plot(spDF[spDF$REG == "05",], main="05") # Valdesia
# plot(spDF[spDF$REG == "06",], main="06") # Enriquillo
# plot(spDF[spDF$REG == "07",], main="07") # El Valle
# plot(spDF[spDF$REG == "08",], main="08") # Yuma
# plot(spDF[spDF$REG == "09",], main="09") # Higuamo
# plot(spDF[spDF$REG == "10",], main="10") # Ozama(Metro)

spDF$strat <- paste0(spDF$REG, "_", spDF$ZONA)

# Somme number checking to the offical documentation looks good
# Recall that the PSUs are the cluster from which data was collected
DF %>%
    select(PSU, Urban) %>%
    unique %>%
    group_by(Urban) %>%
    summarize(n())

DF %>%
    select(PSU, Urban, Region) %>%
    unique %>%
    group_by(Urban, Region) %>%
    summarize(n())

# We want to create large shape files that are urban rural specific.
stratSPDF <- unionSpatialPolygons(spDF, spDF$strat)
shape <- SpatialPolygonsDataFrame(rgeos::gUnaryUnion(spDF), data.frame(ID="1"))

N <- 500

# Example simulation 

ggField(drSim <- simField(sigmaE=.6,
    N=N, shape=shape, beta0=-4,
    offset = c(10000, 20000), 
    max.edge = c(5000,20000),
    betaList = list(),
    rangeE = 120000))

drSim$spdf$strat <- over(drSim$spdf, spDF)$strat

# Simulate the convolution as if we didnt know where the cluster was sampled
M <- 100000
Size <- 36

simDF <- bind_rows(lapply(sort(unique(drSim$spdf$strat)), function(s){
    drSim$spdf@data %>%
        filter(strat == s) %>%
        select(theta) %>%
        unlist %>%
        sample(M, replace=TRUE) %>%
        rbinom(prob=., n=M, size=Size) %>%
        {sapply(0:Size, function(i) sum(i == .))} %>%
        tibble(obs=0:Size, count=., strat=s, prob=./M)}))

simDF %>%
    ggplot(aes(x=obs, y=prob)) +
    geom_line() + 
    facet_wrap(~strat) +
    theme_classic()
