rm(list=ls())

library(DRU5MR)
library(tidyverse)

fieldDR <- DRU5MR::simField(4, rangeE = 1.3)
fieldDR$fqzDF <- left_join(
    fieldDR$fqzDF, 
    select(yearWDF, id, popW, year), 
    by=c("id", "year"))

fieldDR$fqzDF %>%
    ggplot(aes(x=long, y=lat, fill=fqz)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

subDF <- fieldDR$fqzDF %>% filter(strat == "05_0" & year == 2000)
M <- 30

plot(apply(sapply(1:nrow(subDF), function(i){
    dbinom(0:M, size=M, subDF$fqz[i]) * subDF$popW[i]
}), 1, sum), type="l")

fieldDR <- DRU5MR::simField(4, rangeE = .3)
fieldDR$fqzDF <- left_join(
    fieldDR$fqzDF, 
    select(yearWDF, id, popW, year), 
    by=c("id", "year"))

fieldDR$fqzDF %>%
    ggplot(aes(x=long, y=lat, fill=fqz)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

subDF <- fieldDR$fqzDF %>% filter(strat == "10_1" & year == 2000)

plot(apply(sapply(1:nrow(subDF), function(i){
    dbinom(0:M, size=M, subDF$fqz[i]) * subDF$popW[i]
}), 1, sum), type="l")
