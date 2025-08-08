library(terra)
library(tidyverse)
library(tidyterra)
library(lme4)
library(MuMIn)
library(ROCR)
library(car)
library(zoo)
library(tidycensus)

AOI_maui.shp <- vect(file.path('Data', 'AOI_maui.gpkg')) %>% mutate(AOI = 'Maui')
AOI_hawaii.shp <- vect(file.path('Data', 'AOI_hawaii_nobuffer_albers.gpkg')) %>% mutate(AOI = 'Hawaii')

AOI_combined.shp <- rbind(AOI_maui.shp, AOI_hawaii.shp)

GHS.stack <- rast(dir(file.path('Data', 'Landcover'), full.names = T)[grepl('GHS_\\d+\\.tif$', dir(file.path('Data', 'Landcover')))]) %>% project(crs(AOI_combined.shp)) %>% crop(ext(AOI_combined.shp))
plot(GHS.stack)
writeRaster(GHS.stack, file.path('Data', 'Landcover', 'GHS_maui_hawaii.tif'), overwrite=T)

# https://guides.library.manoa.hawaii.edu/c.php?g=105181&p=684171
# https://fred.stlouisfed.org/series/HIMAUI5POP
# https://files.hawaii.gov/dbedt/economic/data_reports/LRF/2050-long-range-projections.pdf

# Hpop.df <- read.csv(file.path('Inputs', 'pop_trend.csv')) %>%
#   mutate(Defacto_pop = X15462+X15461	+X15460	+X15459	+X15458	+X15457	+X15456	+X15455	+X15454	+X15453	+X15452	+X15451	+X15450,
#          AOI = 'Hawaii')
# Mpop.df <- read.csv(file.path('Inputs', 'Maui_historical_pop.csv')) %>%
#   mutate(AOI = 'Maui')
# 
# Pop.df <- bind_rows(Hpop.df, Mpop.df) %>%
#   select(year, AOI, Defacto_pop)

Pop.df <- read.csv(file.path('Inputs', 'Maui_historical_pop.csv')) %>%
  pivot_longer(cols = c(Maui_residentpop, Hawaii_residentpop, Maui_defactopop), names_to = c('AOI', 'poptype'), names_sep = '_', values_to = 'Population') %>% 
  filter(poptype == 'residentpop')

AOI_area.df <- data.frame("AOI" = AOI_combined.shp$AOI,
                          "AOI_area" = expanse(AOI_combined.shp, unit='km'))

class_cutoff <- 100
GHS.df <- ifel(GHS.stack>class_cutoff, 1, 0) %>% 
  zonal(AOI_combined.shp, fun='sum') %>%
  mutate(AOI = AOI_combined.shp$AOI) %>%
  pivot_longer(cols = starts_with('GHS'), names_to ='year', names_pattern = 'GHS_(....)', values_to = 'Dev_px') %>%
  mutate(year = as.numeric(year)) %>%
  full_join(Pop.df) %>%
  left_join(AOI_area.df) %>%
  mutate(Density_pop_px = Population/Dev_px,
         Density_px_area = Dev_px/AOI_area,
         Density_pop_area = Population/AOI_area)

ggplot(data=GHS.df, aes(x=year, y=Density_pop_px, colour=AOI)) +
  geom_line()

p1 <- ggplot(data=GHS.df, aes(x=year, y=Density_pop_area, colour=AOI)) +
  geom_line() + ylab("Population per sq. km.") + ggtitle('Population density: Maui Vs. Hawaii County')+ theme(legend.position="none") + expand_limits(y = 0)

p2 <-ggplot(data=GHS.df, aes(x=year, y=Population, colour=AOI)) +
  geom_line(na.rm=T) + ylab("Population") + ggtitle('Population: Maui Vs. Hawaii County')  + expand_limits(y = 0)

p3 <-ggplot(data=GHS.df, aes(x=year, y=Dev_px, colour=AOI)) +
  geom_line(na.rm=T) + ylab("Developed pixels (100m)") + ggtitle('Developed area: Maui Vs. Hawaii County') + theme(legend.position="none")  + expand_limits(y = 0)

p4 <-ggplot(data=GHS.df, aes(x=year, y=(((Dev_px/100) / AOI_area) * 100), colour=AOI)) +
  geom_line(na.rm=T) + ylab("% Land area developed") + ggtitle('% Developed: Maui Vs. Hawaii County') + theme(legend.position="none")  + expand_limits(y = 0)

library(gridExtra)
grid.arrange(p1, p2, p3, p4)


Maui_stack <- ifel(GHS.stack>class_cutoff, 1, 0) %>%
  crop(AOI_maui.shp, mask = T)

plot(Maui_stack)

r <- ifel(!is.na(Maui_stack[[1]]), 99999, NA)
for (lyr in names(Maui_stack)){
  yr = as.numeric(str_extract(lyr, '\\d\\d\\d\\d'))
  print(yr)
  yrdev <- ifel(Maui_stack[lyr] == 1, Maui_stack[lyr] * yr, 99999)
  r <- min(yrdev, r)
}
plot(r)

mask.r <- ifel(r==99999, 1, NA)
r[r==99999] <- NA

plot(mask.r, col='gray', legend=F)
plot(r, add=T)




