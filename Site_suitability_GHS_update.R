library(terra)
library(tidyverse)
library(tidyterra)
library(lme4)
library(MuMIn)
library(ROCR)
library(car)
library(zoo)
#library(tidycensus)
#census_api_key("95156035001c2244388003f163072c904b153c94", install=T)

AOI <- vect(file.path('Data', 'AOI_albers_buffer.shp'))

#plot(LF01.r); plot(AOI, col='red', add=T)


#LF14.r <- rast(file.path('Data', 'Landcover', 'HI_140_EVC', 'Tif', 'hi_140evc.tif')) |> crop(LF01.r) |> resample(LF01.r, method='near') |> extend(LF01.r)
#LF16.r <- rast(file.path('Data', 'Landcover', 'LF2016_EVC_200_HI', 'Tif', 'LH16_EVC_200.tif')) |> crop(LF01.r) |> resample(LF01.r, method='near') |> extend(LF01.r)
#LF22.r <- rast(file.path('Data', 'Landcover', 'LF2022_EVC_230_HI', 'Tif', 'LH22_EVC_230.tif')) |> crop(LF01.r) |> resample(LF01.r, method='near') |> extend(LF01.r)

GHS80.r <- rast(file.path('Data', 'Landcover', 'GHS_1980.tif')) |> project(crs(AOI)) |> crop(ext(AOI))
GHS90.r <- rast(file.path('Data', 'Landcover', 'GHS_1990.tif')) |> project(crs(AOI)) |> crop(ext(AOI))
GHS00.r <- rast(file.path('Data', 'Landcover', 'GHS_2000.tif')) |> project(crs(AOI)) |> crop(ext(AOI))
GHS10.r <- rast(file.path('Data', 'Landcover', 'GHS_2010.tif')) |> project(crs(AOI)) |> crop(ext(AOI))
GHS15.r <- rast(file.path('Data', 'Landcover', 'GHS_2015.tif')) |> project(crs(AOI)) |> crop(ext(AOI))
GHS20.r <- rast(file.path('Data', 'Landcover', 'GHS_2020.tif')) |> project(crs(AOI)) |> crop(ext(AOI))

LF01.r <- rast(file.path('Data', 'Landcover', 'HI_105_EVC', 'Tif', 'hi_105evc.tif')) |> resample(GHS80.r, method='near') |> crop(ext(AOI))

################################################################################
########### Model osds to pick a cutoff point? #################################
################################################################################

### Look at urban change
# U01.r <- ifel(LF01.r %in% c(21,22,23,24), 1, 0)
# 
# ### check relationship between building area and landfire class
# LFU_resamp.r <- resample(U01.r, GHS00.r, method='near')
# stack_test <- c(GHS00.r, LFU_resamp.r) |> as.data.frame()
# ggplot(data=stack_test, aes(x=GHS_2000, fill=as.factor(VALUE_1))) + geom_density() + xlim(0,200)
# 
# U14.r <- ifel(LF14.r %in% c(21,22,23,24), 1, U01.r)
# U16.r <- ifel(LF16.r %in% c(21,22,23,24), 1, 0)
# U22.r <- ifel(LF22.r %in% c(21,22,23,24), 1, 0)
class_cutoff <- 100  # 10% of cell

U80.r <- as.int(ifel(GHS80.r>class_cutoff, 1, 0)) 
U90.r <- as.int(ifel(GHS90.r>class_cutoff, 1, 0))
U00.r <- as.int(ifel(GHS00.r>class_cutoff, 1, 0))
U10.r <- as.int(ifel(GHS10.r>class_cutoff, 1, 0))
U15.r <- as.int(ifel(GHS15.r>class_cutoff, 1, 0))
U20.r <- as.int(ifel(GHS20.r>class_cutoff, 1, 0))

UC80_20.r <- ifel(U80.r==0&U20.r==1, 1, 0)
UC00_20.r <- ifel(U00.r==0&U20.r==1, 1, 0)
#UC16_22.r <- ifel(U16.r==0&U22.r==1, 1, 0)

plot(UC80_20.r)
plot(UC00_20.r)
#plot(UC14_22.r)
#plot(UC16_22.r)  ## No change!!

sum(U00.r[U00.r == 1])  #517044
sum(U15.r[U15.r == 1])  #477737
sum(U20.r[U20.r == 1])  #570629

#writeRaster(U14.r, file.path('Inputs', 'U14.tif'), overwrite = T)

predictor_stack <- c(LF01.r, U80.r, U90.r, U00.r, U10.r, U15.r, U20.r, UC00_20.r)
names(predictor_stack) <- c('LF01', 'U80', 'U90', 'U00', 'U10', 'U15', 'U20', 'UC00_20')

### Restrictions 
restricted.shp <- vect(file.path('Data', 'Restrictions', 'zones_cantbedeveloped.shp'))
test_zone.r <- restricted.shp |> filter(ZONE == "OPEN") |> project(crs(LF01.r)) |> rasterize(LF01.r, background=0)
names(test_zone.r) <- "Open_restricted"


restricted.r <- restricted.shp |> 
  filter(ZONE != "OPEN") |>  # filter out "Open" zone since it's now a predictor
  project(crs(LF01.r)) |> rasterize(LF01.r, background=0)
restricted.r[LF01.r <=11] <- 1  # restrict water

writeRaster(restricted.r, file.path('Inputs', 'restrictions.tif'), overwrite=T)

roads.r <- vect(file.path('Data', 'Predictors', 'hawaii_OSM', 'gis_osm_roads_free_1.shp')) |> project(crs(LF01.r)) |> crop(LF01.r) |> rasterize(LF01.r, background=0)

Region_rcl.df <- data.frame(
  'from' = c(15450, 15451, 15452, 15453, 15454, 15455, 15456, 15457, 15458, 15459, 15460, 15461, 15462),
  'to' = c(15460, 15460, 15453, 15453, 15453, 15453, 15457, 15457, 15457, 15460, 15460, 15460, 15460)  # 15455 merge with 15453, 15459 and 15461 to 15460 
)
Regions.r <- vect(file.path('Data', 'FAZ.gpkg')) |>
  rasterize(LF01.r, field='FAZID') |> as.int() |>
  classify(Region_rcl.df)
names(Regions.r) <- 'Regions'
predictor_stack <- c(predictor_stack, Regions.r, test_zone.r)

### Load predictors
## Physical
# elevation
if(!file.exists(file.path('Data', 'Predictors', 'elevation.tif'))){
  library(elevatr)
  elev <- elevatr::get_elev_raster(UC00_20.r, z=12)
  elev.r <- rast(elev)
  elev.r <- elev.r |> terra::resample(LF01.r, method='bilinear')
  writeRaster(elev.r, file.path('Data', 'Predictors', 'elevation.tif'), overwrite=T)
}
elev.r <- rast(file.path('Data', 'Predictors', 'elevation.tif')); names(elev.r) <- 'elev'
slope.r <- terrain(elev.r, 'slope'); names(slope.r) <- 'slope'
predictor_stack <- c(predictor_stack, scale(elev.r), scale(slope.r))

# canopy cover
cc01.r <- rast(file.path('Data', 'Landcover', 'HI_105_CC', 'Tif', 'hi_105cc.tif'))  |> crop(LF01.r) |> resample(LF01.r, method='near') |> extend(LF01.r)
cc01.r[cc01.r < 0] <- NA
cc01_sm.r <- cc01.r |> focal(w=5, fun='mean'); names(cc01_sm.r) <- 'cc01_sm'
predictor_stack <- c(predictor_stack, scale(cc01_sm.r))

# lava flow hazard
lava_zones.r <- vect(file.path('Data', 'Predictors', 'Volcano_Lava_Flow_Hazard_Zones', 'Volcano_Lava_Flow_Hazard_Zones.shp')) |> 
  project(crs(LF01.r)) |>
  mutate(hazard = as.numeric(hzone)) |>
  rasterize(LF01.r, field='hazard', background=NA) |>
  focal(w=7, fun='modal', na.policy='only')  # fill in missing cells with a majority filter
lava_zones.r <-   ifel(is.na(lava_zones.r)&!is.na(Regions.r), 9, lava_zones.r)
names(lava_zones.r) <- 'lava_zones'

predictor_stack <- c(predictor_stack, lava_zones.r)

## amenities
# distance to ocean
ocean_dist_km.r <- ifel(LF01.r<=7, 1, 0) |> focal(w=11, fun='modal') |> terra::distance(target = 0, unit='km'); names(ocean_dist_km.r) <- 'ocean_dist_km'
predictor_stack <- c(predictor_stack, scale(ocean_dist_km.r))

# Distance to public access beaches

# Distance to schools

# Distance to parks and nature reserves
parks.r <- vect(file.path('Data', 'Predictors', 'parks_county_haw', 'parks_county_haw.shp')) |> 
  project(crs(LF01.r)) |> 
  rasterize(rast(LF01.r), background=0)
reserves.r <- vect(file.path('Data', 'Predictors', 'parks_state', 'parks_state.shp')) |> 
  project(crs(LF01.r)) |> 
  rasterize(rast(LF01.r), background=0)
park_dist_km.r <- max(parks.r, reserves.r)  |> terra::distance(target = 0, unit='km'); names(park_dist_km.r) <- 'park_dist_km'
predictor_stack <- c(predictor_stack, scale(park_dist_km.r))

## Infrastructure:
# Major road distance
roads_dist_km.r <- vect(file.path('Data', 'Predictors', 'centerlines_haw_maj', 'centerlines_haw_maj.shp')) |> 
  project(crs(LF01.r)) |> 
  rasterize(rast(LF01.r), background=0) |> 
  terra::distance(target = 0, unit='km')
names(roads_dist_km.r) <- 'roads_dist_km'
predictor_stack <- c(predictor_stack, scale(roads_dist_km.r))

# Road network density

# Drive time to airport?
### cost distance needs time taken to traverse a cell, so 1/speed
friction.r <- vect(file.path('Data', 'Predictors', 'hawaii_OSM', 'gis_osm_roads_free_1.shp')) |> project(crs(LF01.r)) |> crop(LF01.r) |> 
  mutate(maxspeed = ifelse(maxspeed==0, 30, maxspeed)) |> 
  mutate(mins_per_px = 3.73/maxspeed) |>
  rasterize(LF01.r, field='mins_per_px', background=1)
airports.r <- vect(file.path('Data', 'Predictors', 'hawaii_OSM', 'gis_osm_transport_a_free_1.shp')) |> project(crs(LF01.r)) |> crop(LF01.r) |> 
  filter(osm_id %in% c(145229278, 146321014)) |> rasterize(LF01.r, background=0)

friction.r[airports.r==1] <- 0

airport_cost_dist.r <- costDist(friction.r, maxiter=200)
airport_cost_dist.r[restricted.r==1] <- NA
plot(airport_cost_dist.r)
names(airport_cost_dist.r) <- 'airport_cost_dist'
predictor_stack <- c(predictor_stack, scale(airport_cost_dist.r))

# Drive time to population centers (Hilo and Kona?)

## Other
# landowners (categorical?)
landowners_conservation.r <- vect(file.path('Data', 'Restrictions', 'landowners_conservationmissions.shp')) |> project(crs(LF01.r)) |> rasterize(LF01.r, background=0) |> as.factor()
names(landowners_conservation.r) <- 'landowners_conservation'
predictor_stack <- c(predictor_stack, landowners_conservation.r)

## Development pressure:
focalmat_grav <- function(size, gamma, scale){
  m <- matrix(nrow=2*size + 1, ncol = 2*size + 1)  #initialize matrix
  center <- size + 1
  for(i in seq(1, 2*size + 1)){  #loop through the rows and columns of the matrix
    for(j in seq(1, 2*size + 1)){
      dist <- sqrt((i - center) * (i - center) + (j - center) * (j - center))  # assign each cell the value of its distance from the center
      if(dist <= size){  # assign distance only to cells in a circle around center
        m[i,j] <- dist
      }
      if(dist > size){  # assign zeroes everywhere else
        m[i,j] <- 0 
      }
    }
  }
  denom <- m^gamma  # raise matrix to power of gamma
  m <- scale/denom  # divide scaling factor by matrix to produce final weights
  m[denom == 0] <- 0  # replace INF value with zero
  return(m)
}

#01 devpr
devpr_sm <- focal(U00.r, w=focalmat_grav(3, 0.5, 1), na.rm=TRUE); names(devpr_sm) <- 'devpr_sm'
devpr_sm2 <- focal(U00.r, w=focalmat_grav(3, 1, 1), na.rm=TRUE); names(devpr_sm2) <- 'devpr_sm2'
devpr_sm3 <- focal(U00.r, w=focalmat_grav(3, 1.5, 1), na.rm=TRUE); names(devpr_sm3) <- 'devpr_sm3'
# devpr_big <- focal(U01.r, w=focalmat_grav(30, 0.5, 1), na.rm=TRUE); names(devpr_big) <- 'devpr_big'
devpr_big <- focal(U00.r, w=focalmat_grav(5, 0.5, 1), na.rm=TRUE); names(devpr_big) <- 'devpr_big'
devpr_big2 <- focal(U00.r, w=focalmat_grav(5, 1, 1), na.rm=TRUE); names(devpr_big2) <- 'devpr_big2'
devpr_big3 <- focal(U00.r, w=focalmat_grav(5, 1.5, 1), na.rm=TRUE); names(devpr_big3) <- 'devpr_big3'

plot(devpr_sm)
plot(devpr_big)

predictor_stack <- c(predictor_stack, scale(devpr_sm), scale(devpr_sm2), scale(devpr_sm3), scale(devpr_big), scale(devpr_big2), scale(devpr_big3))

### Save GRASS input files
writeRaster(predictor_stack, file.path('Data', 'predictor_stack.tif'), overwrite = T)
#predictor_stack <- rast(file.path('Data', 'predictor_stack.tif'))

### Do the model
# mask out non-developable cells
predictor_stack[restricted.r==1|is.na(Regions.r)] <- NA  # exclude restricted cells from modeling
#predictor_stack[U00.r == 1] <- NA  # exclude cells that were already urban from modeling
#predictor_stack[roads.r == 1] <- NA
gc()

# data frame
predictor_stack_train <- predictor_stack; predictor_stack_train[restricted.r==1] <- NA
predictor_stack.df <- as.data.frame(predictor_stack_train)
gc()
# 
# predictor_stack_pca.df <- predictor_stack.df |>
#   select(!c(LF01, U80, U90, U00, U10, U15, U20, UC00_20, devpr_sm, devpr_sm2, devpr_big, devpr_big2, devpr_big3, landowners_conservation, Regions))
# 
# pspca <- princomp(scale(predictor_stack_pca.df))
# summary(pspca)
# loadings(pspca)
# 
# predictor_stack_pca <- predictor_stack |>
#   select(!c(LF01, U80, U90, U00, U10, U15, U20, UC00_20, devpr_sm, devpr_sm2, devpr_big, devpr_big2, devpr_big3, landowners_conservation, Regions))
# pspca_rast <- predict(predictor_stack_pca, pspca)
# plot(pspca_rast)
# 
# library(ggcorrplot)
# 
# ggcorrplot(cor(scale(predictor_stack_pca.df)),  hc.order = TRUE,
#            type = "lower", insig = "blank",
#            lab = TRUE,
#            ggtheme = ggplot2::theme_gray)
# 
# 
# 
# library(caret)
# cors <- caret::findCorrelation(cutoff = 0.7, verbose = T, names = T, cor(scale(predictor_stack_pca.df))); cors

# training and testing sets: all the 1s, and half the zeros
n1s <- nrow(filter(predictor_stack.df, UC00_20 == 1))

sample_df <- bind_rows(
  predictor_stack.df |> filter(UC00_20==1),
  predictor_stack.df |> filter(UC00_20==0) |> group_by(Regions) |> slice_sample(n = 3*n1s) |> ungroup()
) |> mutate(pixel_ID = row_number())

train <- sample_df |> group_by(UC00_20, Regions) |> slice_sample(prop = 0.5) |> ungroup()
test <- sample_df |> anti_join(train)


# kitchen_sink model
names(train)
print(train |> group_by(Regions, UC00_20) |> tally(), n=100)
# ks_mod <- glm(UC01_22 ~ elev + slope + cc01_sm + ocean_dist_km + park_dist_km + roads_dist_km + airport_cost_dist + devpr_sm+devpr_sm2+devpr_sm3+devpr_big+devpr_big2+devpr_big3 + landowners_conservation + lava_zones,
#               data = train,

#               family = binomial(link="logit"),
#               na.action = 'na.fail')
# ks_mod <- glm(UC00_20 ~ elev + slope + cc01_sm + ocean_dist_km + park_dist_km + roads_dist_km + airport_cost_dist + devpr_sm+devpr_sm2+devpr_sm3+devpr_big+devpr_big2+devpr_big3 + landowners_conservation,
#               data = train,
#               family = binomial(link="logit"),
#               na.action = 'na.fail')

ks_mod <- glmer(UC00_20 ~ elev + slope + cc01_sm + ocean_dist_km + park_dist_km + roads_dist_km + airport_cost_dist +  devpr_sm+devpr_sm3+devpr_big+devpr_big3 + landowners_conservation + lava_zones + Open_restricted + (1|Regions),
              data = train,
              family = binomial(link="logit"),
              na.action = 'na.fail')

summary(ks_mod)
vif(ks_mod)


# Set up the cluster
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))

clusterEvalQ(clust, library(lme4))
clusterExport(clust, c("train", 'ks_mod'))

# mod_select_filtered <- dredge(ks_mod, beta = "none", trace = 2,
#                               subset = sum(devpr_sm, devpr_sm2, devpr_sm3, devpr_big, devpr_big2, devpr_big3)==1)
mod_select_filtered <- dredge(ks_mod, beta = "none", trace = 2,
                              subset = sum(devpr_sm, devpr_sm3, devpr_big, devpr_big3)==1, cluster=clust)

print(mod_select_filtered[1,])

top_mod <- glmer(UC00_20 ~ elev + slope + cc01_sm + ocean_dist_km + park_dist_km + roads_dist_km + airport_cost_dist + devpr_big + landowners_conservation + lava_zones + Open_restricted + (1|Regions), 
               data = train,
               nAGQ = 10,
               family = binomial(link="logit"),
               na.action = 'na.fail')

#test <- glm(UC00_20 ~ elev, data=train, family=binomial(link="logit"), na.action='na.fail')

summary(top_mod)
vif(top_mod)
top_devpr <- 'devpr_big'

prediction_test <- predict(top_mod, test)

pred <- prediction(as.numeric(prediction_test), as.numeric(test$UC00_20))
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

# table(observed=test_labels_newUrb, predicted=prediction_test)
auc.perf <- performance(pred, measure="auc")
print(auc.perf@y.values)

#pred_map <- terra::predict(predictor_stack, top_mod)
#plot(pred_map, background='black')


#writeRaster(pred_map, file.path('Outputs', 'Site_suitability.tif'), overwrite=T)

if('glm'%in%coef(top_mod)){  # if a regular glm was used
  coefs <- as.data.frame(t(as.data.frame(coef(top_mod))))
  outdf <- data.frame('ID'=1, 'Intercept'=coefs[,'(Intercept)'], 'devpr_20'=coefs[,top_devpr], coefs |> select(!c(top_devpr, '(Intercept)')))
} else {  # if there is a random effects term
  coefs <- as.data.frame(coef(top_mod)[[1]])
  outdf <- Region_rcl.df |> left_join(data.frame('ID'=as.numeric(rownames(coefs)), 'Intercept'=coefs[,'(Intercept)'], 'devpr_20'=coefs[,top_devpr], coefs |> select(!c(top_devpr, '(Intercept)'))), by=c('to'='ID')) |>
  select(!to) |>
  rename('ID'='from')
}

names(outdf)[grepl('1$', names(outdf))] <- str_replace(names(outdf)[grepl('1$', names(outdf))], "1", "")
write_csv(outdf, file.path('Inputs', 'Hawaii_model.csv'))


#write.table(cbind(rownames(coefs), coefs), file.path('Inputs', 'Hawaii_model.csv'), row.names=FALSE, sep="\t")


# 22 devpr
devpr_20 <- focal(U20.r, w=focalmat_grav(5, 0.5, 1), na.rm=TRUE)
writeRaster(devpr_20, file.path('Inputs', 'devpr_20.tif'), overwrite=T)


# Interpolate population numbers
pop.df <- read.csv(file.path('Data', 'Hawaii_pop_prj_FAZ.csv'))


pop_interp.df <- data.frame('year' = 1980:2045) |>
  left_join(pop.df) |>
  select(!Total) |>
  mutate(across(`X15462`:`X15450`, ~ round(na.approx(.x),0))) 

# Project population numbers using average rate from 2000-2020
yearly_rate <- (pop_interp.df[pop_interp.df$year == 2045,] - pop_interp.df[pop_interp.df$year == 2025,])/20
yearly_rate[yearly_rate<0]<-0
yearly_rate <- yearly_rate * 1; yearly_rate[1] <- 1

for(ext_yr in 2046:2100){
  pop_interp.df[ext_yr-1979, ] <- pop_interp.df[ext_yr-1980, ] + yearly_rate
}
                            
ggplot(data = pop_interp.df |> pivot_longer(cols=`X15462`:`X15450`, names_to = 'FAZ', values_to = 'Population'), aes(x=year, y=Population, colour = FAZ)) + geom_line() + geom_vline(xintercept=c(2020, 2045))

# |>
#   pivot_longer(cols=`X15462`:`X15450`, names_to = 'Region', values_to = 'Defacto_pop') |>
#   mutate(Region = str_replace(Region, 'X', '')) |>
#   left_join(Region_rcl.df)
#   #

#ggplot(data=pop_interp.df, aes(x=year)) + geom_line(aes(y=Pop), color='goldenrod') + geom_line(aes(y=DefactoPop), color='cornflowerblue') + ylim(0, 340000)

# pop_input <- data.frame('year'=pop_interp.df$year, '1' = pop_interp.df$DefactoPop)
# colnames(pop_input)[2]<-1
pop_input <- pop_interp.df
names(pop_input) <- str_replace(names(pop_interp.df), 'X', '')
pop_input_trend <- pop_input |> filter(year %in% c(1980, 1990, 2000, 2010, 2015, 2020))
pop_input_projection <- pop_input |> filter(year > 2020)

write_csv(pop_input_trend, file.path('Inputs', 'pop_trend.csv'))
write_csv(pop_input_projection, file.path('Inputs', 'pop_projection.csv'))




###############################################
future.df <- data.frame(year = 2046:2100,
                        County = 17690+(1:55 * 224))

df <- pop_interp.df |>
  select(year, X15457) |>
  rename(County = X15457) |>
  bind_rows(future.df)

historical.df <- pop_interp.df |> 
  filter(year <= 2021) |> 
  select(year, X15457) |>
  rename(Population = X15457) |>
  mutate(Source = 'Historical')

county.df <- pop_interp.df |> 
  filter(year >= 2020) |> 
  select(year, X15457) |>
  rename(Population = X15457) |>
  bind_rows(data.frame(
    year = 2046:2100,
    Population = 17690 + (1:55 * 224)
  )) |>
  mutate(Source = 'County Planning')

buildout.df <- data.frame(year = 2020:2100,
                          Population = 12087 + (0:80 * 2600),
                          Source = 'Build-out')

plot.df <- bind_rows(historical.df, county.df, buildout.df)

ggplot(data=plot.df, aes(x=year, y=Population, colour = Source)) + geom_line(linewidth=1) + ggtitle('South Kona Villages Population Projections')






