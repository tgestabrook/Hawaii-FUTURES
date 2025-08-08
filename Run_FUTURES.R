library(terra)
library(tidyverse)
library(tidyterra)
library(rgrass)

initGRASS(gisBase = "C:/Program Files/GRASS GIS 8.3",
          gisDbase = "F:/Hawaii",
          home=tempdir(),
          location='GHS',
          mapset='RImport',
          override = TRUE)

execGRASS("g.list", parameters=list(type='all'))
#execGRASS("g.mapset")

repeat_runs <- 100
name_prefix <- '100yr_run'

################################################################################
############### Write all maps to GRASS ########################################
################################################################################

Regions.r <- vect(file.path('Data', 'FAZ.gpkg')) %>%  # overwrite regions to again have all regions in place
  rasterize(LF01.r, field='FAZID') %>% as.int()
Regions.r[is.na(Regions.r)] <- 15450
names(Regions.r) <- 'Regions'
predictor_stack$Regions <- Regions.r

for(map in names(predictor_stack)){
  if(map=="LF01"){next}
  write_RAST(predictor_stack[map], vname=map, overwrite=T)
}


execGRASS("g.list", parameters=list(type='all'))
execGRASS("r.info", parameters=list(map="U90"))

################################################################################
############## Run FUTURES #####################################################
################################################################################
if(!exists("devpr_20")){devpr_20 <- rast(file.path('Inputs', 'devpr_20.tif'))}
write_RAST(devpr_20, vname='devpr_20', overwrite=T)
write_RAST(as.int(U00.r), vname='U00', overwrite=T)
write_RAST(as.int(U20.r), vname='U20', overwrite=T)

execGRASS('r.futures.potsurface', input='Inputs/Hawaii_model.csv', subregions='Regions', output='Site_suitability')
SS.r <- read_RAST('Site_suitability')
plot(SS.r)
writeRaster(SS.r, file.path('Outputs', 'Site_suitability.tif'), overwrite=T)

ss_model.df <- read.csv('Inputs/Hawaii_model.csv')

# do the demand calculation
execGRASS('r.futures.demand', development='U80,U90,U00,U10,U15,U20', observed_population='Inputs/pop_trend.csv', subregions='Regions', 
          projected_population='Inputs/pop_projection.csv', simulation_times=2021:2100, plot='Inputs/plot_demand.png', demand='Inputs/demand.csv', separator='comma', flags = 'overwrite')

# do patch calibration
execGRASS('r.futures.calib', parameters=list(development_start='U00', development_end='U20', subregions='Regions', patch_sizes='Inputs/patches.txt', calibration_results='Inputs/calib.csv', 
          patch_threshold=1800, 'repeat'=10, compactness_mean=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), compactness_range=c(0.1,0.05), discount_factor=c(0.1,0.3,0.5,0.7,0.9),
          predictors=names(ss_model.df)[4:length(names(ss_model.df))], demand='Inputs/demand.csv',
          devpot_params='Inputs/Hawaii_model.csv', num_neighbors=4, seed_search='probability', development_pressure='devpr_20', development_pressure_approach='gravity',
          n_dev_neighbourhood=5, gamma=0.5, scaling_factor=1), flags='overwrite')

calib_results.df <- read.csv(file.path('Inputs', 'calib.csv'))
head(calib_results.df)

# run the model
execGRASS('r.futures.parallelpga', parameters=list('repeat'=repeat_runs, nprocs=4, developed='U20', subregions='Regions', output=name_prefix, patch_sizes='Inputs/patches.txt', 
                                                   predictors=names(ss_model.df)[4:length(names(ss_model.df))],
                                                   demand='Inputs/demand.csv', devpot_params='Inputs/Hawaii_model.csv', num_neighbors=4, seed_search='probability', 
                                                   development_pressure='devpr_20', development_pressure_approach='gravity', n_dev_neighbourhood=5, gamma=0.5, scaling_factor=1,
                                                   compactness_mean=calib_results.df[1, 'compactness_mean'], compactness_range=calib_results.df[1, 'compactness_range'], discount_factor=calib_results.df[1, 'discount_factor']), flags='overwrite')

#out1 <- read_RAST('testrun_final_run1')
#plot(out1)

################################################################################
############## Post-process ####################################################
################################################################################

library(magick)

pop_input_projection<-read.csv(file.path('Inputs', 'pop_projection.csv'))
sim_length <- nrow(pop_input_projection)

elev.r <- rast(file.path('Data', 'Predictors', 'elevation.tif')); names(elev.r) <- 'elev'
hillshade.r <- shade(terrain(elev.r, 'slope', unit='radians'), terrain(elev.r, 'aspect', unit='radians'))  

s_kona_aoi.shp <- vect(file.path('Data', 'South_Kona_limit_V3.shp')) %>% project(elev.r)
#AOI <- ext(s_kona_aoi.shp)
AOI <- ext(elev.r)

### Plot patches
patches <- read_lines(file.path('Inputs', 'patches.txt')) %>% as.numeric()
hist(patches, breaks=200, xlim=c(0,100))

### Loop through and stack rasters
run_stack <- read_RAST(paste0(name_prefix,'_run', 1))
start_urb <- ifel(run_stack[[1]] == 0, 1, NA)

for (i in 2:repeat_runs){
  sim <- read_RAST(paste0(name_prefix,'_run', i))
  run_stack <- c(run_stack, sim)
}

plot(run_stack)

### Misc
#elev.r[is.na(run_stack[[1]])] <- NA
#hillshade.r <- shade(terrain(elev.r, 'slope'), terrain(elev.r, 'aspect')) 

### Make heatmap
binarize <- function(x) {
  ifel(x > 0, 1, 0)
}

final_binary <- sapp(run_stack, binarize)

plot(final_binary)

heatmap <- sum(final_binary)

# png(filename = file.path('Outputs', paste0(name_prefix, '_heatmap.png')), height = 10, width = 6, units = 'in', res = 300)
# plot(heatmap, ext=AOI)
# plot(start_urb, ext=AOI, col='black', add=T,legend=F)
# plot(hillshade.r, col=colorRampPalette(c('black','grey20','white'))(100),add=T,ext=AOI,alpha=0.4,legend=F)
# plot(s_kona_aoi.shp, add=T)
# dev.off()

png(filename = file.path('Outputs', paste0(name_prefix, '_heatmap_full.png')), height = 14, width = 14, units = 'in', res = 300)
plot(heatmap, ext=AOI)
plot(start_urb, ext=AOI, col='black', add=T,legend=F)
plot(hillshade.r, col=colorRampPalette(c('black','grey20','white'))(100),add=T,ext=AOI,alpha=0.4,legend=F)
plot(s_kona_aoi.shp, add=T)
dev.off()

heatmap[start_urb==1] <- -1
writeRaster(heatmap, file.path('Outputs', paste0(name_prefix, '_heatmap.tif')), overwrite=T)

## 2045
clip_to_year <- function(r){
  ifel(r>25, 0, r)
}

heatmap_2045 <- sapp(run_stack, clip_to_year) %>%
  sapp(., binarize) %>% sum()
heatmap_2045[start_urb==1] <- -1
writeRaster(heatmap_2045, file.path('Outputs', paste0(name_prefix, '2045_heatmap.tif')), overwrite=T)

### Make gif
run_to_gif <- 1

for(i in 1:sim_length){
  timestep_map <- ifel(run_stack[[run_to_gif]] > i, -1, run_stack[[run_to_gif]])
  timestep_map <- ifel(timestep_map > 0, 1, timestep_map)
  
  ## Plot frame:
  png(file.path('Outputs', 'Gif_frames',paste0('_GIF_frames_',name_prefix,i,'.png')),width=6,height=8.32,res=500,units='in')
  par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(0,0,0,0),mgp=c(1,0.01,0),tck=-0.002,ps=10,cex=1)
  # plot.new()
  # plot.window(xlim=AOI[1:2], ylim=AOI[3:4],xaxs="i",yaxs="i",asp=1)
  # plot(timestep_map,legend=F,add=T,alpha=0.9)
  plot(timestep_map,col=c('white', 'black', 'green'), ext=AOI, legend=F)
  plot(elev.r,col=colorRampPalette(c('black','white'))(100),add=T,alpha=0.15,ext=AOI,legend=F)
  plot(hillshade.r,col=colorRampPalette(c('black','grey20','white'))(100),add=T,alpha=0.4,ext=AOI,legend=F);box(lwd=3)
  plot(s_kona_aoi.shp, add=T)
  mtext(paste0('year: ',2022+as.numeric(i)),font=2,line=-1,adj=0.95,cex=1)
  dev.off()
}

## Now load images and make gif
cat('\nLoading images...\n')
imgs <- dir(file.path('Outputs', 'Gif_frames'))[grepl(paste0('_GIF_frames_', name_prefix),dir(file.path('Outputs', 'Gif_frames')))]
imgs<-imgs[order(as.numeric(str_match(imgs, '[0-9]+')))]
img_list <- lapply(file.path('Outputs', 'Gif_frames',imgs), image_read)

## join the images together
img_joined <- image_join(img_list)
## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2) 

## view animated image
# img_animated

cat('Writing GIF...')
## save to disk
image_write(image = img_animated,path = file.path('Outputs',paste0(name_prefix, run_to_gif, '_animation.gif')))
cat('Success!\n\n')

################################################################################
############## OSDS-process ####################################################
################################################################################
OSDS_2020.shp <- vect(file.path('Data', 'OSDS', 'OSDS_database_DST_project_11302023.shp'))
OSDS_2020_count <- nrow(OSDS_2020.shp)

demand.df <- read.csv(file.path('Inputs', 'demand.csv')) %>%
  mutate(cumulative_demand = cumsum(X1))

# of cells in South Kona AOI
U20.r <- read_RAST('U20')
sk_2020_dev_clip.r <- U20.r %>% crop(s_kona_aoi.shp, mask=T)
dev_count_2020 <- sum(values(sk_2020_dev_clip.r), na.rm=T)
print(OSDS_2020_count/dev_count_2020)

sk_2100_dev_clip.r <- rast(file.path('Outputs', paste0(name_prefix, '_heatmap.tif')))  %>% crop(s_kona_aoi.shp, mask=T)
sk_2045_dev_clip.r <- rast(file.path('Outputs', paste0(name_prefix, '2045_heatmap.tif')))  %>% crop(s_kona_aoi.shp, mask=T)

plot(sk_2045_dev_clip.r)

dev_count_2045 <- dev_count_2020 + sum(values(sk_2045_dev_clip.r), na.rm=T)/repeat_runs
dev_count_2100 <- dev_count_2020 + sum(values(sk_2100_dev_clip.r), na.rm=T)/repeat_runs

OSDS_2045_count <- (dev_count_2045/dev_count_2020)*nrow(OSDS_2020.shp)
OSDS_2100_count <- (dev_count_2100/dev_count_2020)*nrow(OSDS_2020.shp)

New_OSDS_to_place45 <- as.integer(OSDS_2045_count - OSDS_2020_count)
New_OSDS_to_place100 <- as.integer(OSDS_2100_count - OSDS_2020_count)

Pixels_for_OSDS_2045 <- as.data.frame(sk_2045_dev_clip.r, na.rm=T, xy=T) %>%
  arrange(desc(sum)) %>%
  slice_head(n=New_OSDS_to_place45) %>%
  vect(geom=c('x','y'), crs = crs(sk_2045_dev_clip.r))
writeVector(Pixels_for_OSDS_2045, file.path('Outputs', 'Septic_sites_2045.gpkg'))

Pixels_for_OSDS_2100 <- as.data.frame(sk_2100_dev_clip.r, na.rm=T, xy=T) %>%
  arrange(desc(sum)) %>%
  slice_head(n=New_OSDS_to_place100) %>%
  vect(geom=c('x','y'), crs = crs(sk_2045_dev_clip.r))
writeVector(Pixels_for_OSDS_2100, file.path('Outputs', 'Septic_sites_2100.gpkg'))

################################################
# ID typical run 

runstack_typical <- final_binary * ifel(heatmap==-1, 0, heatmap)

overlaps <- as.data.frame(runstack_typical)

colSums(overlaps)[which.max(colSums(overlaps))]

plot(run_stack['file619824602d3f'])

Pixels_for_OSDS_2100_typical <- as.data.frame(run_stack['file619824602d3f']%>% crop(s_kona_aoi.shp, mask=T), xy=T) %>% 
  filter(file619824602d3f>0)

Extra_OSDS <- as.data.frame(sk_2100_dev_clip.r, na.rm=T, xy=T) %>%
  arrange(desc(sum)) %>%
  anti_join(Pixels_for_OSDS_2100_typical) %>%
  slice_head(n=(New_OSDS_to_place100-nrow(Pixels_for_OSDS_2100_typical)))

Pixels_for_OSDS_2100_typical <- Pixels_for_OSDS_2100_typical %>%
  bind_rows(Extra_OSDS) %>%
  vect(., geom=c('x','y'), crs = crs(sk_2045_dev_clip.r)) %>%
  writeVector(file.path('Outputs', 'Septic_sites_2100_alternate.gpkg'))


par(mfrow=c(1,2))
plot(heatmap, ext=ext(s_kona_aoi.shp), main="Frequency-based placement")
plot(hillshade.r,col=colorRampPalette(c('black','grey20','white'))(100),add=T,alpha=0.25,ext=s_kona_aoi.shp,legend=F)
plot(s_kona_aoi.shp, add=T)
points(Pixels_for_OSDS_2100, alpha=1, cex = 0.1, col='red')

plot(heatmap, ext=ext(s_kona_aoi.shp),  main="Single-simulation placement")
plot(hillshade.r,col=colorRampPalette(c('black','grey20','white'))(100),add=T,alpha=0.25,ext=s_kona_aoi.shp,legend=F)
plot(s_kona_aoi.shp, add=T)
points(Pixels_for_OSDS_2100_typical, cex = 0.1, col='red')













