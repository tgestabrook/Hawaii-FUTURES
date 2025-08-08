library(terra)
library(tidyverse)
library(tidyterra)
library(magick)

name_prefix <- 'test_run'
repeat_runs <- 10
sim_length <- 22

elev.r <- rast(file.path('Data', 'Predictors', 'elevation.tif'))
hillshade.r <- shade(terrain(elev.r, 'slope', unit='radians'), terrain(elev.r, 'aspect', unit='radians'))  

s_kona_aoi.shp <- vect(file.path('Data', 'South_Kona_limit_V3.shp')) %>% project(elev.r)
AOI <- ext(s_kona_aoi.shp)
#AOI <- ext(elev.r)

### Plot patches
patches <- read_lines(file.path('Inputs', 'patches2.txt')) %>% as.numeric()
hist(patches, breaks=200, xlim=c(0,100))

### Loop through and stack rasters
run_stack <- rast(file.path('Outputs', paste0(name_prefix, 1, '.tif')))
start_urb <- ifel(run_stack[[1]] == 0, 1, NA)

for (i in 2:repeat_runs){
  sim <- rast(file.path('Outputs', paste0(name_prefix, i, '.tif')))
  run_stack <- c(run_stack, sim)
}

plot(run_stack)

### Misc
elev.r[is.na(run_stack[[1]])] <- NA
hillshade.r <- shade(terrain(elev.r, 'slope'), terrain(elev.r, 'aspect')) 


### Make heatmap
binarize <- function(x) {
  ifel(x > 0, 1, 0)
}

final_binary <- sapp(run_stack, binarize)

plot(final_binary)

heatmap <- sum(final_binary)

png(filename = file.path('Outputs', paste0(name_prefix, '_heatmap.png')), height = 10, width = 6, units = 'in', res = 300)
plot(heatmap, ext=AOI)
plot(start_urb, ext=AOI, col='black', add=T,legend=F)
plot(hillshade.r, col=colorRampPalette(c('black','grey20','white'))(100),add=T,ext=AOI,alpha=0.4,legend=F)
dev.off()

writeRaster(heatmap, file.path('Outputs', paste0(name_prefix, '_heatmap.tif')))

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






