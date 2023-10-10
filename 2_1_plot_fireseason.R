rm(list = ls())
graphics.off()
gc()

library(RColorBrewer)
library(terra)
library(sf)
library(tidyverse)
library(rmapshaper)
library(rnaturalearth)
library(pals)


## fix parameters
version = 'fireseason_reg'
product='modis' #'modis 


### fix parameters ANDRINA VENTISCA
# dir_out = '/home/andrina/Desktop/climate-fire/out_v2/'
# dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
# dir_shp='/home/andrina/Desktop/climate-fire/datos/shp/'
# # dir_data = 'C:/Users/andri/Desktop/climate/datos/'
# # dir_out = 'C:/Users/andri/Dropbox/model/out_v2/'
#  dir_fwi='/home/andrina/Desktop/climate-fire/scripts/1degree/FWI/'
#  dir_modis='/home/andrina/Dropbox/model/datos/'



dir_shp = '/Users/marco/Documents/virtualBox/temiav/Ecoregions2017/'
dir_data ='/Users/marco/Documents/dati/fire_climate_data/climate_fire/datos/'
dir_fwi='/Users/marco/Documents/dati/obs/ERA5/FWI/'
dir_out = '~/Dropbox/model/out_ecoregions/'
dir_modis ='~/Dropbox/model/datos/'




## load lon lat
# ecoregion
file_shp = file.path(dir_shp, "Ecoregions2017_repaired.shp")
eco <- st_read(file_shp) %>% st_make_valid()


## BA obs_reg
if (product == 'modis') {
  years = 2001:2021
  load(paste0(dir_modis, paste('MCD64CMQ-2001-2021.Rdata')))
}

nreg=dim(obs_reg)[2]
## fire season
# step 0: only cells with annual burned area series (BA>=2) that have at least 2 years were considered for further analysis.”).
BAy = array(0, dim = c(nreg, length(years)))
for (ireg in 1:nreg) {
  for (iyear in 1:length(years)) {
    i1 = (iyear - 1) * 12 + 1
    i2 = (iyear - 1) * 12 + 12
    BAy[ireg, iyear] = sum(obs_reg[ i1:i2,ireg], na.rm = TRUE) #*inout[i,j]
  }
}

mask0 = array(0, dim = c(nreg))
for (ireg in 1:nreg) {
  if (length(which(BAy[ireg,] > 0)) >= 2) {
    mask0[ireg] = 1
  }
}
 if (product == 'modis') {
  load(file.path(dir_out, 'firstmonth_reg_modis.RData'))
  load(file.path(dir_out, 'endmonth_reg_modis.RData'))
}


## aggregate BA over the fire season
BA_esa_y = array(NA, dim = c(nreg, length(years)))
for (ireg in 1:nreg) {
  if (mask0[ireg] == 1) {
    for (iyear in 1:length(years)) {
      if (firstmonth[ireg] >= 1) {
        i1 = (iyear - 1) * 12 + firstmonth[ireg]
        i2 = (iyear - 1) * 12 + endmonth[ireg]
        BA_esa_y[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE)
      } else
        if (iyear == 1) {
          next
        } else {
          i1 = (iyear - 1) * 12 + firstmonth[ireg]
          i2 = (iyear - 1) * 12 + endmonth[ireg]
          BA_esa_y[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE)
        }
    }
  }
}

auxsum=apply(BA_esa_y,c(1),sum,na.rm=T)
auxmean=apply(BA_esa_y,c(1),mean,na.rm=T)

mask = array(0, dim = c(nreg))
for (ireg in 1:nreg) {
  if (length(which(BA_esa_y[ireg,] > 0)) >= 10 && 100*auxsum[ireg]/sum(auxsum)>0.001) {
    mask[ireg] = 1
  }
}

sum(obs_reg,na.rm=T)
sum(obs_reg[,mask0==1],na.rm=T)
100*sum(obs_reg[,mask==1],na.rm=T)/sum(obs_reg,na.rm=T)

dim(eco)
length(which(mask==1))
length(which(mask==0))

if (product == 'modis') {
  save(mask,
       file = paste0(dir_out,
                     "mask_ecoregions_season_ja_2018_modis.RData"))
}


######### mask
wld <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") #https://github.com/ropensci/rnaturalearth
# Exclude Antarctica 
wld <- dplyr::filter(wld, name_long != "Antarctica")

#### end month
  
  # unir a los datos vectoriales
  mask2=mask*NA
  mask2[mask==1]=1
  end_month<-as.vector(endmonth)*mask2
  eco <- mutate(eco, end_month)
  #
  eco_only <- filter(eco, !is.na(end_month))
  
    # Generate the color palette
  br <- seq(1, 12)
  # colors <- stepped(12)
  # colors <- parula(12)
  colors <- ocean.phase(12)
  
  g <- ggplot(eco_only) +
    geom_sf(data = wld, fill = "grey90", colour = "white", linewidth = 0.1) +
    geom_sf(aes(fill = factor(end_month)), colour = "white", linewidth = 0.1) +
    coord_sf(crs = "ESRI:54030") +
    theme_minimal() +
    # theme(legend.position = "bottom") +
    # scale_fill_manual(values = colors)
    # guides(fill = guide_colorbar(title.position = "none", label.position = "bottom")) +
    theme(legend.title=element_blank()) +
    theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 12)) +
    scale_fill_manual(values = colors, labels = month.abb) 
  
  # # Save the map as PDF with improved label and color
  ggsave(paste0(dir_out,"end_month_",tolower(product),".png"), g, height = 10, width = 15, unit = "in",
         bg = "white", dpi = 600)
  
  
  
  #####SEASON LENGTH
  
  # season_length[season_length>=6]=6
  brk_prob <- seq(1, 6, length.out = 6)
  
  colores <-brewer.pal(7, "OrRd")
  
  legend_labels <- c(
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "> 6"
  )
  
  # Classify 
  season_length=endmonth-firstmonth+1
  
  season_length<-as.vector(season_length)*mask2
  eco <- mutate(eco, season_length)
  eco_only <- eco %>%
    filter(!is.na(season_length)) %>%
    mutate(
      discrete_season_length = case_when(
        season_length ==1         ~ legend_labels[1],
        season_length ==2         ~ legend_labels[2],
        season_length ==3         ~ legend_labels[3],
        season_length ==4         ~ legend_labels[4],
        season_length ==5         ~ legend_labels[4],
        season_length ==6         ~ legend_labels[6],
        season_length >= 6              ~ legend_labels[7]
      )
    )
  
  # Creating an empty data frame with all potential classes
  all_classes_df <- data.frame(discrete_season_length = legend_labels)
  
  # Plotting
  g <- ggplot(eco_only) +
    geom_sf(data = wld, fill = "grey90", colour = "white", linewidth = 0.1) +
    geom_sf(aes(fill = discrete_season_length), colour = "white", linewidth = 0.01) +
    # This layer ensures all levels are present but doesn't actually plot anything
    geom_blank(data = all_classes_df, aes(fill = discrete_season_length)) +
    
    # The rest of the plotting code remains the same
    scale_fill_manual(values = colores, breaks = legend_labels, labels = legend_labels, drop = FALSE) +
    labs(fill = "Fire season length") +
    coord_sf(crs = "ESRI:54030") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.key.size = unit(1, "cm"),
          legend.key.width = unit(1.2, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 13, face = "bold")) +
    guides(fill = guide_legend(nrow = 2))
  
  # Save the plot
  ggsave(paste0(dir_out, "lenght_fire_season_",tolower(product),".png"), g, height = 10, width = 15, unit = "in", bg = "white", dpi = 600)


  
