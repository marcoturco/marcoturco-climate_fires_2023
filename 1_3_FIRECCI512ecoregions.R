# Clear the workspace
rm(list = ls())
graphics.off()
gc()

#####################
## LOAD LIBRARIES ##
#####################
library(rgdal)        # To read shapefiles
library(raster)       # To associate .nc files with shapefiles

#############
## SETTING ##
#############

# Output directory
dir_out = '~/Dropbox/model/out_v2/'
# Directories for shapefiles and fire data
dir_shp = '/diskonfire/shapefiles/Ecoregions2017/'
dir_fire= '/diskonfire/FireCCI51'

# Set the range of years and months
years = 2001:2020
months = 1:12

####################
## LOAD SHAPEFILE ##
####################
file_shp = file.path(dir_shp, "Ecoregions2017_repaired.shp")
shp <- readOGR(file_shp)
shp_ll <- spTransform(shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

###############
## PROCESSING ##
###############
obs_nat = matrix(data = NA, nrow = length(years) * length(months), ncol = length(shp_ll))
obs_all = matrix(data = NA, nrow = length(years) * length(months), ncol = length(shp_ll))

k = 0
for (iyear in years) {
  for (imonth in months) {
    k = k + 1
    print(k)
    
    # Natural fires data
    file_nc <- paste0(dir_fire, '/', iyear, '/', iyear, sprintf("%02d", imonth), '01-ESACCI-L4_FIRE-BA-MODIS-fv5.1-natural-fires.nc')
    pr_nc = brick(file_nc)
    obs_nat[k,] = extract(pr_nc, shp_ll, na.rm = TRUE, fun = sum)
    
    # All fires data
    file_nc <- paste0(dir_fire, '/', iyear, '/', iyear, sprintf("%02d", imonth), '01-ESACCI-L4_FIRE-BA-MODIS-fv5.1.nc')
    pr_nc = brick(file_nc)
    obs_all[k,] = extract(pr_nc, shp_ll, na.rm = TRUE, fun = sum)
  }
}

####################
## SAVE PROCESSED ##
####################
save(obs_nat, file = file.path(dir_out, 'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020-natural-fires.Rdata'))
save(obs_all, file = file.path(dir_out, 'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020.Rdata'))
