#####################################################################
## File to compute the fire season as in Abatzoglou et al. (2018)
#####################################################################

# Clear workspace
rm(list = ls())
graphics.off()
gc()

# Directories
dir_data <- '~/Documents/dati/fire_climate_data/climate_fire/datos/'
dir_fwi <- '~/Documents/dati/obs/ERA5/FWI/'
dir_out <- '~/Dropbox/model/out_ecoregions/'
dir_modis ='~/Dropbox/model/datos/'

years = 2001:2021
load(paste0(dir_modis, paste('MCD64CMQ-2001-2021.Rdata')))

# mask0
BAy <- matrix(NA, nrow=length(years), ncol=ncol(obs_reg))
for (i in 1:ncol(obs_reg)) {
  for (iyear in 1:length(years)) {
    i1 <- (iyear - 1) * 12 + 1
    i2 <- (iyear - 1) * 12 + 12
    BAy[iyear, i] <- sum(obs_reg[i1:i2, i], na.rm=TRUE)
  }
}

mask0 <- rep(NA, ncol(obs_reg))
for (i in 1:ncol(obs_reg)) {
  if (length(which(BAy[,i] > 0)) >= 2) {
    mask0[i] <- 1
  }
}

# fire season as in Aabatzoglou et al. (2018)
BAmm <- matrix(NA, nrow=12, ncol=ncol(obs_reg))
for (k in 1:length(mask0)) {
  if (!is.na(mask0[k]) && mask0[k] == 1) {
    BAm <- matrix(NA, nrow=length(years), ncol=12)
    for (im in 1:12) {
      BAm[, im] <- obs_reg[seq(im, by=12, length.out=length(years)), k]
      BAmm[im, k] <- mean(BAm[, im], na.rm=TRUE)
    }
  }
}

# Continue with calculating the peak months
endmonth <- rep(NA, ncol(obs_reg))
lengthmonth <- rep(NA, ncol(obs_reg))
maxperiod <- rep(NA, ncol(obs_reg))
bap <- matrix(NA, nrow=12, ncol=11)

for (k in 1:length(mask0)) {
  if (!is.na(mask0[k])) {
    BAA <- BAmm[,k]
    for (mi in 1:12) {
      for (li in 1:11) {
        mstart <- mi + 1 - li
        mend <- mi
        if (mstart < 1) {
          mstart <- 12 + mstart
          months <- c(1:mend, mstart:12)
        } else {
          months <- mstart:mend
        }
        bap[mi, li] <- sum(BAA[months], na.rm=TRUE) / sum(BAA, na.rm=TRUE)
      }
    }
    
    find80 <- bap > 0.8
    foundit <- FALSE
    ii <- 1
    while (!foundit && ii <= 11) {
      if (any(find80[, ii])) {
        foundit <- TRUE
        f <- which.max(bap[, ii])
        endmonth[k] <- f
        lengthmonth[k] <- ii
        maxperiod[k] <- bap[f, ii]
      } else {
        ii <- ii + 1
      }
      
    }
   
  }
  
  # if (is.na(endmonth[k]) || endmonth[k] == 0) {
  #   stop("Error in calculations")  # dsjdsaf (This was a MATLAB error trigger, adapted for R)
  # }
}

firstmonth <- endmonth - lengthmonth + 1

# Saving the results
save(endmonth, file=paste0(dir_out, 'endmonth_reg_modis.RData'))
save(firstmonth, file=paste0(dir_out, 'firstmonth_reg_modis.RData'))
