rm(list = ls())
graphics.off()
gc()


library(RColorBrewer)
library(terra)
library(sf)
library(tidyverse)
library(dplyr)
library(rmapshaper)
library(rnaturalearth)
library(pals)
library(magrittr)
library(sf)
library(broom)
library(ggplot2)
library(gridExtra)
library(mgcv)


### fix parameters ANDRINA VENTISCA
dir_out = '/home/andrina/Dropbox/model/out_ecoregions/'
dir_data = '/home/andrina/Desktop/climate-fire/datos/'
dir_shp='/home/andrina/Desktop/climate-fire/datos/shp/'

### fix parameters A
#dir_out = '~/Dropbox/model/out_ecoregions/'
#dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
#dir_shp='/Users/marco/Documents/virtualBox/temiav/Ecoregions2017/'

### version
transf="log" #"log" "radq"

if (transf=="log") {
  list_version=c('spei_sfwi_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018','spei_spei_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018')  
} else if (transf=="radq") {
  list_version=c('spei_sfwi_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018','spei_spei_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018')
}
list_version_names=c('spei-sfwi','spei-spei')
list_models=c('spei-sfwi','spei-spei','ANTspei','CONsfwi','CONspei')

load(paste0(dir_out,
                   "mask_ecoregions_season_ja_2018.RData"))

nreg=dim(mask)

all_r2 <- array(NA, dim = c(nreg,length(list_version)))
all_sfwi <- array(NA, dim = c(nreg,length(list_version)))
all_spei <- array(NA, dim = c(nreg,length(list_version)))
all_m_sfwi <- array(NA, dim = c(nreg,length(list_version)))
all_m_spei <- array(NA, dim = c(nreg,length(list_version)))
all_t_sfwi <- array(NA, dim = c(nreg,length(list_version)))
all_t_spei <- array(NA, dim = c(nreg,length(list_version)))

for (k in 1:length(list_version)) {
  version=list_version[k]
  load(paste0(dir_out, "corr_reconstruction_", version, ".RData"))
  load(paste0(dir_out, "sfwi_coef_reconstruction_", version, ".RData"))
  load(paste0(dir_out, "spei_coef_reconstruction_", version, ".RData"))
  
  load(paste0(dir_out, "best_m_SPI_fin_", version, ".RData")) #c(-14,-2)
  load(paste0(dir_out, "best_t_SPI_fin_", version, ".RData"))
  load(paste0(dir_out, "best_m_sfwi_fin_", version, ".RData")) #c(-12, 0)
  load(paste0(dir_out, "best_t_sfwi_fin_", version, ".RData"))
  # best_m_sfwi_fin
  # best_t_sfwi_fin
  # best_m_SPI_fin
  # best_t_SPI_fin
 
  all_r2[,k]=corr*corr
  all_spei[,k]=spei_coef
  all_sfwi[,k]=sfwi_coef
  all_t_spei[,k]=best_t_SPI_fin
  all_t_sfwi[,k]=best_t_sfwi_fin
  all_m_spei[,k]=best_m_SPI_fin
  all_m_sfwi[,k]=best_m_sfwi_fin
}

best_mod <- array(NA, dim = c(nreg))
best_r2 <- array(NA, dim = c(nreg))
best_m_spei <- array(NA, dim = c(nreg))
best_m_sfwi <- array(NA, dim = c(nreg))
best_t_spei <- array(NA, dim = c(nreg))
best_t_sfwi <- array(NA, dim = c(nreg))
best_spei_coef <- array(NA, dim = c(nreg))
best_sfwi_coef <- array(NA, dim = c(nreg))

for (ireg in 1:nreg) {
  if (length(which(is.na(all_r2[ireg, ]))) < length(list_version)) {
    # if (length(which(is.na(all_r2[i,j,])))==1) {cc}
    aux = which(all_r2[ireg, ] == max(all_r2[ireg, ], na.rm = T))
    if (length(aux) == 2) {
      aux = aux[1]
    }
    tmp = gregexpr(pattern = '-', list_version_names[aux])
    best_r2[ireg] = all_r2[ireg, aux]
    # list_models=c('spei-sfwi','spei-spei','ANTspei','CONsfwi','CONspei')
    if (is.na(all_spei[ireg, aux]) & !is.na(all_sfwi[ireg, aux])) {
      aux2 = substr(list_version_names[aux],
                    as.numeric(tmp[[1]]) + 1,
                    nchar(list_version_names[aux]))
      aux2 = paste0('^CON', aux2, '$')
      best_m_sfwi[ireg] <- all_m_sfwi[ireg, aux]
      best_t_sfwi[ireg] <- all_t_sfwi[ireg, aux]
      best_sfwi_coef[ireg] <- all_sfwi[ireg,aux]
    }
    if (!is.na(all_spei[ireg, aux]) & is.na(all_sfwi[ireg, aux])) {
      aux2 = substr(list_version_names[aux], 1, as.numeric(tmp[[1]]) - 1)
      aux2 = paste0('^ANT', aux2, '$')
      best_m_spei[ireg] <- all_m_spei[ireg, aux]
      best_t_spei[ireg] <- all_t_spei[ireg, aux]
      best_spei_coef[ireg] <- all_spei[ireg,aux]
    }
    if (!is.na(all_spei[ireg, aux]) &
        !is.na(all_sfwi[ireg, aux])) {
      aux2 = list_version_names[aux]
      best_m_sfwi[ireg] <- all_m_sfwi[ireg, aux]
      best_t_sfwi[ireg] <- all_t_sfwi[ireg, aux]
      best_m_spei[ireg] <- all_m_spei[ireg, aux]
      best_t_spei[ireg] <- all_t_spei[ireg, aux]
      best_sfwi_coef[ireg] <- all_sfwi[ireg,aux]
      best_spei_coef[ireg] <- all_spei[ireg,aux]
    }
    
    best_mod[ireg] = grep(aux2, list_models)
    if (best_mod[ireg] == 2 ||
        best_mod[ireg] == 5) {
      best_sfwi_coef[ireg] = -best_sfwi_coef[ireg]
    }
  }
}






summary(as.vector(best_r2))

print(paste0(
  "Variance explained (spatial median) = ",
  round(median(best_r2, na.rm = TRUE), digits = 2),
  "\n",
  "Percentage valid models = ",
  round(100 * length(which(!is.na(
    best_r2
  ))) / length(which(mask == 1)), digits = 0),
  '%'))

## climatologies
load(paste0(dir_out, "AI_ERA5_MSWEP_2001_2021_ecoregions.RData")) #AI
load(paste0(dir_out, "PREC_ERA5_MSWEP_2001_2021_ecoregions_annual.RData")) #AI
load(paste0(dir_out, "PET_ERA5_MSWEP_2001_2021_ecoregions_annual.RData")) #AI
# load(paste0(dir_out, "TAS_ERA5_2001_2021_ecoregions.RData")) #AI


# Assuming your data frame is named df
df <- data.frame(
  AI = AI, 
  PREC = PREC,
  PET = PET,
  # TAS = TAS,
  best_r2 = best_r2,
  best_sfwi_coef = abs(best_sfwi_coef),
  best_m_sfwi = best_m_sfwi,
  best_t_sfwi = best_t_sfwi,
  best_spei_coef = best_spei_coef,
  # best_sfwi_coef = best_sfwi_coef,
  best_m_spei = best_m_spei,
  best_t_spei = best_t_spei
)


# A function to extract and format the correlation coefficient based on p-value
get_cor_p <- function(x, y) {
  # test <- cor.test(x, y, method = "spearman")
  test <- cor.test(x, y, method = "pearson")
  cor_val <- round(test$estimate, 2) # rounding to 2 decimal places
  if (test$p.value < 0.01) {
    return(paste0(cor_val, "**"))
  } else if (test$p.value < 0.05) {
    return(paste0(cor_val, "*"))
  } else {
    return(as.character(cor_val))
  }
}

# Variables of interest
best_cols <- grep("^best_", names(df), value = TRUE)
fixed_cols <- c("AI", "PREC", "PET")

# Compute the correlations
results <- sapply(fixed_cols, function(fixed_col) {
  sapply(best_cols, function(best_col) {
    get_cor_p(df[[fixed_col]], df[[best_col]])
  })
}, simplify = "matrix")

# Convert to a data frame for better readability
results_df <- as.data.frame(results)
rownames(results_df) <- best_cols

results_df







#######

df_sfwi <- df[, c("AI", "PREC", "PET", "best_sfwi_coef", "best_m_sfwi", "best_t_sfwi")]
df_sfwi <- df_sfwi[complete.cases(df_sfwi), ]

df_spei <- df[, c("AI", "PREC", "PET", "best_spei_coef", "best_m_spei", "best_t_spei")]
df_spei <- df_spei[complete.cases(df_spei), ]

###### plot

library(ggplot2)
library(gridExtra)
library(dplyr)
library(broom)

# Helper function to format R^2 and p-value
format_labels <- function(r2, p_value){
  asterisks <- ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", ""))
  return(sprintf("%.1f%% %s", r2*100, asterisks))
}


# For p2: best_sfwi_coef vs PET
model_p2 <- lm(best_sfwi_coef ~ PET, data = df_sfwi)
glance_p2 <- glance(model_p2)
p2 <- ggplot(df_sfwi, aes(x = PET, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PET [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), label = "a)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels(glance_p2$r.squared, glance_p2$p.value), hjust = 1.1, vjust = 1.5)

# For p3: best_sfwi_coef vs PREC
model_p3 <- lm(best_sfwi_coef ~ PREC, data = df_sfwi)
glance_p3 <- glance(model_p3)
p3 <- ggplot(df_sfwi, aes(x = PREC, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PREC [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), label = "b)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels(glance_p3$r.squared, glance_p3$p.value), hjust = 1.1, vjust = 1.5)


# For p5: best_spei_coef vs PET
model_p5 <- lm(best_spei_coef ~ PET, data = df_spei)
glance_p5 <- glance(model_p5)
p5 <- ggplot(df_spei, aes(x = PET, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PET [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PET), y = max(df_spei$best_spei_coef), label = "c)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PET), y = max(df_spei$best_spei_coef), 
           label = format_labels(glance_p5$r.squared, glance_p5$p.value), hjust = 1.1, vjust = 1.5)

# For p6: best_spei_coef vs PREC
model_p6 <- lm(best_spei_coef ~ PREC, data = df_spei)
glance_p6 <- glance(model_p6)
p6 <- ggplot(df_spei, aes(x = PREC, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PREC [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PREC), y = max(df_spei$best_spei_coef), label = "d)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PREC), y = max(df_spei$best_spei_coef), 
           label = format_labels(glance_p6$r.squared, glance_p6$p.value), hjust = 1.1, vjust = 1.5)

# Arrange plots
# plot_grid <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
plot_grid <- grid.arrange( p2, p3, p5, p6, ncol = 2)


ggsave(paste0(dir_out,"scatterplots_lm.pdf"), plot_grid, width = 11, height = 8.5, units = "in")





####### GAM



library(mgcv)

# Helper function to format R^2 and p-value for GAM models
format_labels_gam <- function(model){
  r2 <- summary(model)$r.sq
  p_value <- coef(summary(model))[,4][1]  # Assuming you're interested in the p-value for the smooth term
  asterisks <- ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", ""))
  return(sprintf("%.1f%% %s", r2*100, asterisks))
}

# For p2: best_sfwi_coef vs PET
model_p2 <- gam(best_sfwi_coef ~ s(PET), data = df_sfwi)
p2 <- ggplot(df_sfwi, aes(x = PET, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="red") +
  labs(x = "PET [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), label = "a)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels_gam(model_p2), hjust = 1.1, vjust = 1.5)

# For p3: best_sfwi_coef vs PREC
model_p3 <- gam(best_sfwi_coef ~ s(PREC), data = df_sfwi)
p3 <- ggplot(df_sfwi, aes(x = PREC, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="red") +
  labs(x = "PREC [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), label = "b)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels_gam(model_p3), hjust = 1.1, vjust = 1.5)

# For p5: best_spei_coef vs PET
model_p5 <- gam(best_spei_coef ~ s(PET), data = df_spei)
p5 <- ggplot(df_spei, aes(x = PET, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="red") +
  labs(x = "PET [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PET), y = max(df_spei$best_spei_coef), label = "c)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PET), y = max(df_spei$best_spei_coef), 
           label = format_labels_gam(model_p5), hjust = 1.1, vjust = 1.5)

# For p6: best_spei_coef vs PREC
model_p6 <- gam(best_spei_coef ~ s(PREC), data = df_spei)
p6 <- ggplot(df_spei, aes(x = PREC, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="red") +
  labs(x = "PREC [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PREC), y = max(df_spei$best_spei_coef), label = "d)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PREC), y = max(df_spei$best_spei_coef), 
           label = format_labels_gam(model_p6), hjust = 1.1, vjust = 1.5)

# Arrange plots
plot_grid <- grid.arrange(p2, p3, p5, p6, ncol = 2)
ggsave(paste0(dir_out,"scatterplots_gam.pdf"), plot_grid, width = 11, height = 8.5, units = "in")



## boxplot t

########################
##ANTECEDENT CLIMATE
#######################

### BEST M PET
# Calculate the number of samples for each value of best_m_spei
spectral_13_hex <-
  c(
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FFFFBF",
    "#FEE08B",
    "#FDAE61",
    "#F28E2B",
    "#F46D43",
    "#D94801",
    "#D53E4F",
    "#9E0142"
  )
# Reverse the color palette
reversed_palette <- rev(spectral_13_hex)

counts <- df_spei %>%
  group_by(best_m_spei) %>%
  summarise(n = n())

counts$n <- round(100*counts$n/length(which(mask==1)))
# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate and round the 98th percentile of PET for positioning labels and setting y-axis limit
label_position <- round(quantile(df_spei$PET, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f1 <- ggplot(df_spei, aes(x=best_m_spei, y=PET, group=best_m_spei, fill=as.factor(best_m_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_m_spei, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting x-axis breaks and labels to ensure correct positioning and labeling
  scale_x_continuous(breaks = unique(df_spei$best_m_spei),
                     labels = as.character(unique(df_spei$best_m_spei))) +
  
  # Setting reversed color palette for the boxes
  scale_fill_manual(values = reversed_palette) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best m SPEI", y = "PET") + 
  theme_minimal()


print(f1)


### BEST M PPREC
# Calculate the number of samples for each value of best_m_spei
spectral_13_hex <-
  c(
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FFFFBF",
    "#FEE08B",
    "#FDAE61",
    "#F28E2B",
    "#F46D43",
    "#D94801",
    "#D53E4F",
    "#9E0142"
  )
# Reverse the color palette
reversed_palette <- rev(spectral_13_hex)

counts <- df_spei %>%
  group_by(best_m_spei) %>%
  summarise(n = n())

counts$n <- round(100*counts$n/length(which(mask==1)))
# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate and round the 98th percentile of PREC for positioning labels and setting y-axis limit
label_position <- round(quantile(df_spei$PREC, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f2 <- ggplot(df_spei, aes(x=best_m_spei, y=PREC, group=best_m_spei, fill=as.factor(best_m_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_m_spei, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting x-axis breaks and labels to ensure correct positioning and labeling
  scale_x_continuous(breaks = unique(df_spei$best_m_spei),
                     labels = as.character(unique(df_spei$best_m_spei))) +
  
  # Setting reversed color palette for the boxes
  scale_fill_manual(values = reversed_palette) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best m SPEI", y = "PREC") + 
  theme_minimal()

print(f2)


#BEST t PET 
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_spei$best_t_spei))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_spei", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PET for positioning labels
label_position <- round(quantile(df_spei$PET, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f3 <- ggplot(df_spei, aes(x=factor(best_t_spei), y=PET, fill=factor(best_t_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_t_spei, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting color palette for the boxes
  scale_fill_manual(values = selected_Dark2_hex) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best t SPEI", y = "PET") + 
  theme_minimal()

print(f3)



###BEST t PREC
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_spei$best_t_spei))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_spei", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PREC for positioning labels
label_position <- round(quantile(df_spei$PREC, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f4 <- ggplot(df_spei, aes(x=factor(best_t_spei), y=PREC, fill=factor(best_t_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_t_spei, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting color palette for the boxes
  scale_fill_manual(values = selected_Dark2_hex) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best t SPEI", y = "PREC") + 
  theme_minimal()

print(f4)

########################
## CONCURRENT CLIMATE 
#######################

# Create color palette
spectral_13_hex <- c(
  "#4D3E82",
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#ABDDA4",
  "#E6F598",
  "#FFFFBF",
  "#FEE08B",
  "#FDAE61",
  "#F28E2B",
  "#F46D43",
  "#D94801",
  "#D53E4F",
  "#9E0142"
)
reversed_palette <- rev(spectral_13_hex)

# Calculate counts, percentages, and median
counts <- df_sfwi %>%
  group_by(best_m_sfwi) %>%
  summarise(n = n(),
            median_PET = median(PET, na.rm = TRUE)) %>%
  mutate(n = round(100 * n / length(which(mask == 1))),
         label = paste0(n, "%"))


# Create boxplot
f5 <- ggplot(df_sfwi, aes(x = best_m_sfwi, y = PET, group = best_m_sfwi, fill = as.factor(best_m_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(y = 600, label = label),  # Cambiamos el y aquí
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0) +  # Ajustamos vjust para centrar el texto verticalmente
  scale_x_continuous(breaks = unique(df_sfwi$best_m_sfwi),
                     labels = as.character(unique(df_sfwi$best_m_sfwi))) +
  scale_fill_manual(values = spectral_13_hex) +
  guides(fill = FALSE) +
  labs(x = "Best m SFWI", y = "PET") + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 14)
  )

# Print the plot
print(f5)











### BEST M PPREC
# Calculate the number of samples for each value of best_m_sfwi
spectral_13_hex <-
  c(
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FFFFBF",
    "#FEE08B",
    "#FDAE61",
    "#F28E2B",
    "#F46D43",
    "#D94801",
    "#D53E4F",
    "#9E0142"
  )
# Reverse the color palette
reversed_palette <- rev(spectral_13_hex)

counts <- df_sfwi %>%
  group_by(best_m_sfwi) %>%
  summarise(n = n())

counts$n <- round(100*counts$n/length(which(mask==1)))
# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate and round the 98th percentile of PET for positioning labels and setting y-axis limit
label_position <- round(quantile(df_sfwi$PREC, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f6 <- ggplot(df_sfwi, aes(x=best_m_sfwi, y=PREC, group=best_m_sfwi, fill=as.factor(best_m_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_m_sfwi, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting x-axis breaks and labels to ensure correct positioning and labeling
  scale_x_continuous(breaks = unique(df_sfwi$best_m_sfwi),
                     labels = as.character(unique(df_sfwi$best_m_sfwi))) +
  
  # Setting reversed color palette for the boxes
  scale_fill_manual(values = reversed_palette) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best m SFWI", y = "PREC") + 
  theme_minimal()


print(f6)

#BEST t PET 
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_sfwi$best_t_sfwi))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_sfwi", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PET for positioning labels
label_position <- round(quantile(df_sfwi$PET, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f7 <- ggplot(df_sfwi, aes(x=factor(best_t_sfwi), y=PET, fill=factor(best_t_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_t_sfwi, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting color palette for the boxes
  scale_fill_manual(values = selected_Dark2_hex) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best t SFWI", y = "PET") + 
  theme_minimal()

print(f7)



###BEST t PREC
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_sfwi$best_t_sfwi))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_sfwi", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PREC for positioning labels
label_position <- round(quantile(df_sfwi$PREC, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f8 <- ggplot(df_sfwi, aes(x=factor(best_t_sfwi), y=PREC, fill=factor(best_t_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  
  # Add text to show number of samples
  geom_text(data=counts, aes(x=best_t_sfwi, y=label_position, label=label), vjust=-0.5) +
  
  # Set y-axis limit
  coord_cartesian(ylim = c(NA, label_position)) +
  
  # Setting color palette for the boxes
  scale_fill_manual(values = selected_Dark2_hex) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  labs(x = "Best t SFWI", y = "PREC") + 
  theme_minimal()

print(f8)

########GUARDAR FIGURAS- combina ejes automáticamente
#install.packages("patchwork")
library(patchwork)

# Arrange plots
plot_grid <- f5 + f6 + f7 + f8 +f1+f2+f3+f4 + plot_layout(ncol = 2)
ggsave(paste0(dir_out,"boxplot_cc_ac.pdf"), plot_grid, width = 22, height = 28, units = "in")


