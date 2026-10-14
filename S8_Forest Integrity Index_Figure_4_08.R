###########################################################################
###########################################################################
##
##
##  08 Forest Integrity Index 
##
##
##  Author: D.Pollock
##  Affiliation: School of Biological Sciences, University of Bristol
##  Last amended: 2026-09
##
##---------
## Purpose
##---------
## Assessed current forest integrity of the regions identified in this study as priority for
## dung beetle conservation.
##
## For two scenarios:
## 
## High-pressure (BAU)
## Low-pressure (GOV)
##
##---------
## Inputs
##---------
##
## flii_SouthAmerica.tif
##
##---------
## Outputs
##---------
##
## Outcomes for high forest integrity proportion per scenario 
## priority area and corresponding maps for Figure 4. 
##
##---------
## Notes
##---------
##
## priority_sensitivity_bau/gov should carry over from previous script 07. 
##
##
########################
## STEP 1 - LOAD DATA ##
########################


# Load Required Packages and TIFs #

library(terra)
library(raster)
library(sf)
library(ggnewscale)


FoInIx <- rast("pathwayto//flii_SouthAmerica.tif")

# Rescale for calculation & removing storage values 
FoInIx[FoInIx == -9999] <- NA
FoInIx <- FoInIx / 1000

# Convert priority area into shp
Pr_In_Bau <- priority_sensitivity_bau
priority_vec_bau <- as.polygons(Pr_In_Bau, dissolve=TRUE)

Pr_In_Gov <- priority_sensitivity_gov
priority_vec_gov <- as.polygons(Pr_In_Gov, dissolve=TRUE)


crs(FoInIx) <- "EPSG:4326"
crs(Pr_In_Bau)
crs(Pr_In_Gov)

FoInIx_cropped_bau <- crop(FoInIx, priority_vec_bau)

FoInIx_cropped_gov <- crop(FoInIx, priority_vec_gov)


FoInIx_masked_bau <- mask(FoInIx_cropped_bau, priority_vec_bau)
FoInIx.BAU <- FoInIx_masked_bau

FoInIx_masked_gov <- mask(FoInIx_cropped_bau, priority_vec_gov)
FoInIx.GOV <- FoInIx_masked_gov

## Rescale and assign -9999 as NA (nodata) ##

FoInIx_Final_BAU <- FoInIx_masked_bau / 1000
FoInIx_Final_GOV <- FoInIx_masked_gov / 1000

## Export 

writeRaster(FoInIx.BAU, "FoInIx_BAU_Amazonia.tif", filetype = "GTiff", overwrite = TRUE)
writeRaster(FoInIx.GOV, "FoInIx_GOV_Amazonia.tif", filetype = "GTiff", overwrite = TRUE)



##################################
## STEP 2 - PLOT AND CATEGORISE ##
##################################

##Convert to dataframe for plotting
FoInIx.BAU_df <- as.data.frame(FoInIx.BAU, xy = TRUE)

FoInIx.GOV_df <- as.data.frame(FoInIx.GOV, xy = TRUE)

##Assign column names correctly 
colnames(FoInIx.BAU_df) <- c("longitude", "latitude", "value") 
colnames(FoInIx.GOV_df) <- c("longitude", "latitude", "value") 

##Combine BAU and GOV dataframes with scenario label
FoInIx.BAU_df$scenario <- "High-pressure (SSP5-8.5)"
FoInIx.GOV_df$scenario <- "Low-pressure (SSP1-2.6)"

FoInIx_combined <- rbind(FoInIx.BAU_df, FoInIx.GOV_df)
FoInIx_combined$scenario <- factor(FoInIx_combined$scenario,
                                   levels = c("High-pressure (SSP5-8.5)",
                                              "Low-pressure (SSP1-2.6)"))
##Plot 

windows()
ggplot() +
  
  # Base Amazon layer
  geom_sf(data = shp_sf, fill = "#3f5b0b", alpha = 0.1, color = NA) +
  
  # ITPA layers (first fill scale round)
  geom_sf(data = Ing_sf, aes(fill = "Indigenous Areas"), alpha = 0.4, color = NA) +
  geom_sf(data = PA_sf,  aes(fill = "Protected Areas"),  alpha = 0.3, color = NA) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      "Indigenous Areas" = "#7fcdbb",
      "Protected Areas"  = "#2c7fb8"
    )
  ) +
  # Reset fill scale prep for FLII
  new_scale_fill() +
  # FLII layer 
  geom_tile(data = FoInIx.GOV_df[!is.na(FoInIx.GOV_df$value), ],
            aes(x = longitude, y = latitude, fill = value),
            alpha = 0.8) +
  
  scale_fill_gradientn(
    colours  = c("#8B0000", "#E34A33", "#FEE08B", "#66BD63", "#006837"),
    name     = "Forest Landscape\nIntegrity Index",
    limits   = c(0, 10),
    breaks   = c(0, 2.5, 5, 7.5, 10),
    labels   = c("0", "2.5", "5", "7.5", "10"),
    guide    = guide_colorbar(
      ticks.colour = "black",
      frame.colour = "black",
      barwidth     = unit(0.4, "cm"),
      barheight    = unit(3.5, "cm")
    )
  ) +
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    style       = north_arrow_fancy_orienteering,
    height      = unit(1.2, "cm"),
    width       = unit(1.2, "cm")
  ) +
  
  annotation_scale(location = "bl", width_hint = 0.3) +
  coord_sf() +
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    axis.text   = element_blank(),
    axis.title  = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Forest Integrity of Conservation Priority Regions")


########################################################################
## Step 3 - CALCULATE PROPORTION OF HIGHLY INTACT FOREST PER SCENARIO ##
########################################################################

##BAU priority area
vals_bau <- values(FoInIx_masked_bau, na.rm = TRUE)
n_below5_bau <- sum(vals_bau < 5, na.rm = TRUE)
n_total_bau  <- length(vals_bau)
prop_below5_bau <- n_below5_bau / n_total_bau

cat("=== BAU ===\n")
cat("Total cells:       ", n_total_bau, "\n")
cat("Cells below 5:     ", n_below5_bau, "\n")
cat("Proportion below 5:", round(prop_below5_bau, 3), "\n")
cat("Mean FLII:         ", round(mean(vals_bau), 3), "\n")
cat("Median FLII:       ", round(median(vals_bau), 3), "\n")
cat("Min FLII:          ", round(min(vals_bau), 3), "\n")
cat("Max FLII:          ", round(max(vals_bau), 3), "\n")

##GOV priority area
vals_gov <- values(FoInIx_masked_gov, na.rm = TRUE)
n_below5_gov <- sum(vals_gov < 5, na.rm = TRUE)
n_total_gov  <- length(vals_gov)
prop_below5_gov <- n_below5_gov / n_total_gov

cat("\n=== GOV ===\n")
cat("Total cells:       ", n_total_gov, "\n")
cat("Cells below 5:     ", n_below5_gov, "\n")
cat("Proportion below 5:", round(prop_below5_gov, 3), "\n")
cat("Mean FLII:         ", round(mean(vals_gov), 3), "\n")
cat("Median FLII:       ", round(median(vals_gov), 3), "\n")
cat("Min FLII:          ", round(min(vals_gov), 3), "\n")
cat("Max FLII:          ", round(max(vals_gov), 3), "\n")


# cat("=== BAU ===\n")
# === BAU ===
#   > cat("Total cells:       ", n_total_bau, "\n")
# Total cells:        3446157 
# > cat("Cells below 5:     ", n_below5_bau, "\n")
# Cells below 5:      220005 
# > cat("Proportion below 5:", round(prop_below5_bau, 3), "\n")
# Proportion below 5: 0.064 
# > cat("Mean FLII:         ", round(mean(vals_bau), 3), "\n")
# Mean FLII:          8.897 
# > cat("Median FLII:       ", round(median(vals_bau), 3), "\n")
# Median FLII:        9.802 
# > cat("Min FLII:          ", round(min(vals_bau), 3), "\n")
# Min FLII:           0 
# > cat("Max FLII:          ", round(max(vals_bau), 3), "\n")
# Max FLII:           10 

#  
# === GOV ===
#   > cat("Total cells:       ", n_total_gov, "\n")
# Total cells:        736257 
# > cat("Cells below 5:     ", n_below5_gov, "\n")
# Cells below 5:      93700 
# > cat("Proportion below 5:", round(prop_below5_gov, 3), "\n")
# Proportion below 5: 0.127 
# > cat("Mean FLII:         ", round(mean(vals_gov), 3), "\n")
# Mean FLII:          8.191 
# > cat("Median FLII:       ", round(median(vals_gov), 3), "\n")
# Median FLII:        9.422 
# > cat("Min FLII:          ", round(min(vals_gov), 3), "\n")
# Min FLII:           0 
# > cat("Max FLII:          ", round(max(vals_gov), 3), "\n")
# Max FLII:           10 
# 

