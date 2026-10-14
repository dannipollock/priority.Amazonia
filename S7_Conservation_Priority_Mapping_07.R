################################################################################
##
## Conservation Priority Mapping 07
##
## Spatial conservation prioritisation for Amazonian dung beetles
##
## Author:       D. Pollock
## Affiliation:  School of Biological Sciences, University of Bristol
## Last amended: 2026-09
##
## Purpose
## --------------------------------
## Identifies spatial conservation priorities across the study basin by
## combining three gridded pressure and value layers under two scenarios:
##
## Inputs
## --------------------------------
##   limiteteste' - Copy.shp                Study basin boundary
##   shapesIUCN/*_EOO_poly.shp              Per-species EOO polygons
##   CMIP6 TX / T / PR rasters              Historical (1981-2010) and
##                                          mid-century (2041-2060) projections
##                                          under SSP1-2.6 and SSP5-8.5
##   lulc_change_bau.tif + gov              Projected FOREST-LOSS change
##   riv_no_grid.shp                        Rivers
##   ANP_*.shp                              PA
##   Tis_TerritoriosIndigenas.shp           IT
##
## 
## Outputs
## --------------------------------
##   overlap.count_sp_eoos.tif
##   2026_priority.eoo.QScore.bau.95.tif / .gov.95.tif
##   2026_priority.threshold.matching_bau_75.tif / _gov_75.tif
##   2026.priority.sensitivity.bau.tif / .gov.tif
##
## File paths are local to the author's machine and require adjustment
##
##############
## PACKAGES ##
##############

library(terra)
library(raster)
library(ggplot2)
library(sf)
library(ggspatial)
library(patchwork)
library(dplyr)


#############################
## STEP 1 - BASIN AND GRID ##
###########################################

#Load Basin
shp <- vect("F:/MRes (BU)/limiteteste' - Copy.shp")
shp <- project(shp_sf, "EPSG:4326")  # Reproject Indigenous territories to PA CRS
shp_sf <- st_read("F:/MRes (BU)/limiteteste' - Copy.shp")

crs(shp)


# Make template grid

template <- rast(
  xmin = bbox["xmin"],
  xmax = bbox["xmax"],
  ymin = bbox["ymin"],
  ymax = bbox["ymax"],
  resolution = 0.08,
  crs = "EPSG:4326"
)
bbox <- st_bbox(shp)


## Create functions
crop_mask <- function(r) mask(crop(r, shp), shp)

diff_calc <- function(fut, hist) fut - hist


################################################################################
## STEP 2 - POTENTIALLY THREATENED SPECIES RANGE OVERLAP
################################################################################

## Filter for only potentially TH species

#Filter threatened species from ConR output
threatened <- opt %>%                                          #Requires opt to still be in environment
  filter(Category_CriteriaB %in% c("CR", "EN", "VU"))


## Copy their EOO shapefiles to a working folder

# --- Step 2: Define your shapefile folder path ---
shp_folder <- "D:/MRes (BU)/03062026/shapesIUCN"  # update this

output_folder  <- "D:/MRes (BU)/03062026/shapesIUCN_TH"   # update this
dir.create(output_folder, showWarnings = FALSE)

extensions <- c(".shp", ".dbf", ".shx", ".prj")

for (sp in threatened$taxa) {
  base_name <- paste0(gsub(" ", "_", sp), "_EOO_poly")
  
  for (ext in extensions) {
    src  <- file.path(shp_folder, paste0(base_name, ext))
    dest <- file.path(output_folder, paste0(base_name, ext))
    
    if (file.exists(src)) {
      file.copy(src, dest)
    } else {
      cat("Missing:", paste0(base_name, ext), "\n")
    }
  }
}

cat("Done! Copied", nrow(threatened), "species shapefiles to", output_folder, "\n")


## Reading and combining the EOO polygons

shp_folder <- "pathwayto//shapesIUCN_TH/"  # change to your path to rename shp_folder to new TH only
shp_files <- list.files(shp_folder, pattern = "EOO_poly.shp$", full.names = TRUE)
shp_list   <- lapply(shp_files, st_read)
all_shapes <- bind_rows(shp_list) %>%
  st_transform(crs = 4326)

all_shapes$poly <- 1

basename(shp_files)


## Rasterising and counting overlapping ranges

pol_vect <- vect(all_shapes)


r_list <- lapply(1:length(pol_vect), function(i){
  if(i %% 1 == 0) cat("Processing polygon", i, "of", length(pol_vect), "\n")
  rasterize(pol_vect[i], template, field = 1, background = 0)
})

stacked       <- rast(r_list)
overlap_count <- sum(stacked)
overlap_count <- mask(overlap_count, vect(AmazonG))

# Rename layer
names(overlap_count) <- "ovrl___"

# Plot
overlap_df <- as.data.frame(overlap_count, xy = TRUE) %>%
  filter(!is.na(ovrl___) & ovrl___ > 0)# remove NA cells

# Plot Overlaps for Figure 3 and create object for prioritisation
windows()
ggplot() +
  geom_sf(data      = AmazonG,
          fill      = "#EDF2ED", colour = NA) +
  geom_tile(data = overlap_df,
            aes(x = x, y = y, fill = ovrl___)) +
  scale_fill_gradientn(
    colours = c("#EDF2ED", "#F5E6D3", "#F5C9A0", "#F0986A", "#D95F3B", "#A82418", "#6B0000"),
    name    = "Species\noverlap",
    guide   = guide_colorbar(ticks.colour = "black", frame.colour = "black")
  )+
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 8,  colour = "black", face = "bold"),
    legend.text       = element_text(size = 7,  colour = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    plot.margin       = margin(5, 5, 5, 5, "mm"),
    panel.grid        = element_blank()   # <-- add this
  ) +
  labs(title = NULL)


writeRaster(overlap_count,
            "D:/MRes (BU)/Final Species Dataset and Assessment/overlap.count_sp_eoos.tif",
            overwrite = TRUE)

overlap_poly <- as.polygons(overlap_count, dissolve = FALSE)
overlap_poly <- overlap_poly[overlap_poly$ovrl___ > 0, ]


# #check
print(global(overlap_count, range, na.rm = TRUE))


eoo.overlap.data <- mask(overlap_count, shp)
global(eoo.overlap.data, range, na.rm = TRUE)
crs(eoo.overlap.data) # CHECK CRS MATCHES 4326
values(template) <- NA

#crop/mask to exactly match shapefile boundary
template <- mask(template, shp)
crs(template)


####################################
## STEP 3 - CLIMATE CHANGE LAYERS ##
####################################

#############
## Climate ##
#############

#TX
hist.tasmax <- rast("F:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/HIST/Raw Files/CMIP6 - Maximum temperature (TX)  deg C - 1981-2010.tif")

f.tasmax.gov <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/SSP126.F/Raw File/CMIP6 - Maximum temperature (TX)  deg C - Medium Term (2041-2060) SSP1-2.6.tif")

f.tasmax.bau <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/SSP585.F/Raw Files/CMIP6 - Maximum temperature (TX)  deg C - Medium Term (2041-2060) SSP5-8.5.tif")



#Tas
hist.tas <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/HIST/Raw Files/CMIP6 - Mean temperature (T)  deg C - 1981-2010.tif")

f.tas.bau <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/SSP585.F/Raw Files/CMIP6 - Mean temperature (T)  deg C - Medium Term (2041-2060) SSP5-8.5.tif")

f.tas.gov <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/SSP126.F/Raw File/CMIP6 - Mean temperature (T)  deg C - Medium Term (2041-2060) SSP1-2.6.tif")



#PR
hist.pr <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/HIST/Raw Files/CMIP6 - Total precipitation (PR)  mm_day - 1981-2010.tif")

f.pr.gov <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/SSP126.F/Raw File/CMIP6 - Total precipitation (PR)  mm_day - Medium Term (2041-2060) SSP1-2.6.tif")

f.pr.bau <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Climate Change Parameters/Final Clim Data and Outputs/SSP585.F/Raw Files/CMIP6 - Total precipitation (PR)  mm_day - Medium Term (2041-2060) SSP5-8.5.tif")


diff_calc <- function(fut, hist) fut - hist

crop_mask <- function(r) mask(crop(r, shp), shp)

hist.tasmax <- crop_mask(hist.tasmax)
f.tasmax.gov <- crop_mask(f.tasmax.gov)
f.tasmax.bau <- crop_mask(f.tasmax.bau)

hist.tas <- crop_mask(hist.tas)
f.tas.gov <- crop_mask(f.tas.gov)
f.tas.bau <- crop_mask(f.tas.bau)

hist.pr <- crop_mask(hist.pr)
f.pr.gov <- crop_mask(f.pr.gov)
f.pr.bau <- crop_mask(f.pr.bau)


## Change layers

## Calculate change layers
bau <- list(
  tas = diff_calc(f.tas.bau, hist.tas),
  tasmax = diff_calc(f.tasmax.bau, hist.tasmax),
  pr = diff_calc(f.pr.bau, hist.pr)
)

gov <- list(
  tas = diff_calc(f.tas.gov, hist.tas),
  tasmax = diff_calc(f.tasmax.gov, hist.tasmax),
  pr = diff_calc(f.pr.gov, hist.pr)
)


## Composite climate change index

## PR is reversed so that drying contributes positively like the other paramters 

clim_index <- function(lst){
  ( lst$tas +
      lst$tasmax +
      (-lst$pr) ) / 3
}
### Climate index
clim_bau <- clim_index(bau)
clim_gov <- clim_index(gov)


###################################
## STEP 4 - DEFORESTATION LAYERS ##
###################################

### DEFORESTATION

lulc_change_bau <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Final Dataset (IndVal+EOO.AOO)/Climate Change Parameters and Outputs/All Species Climate and LULC Change/Deforestation/lulc_change_bau.tif")
lulc_change_gov <- rast("D:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Final Dataset (IndVal+EOO.AOO)/Climate Change Parameters and Outputs/All Species Climate and LULC Change/Deforestation/lulc_change_gov.tif")

lulc_change_bau[lulc_change_bau==-1] <- 1
lulc_change_gov[lulc_change_gov==-1] <- 1

defor_bau <- lulc_change_bau
defor_gov <- lulc_change_gov


global(defor_bau, range, na.rm=TRUE) # check range


## Aggregating to proportional forest loss

## calculate aggregation factor (fine → coarse)
fact_x <- res(template)[1] / res(defor_bau)[1]
fact_y <- res(template)[2] / res(defor_gov)[2]

fact <- round(mean(c(fact_x, fact_y)))

print(fact) ## should be 9

### aggregate to proportion forest loss
defor_bau <- aggregate(defor_bau, fact=fact, fun=mean, na.rm=TRUE)
defor_gov <- aggregate(defor_gov, fact=fact, fun=mean, na.rm=TRUE)

crs(defor_gov)
global(defor_gov, range, na.rm=TRUE) # check range

# check it is now non-binary
summary(defor_bau)
summary(defor_gov)


## Aligning 

defor_bau <- resample(defor_bau, template, method="bilinear")
defor_gov <- resample(defor_gov, template, method="bilinear")

global(defor_gov, range, na.rm=TRUE) # check range
global(defor_bau, range, na.rm=TRUE) # check range


#################################################
## STEP 5 - STANDARDISING THE COMPONENT LAYERS ##
#################################################

## Standardise (Min/ Max Quantile)

robust_scale <- function(r){
  q <- global(r, quantile, probs=c(0.05,0.95), na.rm=TRUE) #95 to retain higher values for clim
  (r - q[1,1]) / (q[1,2] - q[1,1])
}

clamp01 <- function(r){
  clamp(r, lower=0, upper=1)
}


## Rescaling climate and range overlap

clim_s_bau <- clamp01(robust_scale(clim_bau))
clim_s_gov <- clamp01(robust_scale(clim_gov))
rich_s <- clamp01(robust_scale(overlap_count))

clim_s_bau <- resample(clim_s_bau, rich_s, method="bilinear")
clim_s_gov <- resample(clim_s_gov, rich_s, method="bilinear")


# z-scores measure deviation from mean NOT APPROPRIATE
#
# conservation indices need relative magnitude
#
# min–max scaling preserves spatial gradients

# deforestation is already proportional (0–1), so just rename
defor_s_bau <- defor_bau
defor_s_gov <- defor_gov


#############################################################
## STEP 6 - COMPOSITE INDEX AND 95TH PERCENTILE PRIORITIES ##
#############################################################
## Composite Index Generation 

w_rich <- 0.3
w_def  <- 0.3
w_clim <- 0.3


compound_bau <- w_rich*rich_s + w_def*defor_s_bau + w_clim*clim_s_bau
compound_gov <- w_rich*rich_s + w_def*defor_s_gov + w_clim*clim_s_gov


#calculated by pooling values and taking percentiles from that#

vals_all <- c(values(compound_bau), values(compound_gov))

thresh <- quantile(
  vals_all[!is.na(vals_all) & vals_all > 0],
  0.95
)

priority_bau <- ifel(compound_bau >= thresh,1,NA)
priority_gov <- ifel(compound_gov >= thresh,1,NA)

priority95_bau <- mask(priority_bau, shp)
priority95_gov <- mask(priority_gov, shp)
windows(); plot(priority95_gov)
windows(); plot(priority95_bau)


## Export ##

writeRaster(priority_bau,"2026_priority.eoo.QScore.bau.95.tif",overwrite=TRUE)
writeRaster(priority_gov,"2026_priority.eoo.QScore.gov.95.tif",overwrite=TRUE)

getwd()


##############################################################
## STEP 7 - COINCIDENT HOTSPOTS (75TH PERCENTILE PER LAYER) ##
##############################################################


top25_raster <- function(r){
  vals <- values(r)
  vals <- vals[!is.na(vals) & vals > 0]
  q75  <- quantile(vals, 0.75, na.rm = TRUE)
  ifel(r >= q75, 1, NA)
}

# Now apply top25 to aligned rasters
rich_top      <- top25_raster(overlap_count)
defor_top_bau <- top25_raster(defor_bau)
defor_top_gov <- top25_raster(defor_gov)
clim_top_bau  <- top25_raster(clim_s_bau)
clim_top_gov  <- top25_raster(clim_s_gov)

priority75_bau <- ifel(rich_top & defor_top_bau & clim_top_bau, 1, NA)
priority75_gov <- ifel(rich_top & defor_top_gov & clim_top_gov, 1, NA)

windows(); plot(priority75_bau)
windows(); plot(priority75_gov)

# Confirm
# global(overlap_count, range, na.rm = TRUE)
# windows();plot(rich_top, col = viridis::viridis(50))

# Sanity check
ext(rich_top)
ext(defor_top_bau)
ext(clim_top_bau)

res(rich_top)
res(defor_top_bau)
res(clim_top_bau)

crs(rich_top)
crs(defor_top_bau)
crs(clim_top_bau)



# Export
writeRaster(priority75_bau,"2026_priority.threshold.matching_bau_75.tif",overwrite=TRUE)
writeRaster(priority75_gov,"2026_priority.threshold.matching_gov_75.tif",overwrite=TRUE)


############################################
## STEP 8 - ITPA's AND Method Sensitivity ##
############################################

## Generate the Final Plots + ITPA Calculations ##

## Priority rasters as data frames

priority75_bau_df <- as.data.frame(priority75_bau, xy = TRUE)
priority75_gov_df <- as.data.frame(priority75_gov, xy = TRUE)

priority95_bau_df <- as.data.frame(priority95_bau, xy = TRUE)
priority95_gov_df <- as.data.frame(priority95_gov, xy = TRUE)

# load rvr shp

rvr <- vect("F:/MRes (BU)/riv_no_grid/riv_no_grid.shp")
rvr <- crop(rvr, shp)
rvr <- mask(rvr, shp)
rvr <- project(rvr, shp)

rvr_sf <- st_as_sf(rvr)

rvr_sf <- st_make_valid(rvr_sf)


## Building ITPA layer if needed (not withstanding from 02)

## Create ITPA shp

PA1 <- vect("F:/MRes (BU)/Anps_jun2025/ANP_ReservaFlorestal.shp")
PA2 <- vect("F:/MRes (BU)/Anps_jun2025/ANP_departamental.shp")
PA3 <- vect("F:/MRes (BU)/Anps_jun2025/ANP_nacional.shp")
PA4 <- vect("F:/MRes (BU)/Anps_jun2025/ANP_BosqueProtector.shp")

Ing <- vect("F:/MRes (BU)/Tis_Junho2025/Tis_TerritoriosIndigenas.shp")

PA1 <- project(PA1, shp)
PA2 <- project(PA2, shp)
PA3 <- project(PA3, shp)
PA4 <- project(PA4, shp)

Ing <- project(Ing, shp)

PA <- rbind(PA1, PA2, PA3, PA4)

PA <- crop(PA, shp)
Ing <- crop(Ing, shp)

PA <- mask(PA, shp)
Ing <- mask(Ing, shp)


#Add Label PAs or only Indigenous
PA$type <- rep("PA", nrow(PA))

Ing$type <- rep("Indigenous", nrow(Ing))  # <-- use $ assignment

ITPA <- rbind(PA, Ing)

names(ITPA)

cols <- ifelse(ITPA$type == "PA", "blue", "purple")
windows();plot(ITPA, col = cols, border = NA)


## Reload vectors

writeVector(ITPA, "ITPA_final.shp", overwrite=TRUE)
writeVector(PA, "PA_final.shp", overwrite=TRUE)
writeVector(Ing, "Ing_final.shp", overwrite=TRUE)

PA <- vect("F:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Final Dataset (IndVal+EOO.AOO)/Final Modelling Data and Scripts/Hotspot Generation (Clim + Defor)/2026/PA_final.shp")
Ing <- vect("F:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Final Dataset (IndVal+EOO.AOO)/Final Modelling Data and Scripts/Hotspot Generation (Clim + Defor)/2026/Ing_final.shp")
ITPA <- vect("F:/MRes (BU)/Final Species Dataset and Assessment/Results Biogeographical/Final Dataset (IndVal+EOO.AOO)/Final Modelling Data and Scripts/Hotspot Generation (Clim + Defor)/2026/ITPA_final.shp")


# Convert for plotting
ITPA_sf <- st_as_sf(ITPA)

Ing_sf <- st_as_sf(Ing)

PA_sf <- st_as_sf(PA)


## Four-panel method and scenario comparison

## Plots ##

map_theme <- theme_void(base_size = 9) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(size = 9, face = "bold", colour = "black",
                                margin = margin(b = 3)),
    plot.margin  = margin(2, 2, 2, 2, "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    legend.position  = "bottom",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3, "mm"),
    legend.title     = element_blank()
  )

# Panel A: Coincident Hotspots — Scenario Comparison
pA <- ggplot() +
  geom_sf(data = shp_sf, fill = "#3f5b0b", alpha = 0.1, colour = NA) +
  geom_tile(data = priority75_bau_df, aes(x = x, y = y, fill = "High-pressure (SSP5-8.5)"), alpha = 0.6) +
  geom_tile(data = priority75_gov_df, aes(x = x, y = y, fill = "Low-pressure (SSP1-2.6)"),  alpha = 0.6) +
  scale_fill_manual(values = c("High-pressure (SSP5-8.5)" = "brown3",
                               "Low-pressure (SSP1-2.6)"  = "seagreen4")) +
  coord_sf(expand = FALSE) +
  labs(title = "A) Coincident Hotspots — Scenario Comparison") +
  map_theme

# Panel B: Composite Index — Scenario Comparison
pB <- ggplot() +
  geom_sf(data = shp_sf, fill = "#3f5b0b", alpha = 0.1, colour = NA) +
  geom_tile(data = priority95_bau_df, aes(x = x, y = y, fill = "High-pressure (SSP5-8.5)"), alpha = 0.6) +
  geom_tile(data = priority95_gov_df, aes(x = x, y = y, fill = "Low-pressure (SSP1-2.6)"),  alpha = 0.6) +
  scale_fill_manual(values = c("High-pressure (SSP5-8.5)" = "brown3",
                               "Low-pressure (SSP1-2.6)"  = "seagreen4")) +
  coord_sf(expand = FALSE) +
  labs(title = "B) Composite Index — Scenario Comparison") +
  map_theme

# Panel C: Method Comparison — Low-pressure
pC <- ggplot() +
  geom_sf(data = shp_sf, fill = "#3f5b0b", alpha = 0.1, colour = NA) +
  geom_tile(data = priority95_gov_df, aes(x = x, y = y, fill = "Composite (95th percentile)"),          alpha = 0.6) +
  geom_tile(data = priority75_gov_df, aes(x = x, y = y, fill = "Threshold Matching (75th percentile)"), alpha = 0.6) +
  scale_fill_manual(values = c("Composite (95th percentile)"          = "darkorange",
                               "Threshold Matching (75th percentile)" = "plum3")) +
  coord_sf(expand = FALSE) +
  labs(title = "C) Method Comparison — Low-pressure [SSP1-2.6]") +
  map_theme

# Panel D: Method Comparison — High-pressure
pD <- ggplot() +
  geom_sf(data = shp_sf, fill = "#3f5b0b", alpha = 0.1, colour = NA) +
  geom_tile(data = priority95_bau_df, aes(x = x, y = y, fill = "Composite (95th percentile)"),          alpha = 0.6) +
  geom_tile(data = priority75_bau_df, aes(x = x, y = y, fill = "Threshold Matching (75th percentile)"), alpha = 0.6) +
  scale_fill_manual(values = c("Composite (95th percentile)"          = "darkorange",
                               "Threshold Matching (75th percentile)" = "plum3")) +
  coord_sf(expand = FALSE) +
  labs(title = "D) Method Comparison — High-pressure [SSP5-8.5]") +
  map_theme


final_plot <- (pA | pB) / (pC | pD)

windows();final_plot


################################################################################
## STEP 9 - AGREEMENT BETWEEN PRIORITISATION METHODS
################################################################################

## Generate Final Priority Areas --> Sensitivity

### Calculate overlap between methods first ###
## BAU

priority_sensitivity_bau <- priority95_bau & priority75_bau

n_overlap_bau  <- global(priority_sensitivity_bau, "sum", na.rm = TRUE)
n1_bau         <- global(priority95_bau > 0, "sum", na.rm = TRUE)
n2_bau         <- global(priority75_bau > 0, "sum", na.rm = TRUE)

jaccard_bau    <- n_overlap_bau / (n1_bau + n2_bau - n_overlap_bau)
prop_in_95_bau <- n_overlap_bau / n1_bau
prop_in_75_bau <- n_overlap_bau / n2_bau

cat("=== BAU Sensitivity Overlap ===\n")
cat("Cells in priority95:        ", n1_bau$sum, "\n")

cat("Cells in priority75:        ", n2_bau$sum, "\n")

cat("Overlapping cells:          ", n_overlap_bau$sum, "\n")

cat("Jaccard index:              ", round(jaccard_bau$sum, 3), "\n")

cat("% of priority95 in overlap: ", round(prop_in_95_bau$sum * 100, 1), "%\n")

cat("% of priority75 in overlap: ", round(prop_in_75_bau$sum * 100, 1), "%\n")



## GOV

priority_sensitivity_gov <- priority95_gov & priority75_gov

n_overlap_gov  <- global(priority_sensitivity_gov, "sum", na.rm = TRUE)
n1_gov         <- global(priority95_gov > 0, "sum", na.rm = TRUE)
n2_gov         <- global(priority75_gov > 0, "sum", na.rm = TRUE)

jaccard_gov    <- n_overlap_gov / (n1_gov + n2_gov - n_overlap_gov)
prop_in_95_gov <- n_overlap_gov / n1_gov
prop_in_75_gov <- n_overlap_gov / n2_gov

cat("\n=== GOV Sensitivity Overlap ===\n")
cat("Cells in priority95:        ", n1_gov$sum, "\n")

cat("Cells in priority75:        ", n2_gov$sum, "\n")

cat("Overlapping cells:          ", n_overlap_gov$sum, "\n")


cat("Jaccard index:              ", round(jaccard_gov$sum, 3), "\n")


cat("% of priority95 in overlap: ", round(prop_in_95_gov$sum * 100, 1), "%\n")


cat("% of priority75 in overlap: ", round(prop_in_75_gov$sum * 100, 1), "%\n")



## Cells agreed by both methods and both scenarios

priority_sensitivity_bau.gov <- priority_sensitivity_bau & priority_sensitivity_gov


## GENERATE FINAL RASTERS ##

priority_sensitivity_gov <- ifel(priority_sensitivity_gov, 1, NA)
priority_sensitivity_bau <- ifel(priority_sensitivity_bau, 1, NA)
priority_sensitivity_bau.gov <- ifel(priority_sensitivity_bau.gov, 1, NA)


writeRaster(priority_sensitivity_bau,"2026.priority.sensitivity.bau.tif",overwrite=TRUE)
writeRaster(priority_sensitivity_gov,"2026.priority.sensitivity.gov.tif",overwrite=TRUE)

# Count priority cells
n_cells_bau <- global(priority_sensitivity_bau, "sum", na.rm = TRUE)
n_cells_gov <- global(priority_sensitivity_gov, "sum", na.rm = TRUE)

cat("Priority cells BAU:", n_cells_bau$sum, "\n")


cat("Priority cells GOV:", n_cells_gov$sum, "\n")



#########################################################################
## STEP 11 - PRIORITY AREA OVERLAP WITH PROTECTED AND INDIGENOUS LANDS ##
#########################################################################

####################################################
### Analysis of % in PA's and Indigenous Regions ###
####################################################

# Rasterise 
pa_r  <- rasterize(PA, priority_sensitivity_bau, field=1)
ind_r <- rasterize(Ing, priority_sensitivity_bau, field=1)



#--------------------------
# PA overlap GOV
#--------------------------
priority_pa_mask <- mask(priority_sensitivity_gov, pa_r)
# Percent of priority inside PA
percent_pa <- global(priority_pa_mask, "sum", na.rm=TRUE) /
  global(priority_sensitivity_gov, "sum", na.rm=TRUE) * 100
print(percent_pa)




#--------------------------
# Indigenous overlap GOV
#--------------------------
priority_ind_mask <- mask(priority_sensitivity_gov, ind_r)
# Percent of priority inside PA
percent_ind <- global(priority_ind_mask, "sum", na.rm=TRUE) /
  global(priority_sensitivity_gov, "sum", na.rm=TRUE) * 100
percent_ind



#--------------------------
# PA overlap BAU
#--------------------------
priority_pa_mask_bau <- mask(priority_sensitivity_bau, pa_r)
# Percent of priority inside PA
percent_pa_bau <- global(priority_pa_mask_bau, "sum", na.rm=TRUE) /
  global(priority_sensitivity_bau, "sum", na.rm=TRUE) * 100
percent_pa_bau


#--------------------------
# Indigenous overlap BAU
#--------------------------
priority_ind_mask_bau <- mask(priority_sensitivity_bau, ind_r)
# Percent of priority inside PA
percent_ind_bau <- global(priority_ind_mask_bau, "sum", na.rm=TRUE) /
  global(priority_sensitivity_bau, "sum", na.rm=TRUE) * 100
percent_ind_bau


#--------------------------
# PA + Indigenous overlap GOV
#--------------------------
pa_ind_mask <- app(c(pa_r, ind_r), fun=function(x) ifelse(x[1]==1 & x[2]==1, 1, NA))
priority_pa_ind_mask <- mask(priority_sensitivity_gov, pa_ind_mask)

percent_pa_ind <- global(priority_pa_ind_mask, "sum", na.rm=TRUE) /
  global(priority_sensitivity_gov, "sum", na.rm=TRUE) * 100
percent_pa_ind

#--------------------------
# PA + Indigenous overlap BAU
#--------------------------
pa_ind_mask_bau <- app(c(pa_r, ind_r), fun=function(x) ifelse(x[1]==1 & x[2]==1, 1, NA))
priority_pa_ind_mask_bau <- mask(priority_sensitivity_bau, pa_ind_mask_bau)

percent_pa_ind_bau <- global(priority_pa_ind_mask_bau, "sum", na.rm=TRUE) /
  global(priority_sensitivity_bau, "sum", na.rm=TRUE) * 100
percent_pa_ind_bau


#######################################
## STEP 10 - FINAL PRIORITY AREA MAP ##
#######################################

priority_sensitivity_bau <- rast("F:/MRes (BU)/03062026/2026.priority.sensitivity.bau.tif")
priority_sensitivity_gov <- rast("F:/MRes (BU)/03062026/2026.priority.sensitivity.gov.tif")


# Convert raster to dataframe
priority_gov_df <- as.data.frame(priority_sensitivity_gov, xy = TRUE)
priority_bau_df <- as.data.frame(priority_sensitivity_bau, xy = TRUE)


# Ensure column name for fill is consistent
names(priority_gov_df)[3] <- "value"
names(priority_bau_df)[3] <- "value"
names(priority_gov_df)


## PLOT

pmain <- ggplot() +
  geom_sf(
    data = shp_sf,
    fill = "#3f5b0b",
    alpha = 0.1,
    color = NA
  ) +
  
  geom_sf(data = rvr_sf, color = "darkblue", linewidth = 0.1, alpha = 0.5) +
  
  
  geom_sf(data = Ing_sf,
          aes(fill = "Indigenous Areas"),
          alpha = 0.4, color = NA) +
  
  geom_sf(data = PA_sf,
          aes(fill = "Protected Areas"),
          alpha = 0.3, color = NA) +

  geom_tile(data = priority_gov_df,
            aes(x = x, y = y, fill = "Priority Areas (Low-pressure [SSP1-2.6])"),
            alpha = 0.3) +
  
  geom_tile(data = priority_bau_df,
            aes(x = x, y = y, fill = "Priority Areas (High-pressure [SSP5-8.5])"),
            alpha = 0.3) +
  
  scale_fill_manual(
    values = c(
      "Indigenous Areas" = "#7fcdbb",
      "Protected Areas" = "#2c7fb8",
      "Priority Areas (Low-pressure [SSP1-2.6])" = "firebrick3",
      "Priority Areas (High-pressure [SSP5-8.5])" = "firebrick3"
    )
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm")
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.3
  ) +
  
  coord_sf() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
  ) +
  labs(title = "Priority Regions for Dung Beetle Conservation")

