################################################################################
##
## Climate + Deforestation Clipping, Change, and Overlap Calculations 
##
## Species-level climatic and forest-loss exposure for Amazonian dung beetles
##
## Author:       D. Pollock
## Affiliation:  School of Biological Sciences, University of Bristol
## Last revised: 2026-09
##
## --------
## Purpose
## --------
## Species-specific exposure to projected climatic and forest-loss within its extent of occurrence:
##
##   1  EOO shapefiles
##   2  Clipping climate surfaces to each species EOO
##   3  Kernel density overlap between baseline and projected climate
##   4  Magnitude of projected climatic change per species
##   5  Percentage of each EOO subject to projected deforestation
##
## Produces the necessary components to calculate risk index.
##
## ---------------------------------------------------------------------------
## Inputs
## ---------------------------------------------------------------------------
##   Data carries over (as 'opt') from S3_ConR_Eval_DungBeetles_All_Merged_03
##
##   shapesIUCN/         -  per-species EOO polygons from ConR::IUCN.eval()
##   Climate_FINAL_DATA/ -  IPCC climate surfaces (TAS, TASMAX, PR) per
##                          scenario
##   amazonbau2002.tif, amazonbau2050.tif,
##   amazongov2002.tif, amazongov2050.tif (Soares-Filho et al., 2006)
##                         
##   limiteteste' - Copy.shp   Amazon basin boundary
##
## ---------------------------------------------------------------------------
## Outputs
## ---------------------------------------------------------------------------
##   Per-species clipped climate rasters (one directory per parameter/scenario)
##   EOO_species_areas_km2.csv
##   Per-species density-curve plots
##   results_{tas,tasmax,pr}_{BAU,GOV}.csv          climate overlap
##   raw_change_{TAS,TASMAX,PR}_{BAU,GOV}.csv       climate change magnitude
##   top_change_BAU.csv, top_change_GOV.csv         combined change, ranked
##   lulc_change_bau.tif, lulc_change_gov.tif
##   individual_species_deforestation_summary.csv
##
## ---------------------------------------------------------------------------
## IMPORTANT - manual configuration
## ---------------------------------------------------------------------------
## Section 2 is run once per climate parameter (TAS, TASMAX, PR) and per
## scenario (HIST, BAU, GOV). Two lines must be edited between runs:
##
##   the input raster path   (section 2.2)
##   output_raster_dir       (section 2.4)
##
## These are not parameterised. Failing to update both together will write one
## scenario's values into another scenario's output directory.
##
################################################################################
## PACKAGES AND TERRA CONFIGURATION
################################################################################

library(raster)
library(terra)
library(tools)
library(sf)
library(readxl)

################################################################################
## STEP 1 - SELECTING AND COPYING PER-SPECIES EOO SHAPEFILES
################################################################################

## Selecting species from the ConR assessment

threatened <- opt %>%
  filter(Category_CriteriaB %in% c("CR", "EN", "VU", "LC or NT"))


## Source and destination directories

shp_folder <- "D:/MRes (BU)/03062026/shapesIUCN"

output_folder <- "D:/MRes (BU)/03062026/shapesIUCN_all"
dir.create(output_folder, showWarnings = FALSE)

extensions <- c(".shp", ".dbf", ".shx", ".prj")


## Copying each shapefile component per species

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


################################################################################
## STEP 2 - CLIPPING CLIMATE SURFACES TO EACH SPECIES EOO
################################################################################

## Run once per climate parameter and scenario; see the note in the header
## regarding the two paths that must be edited between runs.


## ---- 2.1 Reference grid ---------------------------------------------------

## The land-use change layer defines the target resolution to which climate
## surfaces are resampled.

defor <- rast("F:/MRes (BU)/03062026/Species Climate and LULC/lulc_change_bau.tif")

ref_raster <- rast(ext(defor), resolution = 0.0089278, crs = crs(defor))


## ---- 2.2 Study-area boundary and climate surface --------------------------

boundary <- vect("F:/MRes (BU)/limiteteste' - Copy.shp")
boundary <- project(boundary, "EPSG:4326")

## EDIT PER PARAMETER AND SCENARIO
climate_parameter <- terra::rast("F:/climate/CHELSA_FINAL_DATA/tasmax_baseline_annmean_1981-2010.tif")


## Resampling and masking to the study area ##

climate_parameter <- terra::resample(climate_parameter, ref_raster,
                                     method = "bilinear")
climate_parameter <- terra::mask(climate_parameter, boundary)


## Output ##

shapefile_dir <- "F:/MRes (BU)/03062026/shapesIUCN_all"

eoo_files <- list.files(shapefile_dir, pattern = "_EOO_.*\\.shp$",
                        full.names = TRUE)

## ! EDIT PER PARAMETER AND SCENARIO ! ##
output_raster_dir <- "F:/climate/CHELSA_FINAL_DATA/HIST/TASMAX"
dir.create(output_raster_dir, showWarnings = FALSE)


## Clipping loop

for (shp in eoo_files) {
  
  message("Processing: ", basename(shp))
  
  ## Read and repair polygon geometry
  eoo_poly <- st_read(shp, quiet = TRUE)
  eoo_poly <- st_make_valid(eoo_poly)
  
  ## Reproject with sf before converting to a terra object
  target_crs <- st_crs(crs(climate_parameter))
  if (st_crs(eoo_poly) != target_crs) {
    eoo_poly <- st_transform(eoo_poly, target_crs)
  }
  
  eoo_vect <- vect(eoo_poly)
  
  ## Clip the climate surface to the polygon
  clipped <- mask(crop(climate_parameter, eoo_vect), eoo_vect)
  
  ## Write output
  species_name <- gsub("\\.shp$", "", basename(shp))
  output_file  <- file.path(output_raster_dir,
                            paste0(species_name, "_HIST_TASMAX.tif"))
  
  writeRaster(clipped, output_file, overwrite = TRUE)
  
  message("Saved: ", output_file)
}

## Note: this loop may write corrupt TIFFs. Outputs must be
## checked and the loop re-run for any affected species.


########################################
## STEP 3 - EXTENT OF OCCURRENCE AREA ##
########################################

library(sf)
library(tools)

eoo_files <- list.files(
  "D:/MRes (BU)/03062026/shapesIUCN_all/",
  pattern    = "\\.shp$",
  full.names = TRUE
)

results_area <- data.frame(
  species          = character(),
  area_km2         = numeric(),
  stringsAsFactors = FALSE
)

for (shp in eoo_files) {
  
  message("Calculating area: ", basename(shp))
  
  eoo_poly <- st_read(shp, quiet = TRUE)
  eoo_poly <- st_make_valid(eoo_poly)
  
  ## Geodesic area; sf computes spherical area for geographic coordinates
  area_m2  <- st_area(eoo_poly)
  area_km2 <- as.numeric(area_m2) / 1e6
  
  species_name <- gsub("\\.shp$", "", basename(shp))
  
  results_area <- rbind(
    results_area,
    data.frame(
      species  = species_name,
      area_km2 = area_km2
    )
  )
}

write.csv(results_area, "EOO_species_areas_km2.csv", row.names = FALSE)


############################################################################
## STEP 4 - CLIMATE RANGE OVERLAP BETWEEN BASELINE AND PROJECTED CLIMATE ##
############################################################################

library(terra)
library(ggplot2)
library(dplyr)


## Kernel density overlap

## Overlap is the kernel density estimates on a common grid. Values range from 0 (no climate overlap) 
#  to 1 (identical distributions).

calc_density_overlap <- function(vals1, vals2, n = 2000) {
  d1 <- density(vals1, na.rm = TRUE)
  d2 <- density(vals2, na.rm = TRUE)
  
  grid_x <- seq(min(d1$x, d2$x), max(d1$x, d2$x), length.out = n)
  y1 <- approx(d1$x, d1$y, xout = grid_x)$y
  y2 <- approx(d2$x, d2$y, xout = grid_x)$y
  y1[is.na(y1)] <- 0
  y2[is.na(y2)] <- 0
  sum(pmin(y1, y2)) * (grid_x[2] - grid_x[1])
}


## Density-curve plots

outdir <- "F:/climate/CHELSA_FINAL_DATA/desnity curve"

plot_smooth_density <- function(vals_hist, vals_fut, species, variable,
                                scenario, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  df <- data.frame(
    value    = c(vals_hist, vals_fut),
    scenario = rep(c("HIST", scenario),
                   c(length(vals_hist), length(vals_fut)))
  )
  
  medians <- df %>%
    group_by(scenario) %>%
    summarize(median_val = median(value, na.rm = TRUE))
  
  outfile <- file.path(outdir, paste0(species, "_", variable, "_", scenario,
                                      "_smooth_density.png"))
  
  p <- ggplot(df, aes(x = value, color = scenario, fill = scenario)) +
    geom_density(alpha = 0.3, size = 1.2) +
    geom_vline(data = medians, aes(xintercept = median_val, color = scenario),
               linetype = "dashed", size = 1) +
    scale_color_manual(values = c("HIST" = "#1f78b4", scenario = "#ff7f00")) +
    scale_fill_manual(values = c("HIST" = "#a6cee3", scenario = "#fdbf6f")) +
    labs(title    = species,
         subtitle = paste(variable, "- HIST vs", scenario),
         x        = variable,
         y        = "Density") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top",
          plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(face = "italic")) +
    guides(color = guide_legend(title = NULL), fill = "none")
  
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 200)
}


## Overlap across species

run_overlap <- function(hist_files, fut_files, sp_hist, sp_fut,
                        variable, scenario, plot_dir) {
  common_sp <- intersect(sp_hist, sp_fut)
  out <- list()
  
  for (i in seq_along(common_sp)) {
    sp <- common_sp[i]
    
    hist_r <- rast(hist_files[sp_hist == sp])
    fut_r  <- rast(fut_files[sp_fut == sp])
    
    v1 <- values(hist_r, na.rm = TRUE)
    v2 <- values(fut_r,  na.rm = TRUE)
    
    ov <- calc_density_overlap(v1, v2)
    
    plot_smooth_density(v1, v2, sp, variable, scenario, plot_dir)
    
    out[[i]] <- data.frame(species  = sp,
                           variable = variable,
                           scenario = scenario,
                           overlap  = ov)
    
    cat("Completed:", sp, variable, scenario, "- overlap:", ov, "/n")
  }
  
  do.call(rbind, out)
}


## Input directories of species clipped cliamte rasters

hist_tas_dir    <- "pathwayto//HIST/TAS/"
hist_tasmax_dir <- "pathwayto//HIST/TASMAX/"
hist_pr_dir     <- "pathwayto//HIST/PR/"

fut_tas_BAU_dir    <- "pathwayto//BAU/TAS/"
fut_tasmax_BAU_dir <- "pathwayto//BAU/TASMAX/"
fut_pr_BAU_dir     <- "F:/climate/CHELSA_FINAL_DATA/BAU/PR/"

fut_tas_GOV_dir    <- "pathwayto//GOV/TAS/"
fut_tasmax_GOV_dir <- "pathwayto//GOV/TASMAX/"
fut_pr_GOV_dir     <- "pathwayto//GOV/PR/"


## Listing per-species rasters

hist_tas_files    <- list.files(hist_tas_dir,    pattern = "_HIST_TAS.tif$",    full.names = TRUE)
hist_tasmax_files <- list.files(hist_tasmax_dir, pattern = "_HIST_TASMAX.tif$", full.names = TRUE)
hist_pr_files     <- list.files(hist_pr_dir,     pattern = "_HIST_PR.tif$",     full.names = TRUE)

fut_tas_bau_files    <- list.files(fut_tas_BAU_dir,    pattern = "_BAU_TAS.tif$",    full.names = TRUE)
fut_tasmax_bau_files <- list.files(fut_tasmax_BAU_dir, pattern = "_BAU_TASMAX.tif$", full.names = TRUE)
fut_pr_bau_files     <- list.files(fut_pr_BAU_dir,     pattern = "_BAU_PR.tif$",     full.names = TRUE)

fut_tas_GOV_files    <- list.files(fut_tas_GOV_dir,    pattern = "_GOV_TAS.tif$",    full.names = TRUE)
fut_tasmax_GOV_files <- list.files(fut_tasmax_GOV_dir, pattern = "_GOV_TASMAX.tif$", full.names = TRUE)
fut_pr_GOV_files     <- list.files(fut_pr_GOV_dir,     pattern = "_GOV_PR.tif$",     full.names = TRUE)


## Deriving species names from filenames

extract_species <- function(x) {
  name <- basename(x)
  name <- tolower(name)
  name <- gsub(" ", "_", name)
  name <- gsub("_(hist|bau|gov)_(tas|tasmax|pr)\\.tif$", "", name)
  return(name)
}

sp_hist_tas    <- extract_species(hist_tas_files)
sp_hist_tasmax <- extract_species(hist_tasmax_files)
sp_hist_pr     <- extract_species(hist_pr_files)

sp_bau_tas    <- extract_species(fut_tas_bau_files)
sp_bau_tasmax <- extract_species(fut_tasmax_bau_files)
sp_bau_pr     <- extract_species(fut_pr_bau_files)

sp_gov_tas    <- extract_species(fut_tas_GOV_files)
sp_gov_tasmax <- extract_species(fut_tasmax_GOV_files)
sp_gov_pr     <- extract_species(fut_pr_GOV_files)


## Plot output directory

plot_base_dir <- "pathwayto//outputdirectory"


## Running overlaps for each parameter and scenario

## High-pressure

results_tas_BAU <- run_overlap(hist_tas_files, fut_tas_bau_files,
                               sp_hist_tas, sp_bau_tas, "TAS", "BAU",
                               file.path(plot_base_dir, "TAS_BAU"))
write.csv(results_tas_BAU, "results_tas_BAU.csv", row.names = FALSE)

results_tx_BAU <- run_overlap(hist_tasmax_files, fut_tasmax_bau_files,
                              sp_hist_tasmax, sp_bau_tasmax, "TASMAX", "BAU",
                              file.path(plot_base_dir, "TX_BAU"))
write.csv(results_tx_BAU, "results_tasmax_BAU.csv", row.names = FALSE)

results_pr_BAU <- run_overlap(hist_pr_files, fut_pr_bau_files,
                              sp_hist_pr, sp_bau_pr, "PR", "BAU",
                              file.path(plot_base_dir, "PR_BAU"))
write.csv(results_pr_BAU, "results_pr_BAU.csv", row.names = FALSE)


## Low-pressure

results_tas_GOV <- run_overlap(hist_tas_files, fut_tas_GOV_files,
                               sp_hist_tas, sp_gov_tas, "TAS", "GOV",
                               file.path(plot_base_dir, "TAS_GOV"))
write.csv(results_tas_GOV, "results_tas_GOV.csv", row.names = FALSE)

results_tx_GOV <- run_overlap(hist_tasmax_files, fut_tasmax_GOV_files,
                              sp_hist_tasmax, sp_gov_tasmax, "TASMAX", "GOV",
                              file.path(plot_base_dir, "TX_GOV"))
write.csv(results_tx_GOV, "D:/MRes (BU)/03062026/results_tasmax_GOV.csv", row.names = FALSE)

results_pr_GOV <- run_overlap(hist_pr_files, fut_pr_GOV_files,
                              sp_hist_pr, sp_gov_pr, "PR", "GOV",
                              file.path(plot_base_dir, "PR_GOV"))
write.csv(results_pr_GOV, "results_pr_GOV.csv", row.names = FALSE)


#####################################################
## STEP 5 - MAGNITUDE OF PROJECTED CLIMATIC CHANGE ##
#####################################################

## Shapefile lookup by species name

shp_files <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)

sp_shp <- tolower(gsub(" ", "_",
                       tools::file_path_sans_ext(basename(shp_files))))

shp_lookup <- setNames(shp_files, sp_shp)


## Per-species change statistics

calc_change <- function(hist_file, fut_file, poly_file) {
  
  poly <- st_read(poly_file, quiet = TRUE)
  poly <- st_make_valid(poly)
  poly <- terra::vect(poly)
  
  r_hist <- terra::rast(hist_file)
  r_fut  <- terra::rast(fut_file)
  
  if (!identical(terra::crs(r_hist), terra::crs(poly))) {
    poly <- terra::project(poly, terra::crs(r_hist))
  }
  
  h <- terra::mask(terra::crop(r_hist, poly), poly)
  f <- terra::mask(terra::crop(r_fut,  poly), poly)
  
  hv <- terra::values(h, na.rm = TRUE)
  fv <- terra::values(f, na.rm = TRUE)
  
  if (length(hv) == 0 || length(fv) == 0) return(NULL)
  
  diff <- fv - hv
  
  data.frame(
    mean_change   = mean(diff, na.rm = TRUE),
    median_change = median(diff, na.rm = TRUE),
    max_change    = max(diff, na.rm = TRUE),
    p95_change    = as.numeric(quantile(diff, 0.95, na.rm = TRUE)),
    n_cells       = length(diff)
  )
}


## Change across species

run_change <- function(hist_files, fut_files, sp_hist, sp_fut, label) {
  
  common_sp <- intersect(sp_hist, sp_fut)
  results <- list()
  
  for (sp in common_sp) {
    
    hist_file <- hist_files[sp_hist == sp]
    fut_file  <- fut_files[sp_fut  == sp]
    poly_file <- shp_lookup[[sp]]
    
    ## Diagnostics for unmatched species
    if (is.null(poly_file)) {
      message("No shapefile: ", sp)
    }
    
    if (length(hist_file) != 1) {
      message("Hist match problem: ", sp, " (", length(hist_file), ")")
    }
    
    if (length(fut_file) != 1) {
      message("Fut match problem: ", sp, " (", length(fut_file), ")")
    }
    
    if (!is.null(poly_file) &&
        length(hist_file) == 1 &&
        length(fut_file)  == 1) {
      
      out <- calc_change(hist_file, fut_file, poly_file)
      
      if (!is.null(out)) {
        out$species   <- sp
        results[[sp]] <- out
        message("Done ", label, ": ", sp)
      }
    }
  }
  
  dplyr::bind_rows(results)
}


## High-pressure scenario

results_tas_change    <- run_change(hist_tas_files, fut_tas_bau_files,
                                    sp_hist_tas, sp_bau_tas, "TAS BAU")

results_tasmax_change <- run_change(hist_tasmax_files, fut_tasmax_bau_files,
                                    sp_hist_tasmax, sp_bau_tasmax, "TASMAX BAU")

results_pr_change     <- run_change(hist_pr_files, fut_pr_bau_files,
                                    sp_hist_pr, sp_bau_pr, "PR BAU")
getwd()

write.csv(results_tas_change,    "raw_change_TAS_BAU.csv",    row.names = FALSE)
write.csv(results_tasmax_change, "raw_change_TASMAX_BAU.csv", row.names = FALSE)
write.csv(results_pr_change,     "raw_change_PR_BAU.csv",     row.names = FALSE)


## Combining parameters

combined_change <- results_tas_change %>%
  rename_with(~paste0("tas_", .), -species) %>%
  left_join(
    results_tasmax_change %>%
      rename_with(~paste0("tasmax_", .), -species),
    by = "species"
  ) %>%
  left_join(
    results_pr_change %>%
      rename_with(~paste0("pr_", .), -species),
    by = "species"
  ) %>%
  mutate(total_mean_change = tas_mean_change +
           tasmax_mean_change +
           pr_mean_change)

top_change <- combined_change %>%
  arrange(desc(total_mean_change))

print(top_change)
write.csv(top_change, "top_change_BAU.csv", row.names = FALSE)


## Low-pressure scenario

results_tas_change_GOV    <- run_change(hist_tas_files, fut_tas_GOV_files,
                                        sp_hist_tas, sp_gov_tas, "TAS GOV")

results_tasmax_change_GOV <- run_change(hist_tasmax_files, fut_tasmax_GOV_files,
                                        sp_hist_tasmax, sp_gov_tasmax, "TASMAX GOV")

results_pr_change_GOV     <- run_change(hist_pr_files, fut_pr_GOV_files,
                                        sp_hist_pr, sp_gov_pr, "PR GOV")

write.csv(results_tas_change_GOV,    "raw_change_TAS_GOV.csv",    row.names = FALSE)
write.csv(results_tasmax_change_GOV, "raw_change_TASMAX_GOV.csv", row.names = FALSE)
write.csv(results_pr_change_GOV,     "raw_change_PR_GOV.csv",     row.names = FALSE)


## Combining parameters, governance

combined_change_GOV <- results_tas_change_GOV %>%
  rename_with(~paste0("tas_", .), -species) %>%
  left_join(
    results_tasmax_change_GOV %>%
      rename_with(~paste0("tasmax_", .), -species),
    by = "species"
  ) %>%
  left_join(
    results_pr_change_GOV %>%
      rename_with(~paste0("pr_", .), -species),
    by = "species"
  ) %>%
  mutate(total_mean_change = tas_mean_change +
           tasmax_mean_change +
           pr_mean_change)

top_change_GOV <- combined_change_GOV %>%
  arrange(desc(total_mean_change))

print(top_change_GOV)
write.csv(top_change_GOV, "top_change_GOV.csv", row.names = FALSE)


##############################################################
## STEP 6 - PROJECTED DEFORESTATION WITHIN EACH SPECIES EOO ##
##############################################################

library(terra)
library(sf)
library(tools)


## Forest-loss change layers

# Source: Soares-Filho et al.. SimAmazonia

amazon_bau_f <- rast("pathwayto//amazonbau2050.tif")

amazon_bau_h <- rast("pathwayto//amazonbau2002.tif")

amazon_gov_f <- rast("pathwayto//amazongov2050.tif")

amazon_gov_h <- rast("pathwayto//amazongov2002.tif")

lulc_change_bau <- amazon_bau_f - amazon_bau_h

plot(lulc_change_gov)

lulc_change_gov <- amazon_gov_f - amazon_gov_h


## Export

writeRaster(lulc_change_bau, "lulc_change_bau.tif", overwrite = TRUE)
writeRaster(lulc_change_gov, "lulc_change_gov.tif", overwrite = TRUE)


## Recoding change values

lulc_change_bau[lulc_change_bau == -1] <- 1
lulc_change_gov[lulc_change_gov == -1] <- 1

lulc_change_bau[is.nan(lulc_change_bau)] <- NA
lulc_change_gov[is.nan(lulc_change_gov)] <- NA


## Per-species deforestation percentage

shp_files <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)

results <- data.frame(
  Species          = character(),
  BAU_defor_cells  = numeric(),
  BAU_total_cells  = numeric(),
  BAU_percent      = numeric(),
  GOV_defor_cells  = numeric(),
  GOV_total_cells  = numeric(),
  GOV_percent      = numeric(),
  stringsAsFactors = FALSE
)

for (shp in shp_files) {
  
  species_name <- tools::file_path_sans_ext(basename(shp))
  species_vect <- vect(shp)
  
  ## BAU
  bau_vals  <- extract(lulc_change_bau, species_vect)[, "amazonbau2050"]
  bau_defor <- sum(bau_vals == 1, na.rm = TRUE)
  bau_total <- length(bau_vals)
  bau_pct   <- (bau_defor / bau_total) * 100
  
  ## GOV
  gov_vals  <- extract(lulc_change_gov, species_vect)[, "amazongov2050"]
  gov_defor <- sum(gov_vals == 1, na.rm = TRUE)
  gov_total <- length(gov_vals)
  gov_pct   <- (gov_defor / gov_total) * 100
  
  results <- rbind(results, data.frame(
    Species         = species_name,
    BAU_defor_cells = bau_defor,
    BAU_total_cells = bau_total,
    BAU_percent     = bau_pct,
    GOV_defor_cells = gov_defor,
    GOV_total_cells = gov_total,
    GOV_percent     = gov_pct
  ))
  
  cat("Completed:", species_name, "\n")
}

print(results)

write.csv(results,
          "individual_species_deforestation_summary.csv",
          row.names = FALSE)

#################################################################################
## Necessary extractions to build final dataset ##
# ---------------------------------------------- #
#  1.Species list (LC/NT, VU, EN, CR [without DD])
#  2.Corresponding potential threat status 
#  3.Climate parameters overlap metric (from results_ outputs)
#  4.Climate parameters median change (from raw_change_ outputs)
#  5.Forest-loss parameters from (individual_species_deforestation_summary.csv)
# ---------------------------------------------- #
#  An excel file containing scenario-specific sheets containing this information is necessary for the modelling steps. 



