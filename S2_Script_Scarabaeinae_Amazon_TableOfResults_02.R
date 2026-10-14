library(dplyr)
library(sf)
library(raster)
library(redlistr)
library(ggplot2)

# Initialize an empty list to store results
results <- list()

for(i in sp){         
  
  occ.sp <- occ %>% 
    dplyr::filter(Species == i) %>% 
    dplyr::select(Longitude, Latitude)
  
  r.occ <- raster::rasterize(occ.sp, var.SA)
  
  EOO.polygon <- redlistr::makeEOO(r.occ)
  EOO.sf <- st_as_sf(EOO.polygon)
  
  # Convert to spatial objects
  EOO.sf <- st_as_sf(EOO.polygon)
  
  # Area of the EOO within Amazon
  EOO_in_Amazon <- st_intersection(EOO.sf, Amazon)
  area_EOO_in_Amazon <- sum(st_area(EOO_in_Amazon))
  total_area_EOO <- sum(st_area(EOO.sf))
  perc_area_EOO_in_Amazon <- as.numeric(area_EOO_in_Amazon / total_area_EOO) * 100
  
  # Occurrences within Amazon
  occ.sf <- st_as_sf(occ.sp, coords = c("Longitude", "Latitude"), crs = st_crs(Amazon))
  occ_in_Amazon <- st_intersection(occ.sf, Amazon)
  num_occ_in_Amazon <- nrow(occ_in_Amazon)
  total_occ <- nrow(occ.sf)
  perc_occ_in_Amazon <- (num_occ_in_Amazon / total_occ) * 100
  
  # Store results
  results[[i]] <- data.frame(
    Species = i,
    Percentage_Area_EOO_in_Amazon = perc_area_EOO_in_Amazon,
    Percentage_Occurrences_in_Amazon = perc_occ_in_Amazon
  )
  
  # Plotting
  Amazon_all_occ <- ggplot() +
    geom_sf(data = America, fill = NA, color = "gray", size = 0.05) +
    geom_sf(data = Amazon, fill = NA, color = "green4", size = 0.05) +
    geom_point(data = occ.sp, aes(Longitude, Latitude), fill = NA, size = 0.4, color = "red", alpha = 1) +
    geom_sf(data = EOO.sf, color = "red", fill = "gray", alpha = 0.1) +
    labs(x = "Longitude", y = "Latitude", 
         fill = "",
         title = paste0(stringr::str_to_sentence(i))) +
    theme_classic() + 
    theme(legend.position = "none",
          title = element_text(face = "italic", size = 12),
          axis.title = element_text(face = "plain", size = 12),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 12))
  
  print(Amazon_all_occ)
}

# Combine results into a single data frame
results_df <- do.call(rbind, results)

# Save results to a CSV file
write.csv(results_df, "species_amazon_analysis.csv", row.names = FALSE)