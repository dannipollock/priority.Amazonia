
# memory
rm(list = ls())


library(raster)
library(ggplot2)
library(viridis)
library(rainbow)
library(maptools)
library(sf)
library(redlistr)
library(ggplot2)
library(ColorBrewer)
library(ggsn)
library(rgdal)
library(wesanderson)
library(tidyverse)
library(terra)

pkg <- c("raster",
         "ggplot2",
         "viridis",
         "rainbow",
         "maptools",
         "sf",
         "redlistr",
         "ggplot2",
         "RColorBrewer",
         "ggsn",
         "rgdal",
         "wesanderson",
         "tidyverse",
         "terra"
)

sapply(pkg, library, character.only=TRUE, logical.return=TRUE)

setwd("G:/Meu Drive/Artigo_Scarabaeinae_Amazonia_Pollock")
dir()


var.SA<- raster::raster("var.SA.tif")### I use a raster to do the rasterization 
plot(var.SA)

SA_limpo<- sf::st_read("SA_limpo.shp") #### use the south american mask to make the cut
 
SA_limpo
plot(SA_limpo)

America<- sf::st_read("America do Sul_Continente.shp")  ### shape of south america with the countries to look beautiful on the map 
America
plot(America)

Amazon<- sf::st_read("Limites_Amazon.shp")  ### limits of Amazonia 
Amazon
plot(Amazon)



########################## occ 


occ<- readr::read_csv("Tab_604_spp_All_occ_Amazon.csv")
occ

occ$Longitude <- as.numeric(occ$Longitude)
occ$Latitude <- as.numeric(occ$Latitude)
nrow(occ)

# 35654

sp <- occ %>% 
  dplyr::select(Species) %>% 
  dplyr::distinct() %>% 
  dplyr::pull()
sp




################################
#Initialize an empty list to store results
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
  perc_area_EOO_in_Amazon  <- as.numeric(area_EOO_in_Amazon / total_area_EOO) * 100
  
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
    geom_point(data = occ.sp, aes(Longitude, Latitude), fill = NA, size = 0.4, color = "orange", alpha = 1) +
    geom_sf(data =  occ_in_Amazon, fill = NA, color = "red", size = 0.05) +
   
    
     geom_sf(data = EOO.sf, color = "orange", fill = "gray95", alpha = 0.1) +
    geom_sf(data =  EOO_in_Amazon, color = "red", fill = "gray1", alpha = 0.1) +
   
    
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
 setwd("F:/Mapas_Scarabaeinae_Pollock")
  
  ggsave(filename = paste0(i,".jpg"), plot =Amazon_all_occ, wi = 10, he = 12, un = "cm", dpi = 200)
  
  
}

# Combine results into a single data frame
results_df <- do.call(rbind, results)

# Save results to a CSV file
setwd("F:/Mapas_Scarabaeinae_Pollock")

write.csv(results_df, "species_amazon_analysis.csv", row.names = FALSE)
