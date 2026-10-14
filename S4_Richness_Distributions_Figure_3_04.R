################################################################################
##
## Supplementary Script S4
##
## Spatial distribution of dung beetle species richness, potential threat
## status, and range overlap across Amazonia
##
## Author:       D. Pollock
## Affiliation:  TropEcoLab, School of Biological Sciences, University of Bristol
## Last revised: 2026-09
##
## ---------------------------------------------------------------------------
## Purpose
## ---------------------------------------------------------------------------
## Produces a four-panel figure:
##
##   A  Total dung beetle species richness
##   B  Richness of potentially threatened species (Criterion B: CR, EN, VU)
##   C  Proportion of species that are potentially threatened
##   D  Potspots of extent-of-occurrence (EOO) polygon overlap among
##      potentially threatened species (25% upper quartile of cells)
##
## ---------------------------------------------------------------------------
## Inputs
## ---------------------------------------------------------------------------
##   IUCN_Final.Dataset.biogeo.with.DD.Coords.2026.csv.xlsx
##     Species occurrence records with Criterion B assessments and coordinates
##
##   shapesIUCN_TH/*EOO_poly.shp
##     Per-species EOO polygons exported by ConR::IUCN.eval() for potentially
##     threatened species
##
##   Amazon
##     Amazon basin boundary (sf object; created in Script S2)
##
## ---------------------------------------------------------------------------
## Output
## ---------------------------------------------------------------------------
##  4 Panel plot of richness parameters. 
##
############################################
##
## PACKAGES
##
############################################

library(readr)
library(sf)
library(patchwork)
library(scales)
library(rgdal)
library(rgeos)
library(ggplot2)
library(rworldmap)
library(dplyr)
library(readxl)


#########################
## STEP 1 - DATA ENTRY ##
#########################

data <- read_xlsx("F:/MRes (BU)/03062026/IUCN_Final.Dataset.biogeo.with.DD.Coords.2026.csv.xlsx")


##############################################################
## STEP 2 - GRIDDED RICHNESS AND THREAT STATUS (PANELS A-C) ##
##############################################################

## Species assessed as CR, EN or VU treated as potentially threatened.

coords_dat <- data %>%
  mutate(
    threatened = Category_CriteriaB %in% c("CR", "EN", "VU"),
    cell_lon   = round(ddlon / 1.5) * 1.5,
    cell_lat   = round(ddlat / 1.5) * 1.5
  )


## Per-cell richness and proportion threatened

cell_richness <- coords_dat %>%
  group_by(cell_lon, cell_lat) %>%
  summarise(
    total_spp    = n_distinct(taxa),
    threatened_n = n_distinct(taxa[threatened]),
    .groups      = "drop"
  ) %>%
  mutate(proportion = threatened_n / total_spp)


## Clipping cells to the Amazon basin

cell_sf <- cell_richness %>%
  st_as_sf(coords = c("cell_lon", "cell_lat"), crs = 4326) %>%
  st_intersection(st_union(Amazon))

cell_richness_clipped <- cell_richness %>%
  filter(paste(cell_lon, cell_lat) %in%
           paste(st_coordinates(cell_sf)[, 1],
                 st_coordinates(cell_sf)[, 2]))


############################################
## STEP 3 - GRIDDED EOO OVERLAP (PANEL D) ##
############################################

## Reading and combining per-species EOO polygons

shp_folder <- "F:/MRes (BU)/03062026/shapesIUCN_TH/"
shp_files  <- list.files(shp_folder, pattern = "EOO_poly.shp$", full.names = TRUE)
shp_list   <- lapply(shp_files, st_read)
all_shapes <- bind_rows(shp_list)
all_shapes <- st_transform(all_shapes, crs = 4326)


## Pairwise range overlap

## The number of other species polygons intersected excluding itself

all_shapes$poly <- 1
sf_use_s2(TRUE)
overlaps <- st_intersects(all_shapes)
all_shapes$overlap_count_no_self <- lengths(overlaps) - 1


## Constructing the analysis grid

amazon_extent <- st_bbox(c(xmin = -80, ymin = -20,
                           xmax = -43, ymax =  9),
                         crs = st_crs(4326)) %>%
  st_as_sfc()

grid_cells <- st_make_grid(amazon_extent,
                           cellsize = 1.5,
                           square   = TRUE) %>%
  st_as_sf() %>%
  mutate(
    cell_id  = row_number(),
    cell_lon = st_coordinates(st_centroid(x))[, 1],
    cell_lat = st_coordinates(st_centroid(x))[, 2]
  )


## Counting EOO polygons per grid cell

eoo_grid <- st_join(grid_cells, all_shapes, join = st_intersects)

cell_overlap <- eoo_grid %>%
  st_drop_geometry() %>%
  group_by(cell_id, cell_lon, cell_lat) %>%
  summarise(
    eoo_count = sum(!is.na(overlap_count_no_self)),
    .groups   = "drop"
  ) %>%
  filter(eoo_count > 0)


## Clipping cells to the Amazon basin

cell_overlap_clipped <- cell_overlap %>%
  left_join(grid_cells %>% st_drop_geometry(),
            by = c("cell_id", "cell_lon", "cell_lat")) %>%
  st_as_sf(coords = c("cell_lon", "cell_lat"), crs = 4326) %>%
  st_filter(Amazon) %>%
  mutate(
    cell_lon = st_coordinates(.)[, 1],
    cell_lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()


## Restricting to upper-quartile cells

overlap_threshold <- quantile(cell_overlap_clipped$eoo_count, 0.75, na.rm = TRUE)

cell_overlap_top25 <- cell_overlap_clipped %>%
  filter(eoo_count >= overlap_threshold)


###################################
## STEP 4 - SHARED PLOT ELEMENTS ##
###################################

map_theme <- theme_void(base_size = 10) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 8, colour = "black",
                                     face = "bold"),
    legend.text       = element_text(size = 7, colour = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    plot.title        = element_text(size = 9, face = "bold",
                                     colour = "black",
                                     margin = margin(b = 4)),
    plot.margin       = margin(2, 2, 2, 2, "mm"),
    panel.background  = element_rect(fill = "white", colour = NA)
  )

richness_colours <- c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000")

amazon_bg <- list(
  geom_sf(data = Amazon, fill = "#F0F4F0",
          colour = "grey40", linewidth = 0.3),
  coord_sf(xlim = c(-80, -43), ylim = c(-20, 9), expand = FALSE)
)


################################################################################
## STEP 5 - FIGURE PANELS
################################################################################

## ---- Panel A: total species richness
pA <- ggplot() +
  amazon_bg +
  geom_tile(data   = cell_richness_clipped,
            aes(x = cell_lon, y = cell_lat, fill = total_spp),
            width = 1.5, height = 1.5,
            colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = richness_colours,
    name    = "Richness",
    guide   = guide_colorbar(ticks.colour = "black",
                             frame.colour = "black")
  ) +
  labs(title = "A) Species richness") +
  map_theme


## ---- Panel B: richness of potentially threatened species

pB <- ggplot() +
  amazon_bg +
  geom_tile(data   = cell_richness_clipped,
            aes(x = cell_lon, y = cell_lat, fill = threatened_n),
            width = 1.5, height = 1.5,
            colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = richness_colours,
    name    = "Threatened\nspp.",
    guide   = guide_colorbar(ticks.colour = "black",
                             frame.colour = "black")
  ) +
  labs(title = "B) Potentially threatened species richness") +
  map_theme


## ---- Panel C: proportion of species potentially threatened

pC <- ggplot() +
  amazon_bg +
  geom_tile(data   = cell_richness_clipped,
            aes(x = cell_lon, y = cell_lat, fill = proportion),
            width = 1.5, height = 1.5,
            colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = richness_colours,
    limits  = c(0, 1),
    labels  = scales::percent_format(accuracy = 1),
    name    = "Proportion\nthreatened",
    guide   = guide_colorbar(ticks.colour = "black",
                             frame.colour = "black")
  ) +
  labs(title = "C) Proportion of potentially threatened species") +
  map_theme


pD <- ggplot() +
  amazon_bg +
  geom_tile(data   = cell_overlap_top25,
            aes(x = cell_lon, y = cell_lat, fill = eoo_count),
            width = 1.5, height = 1.5,
            colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = richness_colours,
    name    = "EOO count\n(top 25%)",
    guide   = guide_colorbar(ticks.colour = "black",
                             frame.colour = "black")
  ) +
  labs(title = "D) Potentially threatened species range hotspots") +
  map_theme


#################################################
## STEP 6 - COMBINING AND EXPORTING THE FIGURE ##
#################################################

final_map <- (pA | pB) / (pC | pD)

windows(width = 14, height = 10)
final_map
