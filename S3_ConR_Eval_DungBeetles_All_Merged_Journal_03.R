################################################################################
##
## Script 03
##
## Extent of occurrence, area of occupancy and preliminary IUCN Criterion B
## assessment for Amazonian dung beetles (Coleoptera: Scarabaeidae)
##
## Author:       D. Pollock
## Affiliation:  School of Biological Sciences, University of Bristol
## Last revised: 2026-06
## --------
## Purpose
## --------
## Clean, merged and final species list and records
## Creates deduplicated datasets
## Estimates EOO and AOO 
## Provides preliminary conservation assessments under IUCN Red List Criterion B,
## incorporating protected areas and Indigenous territories
## ------
## Inputs
## ------
##   Occurrence data (final, merged)
##
##   Spatial layers
##     ANP_ReservaFlorestal.shp      protected areas (forest reserves)
##     ANP_departamental.shp         protected areas (departmental)
##     ANP_nacional.shp              protected areas (national)
##     ANP_BosqueProtector.shp       protected areas (protective forest)
##     Tis_TerritoriosIndigenas.shp  Indigenous territories
##     America_Sul.shp               South America country boundaries
##     limiteteste' - Copy.shp       Amazon basin boundary
##     neotropics.shp                Neotropical region boundary
##
##################
## SESSION SETUP #
##################

rm(list = ls())

options(scipen = 999)          # suppress scientific notation
options("max.print" = 100000000)
set.seed(29041988)             # fixed seed for reproducibility

##############
## PACKAGES ##
##############
## ---- Installation ---------------------------------------------------------
## Alternative installation routes retained for reference:
##   install.packages("ConR")
##   remotes::install_github("rspatial/rspatial")
##   install.packages("devtools")
##   devtools::install_github("gdauby/ConR")

install.packages('plantR')
library(plantR)

## ---- Loading --------------------------------------------------------------

pkg <- c("dplyr",
         "readr",
         "rgbif",
         "scrubr",
         "maps",
         "plantR",
         "remotes",
         "rotl",
         "tidyr",
         "taxize",
         "RColorBrewer",
         "raster",
         "rgdal",
         "rspatial",
         "sf",
         "colorRamps",
         "rasterVis",
         "ggplot2",
         "redlistr",
         "adehabitatHR",
         "sp",
         "ConR",
         "sf",
         "tidyverse",
         "rnaturalearth",
         "magrittr",
         "tidyterra",
         "terra",
         "rasterVis",
         "scrubr"
)

sapply(pkg, library, character.only = TRUE, logical.return = TRUE)


## STEP 1 - DATA ENTRY

beetles <- readxl::read_xlsx("FilePath// to Merged Dataset.xlsx",
                             sheet = "IndVal.amazonian.sp.only",
                             col_names = TRUE, na = "NA")

beetles <- as.data.frame(beetles)

nrow(beetles)  
str(beetles)

## Record counts by dataset variant:
##   all species                    80,498 records, 4 fields
##   CLAM only                      39,772 records, 4 fields
##   EOO / occurrences > 70         51,572 records, 4 fields
##   EOO + occurrences > 70         36,545 records, 4 fields
##   EOO + occurrences > 70, CLAM   41,062 records, 4 fields



## STEP 2 - TAXONOMIC DATA CLEANING

names(beetles)


## Records per species 

beetle_list <- beetles %>% dplyr::count(sp)
beetle_list

## Species richness by dataset variant:
##   all species                    500 species
##   CLAM only                      208 species
##   EOO + AOO > 70                 281 species
##   EOO / occurrences > 70         350 species
##   EOO / occurrences and CLAM     293 species


## Matching against the Open Tree Taxonomy 

beetle_OTTmatch <- rotl::tnrs_match_names(                # will produce list of flagged species
  names                   = beetle_list$sp,
  context_name            = "All life",
  do_approximate_matching = TRUE,
  ids                     = NULL,
  include_suppressed      = FALSE)



head(beetle_OTTmatch)
write.csv(beetle_OTTmatch,
          'GBIF_DungBeetle_species_OTTmatching.csv',
          row.names = F)


## Matching scores from the Global Names Verifier 

remotes::install_github("ropensci/taxize")

beetle_names <- beetle_list %>%
  filter(!is.na(sp), sp != "") %>%
  distinct(sp) %>%
  pull(sp)

chunks <- split(beetle_names, ceiling(seq_along(beetle_names) / 40))         # names are submitted in chunks of 40 to avoid limits

beetle_matching_list <- map(chunks, function(chunk) {
  tryCatch({
    taxize::gna_verifier(
      names = chunk,
      data_source_ids = c(3, 4),   # ITIS and NCBI
      resolve_once = TRUE,
      with_context = FALSE,
      canonical = FALSE,
      highestscore = TRUE,
      best_match_only = TRUE,
      ttp = "get",
      cap_first = TRUE,
      fields = "minimal"
    )
  }, error = function(e) {
    message("Chunk failed: ", e$message)
    NULL
  })
})

## Combine all successful chunks into one data frame
beetle_matching <- bind_rows(beetle_matching_list)

## Check how many names were successfully matched
matched_species <- unique(beetle_matching$submittedName)
missing_species <- setdiff(beetle_names, matched_species)

cat("Successfully matched:", length(matched_species), "species\n")
cat("Missing or dropped:", length(missing_species), "species\n")

str(beetle_matching)
write.csv(beetle_matching,
          'GBIF_DungBeetle_species_matching-scores.csv',
          row.names = F)


## Mean matching score per species 

beetle_match_sum <- aggregate(beetle_matching$sortScore,
                              list(beetle_matching$submittedName), FUN = mean)

names(beetle_match_sum)
colnames(beetle_match_sum) <- c("Species", "Averaged scores")
names(beetle_match_sum)

write.csv(beetle_match_sum,
          'FINAL_merged_DungBeetle_species_AVERAGE_scores_INDVAL.FINAL.amazonian.sp.only.amazonian.sp.only.amazonian.sp.only.csv',
          row.names = F)


##############################################
## STEP 3 - DUPLICATE DETECTION AND REMOVAL ##
##############################################

## Records are screened for duplication 

str(beetles)

## Record counts by dataset variant:
##   all species                    80,498 records, 4 fields
##   CLAM only                      39,792 records, 4 fields
##   EOO or occurrences             51,611 records, 4 fields
##   EOO + occurrences > 70         36,545 records, 4 fields
##   EOO + occurrences > 70, CLAM   41,062 records, 4 fields

dups <- duplicated(beetles[c("sp", "lat", "long")])

sum(dups)

## Duplicates flagged by dataset variant:
##   all species                    38,547
##   CLAM only                      18,268
##   EOO or occurrences             25,345
##   EOO + occurrences > 70         17,383
##   EOO + occurrences + CLAM       19,240

str(dups)   


## Flagging duplicates in the main data frame 

beetles_NEW <- cbind(beetles, dups)

## Cleaned dataset for the Americas
write.csv(beetles_NEW,
          'FINAL_merged_DungBeetles.csv',
          row.names = F)


## Removing duplicates 

beetles_NEW <- beetles_NEW[!beetles_NEW$dups == 'TRUE', ]

head(beetles_NEW)
str(beetles_NEW)

## Records retained by dataset variant:
##   all species                    41,998 records, 5 fields
##   CLAM only                      21,549 records, 5 fields
##   EOO or occurrences             26,266 records, 5 fields
##   EOO + occurrences > 70         19,162 records, 5 fields
##   EOO + occurrences + CLAM       21,822 records, 5 fields

write.csv(beetles_NEW, file = 'FINAL.amazonian.csv')


###################################################################
## STEP 4 - DISTRIBUTION MAPPING AND IUCN CRITERION B ASSESSMENT ##
###################################################################

## Checking for missing coordinates and identifications

sum(is.na(beetles_NEW$sp))     # 0
sum(is.na(beetles_NEW$lat))    # 0
sum(is.na(beetles_NEW$long))   # 0


## Retaining identification, coordinate and source fields 

names(beetles_NEW)
occ <- beetles_NEW[, c("sp", "lat", "long", "source")]
head(occ)


## Data source list 

Source.List <- occ %>%
  dplyr::select(source) %>%
  dplyr::distinct() %>%
  dplyr::pull()
Source.List

## Final occurrence dataset used in subsequent analyses
write.csv(occ,
          'FINAL_merged_DungBeetles_Amazonia.csv',
          row.names = F)


## Species list with matching scores 

sp <- occ %>%
  dplyr::select(sp) %>%
  dplyr::distinct() %>%
  dplyr::pull()

head(sp)

names(beetle_match_sum)
sp_final <- beetle_match_sum %>% filter(Species %in% sp)

## Final dung beetle species list for Amazonia
write.csv(sp_final,
          'FINAL_merged_DungBeetles_SpeciesList.csv',
          row.names = F)


## Spatial layers 

## Protected areas and Indigenous territories
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

PA  <- crop(PA, shp)
Ing <- crop(Ing, shp)
PA  <- mask(PA, shp)
Ing <- mask(Ing, shp)

## Label protected areas and Indigenous territories, then combine
PA$type  <- rep("PA", nrow(PA))
Ing$type <- rep("Indigenous", nrow(Ing))
ITPA     <- rbind(PA, Ing)
ITPA_sp  <- as(st_as_sf(ITPA), "Spatial")
ITPA$PA_id <- as.character(seq_len(nrow(ITPA)))
ITPA_sp  <- as(st_as_sf(ITPA), "Spatial")

names(ITPA_sp)

## South America
America  <- sf::st_read("America_Sul.shp")
America2 <- readOGR("America_Sul.shp")

## Amazonia
Amazon  <- sf::st_read("F:/MRes (BU)/limiteteste' - Copy.shp")
Amazon  <- st_make_valid(Amazon)
Amazon  <- st_transform(Amazon, 4326)
Amazon2 <- readOGR("F:/MRes (BU)/limiteteste' - Copy.shp")

## Neotropics
nt_shp  <- readOGR("neotropics.shp")
nt_shp2 <- sf::st_read("neotropics.shp")


## Renaming coordinate and taxon fields 

names(occ)
occ_final <- occ %>%
  dplyr::select(sp, lat, long, source)

names(occ_final)
names(occ_final)[1] <- paste("taxa")
names(occ_final)[2] <- paste("ddlat")
names(occ_final)[3] <- paste("ddlon")
names(occ_final)


## Coercing coordinates to numeric 

occ_final$ddlon <- as.numeric(occ_final$ddlon)
occ_final$ddlat <- as.numeric(occ_final$ddlat)
head(occ_final)

summary(occ_final)


## Restricting records to the Amazon basin 

points_sf <- st_as_sf(occ_final, coords = c("ddlon", "ddlat"),
                      crs = st_crs(Amazon))

sf_use_s2(FALSE)   # required where spherical geometry errors occur

points_inside_polygon <- st_filter(points_sf, Amazon)

coordinates_df <- as.data.frame(st_coordinates(points_inside_polygon))

names(coordinates_df) <- c("ddlon", "ddlat")

attributes_df <- st_drop_geometry(points_inside_polygon)

points_inside_df    <- cbind(coordinates_df, attributes_df)

points_inside_final <- cbind(coordinates_df, attributes_df)

## Check for missing coordinates
sum(is.na(points_inside_final$ddlon))
sum(is.na(points_inside_final$ddlat))

## Bounding-box filter approach:
##   occ_final <- occ_final %>%
##     filter(ddlat <= 5.283333 & ddlon < -41.81667)
##   summary(occ_final)


## Visual check  
windows(width = 14, height = 14)
ggplot() +
  geom_sf(data = America) +
  geom_sf(data = points_inside_final %>%
            st_as_sf(coords = c("ddlon", "ddlat"),
                     crs = 4326)) +
  theme_minimal()

## Removal of erroneous records manually here, if still present
.
## Attributing records to Amazonian countries 

## Approach after:
## https://www.r-bloggers.com/2022/04/robservations-29-classifying-and-filtering-coordinates-by-using-the-sf-library/

## Isolating country polygons
BoliviaPolygons      <- America %>% filter(COUNTRY == "Bolivia") %$% geometry
BrazilPolygons       <- America %>% filter(COUNTRY == "Brazil") %$% geometry
ColombiaPolygons     <- America %>% filter(COUNTRY == "Colombia") %$% geometry
EcuadorPolygons      <- America %>% filter(COUNTRY == "Ecuador") %$% geometry
FrenchGuianaPolygons <- America %>% filter(COUNTRY == "French Guiana (France)") %$% geometry
GuyanaPolygons       <- America %>% filter(COUNTRY == "Guyana") %$% geometry
PeruPolygons         <- America %>% filter(COUNTRY == "Peru") %$% geometry
SurinamePolygons     <- America %>% filter(COUNTRY == "Suriname") %$% geometry
VenezuelaPolygons    <- America %>% filter(COUNTRY == "Venezuela") %$% geometry

## Intersecting records with each country polygon

BoliviaIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(BoliviaPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

BrazilIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(BrazilPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

ColombiaIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(ColombiaPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

EcuadorIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(EcuadorPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

FrenchGuianaIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(FrenchGuianaPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

GuyanaIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(GuyanaPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

PeruIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(PeruPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

SurinameIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(SurinamePolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

VenezuelaIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(VenezuelaPolygons) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()

AmazonIndex <- points_inside_final %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326) %$% geometry %>%
  st_intersects(Amazon) %>%
  lapply(function(x) ifelse(length(x) == 0, FALSE, TRUE)) %>%
  unlist()


## Assigning country labels

points_inside_final$Country <- c()
points_inside_final$Country[BoliviaIndex]      <- "Bolivia"
points_inside_final$Country[BrazilIndex]       <- "Brazil"
points_inside_final$Country[ColombiaIndex]     <- "Colombia"
points_inside_final$Country[EcuadorIndex]      <- "Ecuador"
points_inside_final$Country[FrenchGuianaIndex] <- "French Guiana"
points_inside_final$Country[GuyanaIndex]       <- "Guyana"
points_inside_final$Country[PeruIndex]         <- "Peru"
points_inside_final$Country[SurinameIndex]     <- "Suriname"
points_inside_final$Country[VenezuelaIndex]    <- "Venezuela"


## Record distribution by country 

windows(width = 14, height = 14)
ggplot() +
  geom_sf(data = America) +
  geom_sf(data = points_inside_final %>%
            filter(!is.na(Country)) %>%
            st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326),
          mapping = aes(color = Country)) +
  theme_minimal()

## Equivalent plot using occ_final, retained for reference:
##   ggplot() +
##     geom_sf(data = America) +
##     geom_sf(data = occ_final %>%
##               filter(!is.na(Country)) %>%
##               st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326),
##             mapping = aes(color = Country)) +
##     theme_minimal()

ggsave("Country-level data distribution_INDVAL.FINAL.amazonian.sp.only.pdf",
       width = 14, height = 14)

## Final dung beetle dataset with country attribution
write.csv(points_inside_df,
          'FINAL_merged_clean_Amazon_DungBeetle_dataset.csv',
          row.names = F)


## Extent of occurrence 

names(points_inside_final)
occ_final <- points_inside_final %>% dplyr::select(ddlat, ddlon, taxa)

EOO <- EOO.computing(occ_final,
                     exclude.area      = TRUE,
                     country_map       = nt_shp,
                     parallel          = TRUE,
                     export_shp        = TRUE,
                     write_shp         = TRUE,
                     method.less.than3 = 'arbitrary')

str(EOO)


## Preliminary assessment under IUCN Criterion B 

names(occ_final)

opt <- ConR::IUCN.eval(occ_final,                # cleaned occurrence data
                       
                       country_map   = nt_shp,   # used to crop spatial polygons
                       # and compute EOO
                       
                       exclude.area  = TRUE,     # remove areas outside
                       # country_map
                       
                       protec.areas  = ITPA_sp,  # protected areas and
                       # Indigenous territories
                       
                       ID_shape_PA   = "PA_id",
                       
                       write_shp     = TRUE,     # export EOO polygons as ESRI
                       # shapefiles
                       
                       DrawMap       = TRUE,     # export a map per species
                       
                       Cell_size_AOO = 2,        # IUCN-standard AOO grid (km)
                       
                       showWarnings  = FALSE,
                       
                       SubPop        = TRUE,     # estimate subpopulations
                       
                       parallel      = TRUE,
                       
                       NbeCores      = 4)        # number of cores for parallel processing


## Assigning Data Deficient status 

## Species represented by < four records are reclassified as Data Deficient

occurrences_count <- occ_final %>%
  group_by(taxa) %>%
  summarise(count = n())

opt$Category_CriteriaB <- ifelse(
  opt$taxa %in% occurrences_count$taxa[occurrences_count$count < 4],
  "DD",
  opt$Category_CriteriaB)

head(opt)
class(opt)

write.csv(opt, file = "IUCN_eval_output.DD.csv")

# Manually remove these in excel for further prelinimary threatened species-based analyses

## Joining assessments to occurrence records 

occ_optmerged <- merge(occ_final, opt, by = "taxa")

typeof(occ_optmerged)
colnames(occ_optmerged)

write.csv(occ_optmerged, file = "IUCN_eval_output_and_coords.csv")


## Mapping assessment results in grid cells 

windows(width = 20, height = 20)
map_DB_all <- map.res(Results     = opt,
                      Occurrences = occ_final,
                      country_map = nt_shp,
                      threshold   = 5,   # minimum records per grid cell shown;
                      # set to 5 for consistency with the
                      # companion plant analysis
                      export_data = TRUE,
                      file_name   = "map_DB_all.CLAM.EOO+.FINAL.DD")

write.csv(map_DB_all, file = "map_DB_all_EOO.or.AOO.DD.csv")


