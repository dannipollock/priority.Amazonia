################################################################################
##
## Script S1
##
## Compilation and cleaning of Scarabaeidae (dung beetle) occurrence records
## from the Global Biodiversity Information Facility (GBIF)
##
## Author:       D. Pollock
## Affiliation:  School of Biological Sciences, University of Bristol
## Date Amended: 2026-06
##
## ---------------------------------------------------------------------------
## Purpose
## ---------------------------------------------------------------------------
## This script documents the full cleaning pipeline applied to a GBIF
## occurrence download for the family Scarabaeidae, from raw Darwin Core
## archive through to the final validated dataset used in subsequent analyses.
##
## The workflow follows the protocol of Lima et al. (2021) as implemented in
## the plantR package:
##   Lima, R.A.F. et al. (2021) plantR: An R package and workflow for managing
##   species records from biological collections. Methods in Ecology and
##   Evolution, 12, 1353-1359. https://doi.org/10.1111/2041-210X.13779
##
## ---------------------------------------------------------------------------
## Pipeline overview
## ---------------------------------------------------------------------------
##   STEP 1  Data entry
##   STEP 2  Standardising collector names and collection information
##   STEP 3  Taxonomic data cleaning and name reconciliation
##   STEP 4  Standardising locality and spatial information
##   STEP 5  Validating locality and spatial information
##   STEP 6  Data summary
##   STEP 7  Duplicate detection and export
##
## ---------------------------------------------------------------------------
## Input
## ---------------------------------------------------------------------------
##   GBIF occurrence download 0082819-240506114902167 (DOI in manuscript)
##   FixLoc2_function.R  -- local modification of plantR::fixLoc() providing
##                          character-encoding correction for Latin American
##                          locality strings
##
##
################################################################################
## SESSION SETUP
################################################################################

rm(list = ls())

options(scipen = 999)        # suppress scientific notation
options(max.print = 1e8)
set.seed(29041988)           # fixed seed for reproducibility


################################################################################
## PACKAGES
################################################################################

## ---- Installation ---------------------------------------------------------
## Run once. Development versions are required for plantR, scrubr, Rocc and
## bold; CRAN versions of the remaining packages are sufficient.

install.packages("remotes")
library("remotes")

remotes::install_github("LimaRAF/plantR")
library("plantR")

remotes::install_github("ropensci/scrubr")
library("scrubr")

install.packages("ggplot2")
install.packages("cli")
install.packages("countrycode")

devtools::install_github("liibre/Rocc")
remotes::install_github("ropensci/bold")

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
         "Rocc",
         "countrycode",
         "jsonlite")

sapply(pkg, library, character.only = TRUE, logical.return = TRUE)

## Modified plantR::fixLoc() with character-encoding correction, called
## internally by formatLoc() in STEP 4.
source("FixLoc2_function.R")


## ---------------------------------------------------------------------------
## gnr_resolve2()
##
## Replacement for taxize::gnr_resolve(), which is no longer functional. Queries
## the Global Names Verifier API (https://verifier.globalnames.org) and returns
## the same five columns as gnr_resolve(fields = "minimal"):
## user_supplied_name, submitted_name, matched_name, data_source_title, score.
## ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

gnr_resolve2 <- function(sci,
                         data_source_ids = c(1, 11, 3),
                         chunk = 50,
                         pause = 1) {
  
  sci    <- unique(trimws(sci[!is.na(sci) & trimws(sci) != ""]))
  chunks <- split(sci, ceiling(seq_along(sci) / chunk))
  
  out <- lapply(seq_along(chunks), function(i) {
    
    url <- paste0(
      "https://verifier.globalnames.org/api/v1/verifications/",
      paste(vapply(chunks[[i]], utils::URLencode, character(1), reserved = TRUE),
            collapse = "%7C"),
      "?data_sources=", paste(data_source_ids, collapse = "%7C"),
      "&all_matches=true")
    
    r <- tryCatch(fromJSON(url, simplifyDataFrame = FALSE)$names,
                  error = function(e) {
                    message("chunk ", i, " failed")
                    NULL
                  })
    
    cat(sprintf("\r%d / %d", i, length(chunks)))
    flush.console()
    Sys.sleep(pause)
    r
  })
  cat("\n")
  
  nms <- unlist(out[!vapply(out, is.null, logical(1))], recursive = FALSE)
  
  res <- bind_rows(lapply(nms, function(n) {
    if (length(n$results %||% list()) == 0) return(NULL)
    bind_rows(lapply(n$results, function(r) tibble(
      user_supplied_name = n$name,
      submitted_name     = n$name,
      matched_name       = r$matchedName %||% NA_character_,
      data_source_title  = r$dataSourceTitleShort,
      score              = r$sortScore %||% NA_real_
    )))
  }))
  
  message(length(unique(res$user_supplied_name)), " of ", length(sci),
          " names returned matches")
  res
}

################################################################################
## STEP 1 - DATA ENTRY
################################################################################

beetles <- readData(file = "0082819-240506114902167.zip",
                    path = "https://api.gbif.org/v1/occurrence/download/request/")$occurrence

nrow(beetles)   # 506,522 records, 234 fields
str(beetles)

################################################################################
## STEP 2 - STANDARDISING COLLECTOR NAMES AND COLLECTION INFORMATION
################################################################################

## This step precedes taxonomic cleaning because formatOcc() generates columns
## required by the subsequent filtering steps (notably year.new and
## yearIdentified.new).

## formatOcc() edits collector and identifier names, collector numbers and
## collection years.

beetles1 <- formatOcc(beetles)

str(beetles1)   # 506,522 records, 244 fields (10 columns added)


new.collect.cols <- data.frame(beetles1$collectionCode.new,
                               beetles1$collectionObs,
                               beetles1$recordedBy.new,
                               beetles1$recordNumber.new,
                               beetles1$year.new,
                               beetles1$identifiedBy.new,
                               beetles1$yearIdentified.new,
                               beetles1$recordedBy.aux,
                               beetles1$identifiedBy.aux,
                               beetles1$last.name)


################################################################################
## STEP 3 - TAXONOMIC DATA CLEANING
################################################################################

names(beetles1)


## ---- 3.1 Remove records lacking a specific epithet -------------------------

beetles2 <- beetles1[beetles1$specificEpithet != "", ]

str(beetles2)   # 354,478 records retained; 152,044 removed

write.csv(beetles2, "GBIF_DungBeetle_specificEpithet_CLEAN.csv")


## Remove records lacking identification or collection year

beetles3 <- beetles2[beetles2$yearIdentified.new != "n.d.", ]
str(beetles3)   # 241,024 records retained; 113,454 removed

beetles4 <- beetles3[beetles3$year.new != "n.d.", ]
str(beetles4)   # 241,024 records retained; 43,672 removed


beetles5 <- filter(beetles4, yearIdentified.new >= year.new)

str(beetles5)   # 240,696 records retained; 328 removed


beetles6 <- beetles5 %>%
  tidyr::unite("scientificName.new", c(genus, specificEpithet), sep = " ") %>%
  mutate(across(where(is.character), ~ na_if(., "")))

str(beetles6)    # 240,696 records, 243 fields
names(beetles6)  # 'genus' and 'specificEpithet' consumed by unite()


## Records per species 

beetle_list <- beetles6 %>% dplyr::count(scientificName.new)
beetle_list


## Identification fields for expert verification 

ID.cols <- dplyr::select(beetles6, c("materialSampleID",
                                     "eventID",
                                     "parentEventID",
                                     "eventDate",
                                     "eventTime",
                                     "verbatimEventDate",
                                     "samplingProtocol",
                                     "locationID",
                                     "nameAccordingToID",
                                     "order",
                                     "family",
                                     "subfamily",
                                     "scientificName.new",
                                     "iucnRedListCategory",
                                     "collectionCode.new",
                                     "collectionObs",
                                     "recordedBy.new",
                                     "recordNumber.new",
                                     "year.new",
                                     "identifiedBy.new",
                                     "yearIdentified.new",
                                     "recordedBy.aux",
                                     "identifiedBy.aux",
                                     "last.name"))

write.csv(ID.cols, "GBIF_DungBeetle_CLEAN_ID_cols.csv")


## Matching against the Open Tree Taxonomy 

beetle_OTTmatch <- rotl::tnrs_match_names(
  names                   = beetle_list$scientificName.new,
  context_name            = "All life",
  do_approximate_matching = TRUE,
  ids                     = NULL,
  include_suppressed      = FALSE)

## Forty-three names returned no match and were flagged for manual review as
## potential synonyms or misidentifications (listed in Table SX).

head(beetle_OTTmatch)
write.csv(beetle_OTTmatch, "GBIF_DungBeetle_species_OTTmatching.csv")


## Matching scores from the Global Names Verifier -------------------

## Each name in the species list is matched against Catalogue of Life, the
## GBIF Backbone Taxonomy and ITIS, retaining all matches per name.

beetle_matching <- gnr_resolve2(beetle_list$scientificName.new)

## 2,136 of 2,136 names returned matches

length(unique(beetle_matching$user_supplied_name))   # 2,136
str(beetle_matching)

write.csv(beetle_matching, "GBIF_DungBeetle_species_matching-scores.csv")


## ---- 3.9 Mean matching score per species ----------------------------------

beetle_match_sum <- aggregate(beetle_matching$score,
                              list(beetle_matching$user_supplied_name),
                              FUN = mean)

beetle_match_sum
write.csv(beetle_match_sum, "GBIF_DungBeetle_species_AVERAGE_scores.csv")


################################################################################
## STEP 4 - STANDARDISING LOCALITY AND SPATIAL INFORMATION
################################################################################

names(beetles6)


## Overall locality standardisation

## formatLoc() standardises the Darwin Core fields 'country', 'stateProvince',
## 'municipality' and 'locality': names are transliterated to English and
## notation standardised. It calls fixLoc(), prepLoc(), strLoc() and getLoc()
## internally.

beetles7 <- formatLoc(beetles6)

str(beetles7)    # 254 fields (11 columns added)
names(beetles7)


## Inspection of new columns 

new.loc.cols <- data.frame(beetles7$country.new,
                           beetles7$stateProvince.new,
                           beetles7$municipality.new,
                           beetles7$locality.new,
                           beetles7$locality.scrap,
                           beetles7$resol.orig,
                           beetles7$loc,
                           beetles7$loc.correct,
                           beetles7$latitude.gazetteer,
                           beetles7$longitude.gazetteer,
                           beetles7$resolution.gazetteer)

head(new.loc.cols)


## Geographical coordinates 

beetles8 <- formatCoord(beetles7)

str(beetles8)   # 24,222 records, 259 fields


################################################################################
## STEP 5 - VALIDATING LOCALITY AND SPATIAL INFORMATION
################################################################################

## Locality validation 

beetles_LocVal <- validateLoc(beetles8)

str(beetles_LocVal)


## Coordinate validation 

## validateCoord()

beetles9 <- validateCoord(x        = beetles8,
                          tax.name = "scientificName.new")

str(beetles9)   # 24,222 records, 262 fields


## Joining validation output 

beetles10 <- left_join(beetles8, beetles_LocVal)

beetles10_prepdup <- prepDup(
  beetles10,
  col.names = c(family        = "family.new",
                species       = "scientificName.new",
                col.name      = "recordedBy.new",
                col.last.name = "last.name",
                col.number    = "recordNumber.new",
                col.year      = "year.new",
                col.loc       = "municipality.new",
                loc.str       = "loc.correct"),
  rec.ID    = "gbifID")

beetles10_prepdup[, "numTombo"] <- factor(beetles10_prepdup[, "numTombo"])

beetles10_getdup <- getDup(beetles10_prepdup)

beetles10[, "gbifID"] <- factor(beetles10[, "gbifID"])

beetles11 <- beetles10 %>%
  left_join(beetles10_getdup, join_by(gbifID == numTombo))


###########################
## STEP 6 - DATA SUMMARY ##
###########################

## Summary of the compiled dataset 

beetles_summ <- summaryData(beetles11)
beetles_summ



## 6.2 Validation flags 

beetles_flags <- summaryFlags(beetles11)
beetles_flags


###########################################
## STEP 7 - DUPLICATE REMOVAL AND EXPORT ##
###########################################

## Records are additionally screened for duplication at the level of species
## and coordinates, irrespective of collector or collection year.

str(beetles11)   # 24,222 records, 270 fields

dups <- duplicated(beetles11[c("scientificName.new",
                               "decimalLatitude.new",
                               "decimalLongitude.new")])

sum(dups)   # 19,448 duplicate records flagged
str(dups)


## FINAL DATASET FOR MERGING ##

beetles_FINAL <- cbind(beetles11, dups)

write.csv(beetles_FINAL, "GBIF_DungBeetles_FINAL_CLEAN-V1.csv")


