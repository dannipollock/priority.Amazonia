#########################################################
## IUCN / Study Assessment Comparison ##
########################################
##
## Author: Danni Pollock
## Affiliation: School of Health and Biological Sciences, University of Bristol
## Last amended: 2026-06
##
##---------
## Purpose
##---------
## To perform comparisons of current IUCN species assessments and assessments from this study. 
##
##---------
## Data
##---------
## Data downloaded from IUCN (https://www.iucnredlist.org/search?dl=true&permalink=4eef9087-f64b-4090-93dc-6e5d0359bec6)
## Species and assessments from this study. 
##
##----------
## Packages
##----------
library(readxl)
library(dplyr)
library(ggplot2)

####################################
## STEP 1 - DATA ENTRY AND FILTER ##
####################################

data_sheet   <- read_excel("pathto// IUCN DATA")
filter_sheet <- read_excel("pathto// STUDY SPECIES LIST")

## Filter IUCN data by Amazonian species in this study

filtered_data <- data_sheet %>%
  filter(data_sheet$scientificName %in% filter_sheet$taxa)

## Export
getwd()
write.csv(filtered_data,"NAME OF FILTERED IUCN DATA")


###############################################
## STEP 2 - CALCULATE ASSESSMENT PROPORTIONS ##
###############################################

props <- filtered_data %>%
  count(data_sheet$redlistCategory, name = "n") %>%
  mutate(proportion = n / sum(n),
         percent    = round(100 * proportion, 1)) %>%
  arrange(desc(n))

props

## For IUCN Amazonian Species 

props <- filter_sheet %>%
  count(filter_sheet$Category_CriteriaB, name = "n") %>%
  mutate(proportion = n / sum(n),
         percent    = round(100 * proportion, 1)) %>%
  arrange(desc(n))

props

####################################
## STEP 3 - PLOT COMPARATIVE PLOT ##
####################################

bar_data <- tibble::tribble(
  ~source,      ~category,                     ~percent,
  "This study", "Data Deficient",                  24.9,
  "This study", "Least Concern/ Near Threatened",  42.3,
  "This study", "Vulnerable",                      19.6,
  "This study", "Endangered",                      12.6,
  "This study", "Critically Endangered",            0.6,
  "IUCN 2026",  "Data Deficient",                  37.7,
  "IUCN 2026",  "Least Concern/ Near Threatened",  60.6,
  "IUCN 2026",  "Vulnerable",                       1.6,
  "IUCN 2026",  "Endangered",                       0.0,
  "IUCN 2026",  "Critically Endangered",            0.0
)

cat_levels <- c("Data Deficient", "Least Concern/ Near Threatened",
                "Vulnerable", "Endangered", "Critically Endangered")

cat_cols <- c("Critically Endangered"          = "#D81E07",
              "Endangered"                     = "#FB763F",
              "Vulnerable"                     = "#F5A643",
              "Least Concern/ Near Threatened" = "#7FC84D",
              "Data Deficient"                 = "#CDCDC4")

bar_data <- bar_data %>%
  mutate(category = factor(category, levels = cat_levels),
         source   = factor(source, levels = c("This study", "IUCN 2026")))

names(bar_data)

windows()
ggplot(bar_data, aes(x = source, y = percent, fill = category)) +
  geom_col(width = 0.45) +
  scale_fill_manual(values = cat_cols,
                    breaks = rev(cat_levels),   # change depending on desired direction of colours
                    name   = NULL) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 10),
                     labels = function(x) paste0(x, "%"),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of species (%)") +
  theme_classic(base_size = 11) +
  theme(axis.line.x   = element_line(colour = "grey40"),
        axis.line.y   = element_blank(),
        axis.ticks.y  = element_line(colour = "grey60"),
        legend.position = "right",
        legend.key.size = unit(0.8, "lines"))
