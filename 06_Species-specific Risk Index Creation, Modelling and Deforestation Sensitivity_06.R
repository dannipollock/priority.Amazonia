################################################################################
##
## Species-specific Risk Index Creation, Modelling and Deforestation Sensitivity 06
##
## Individual risk index, ordinal regression (and Spearman's rank) testing with brant tests for association patterns 
## between risk index and preliminary IUCN threat status 
##
## Author:       D. Pollock
## Affiliation:  School of Biological Sciences, University of Bristol
## Last amended: 2026-09
##
## ---------
##  Purpose
## ---------
## Combines climate exposure and deforestation exposure into an
## environmental risk index for two scenarios:
##
##   BAU  high-pressure  
##   GOV  low-pressure     
##
## --------
## Input
## --------
##     Sheets "GOV" and "BAU"; one row per species, with niche-overlap values,
##     median climate changes, deforestation exposure, and Criterion B category
##
##############
## PACKAGES ##
##############


library(readxl)
library(brant)
library(MASS)
library(ggplot2)
library(patchwork)
library(dplyr)
library(writexl)


#########################
## STEP 1 - DATA ENTRY ##
#########################
# Load Data

GOV <- readxl::read_xlsx("pathwayto//dataset",
                         sheet = "GOV",
                         col_names = TRUE, na= "NA")

BAU <- readxl::read_xlsx("pathwayto//dataset",
                         sheet = "BAU",
                         col_names = TRUE, na= "NA")


###############################
## STEP 2 - CLIMATE EXPOSURE ##
###############################

### Climate parameters
#convert to climate change exposure (higher value = more severe)

## Low-pressure scenario 

GOV$PR_overlap <- as.numeric(GOV$PR_overlap) #IF NEEDED
GOV$PR_exposure     <- 1 - GOV$PR_overlap

GOV$TAS_overlap <- as.numeric(GOV$TAS_overlap) #IF NEEDED
GOV$TAS_exposure    <- 1 - GOV$TAS_overlap

GOV$TASMAX_overlap <- as.numeric(GOV$TASMAX_overlap) #IF NEEDED
GOV$TASMAX_exposure <- 1 - GOV$TASMAX_overlap


## High-pressure scenario 

BAU$PR_overlap      <- as.numeric(BAU$PR_overlap) #IF NEEDED
BAU$PR_exposure     <- 1 - BAU$PR_overlap

BAU$TAS_overlap <- as.numeric(BAU$TAS_overlap) #IF NEEDED

BAU$TAS_exposure    <- 1 - BAU$TAS_overlap

BAU$TASMAX_overlap <- as.numeric(BAU$TASMAX_overlap) #IF NEEDED

BAU$TASMAX_exposure <- 1 - BAU$TASMAX_overlap


###############################################
## STEP 3 - STANDARDISATION ACROSS SCENARIOS ##
###############################################

##Standardisation function

z <- function(x, ref){
  (x - mean(ref, na.rm=TRUE)) / sd(ref, na.rm=TRUE)
}


## Pooled reference distributions

## Each variable is centred and scaled on the two scenarios combined
## z-scores have a common scale 

# pooled references
ref_PR        <- c(GOV$PR_exposure, BAU$PR_exposure)
ref_TAS       <- c(GOV$TAS_exposure, BAU$TAS_exposure)
ref_TASMAX    <- c(GOV$TASMAX_exposure, BAU$TASMAX_exposure)

ref_PRchg     <- c(GOV$PR_median_change, BAU$PR_median_change)
ref_TASchg    <- c(GOV$TAS_median_change, BAU$TAS_median_change)
ref_TASMAXchg <- c(GOV$TASMAX_median_change, BAU$TASMAX_median_change)

ref_defor     <- c(GOV$GOV_percent, BAU$BAU_percent)
ref_cells     <- c(GOV$GOV_total_cells, BAU$BAU_total_cells)


## Standardised variables: low-pressure scenario

## Precipitation change and range size are negated before standardisation
## making larger values of every standardised term equal worsening 

## GOV
GOV$PR_z        <- z(GOV$PR_exposure, ref_PR)
GOV$TAS_z       <- z(GOV$TAS_exposure, ref_TAS)
GOV$TASMAX_z    <- z(GOV$TASMAX_exposure, ref_TASMAX)

GOV$PR_change_z     <- z(-GOV$PR_median_change, -ref_PRchg)
GOV$TAS_change_z    <- z(GOV$TAS_median_change, ref_TASchg)
GOV$TASMAX_change_z <- z(GOV$TASMAX_median_change, ref_TASMAXchg)

GOV$defor_percent_z <- z(GOV$GOV_percent, ref_defor)
GOV$total_cells_z   <- z(-GOV$GOV_total_cells, -ref_cells)


## Standardised variables: high-pressure scenario

## BAU 
BAU$PR_z        <- z(BAU$PR_exposure, ref_PR)
BAU$TAS_z       <- z(BAU$TAS_exposure, ref_TAS)
BAU$TASMAX_z    <- z(BAU$TASMAX_exposure, ref_TASMAX)

BAU$PR_change_z     <- z(-BAU$PR_median_change, -ref_PRchg)
BAU$TAS_change_z    <- z(BAU$TAS_median_change, ref_TASchg)
BAU$TASMAX_change_z <- z(BAU$TASMAX_median_change, ref_TASMAXchg)

BAU$defor_percent_z <- z(BAU$BAU_percent, ref_defor)
BAU$total_cells_z   <- z(-BAU$BAU_total_cells, -ref_cells)


#################################################
## STEP 4 - COMPOSITE ENVIRONMENTAL RISK INDEX ##
#################################################

##GOV

GOV$threat_score_no_eoo <- rowMeans(cbind(
  GOV$PR_z,
  GOV$TAS_z,
  GOV$TASMAX_z,
  GOV$PR_change_z,
  GOV$TAS_change_z,
  GOV$TASMAX_change_z,
  GOV$defor_percent_z,
  GOV$defor_percent_z,
  GOV$defor_percent_z),na.rm = TRUE)



##BAU

BAU$threat_score_no_eoo <- rowMeans(cbind(
  BAU$PR_z,
  BAU$TAS_z,
  BAU$TASMAX_z,
  BAU$PR_change_z,
  BAU$TAS_change_z,
  BAU$TASMAX_change_z,
  BAU$defor_percent_z,
  BAU$defor_percent_z,
  BAU$defor_percent_z),na.rm = TRUE)


########################################
## STEP 5 - COMPONENT-SPECIFIC SCORES ##
########################################

##GOV

# Climate threat score = mean of climate parameters
GOV$climate_threat_score <- rowMeans(data.frame(GOV$PR_z, GOV$TAS_z, GOV$TASMAX_z), na.rm=TRUE)
# Deforestation threat score
GOV$deforestation_threat_score <- GOV$defor_percent_z

##BAU

# Climate threat score = mean of climate parameters
BAU$climate_threat_score <- rowMeans(data.frame(BAU$PR_z, BAU$TAS_z, BAU$TASMAX_z), na.rm=TRUE)

# Deforestation threat score
BAU$deforestation_threat_score <- BAU$defor_percent_z


###############################
## STEP 6 - SPECIES RANKINGS ##
###############################

# GOV ranking
GOV_ranked <- GOV[order(-GOV$threat_score_no_eoo), ]

# BAU ranking
BAU_ranked <- BAU[order(-BAU$threat_score_no_eoo), ]

# View top 20 most threatened species
head(GOV_ranked[, c("taxa", "threat_score_no_eoo")], 20)
head(BAU_ranked[, c("taxa", "threat_score_no_eoo")], 20)

## Output
##
## > head(GOV_ranked[, c("taxa", "threat_score_no_eoo")], 20)
## # A tibble: 20 x 2
## taxa                     threat_score_no_eoo
## <chr>                                  <dbl>
##   1 Hansreia grossii                      1.40
## 2 Hansreia krinskii                     1.36
## 3 Canthon obscuriellus                  0.733
## 4 Cryptocanthon genieri                 0.629
## 5 Scybalocanthon kaestneri              0.280
## 6 Canthon arriagadai                    0.0750
## 7 Ateuchus vitulus                      0.0745
## 8 Bdelyrus paraensis                    0.0598
## 9 Coprophanaeus terrali                 0.0342
## 10 Besourenga vejdovskyi                -0.0174
## 11 Ateuchus laevicollis                 -0.0393
## 12 Ateuchus carinifer                   -0.0546
## 13 Hansreia oxygona                     -0.0579
## 14 Canthidium sladeni                   -0.0740
## 15 Canthidium rutilum                   -0.103
## 16 Dichotomius zikani                   -0.104
## 17 Onthophagus xanthomerus              -0.109
## 18 Ateuchus marquesi                    -0.113
## 19 Ateuchus imperator                   -0.117
## 20 Hansreia peugeoti                    -0.117
## > head(BAU_ranked[, c("taxa", "threat_score_no_eoo")], 20)
## # A tibble: 20 x 2
## taxa                      threat_score_no_eoo
## <chr>                                   <dbl>
##   1 Hansreia grossii                         2.82
## 2 Hansreia krinskii                        2.40
## 3 Bdelyrus paraensis                       1.76
## 4 Ateuchus caaetepora                      1.62
## 5 Onthophagus xanthomerus                  1.60
## 6 Ateuchus laevicollis                     1.59
## 7 Canthon arriagadai                       1.51
## 8 Ateuchus marquesi                        1.50
## 9 Deltochilum schefflerorum                1.48
## 10 Eurysternus arnaudi                      1.48
## 11 Anomiopus myrmidon                       1.42
## 12 Canthon xanthopus                        1.40
## 13 Eutrichillum racliffei                   1.39
## 14 Canthon quadrimaculatus                  1.38
## 15 Canthidium dohrni                        1.37
## 16 Hansreia peugeoti                        1.36
## 17 Ateuchus carinifer                       1.34
## 18 Dichotomius zikani                       1.34
## 19 Phanaeus sororibispinus                  1.33
## 20 Dichotomius leopoldinae                  1.32


## ---- 6.2 Centring check ---------------------------------------------------

## Pooled standardisation implies a mean of zero across the two scenarios
## combined; a departure from zero indicates a malformed index.

mean(c(BAU$threat_score_no_eoo, GOV$threat_score_no_eoo))


## ---- 6.3 Export -----------------------------------------------------------

write_xlsx(
  list(GOV = GOV,
       BAU = BAU),
  "Final_Threat_Scores_no_EOO.SSP.matched.All.Species_2026.xlsx"
)


#################################
### STEP 7. DESCRIPTIVE STATS ###
#################################


BAU %>%
  summarise(
    mean_score = mean(threat_score_no_eoo, na.rm = TRUE),
    sd_score   = sd(threat_score_no_eoo,   na.rm = TRUE),
    min_score  = min(threat_score_no_eoo,  na.rm = TRUE),
    max_score  = max(threat_score_no_eoo,  na.rm = TRUE)
  ) %>% mutate(scenario = "BAU")

## Output
##
## # A tibble: 1 x 5
## mean_score sd_score min_score max_score scenario
## <dbl>    <dbl>     <dbl>     <dbl> <chr>
##   1      0.569    0.595    -0.913      2.82 BAU


GOV %>%
  summarise(
    mean_score = mean(threat_score_no_eoo, na.rm = TRUE),
    sd_score   = sd(threat_score_no_eoo,   na.rm = TRUE),
    min_score  = min(threat_score_no_eoo,  na.rm = TRUE),
    max_score  = max(threat_score_no_eoo,  na.rm = TRUE)
  ) %>% mutate(scenario = "GOV")

## Output
##
## # A tibble: 1 x 5
## mean_score sd_score min_score max_score scenario
## <dbl>    <dbl>     <dbl>     <dbl> <chr>
##   1     -0.569    0.381     -1.40      1.40 GOV



#Per threat category
BAU %>%
  group_by(Category_CriteriaB) %>%
  summarise(
    n          = n(),
    mean_score = round(mean(threat_score_no_eoo, na.rm = TRUE), 3),
    sd_score   = round(sd(threat_score_no_eoo,   na.rm = TRUE), 3),
    min_score  = round(min(threat_score_no_eoo,  na.rm = TRUE), 3),
    max_score  = round(max(threat_score_no_eoo,  na.rm = TRUE), 3)
  ) %>% mutate(scenario = "BAU")

## Output
##
## A tibble: 4 x 7
## Category_CriteriaB     n mean_score sd_score min_score max_score scenario
## <chr>              <int>      <dbl>    <dbl>     <dbl>     <dbl> <chr>
## 1 CR                     1      1.32    NA         1.32       1.32 BAU
## 2 EN                    45      0.584    0.719    -0.49       2.82 BAU
## 3 LC or NT             151      0.577    0.499    -0.913      1.50 BAU
## 4 VU                    70      0.53     0.697    -0.877      1.76 BAU

GOV %>%
  group_by(Category_CriteriaB) %>%
  summarise(
    n          = n(),
    mean_score = round(mean(threat_score_no_eoo, na.rm = TRUE), 3),
    sd_score   = round(sd(threat_score_no_eoo,   na.rm = TRUE), 3),
    min_score  = round(min(threat_score_no_eoo,  na.rm = TRUE), 3),
    max_score  = round(max(threat_score_no_eoo,  na.rm = TRUE), 3)
  ) %>% mutate(scenario = "GOV")

## Output
##
## Category_CriteriaB     n mean_score sd_score min_score max_score scenario
## <chr>              <int>      <dbl>    <dbl>     <dbl>     <dbl> <chr>
## 1 CR                     1      0.733   NA         0.733     0.733 GOV
## 2 EN                    45     -0.404    0.53     -1.23      1.40  GOV
## 3 LC or NT             151     -0.638    0.28     -1.40     -0.074 GOV
## 4 VU                    70     -0.543    0.395    -1.27      0.28  GOV



#Scenario Difference Wilcox

wilcox.test(BAU$threat_score_no_eoo, GOV$threat_score_no_eoo, paired = TRUE)

## Output
##
## Wilcoxon signed rank test with continuity correction
##
## data:  BAU$threat_score_no_eoo and GOV$threat_score_no_eoo
## V = 35777, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0


median(BAU$threat_score_no_eoo - GOV$threat_score_no_eoo)

## 1.21912


#########################################################################
#### Ordinal regression Sensitivity Test for Deforestation Weighting ####
#########################################################################

##Ordering
BAU$Category_CriteriaB <- factor(BAU$Category_CriteriaB,
                                 levels  = c("LC or NT", "VU", "EN", "CR"),
                                 ordered = TRUE)
GOV$Category_CriteriaB <- factor(GOV$Category_CriteriaB,
                                 levels  = c("LC or NT", "VU", "EN", "CR"),
                                 ordered = TRUE)



## BRANT ##
#Proportional Odds Assessments

brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = BAU, Hess = TRUE))
brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = GOV, Hess = TRUE))

## Output
##
## > brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = BAU, Hess = TRUE))
## -----------------------------------------------
##   Test for             X2      df      probability
## -----------------------------------------------
##   Omnibus              2.21    2       0.33
## threat_score_no_eoo  2.21    2       0.33
## -----------------------------------------------
##
##   H0: Parallel Regression Assumption holds
## > brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = GOV, Hess = TRUE))
## -----------------------------------------------
##   Test for             X2      df      probability
## -----------------------------------------------
##   Omnibus              2.22    2       0.33
## threat_score_no_eoo  2.22    2       0.33
## -----------------------------------------------
##
##   H0: Parallel Regression Assumption holds


##Coefficients

ctableBAU <- coef(summary(modelBAU <- MASS::polr(Category_CriteriaB ~ threat_score_no_eoo, data = BAU, Hess=TRUE)))

ctableGOV <- coef(summary(modelGOV <- MASS::polr(Category_CriteriaB ~ threat_score_no_eoo, data = GOV, Hess=TRUE)))

#
pvalsBAU <- pnorm(abs(ctableBAU[, "t value"]), lower.tail = FALSE) * 2

pvalsGOV <- pnorm(abs(ctableGOV[, "t value"]), lower.tail = FALSE) * 2
#
cbind(ctableGOV, "p value" = pvalsGOV)
cbind(ctableBAU, "p value" = pvalsBAU)

## Output
##
## > cbind(ctableGOV, "p value" = pvalsGOV)
## Value Std. Error   t value      p value
## threat_score_no_eoo  1.3228357  0.3309403  3.997204 6.409515e-05
## LC or NT|VU         -0.4848183  0.2209201 -2.194542 2.819649e-02
## VU|EN                0.8895947  0.2297343  3.872276 1.078237e-04
## EN|CR                5.0183565  1.0083222  4.976937 6.459818e-07
## > cbind(ctableBAU, "p value" = pvalsBAU)
## Value Std. Error     t value      p value
## threat_score_no_eoo -0.01244172  0.1978059 -0.06289862 9.498472e-01
## LC or NT|VU          0.25671875  0.1657733  1.54861330 1.214747e-01
## VU|EN                1.56260094  0.1958024  7.98049860 1.457434e-15
## EN|CR                5.57647177  1.0080178  5.53211631 3.163900e-08



# DEFORESTATION ONLY

brant(polr(Category_CriteriaB ~ defor_percent_z, data = BAU, Hess = TRUE))
brant(polr(Category_CriteriaB ~ defor_percent_z, data = GOV, Hess = TRUE))

## Output
##
## > brant(polr(Category_CriteriaB ~ defor_percent_z, data = BAU, Hess = TRUE))
## -------------------------------------------
##   Test for         X2      df      probability
## -------------------------------------------
##   Omnibus          0.06    2       0.97
## defor_percent_z  0.06    2       0.97
## -------------------------------------------
##
##   H0: Parallel Regression Assumption holds
## > brant(polr(Category_CriteriaB ~ defor_percent_z, data = GOV, Hess = TRUE))
## -------------------------------------------
##   Test for         X2      df      probability
## -------------------------------------------
##   Omnibus          0.87    2       0.65
## defor_percent_z  0.87    2       0.65
## -------------------------------------------
##
##   H0: Parallel Regression Assumption holds

ctableBAU <- coef(summary(modelBAU <- MASS::polr(Category_CriteriaB ~ defor_percent_z, data = BAU, Hess=TRUE)))

ctableGOV <- coef(summary(modelGOV <- MASS::polr(Category_CriteriaB ~ defor_percent_z, data = GOV, Hess=TRUE)))

#
pvalsBAU <- pnorm(abs(ctableBAU[, "t value"]), lower.tail = FALSE) * 2

pvalsGOV <- pnorm(abs(ctableGOV[, "t value"]), lower.tail = FALSE) * 2
#
cbind(ctableGOV, "p value" = pvalsGOV)
cbind(ctableBAU, "p value" = pvalsBAU)

## Output
##
## Value Std. Error    t value      p value
## defor_percent_z  0.57067687  0.1955769  2.9179149 3.523805e-03
## LC or NT|VU     -0.05941426  0.1645542 -0.3610619 7.180532e-01
## VU|EN            1.28687380  0.1865428  6.8985454 5.253772e-12
## EN|CR            5.36957847  1.0032743  5.3520544 8.696123e-08
## Value Std. Error   t value      p value
## defor_percent_z -0.02909779  0.1279384 -0.227436 8.200847e-01
## LC or NT|VU      0.24694391  0.1434587  1.721359 8.518564e-02
## VU|EN            1.55311642  0.1771210  8.768677 1.807778e-18
## EN|CR            5.56692938  1.0041216  5.544079 2.955051e-08

#---------------------------------------------------------------#


#####################################################
## STEP 9 - SENSITIVITY TO DEFORESTATION WEIGHTING ##
#####################################################

### Sensitivity Checks for Deforestation Weight  ###


calc_score <- function(data, w_clim, w_defor){
  
  climate <- rowMeans(cbind(data$PR_z,
                            data$TAS_z,
                            data$TASMAX_z,
                            data$PR_change_z,
                            data$TAS_change_z,
                            data$TASMAX_change_z),
                      na.rm = TRUE)
  
  defor <- data$defor_percent_z
  
  score <- (w_clim * climate) + (w_defor * defor)
  
  return(score)
}


##BAU
#ordered response variable
BAU$Category_ord <- factor(
  BAU$Category_CriteriaB,
  levels = c("LC or NT", "VU", "EN", "CR"),
  ordered = TRUE
)

weights <- seq(0.1, 0.9, by = 0.1)

resultsBAU <- data.frame()
for(w in weights){
  BAU$score <- calc_score(BAU, 1 - w, w)
  temp <- na.omit(data.frame(Category_ord = BAU$Category_ord, score = BAU$score))
  model      <- polr(Category_ord ~ score, data = temp, Hess = TRUE)
  coef_table <- coef(summary(model))
  beta       <- coef_table["score", "Value"]
  se         <- coef_table["score", "Std. Error"]
  z_value    <- beta / se
  p_value    <- 2 * (1 - pnorm(abs(z_value)))
  resultsBAU <- rbind(resultsBAU,
                      data.frame(defor_weight = w, beta = beta,
                                 se = se, z = z_value, p = p_value))
}


resultsBAU

## Output
##
## defor_weight         beta        se           z         p
## 1          0.1  0.008211416 0.2049302  0.04006934 0.9680378
## 2          0.2 -0.001028319 0.2046304 -0.00502525 0.9959904
## 3          0.3 -0.009799896 0.2001220 -0.04896961 0.9609435
## 4          0.4 -0.017172809 0.1921942 -0.08935132 0.9288027
## 5          0.5 -0.022680562 0.1820027 -0.12461660 0.9008271
## 6          0.6 -0.026341559 0.1706909 -0.15432312 0.8773550
## 7          0.7 -0.028454504 0.1591577 -0.17878182 0.8581090
## 8          0.8 -0.029392955 0.1479969 -0.19860516 0.8425716
## 9          0.9 -0.029489785 0.1375413 -0.21440671 0.8302299


##GOV
# ordered response variable
GOV$Category_ord <- factor(
  GOV$Category_CriteriaB,
  levels = c("LC or NT", "VU", "EN", "CR"),
  ordered = TRUE
)

weights <- seq(0.1, 0.9, by = 0.1)

resultsGOV <- data.frame()

resultsGOV <- data.frame()
for(w in weights){
  GOV$score <- calc_score(GOV, 1 - w, w)
  temp <- na.omit(data.frame(Category_ord = GOV$Category_ord, score = GOV$score))
  model      <- polr(Category_ord ~ score, data = temp, Hess = TRUE)
  coef_table <- coef(summary(model))
  beta       <- coef_table["score", "Value"]
  se         <- coef_table["score", "Std. Error"]
  z_value    <- beta / se
  p_value    <- 2 * (1 - pnorm(abs(z_value)))
  resultsGOV <- rbind(resultsGOV,
                      data.frame(defor_weight = w, beta = beta,
                                 se = se, z = z_value, p = p_value))
}

resultsGOV

## Output
##
## defor_weight      beta        se        z            p
## 1          0.1 1.0209568 0.2957337 3.452284 5.558619e-04
## 2          0.2 1.1883828 0.3162262 3.758015 1.712667e-04
## 3          0.3 1.3035255 0.3291186 3.960656 7.474400e-05
## 4          0.4 1.3272957 0.3299636 4.022552 5.757087e-05
## 5          0.5 1.2533587 0.3170013 3.953797 7.692078e-05
## 6          0.6 1.1148213 0.2940317 3.791501 1.497398e-04
## 7          0.7 0.9563032 0.2669403 3.582461 3.403722e-04
## 8          0.8 0.8060357 0.2402234 3.355359 7.926204e-04
## 9          0.9 0.6770427 0.2162381 3.131006 1.742085e-03


##########################################
## STEP 10 - SPEARMAN RANK CORRELATIONS ##
##########################################

#Spearman Robustness
##GOV
severity <- as.numeric(factor(GOV$Category_CriteriaB,
                              levels = c("LC or NT", "VU","EN","CR")))
cor.test(GOV$threat_score_no_eoo, severity, method="spearman")

## Output
##
## Spearman's rank correlation rho
##
## data:  GOV$threat_score_no_eoo and severity
## S = 2507115, p-value = 0.0005633
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho
## 0.2096893


##BAU
severity <- as.numeric(factor(BAU$Category_CriteriaB,
                              levels = c("LC or NT", "VU","EN","CR")))
cor.test(BAU$threat_score_no_eoo, severity, method="spearman")

## Output
##
## Spearman's rank correlation rho
##
## data:  BAU$threat_score_no_eoo and severity
## S = 3221933, p-value = 0.7992
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##         rho
## -0.01564076


#######################################################################
## STEP 11 - VERSION EXCLUDING THE CRITICALLY ENDANGERED SPECIES ##
#######################################################################


GOV_noCR <- droplevels(subset(GOV, Category_CriteriaB != "CR"))
BAU_noCR <- droplevels(subset(BAU, Category_CriteriaB != "CR"))

levels(GOV_noCR$Category_CriteriaB)   # should be 3
nrow(GOV_noCR)

levels(BAU_noCR$Category_CriteriaB)   # shoould be 3
nrow(BAU_noCR)


## Coefficients 

ctableBAU_noCR <- coef(summary(modelBAU_noCR <- MASS::polr(Category_CriteriaB ~ threat_score_no_eoo, data = BAU_noCR, Hess=TRUE)))

ctableGOV_noCR <- coef(summary(modelGOV_noCR <- MASS::polr(Category_CriteriaB ~ threat_score_no_eoo, data = GOV_noCR, Hess=TRUE)))


pvalsBAU_noCR <- pnorm(abs(ctableBAU_noCR[, "t value"]), lower.tail = FALSE) * 2

pvalsGOV_noCR <- pnorm(abs(ctableGOV_noCR[, "t value"]), lower.tail = FALSE) * 2

cbind(ctableGOV_noCR, "p value" = pvalsGOV_noCR)
cbind(ctableBAU_noCR, "p value" = pvalsBAU_noCR)

## Output
##
## > cbind(ctableGOV_noCR, "p value" = pvalsGOV_noCR)
## Value Std. Error   t value      p value
## threat_score_no_eoo  1.2088750  0.3321351  3.639708 2.729469e-04
## LC or NT|VU         -0.4144097  0.2223462 -1.863804 6.234924e-02
## VU|EN                0.9653711  0.2322273  4.157010 3.224397e-05
## > cbind(ctableBAU_noCR, "p value" = pvalsBAU_noCR)
## Value Std. Error    t value      p value
## threat_score_no_eoo -0.04218651  0.1986342 -0.2123829 8.318083e-01
## LC or NT|VU          0.24878930  0.1659545  1.4991412 1.338370e-01
## VU|EN                1.56813189  0.1967354  7.9707644 1.576959e-15


##Proportional odds

brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = GOV_noCR, Hess = TRUE))
brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = BAU_noCR, Hess = TRUE))

## Output
##
## brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = GOV_noCR, Hess = TRUE))
## -----------------------------------------------
##   Test for             X2      df      probability
## -----------------------------------------------
##   Omnibus              0.29    1       0.59
## threat_score_no_eoo  0.29    1       0.59
## -----------------------------------------------
##
##   H0: Parallel Regression Assumption holds
## > brant(polr(Category_CriteriaB ~ threat_score_no_eoo, data = BAU_noCR, Hess = TRUE))
## -----------------------------------------------
##   Test for             X2      df      probability
## -----------------------------------------------
##   Omnibus              0.3     1       0.58
## threat_score_no_eoo  0.3     1       0.58
## -----------------------------------------------
##
##   H0: Parallel Regression Assumption holds


################
## Lets Plot! ##
################

resultsBAU$scenario <- "Business-as-Usual"
resultsGOV$scenario <- "Governance"

combined_results <- bind_rows(resultsBAU, resultsGOV) %>%
  mutate(
    ci_lo    = beta - 1.96 * se,
    ci_hi    = beta + 1.96 * se,
    scenario = factor(scenario, levels = c("Business-as-Usual", "Governance"))
  )

scenario_colors <- c("Business-as-Usual" = "darkred",
                     "Governance"         = "orange2")

#Plot
windows()
ggplot(combined_results,
       aes(x = defor_weight, y = beta,
           colour = scenario, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_colour_manual(values = scenario_colors) +
  scale_fill_manual(values = scenario_colors) +
  labs(
    x      = "Deforestation weight",
    y      = expression(paste("Regression coefficient (", beta, ")")),
    colour = NULL,
    fill   = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text      = element_text(size = 9,  colour = "black"),
    axis.title     = element_text(size = 10, colour = "black"),
    axis.line      = element_line(colour = "black", linewidth = 0.4),
    axis.ticks     = element_line(colour = "black", linewidth = 0.4),
    legend.text    = element_text(size = 9),
    legend.position = "top",
    plot.margin    = margin(5, 5, 5, 5, "mm")
  )



#Boxplot

scenario_colors <- c("BAU" = "darkred",
                     "GOV"         = "orange2")

GOV$Scenario <- "GOV"
BAU$Scenario <- "BAU"

plot_data <- bind_rows(GOV, BAU)

plot_data$Category_CriteriaB <- factor(plot_data$Category_CriteriaB,
                                       levels = c("LC or NT","VU","EN", "CR"),
                                       ordered = TRUE)
ord_labels <- data.frame(
  Scenario = c("BAU", "GOV"),
  label = c("High-pressure",
            "Low-pressure"),
  x = c(1, 1),
  y = c(max(plot_data$threat_score_no_eoo, na.rm = TRUE) * 1.15,
        max(plot_data$threat_score_no_eoo, na.rm = TRUE) * 1.05)
)

windows()
ggplot(plot_data, aes(x = Category_CriteriaB,
                      y = threat_score_no_eoo,
                      fill = Scenario)) +
  geom_boxplot(aes(colour = Scenario),
               outlier.shape = NA,
               alpha         = 0.3,
               width         = 0.6,
               linewidth     = 0.1,
               position      = position_dodge(width = 0.75)
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size     = 0.8,
    alpha    = 0.4,
    aes(colour = Scenario),
    show.legend = FALSE
  ) +
  stat_summary(
    aes(x = as.numeric(Category_CriteriaB), colour = Scenario),
    fun = mean, geom = "line",
    position = position_dodge(width = 0.75),
    linewidth = 1.2, alpha = 0.6, lineend = "round", linejoin = "round",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = scenario_colors,
                    labels = c("BAU" = "High-pressure (SSP5-8.5)",
                               "GOV" = "Low-pressure (SSP1-2.6)")) +
  scale_colour_manual(values = scenario_colors, guide = "none") +
  labs(
    x = "Potential Threat Category (IUCN Criterion B)",
    y = "Species-specific Environmental Risk Index",
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title      = element_text(size = 9.5, colour = "black"),
    axis.text       = element_text(size = 8.5, colour = "black"),
    axis.line       = element_line(colour = "black", linewidth = 0.4),
    axis.ticks      = element_line(colour = "black", linewidth = 0.4),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(3.5, "mm"),
    plot.margin     = margin(4, 6, 4, 4, "mm")
  )



