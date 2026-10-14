

###############

setwd("G:/Meu Drive/Artigo_Scarabaeinae_Amazonia_Pollock/All pontos")
dir()

Comp_Amazon_Clamm_all<- readr::read_csv("Comp_Amazon_ID_xy_all_pontos.csv")
Comp_Amazon_Clamm_all
nrow(Comp_Amazon_Clamm_all)
ncol(Comp_Amazon_Clamm_all)

comp<-Comp_Amazon_Clamm_all[,5:607]
ncol(comp)



Site<-Comp_Amazon_Clamm_all$Site


Clamm_AM_all <- clamtest(comp, groups=Site,alpha=0.05)

summary(Clamm_AM_all)
head(Clamm_AM_all)
plot(Clamm_AM_all)

setwd("G:/Meu Drive/Artigo_Scarabaeinae_Amazonia_Pollock/All pontos")
write.csv(Clamm_AM_all, "Clamm_AM_all_pontos_603_spp_correto.csv", row.names = FALSE)

##################################################################################


indval_Amazon_all_pontos = multipatt(comp, Site)

summary(indval_Amazon_all_pontos,alpha=0.05) 
