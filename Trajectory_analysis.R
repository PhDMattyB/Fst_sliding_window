##################################################################
## Phenotypic trajectory analysis
##
## Matt Brachmann (PhDMattyB)
##
## 2022-04-04
##
##################################################################

setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

## Load the data manipulation work horse
library(tidyverse)
library(geomorph)
library(RRPP)

PWS_data = read_csv('GSTV_PartialWarpScores_AllometryExcluded.csv') 

Geno_indv = read_csv('GSTVMF_Morph_Eco_Geno.csv') %>% 
  filter(Lake %in% c('Galtabol', 
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  na.omit() %>% 
  rename(SpecimenID = Specimen.ID)

PWS_data_clean = inner_join(PWS_data, 
           Geno_indv, 
           by = 'SpecimenID') %>% 
  filter(LaMorph.x != 'S.LGB2', 
         LaMorph.x != 'T.PL2') %>% 
  rename(Lake = Lake.x, 
         Morph = Morph.x, 
         BP = BP.x, 
         LaMorph = LaMorph.x, 
         Sex = Sex.x)


rrpp = rrpp.data.frame(PWS_data_clean)

Bodyshape_PWS = cbind(rrpp$PW19X, 
                      rrpp$PW19Y, 
                      rrpp$PW18X, 
                      rrpp$PW18Y, 
                      rrpp$PW17X,
                      rrpp$PW17Y, 
                      rrpp$PW16X,
                      rrpp$PW16Y, 
                      rrpp$PW15X, 
                      rrpp$PW15Y, 
                      rrpp$PW14X, 
                      rrpp$PW14Y, 
                      rrpp$PW13X, 
                      rrpp$PW13Y,
                      rrpp$PW12X, 
                      rrpp$PW12Y, 
                      rrpp$PW11X, 
                      rrpp$PW11Y, 
                      rrpp$PW10X, 
                      rrpp$PW10Y, 
                      rrpp$PW9X, 
                      rrpp$PW9Y, 
                      rrpp$PW8X, 
                      rrpp$PW8Y, 
                      rrpp$PW7X, 
                      rrpp$PW7Y, 
                      rrpp$PW6X, 
                      rrpp$PW6Y, 
                      rrpp$PW5X, 
                      rrpp$PW5Y, 
                      rrpp$PW4X, 
                      rrpp$PW4Y, 
                      rrpp$PW3X, 
                      rrpp$PW3Y, 
                      rrpp$PW2X, 
                      rrpp$PW2Y, 
                      rrpp$PW1X, 
                      rrpp$PW1Y, 
                      rrpp$UNIX, 
                      rrpp$UNIY, 
                      rrpp$CS)
BP_Morph_Pair = rrpp$BP
Lake = rrpp$Lake


