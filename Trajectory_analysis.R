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
