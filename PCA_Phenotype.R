##################################################################
## Body shape PCA
##
## Matt Brachmann (PhDMattyB)
##
## 2022-04-04
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)


setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

Phenotype = read_csv('RDA_Morphological_variables_PCA.csv')
Phenotype = read_csv('RDA_Morphological_Variables_AllPops.csv')
Phenotype = read_csv('GSTVMF_Morph_Eco_Geno.csv')

Phenotype %>% 
  select(Lake) %>% 
  distinct()

Phenotype_clean = Phenotype %>% 
  filter(Lake %in% c('Galtabol', 
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  select(1:4, 
         LaMorph, 
         Fork.length, 
         CO_NA_PCA_1, 
         CO_NA_PCA_2, 
         NA_PCA_1, 
         NA_PCA_2)
# Phenotype_clean = mutate(.data = Phenotype, 
#                          Population_code = as.factor(case_when(
#                            POP_name == 'G.SB' ~ 'Galtabol', 
#                            POP_name == 'G.PI' ~ 'Galtabol', 
#                            POP_name == 'S.LGB' ~ 'Svinavatn', 
#                            POP_name == 'S.PL' ~ 'Svinavatn', 
#                            POP_name == 'S.PI' ~ 'Svinavatn', 
#                            POP_name == 'T.LGB' ~ 'Thingvallavatn', 
#                            POP_name == 'T.SB' ~ 'Thingvallavatn', 
#                            POP_name == 'T.PL' ~ 'Thingvallavatn', 
#                            POP_name == 'V.BR' ~ 'Vatnshlidarvatn', 
#                            POP_name == 'V.SIL' ~ 'Vatnshlidarvatn')))

# Phenotype_clean = mutate(.data = Phenotype, 
#                          Population_code = as.factor(case_when(
#                            POP_name == 'G.SB' ~ 'Galtabol', 
#                            POP_name == 'G.PI' ~ 'Galtabol', 
#                            POP_name == 'S.LGB' ~ 'Svinavatn', 
#                            POP_name == 'S.PEL' ~ 'Svinavatn', 
#                            POP_name == 'S.PEL' ~ 'Svinavatn', 
#                            POP_name == 'T.LGB' ~ 'Thingvallavatn', 
#                            POP_name == 'T.SB' ~ 'Thingvallavatn', 
#                            POP_name == 'T.PL' ~ 'Thingvallavatn', 
#                            POP_name == 'V.BR' ~ 'Vatnshlidarvatn', 
#                            POP_name == 'V.SIL' ~ 'Vatnshlidarvatn')))

PCA_colours = c('#467BB3', 
                '#FF1E0C', 
                '#108565', 
                '#D1BA0A')
theme_set(theme_bw())


ggplot(Phenotype_clean, 
       aes(x = CO_NA_PCA_1, 
           y = CO_NA_PCA_2))+
  geom_point(aes(col = Lake))+
  scale_color_manual(values = PCA_colours)
