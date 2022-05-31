##################################################################
## Fixing data frames for submission to dryad
##
## Matt Brachmann (PhDMattyB)
##
## 2022-05-30
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)

setwd('~/PhD_Genomics_Chapter3/Fst_Iceland_pops/')

GSBPI = read_csv('GSBPI_RDA_outlier_clean.csv') %>% 
  mutate(BP_Morph_Pair = 'Galtabol_BP')
SLGBPL = read_csv('SLGBPL_RDA_outlier_clean.csv') %>% 
  mutate(BP_Morph_Pair = 'Svinavatn_BP')
TSBPL = read_csv('TSBPL_RDA_outlier_clean.csv') %>% 
  mutate(BP_Morph_Pair = 'Thingvallavatn_BP2')
TLGBPL = read_csv('TLGBPL_RDA_outlier_clean.csv') %>% 
  mutate(BP_Morph_Pair = 'Thingvallavatn_BP1')

RDA_outliers_clean = bind_rows(GSBPI, 
                               SLGBPL, 
                               TSBPL, 
                               TLGBPL)

##test
RDA_outliers_clean %>% 
  filter(BP_Morph_Pair == 'Thingvallavatn_BP1')

RDA_outliers_clean %>% 
  write_csv('RDA_Outliers_All_Pops_Combined.csv')
