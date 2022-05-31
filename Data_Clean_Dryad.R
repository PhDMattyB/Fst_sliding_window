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


# RDA outlier clean up ----------------------------------------------------


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

# FST clean up ------------------------------------------------------------

GSBPI_data = read_tsv('Galtabol_chr_fix.fst')%>% 
  mutate(BP_Morph_Pair = 'Galtabol_BP')
TLGBPL_data = read_tsv('TLGBPL_fst.fst')%>% 
  mutate(BP_Morph_Pair = 'Thingvallavatn_BP1')
TSBPL_data = read_tsv('TSBPL_fst.fst')%>% 
  mutate(BP_Morph_Pair = 'Thingvallavatn_BP2')
VBRSIL_data = read_tsv('VBRSIL_fst.fst')%>% 
  mutate(BP_Morph_Pair = 'Vatnshlidarvatn_BP')
SLGBPEL_data = read_tsv('SLGBPEL_fst.fst')%>% 
  mutate(BP_Morph_Pair = 'Svinavatn_BP')

Fst_data_clean = bind_rows(GSBPI_data, 
                           SLGBPEL_data, 
                           TSBPL_data, 
                           TLGBPL_data,
                           VBRSIL_data)

##test
Fst_data_clean %>% 
  filter(BP_Morph_Pair == 'Thingvallavatn_BP1')

Fst_data_clean %>% 
  write_csv('FST_data_All_Pops_Combined.csv')

