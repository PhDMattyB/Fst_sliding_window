##############################
## Fst Chromosome fix
##
## Matt Brachmann (PhDMattyB)
##
## 2020-04-09
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory/')

library(tidyverse)

MAP = read_tsv('Feb202019_Poly_Plink_input.map') %>% 
  arrange(`Marker ID`)

MAP3 = mutate(.data = MAP,
              chromosome = as.factor(case_when(
                `#Chromosome` == '1' ~ '1',
                `#Chromosome` == '2' ~ '2',
                `#Chromosome` == '3' ~ '3',
                `#Chromosome` == '4' ~ '4',
                `#Chromosome` == '5' ~ '5',
                `#Chromosome` == '6' ~ '6',
                `#Chromosome` == '7' ~ '7',
                `#Chromosome` == '8' ~ '8',
                `#Chromosome` == '9' ~ '9',
                `#Chromosome` == '10' ~ '10',
                `#Chromosome` == '11' ~ '11',
                `#Chromosome` == '12' ~ '12',
                `#Chromosome` == '13' ~ '13',
                `#Chromosome` == '14' ~ '14',
                `#Chromosome` == '15' ~ '15',
                `#Chromosome` == '16' ~ '16',
                `#Chromosome` == '17' ~ '17',
                `#Chromosome` == '18' ~ '18',
                `#Chromosome` == '19' ~ '19',
                `#Chromosome` == '20' ~ '20',
                `#Chromosome` == '21' ~ '21',
                `#Chromosome` == '22' ~ '22',
                `#Chromosome` == '23' ~ '23',
                `#Chromosome` == '24' ~ '24',
                `#Chromosome` == '25' ~ '25',
                `#Chromosome` == '26' ~ '26',
                `#Chromosome` == '27' ~ '27',
                `#Chromosome` == '28' ~ '28',
                `#Chromosome` == '29' ~ '29',
                `#Chromosome` == '30' ~ '30',
                `#Chromosome` == '31' ~ '31',
                `#Chromosome` == '32' ~ '32',
                `#Chromosome` == '33' ~ '33',
                `#Chromosome` == '34' ~ '34',
                `#Chromosome` == '35' ~ '35',
                `#Chromosome` == '36' ~ '36',
                `#Chromosome` == '37' ~ '37',
                `#Chromosome` == '38' ~ '38',
                `#Chromosome` == '39' ~ '39',
                `#Chromosome` > '39' ~ '0')))

MAP3$chromosome[is.na(MAP3$chromosome)] = '0'

MAP_fixed = MAP3 %>% 
  group_by(chromosome) %>% 
  select(-`#Chromosome`) %>% 
  rename(`#Chromosome` = chromosome) %>% 
  select(`#Chromosome`, everything())

MAP_fixed %>% 
  select(`Physical position`) 

  
write_tsv(MAP_fixed, 
          'MAP_Fixed_Chromosome.MAP')


# nonstandard_chromes = mutate(.data = MAP_fixed,
#               ns_chromosomes = as.factor(case_when(
#                 `#Chromosome` == '1' ~ 'AC01',
#                 `#Chromosome` == '2' ~ 'AC02',
#                 `#Chromosome` == '3' ~ 'AC03',
#                 `#Chromosome` == '4' ~ 'AC04p',
#                 `#Chromosome` == '5' ~ 'AC04q.1:29',
#                 `#Chromosome` == '6' ~ 'AC04q.2',
#                 `#Chromosome` == '7' ~ 'AC05',
#                 `#Chromosome` == '8' ~ 'AC06',
#                 `#Chromosome` == '9' ~ 'AC06',
#                 `#Chromosome` == '10' ~ 'AC07',
#                 `#Chromosome` == '11' ~ 'AC08',
#                 `#Chromosome` == '12' ~ 'AC09',
#                 `#Chromosome` == '13' ~ 'AC10',
#                 `#Chromosome` == '14' ~ 'AC11',
#                 `#Chromosome` == '15' ~ 'AC12',
#                 `#Chromosome` == '16' ~ 'AC13',
#                 `#Chromosome` == '17' ~ 'AC14',
#                 `#Chromosome` == '18' ~ 'AC15',
#                 `#Chromosome` == '19' ~ 'AC16',
#                 `#Chromosome` == '20' ~ 'AC17',
#                 `#Chromosome` == '21' ~ 'AC18',
#                 `#Chromosome` == '22' ~ 'AC19',
#                 `#Chromosome` == '23' ~ 'AC20',
#                 `#Chromosome` == '24' ~ 'AC21',
#                 `#Chromosome` == '25' ~ 'AC22',
#                 `#Chromosome` == '26' ~ 'AC23',
#                 `#Chromosome` == '27' ~ 'AC24',
#                 `#Chromosome` == '28' ~ 'AC25',
#                 `#Chromosome` == '29' ~ 'AC26',
#                 `#Chromosome` == '30' ~ 'AC27',
#                 `#Chromosome` == '31' ~ 'AC28',
#                 `#Chromosome` == '32' ~ 'AC30',
#                 `#Chromosome` == '33' ~ 'AC31',
#                 `#Chromosome` == '34' ~ 'AC32',
#                 `#Chromosome` == '35' ~ 'AC33',
#                 `#Chromosome` == '36' ~ 'AC34',
#                 `#Chromosome` == '37' ~ 'AC35',
#                 `#Chromosome` == '38' ~ 'AC36',
#                 `#Chromosome` == '39' ~ 'AC37',
#                 `#Chromosome` == '40' ~ 'Contigs')))
# 
# nonstandard_chromes %>% 
#   select(ns_chromosomes)
