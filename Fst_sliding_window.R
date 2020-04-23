##############################
## Fst sliding window
##
## Matt Brachmann (PhDMattyB)
##
## 2020-04-01
##
##############################
# dir.create('~/Fst_sliding_window/Galtabol_fst/')

library(devtools)
library(windowscanr)
library(patchwork)
library(tidyverse)

setwd('~/Fst_sliding_window/Galtabol_fst/')
setwd('~/Fst_sliding_window/TLGBPL/')
setwd('~/Fst_sliding_window/TSBPL/')
setwd('~/Fst_sliding_window/VBRSIL/')
setwd('~/Fst_sliding_window/SLGBPL/')
setwd('~/Fst_sliding_window/SLGBPI/')

# data --------------------------------------------------------------------
## This is the data needed for the sliding window analysis. 
## We get this data from PLINK 1.9 --fst --within flags
## Need to make sure we specify the --chr-set to 39
## and make sure the contigs or unplaced regions are 
## set to 0
data = read_tsv('Galtabol_chr_fix.fst')
data = read_tsv('TLGBPL_fst.fst')
data = read_tsv('TSBPL_fst.fst')
data = read_tsv('VBRSIL_fst.fst')
data = read_tsv('SLGBPL_fst.fst')
data = read_tsv('SLGBPI_fst.fst')

# Sliding window analysis -------------------------------------------------

## This is performed on the lab computer
## it makes mine get emotional. 
## Use the data above, from PLINK, to perform the analysis

fst_position = winScan(x = data, 
                       groups = 'CHR', 
                       position = 'POS',
                       values = 'FST', 
                       win_size = 200000, 
                       win_step = 1000, 
                       funs = c('mean', 'sd'))

# write avg fst per window ------------------------------------------------
## Write the txt file for each window size. 
## Need to compare the different window sizes to see which one
## is the most appropriate. 
## small window size == greater chance for false positives
## large window size == less chance to find differences

write_tsv(fst_position, 
          'SLGBPI_Fst_200Kb_window.txt')

# Plot data ---------------------------------------------------------------

setwd('~/Fst_sliding_window/Galtabol_fst/')
setwd('~/Fst_sliding_window/TLGBPL/') 
setwd('~/Fst_sliding_window/TSBPL/')
setwd('~/Fst_sliding_window/VBRSIL/')
setwd('~/Fst_sliding_window/SLGBPL/')
setwd('~/Fst_sliding_window/SLGBPI/')

data = read_tsv('VBRSIL_Fst_200Kb_window.txt') %>% 
  mutate(AC_CHR = as.factor(case_when(
           CHR == '1' ~ 'AC01',
           CHR == '2' ~ 'AC02',
           CHR == '3' ~ 'AC03',
           CHR == '4' ~ 'AC04p',
           CHR == '5' ~ 'AC04q.1:29',
           CHR == '6' ~ 'AC04q.2',
           CHR == '7' ~ 'AC05',
           CHR == '8' ~ 'AC06',
           CHR == '9' ~ 'AC06',
           CHR == '10' ~ 'AC07',
           CHR == '11' ~ 'AC08',
           CHR == '12' ~ 'AC09',
           CHR == '13' ~ 'AC10',
           CHR == '14' ~ 'AC11',
           CHR == '15' ~ 'AC12',
           CHR == '16' ~ 'AC13',
           CHR == '17' ~ 'AC14',
           CHR == '18' ~ 'AC15',
           CHR == '19' ~ 'AC16',
           CHR == '20' ~ 'AC17',
           CHR == '21' ~ 'AC18',
           CHR == '22' ~ 'AC19',
           CHR == '23' ~ 'AC20',
           CHR == '24' ~ 'AC21',
           CHR == '25' ~ 'AC22',
           CHR == '26' ~ 'AC23',
           CHR == '27' ~ 'AC24',
           CHR == '28' ~ 'AC25',
           CHR == '29' ~ 'AC26',
           CHR == '30' ~ 'AC27',
           CHR == '31' ~ 'AC28',
           CHR == '32' ~ 'AC30',
           CHR == '33' ~ 'AC31',
           CHR == '34' ~ 'AC32',
           CHR == '35' ~ 'AC33',
           CHR == '36' ~ 'AC34',
           CHR == '37' ~ 'AC35',
           CHR == '38' ~ 'AC36',
           CHR == '39' ~ 'AC37',
           CHR == '0' ~ 'Contigs')))

# Split AC04q.1:29 --------------------------------------------------------

AC04q.1_29_split = function(data){
  AC04q.1_29 = data %>% 
    filter(AC_CHR == 'AC04q.1:29') %>% 
    mutate(AC_CHR = as.factor(case_when(
      win_mid > '40000000' ~ 'AC04q.1',
      win_mid < '40000000' ~ 'AC29')))
  
  data = data %>% 
    filter(AC_CHR != 'AC04q.1:29')
  
  data = bind_rows(data, 
                   AC04q.1_29)
}

data = AC04q.1_29_split(data = data)

# Conversion --------------------------------------------------------------

Mb_Conversion = function(data){
  data %>% 
    group_by(AC_CHR) %>% 
    mutate(win_mid_mb = win_mid/1000000)
}

data = Mb_Conversion(data)

#
# SNP density plot -------------------------------------------------------------

SNP_density_plot = data %>% 
  group_by(AC_CHR) %>% 
  mutate(window = row_number()) %>% 
  # filter(AC_CHR == 'AC29') %>%
  filter(AC_CHR != 'Contigs') %>%
  ggplot(aes(x = window,
             y = FST_n))+
  geom_step(stat = 'identity', 
            size = 1)+
  scale_y_continuous(breaks = seq(0, 5000, by = 1))+
  # scale_x_continuous(breaks = seq(0, 1000000, by = 1000))+
  geom_hline(yintercept = 3, 
             col = 'red',
             size = 3)+
  # geom_density(stat = 'identity')
  facet_grid(~AC_CHR, 
             scales = 'free')+
  labs(x = 'Window number per chromosome',
       y = 'Number of SNPS',
       title = 'Number of SNPS per window')+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.ticks = element_line(size = 1), 
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'),
        plot.title = element_text(size = 18,
                                  face = 'bold'))

setwd('~/Fst_sliding_window/')
ggsave('SNP_density_per_window.tiff',
       plot = SNP_density_plot, 
       width = 30,
       height = 10,
       limitsize = FALSE)


# Plot the results --------------------------------------------------------

## set the theme for all of the plots
theme_set(theme_bw())

## Drop windows with less than 3 SNPs

plot_fst = data %>% 
  filter(FST_n > 3) %>% 
  ggplot(aes(x = win_mid_mb,
           y = FST_mean, 
           group = AC_CHR)) +
  geom_point(col = 'grey49')+
  # geom_step(aes(x = win_mid,
  #               y = FST_mean),
  #           direction = 'mid',
  #           size = 2,
  #           col = '#1F2440')+
    geom_smooth(col = '#1F2440',
      # col = '#02E084',
                size = 2)+
  facet_grid(~ AC_CHR, 
             scales = 'free')+
  labs(title = 'VBRSIL Fst 200kb window', 
       x = 'Chromosomal position (Mb)', 
       y = 'Fst')+
  ylim(0.00, 1.00)+
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = 'bold',
                                  size = 10),
        strip.background = element_rect(fill = 'white',
                                        colour = 'black'), 
        plot.title = element_text(size = 12,
                                  face = 'bold'))

plot_fst
# final = plot_10kb/plot_50kb/plot_100kb/plot_150kb/plot_200kb/plot_250kb

ggsave('VBRSIL_Fst_200kb_window.tiff',
       plot = plot_fst, 
       width = 25,
       height = 10,
       limitsize = FALSE)
