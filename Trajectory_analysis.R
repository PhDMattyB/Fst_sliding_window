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
  filter(LaMorph.x != 'S.LGB2') %>% 
  rename(Lake = Lake.x, 
         Morph = Morph.x, 
         BP = BP.x, 
         LaMorph = LaMorph.x, 
         Sex = Sex.x)
# View(PWS_data_clean)

# PWS_data_clean = mutate(.data = PWS_data_clean,
#                  BP2 = as.factor(case_when(
#                    LaMorph == 'G.SB' ~ 'Benthic',
#                    LaMorph == 'G.PI' ~ 'Pelagic',
#                    LaMorph == 'S.LGB' ~ 'Benthic',
#                    LaMorph == 'S.PL' ~ 'Pelagic',
#                    LaMorph == 'S.PI' ~ 'Pelagic',
#                    LaMorph == 'T.LGB' ~ 'Benthic',
#                    LaMorph == 'T.SB' ~ 'Benthic2',
#                    LaMorph == 'T.PL' ~ 'Pelagic',
#                    LaMorph == 'V.BR' ~ 'Pelagic',
#                    LaMorph == 'V.SIL' ~ 'Benthic')))


# PWS_data_clean = mutate(.data = PWS_data_clean,
#                         BP2 = as.factor(case_when(
#                           LaMorph == 'G.SB' ~ 'G.Benthic',
#                           LaMorph == 'G.PI' ~ 'G.Pelagic',
#                           LaMorph == 'S.LGB' ~ 'S.Benthic',
#                           LaMorph == 'S.PL' ~ 'S.Pelagic',
#                           LaMorph == 'S.PI' ~ 'S.Pelagic',
#                           LaMorph == 'T.LGB' ~ 'T.Benthic',
#                           LaMorph == 'T.SB' ~ 'T.Benthic2',
#                           LaMorph == 'T.PL' ~ 'T.Pelagic',
#                           LaMorph == 'V.BR' ~ 'V.Pelagic',
#                           LaMorph == 'V.SIL' ~ 'V.Benthic')))

PWS_data_clean = mutate(.data = PWS_data_clean,
                        Vector = as.factor(case_when(
                          LaMorph == 'G.SB' ~ 'GBP',
                          LaMorph == 'G.PI' ~ 'GBP',
                          LaMorph == 'S.LGB' ~ 'SBP',
                          LaMorph == 'S.PL' ~ 'SBP',
                          LaMorph == 'S.PI' ~ 'SBP',
                          LaMorph == 'T.LGB' ~ 'TBP',
                          LaMorph == 'T.SB' ~ 'TBP2',
                          LaMorph == 'T.PL' ~ 'TBP',
                          LaMorph == 'T.PL2' ~ 'TBP2',
                          LaMorph == 'V.BR' ~ 'VBP',
                          LaMorph == 'V.SIL' ~ 'VBP')))
# RRPP --------------------------------------------------------------------


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
Vector = rrpp$Vector
Lake = rrpp$Lake

BP_test = lm.rrpp(Bodyshape_PWS ~ Lake*BP_Morph_Pair, 
                  data = rrpp, 
                  iter = 999)

summary(BP_test)
anova(BP_test)


# Trajectory Analysis -----------------------------------------------------

Traj_Analysis = trajectory.analysis(BP_test, 
                                    groups = BP_Morph_Pair, 
                                    traj.pts = Lake, 
                                    print.progress = F)

# Traj_Analysis = trajectory.analysis(BP_test, 
#                                     groups = Vector, 
#                                     traj.pts = Lake, 
#                                     print.progress = F)


summary(Traj_Analysis, 
        attribute = 'MD') ## Magnitude of differences is significant
summary(Traj_Analysis, 
        attribute = 'TC', 
        angle.type = 'deg') ## Differences in angles
summary(Traj_Analysis, 
        attribute = 'SD') ##Shape differences among trajectories




# Partial warps based trajectory analysis ---------------------------------

library(remotes)
install_github("fruciano/GeometricMorphometricsMix")
library(GeometricMorphometricsMix)
library(Morpho)

PWS_PCA = prcomp(Bodyshape_PWS)

qplot(PWS_PCA$x[,1], 
      PWS_PCA$x[,2], 
      col = rrpp$BP, 
      shape = rrpp$Vector, 
      size = 4) +
  coord_equal()


qplot(PCA_general$x[,1], PCA_general$x[,2], 
      col=allfishPW$pair, shape=allfishPW$temp, size=4)+
  coord_equal()+theme_classic()


Coords_SPLIT = split(as.data.frame(Bodyshape_PWS), 
                     list(rrpp$BP, 
                          rrpp$Vector), 
                     drop = T)
Coords = lapply(Coords_SPLIT, 
                colMeans)

