#####################Genetic diversity for Epidendrum fulgens#####################
####################################2022##########################################

##packages
install.packages("gdm")
install.packages("ggplot2")

library(ggplot2)
library(stats)

##setting working directory
setwd("")
##reading main table with population parameters and distances from centroids
pop = read.table("", sep=",", header = T)
pop


##separating variables

distance = pop$distance_centroid_km
eco_distance = pop$ecological_distance
allelic_richness = pop$riqueza_alelica
expected_heter = pop$heterozigozidade_esperada
fis = pop$fis
pop_specific_fst = pop$pop_specific_Fst

##transforming in numeric data
dist = as.numeric(distance)
allele = as.numeric(allelic_richness)
exp_het = as.numeric(expected_heter)

##transforming and leaving on the data frame

str(pop) #checking the levels of the variables
pop$distance_centroid_km = as.numeric(as.character(pop$distance_centroid_km))
pop$riqueza_alelica = as.numeric(as.character(pop$riqueza_alelica))
pop$heterozigozidade_esperada = as.numeric(as.character(pop$heterozigozidade_esperada))
str(pop)
par(mfrow = c(2,2))


###############################################################
#######################Linear Models###########################
###############################################################


##Using the log 10 of the geographical distance##

dist_geo_log = log10(dist)

##allelic richness##
reg.AR = lm(allele~eco_distance+dist_geo_log, data=pop[,1:7])
summary(reg.AR)
coef(lm(reg.AR))

##expected heterozygosity##
reg.EH = lm(exp_het~eco_distance+dist_geo_log, data=pop[,1:7])
summary(reg.EH)
coef(lm(reg.EH))

##Fis##
reg.Fis = lm(fis~eco_distance+dist_geo_log, data=pop_with_fis)
summary(reg.Fis)
coef(lm(reg.Fis))


##Pop specific Fst##
reg.Fst_popspecific = lm(pop$pop_specific_fst~eco_distance+dist_geo_log, data=pop[,1:7])
summary(reg.Fst_popspecific)
coef(lm(reg.Fst_popspecific))

#######################################################
####variance partitioning with the package relaimpo####
#######################################################

install.packages("relaimpo")
install.packages("RColorBrewer")
library(relaimpo)
library(RColorBrewer)

relaimpo::calc.relimp(reg.AR)

relaimpo::calc.relimp(reg.EH)

relaimpo::calc.relimp(reg.Fst)

relaimpo::calc.relimp(reg.Fis)

relaimpo::calc.relimp(reg.Fst_amova)

relaimpo::calc.relimp(reg.Fst_popspecific)


###############################
######plotting the models######
###############################

##allelic richness ## filtering by colour

ggplot(data=pop, mapping=aes(x=dist_geo_log, y=allele))+
  geom_point(aes(color=eco_distance), size=3.5)+
  stat_smooth(method = "lm", color="grey") +
  xlab("Log10 of Distance (km)") +
  ylab("Allelic Richness") +
  theme_bw() +
  theme(axis.title = element_text(size=10,face="bold"), axis.text = element_text(size=8))+
  scale_color_continuous(name="Environmental Distance", breaks = c(250, 750), 
                        labels = c("Low", "High"))


##expected heterozygosity ## filtering by colour

ggplot(data=pop, mapping=aes(x=dist_geo_log, y=exp_het))+
  geom_point(aes(color=eco_distance), size = 3.5)+
  stat_smooth(method = "lm", color="grey") +
  xlab("Log10 of Distance (km)") +
  ylab("Expected Heterozygosity") +
  theme_bw() +
  theme(axis.title = element_text(size=10,face="bold"), axis.text = element_text(size=8))+
  scale_color_continuous(name="Environmental Distance", breaks = c(250, 750), 
                        labels = c("Low", "High"))


##Fis ## filtering by colour

ggplot(data=pop, mapping=aes(x=dist_geo_log, y=pop$fis))+
  geom_point(aes(color=eco_distance), size= 3.5) +
  stat_smooth(method = "lm", color="grey") +
  xlab("Log10 of Distance (km)") +
  ylab("Fis") +
  theme_bw() +
  theme(axis.title = element_text(size=10,face="bold"), axis.text = element_text(size=8))+
  scale_color_continuous(name="Environmental Distance", breaks = c(250, 750), 
                         labels = c("Low", "High"))

##Pop specific Fst ## filtering by colour

ggplot(data=pop, mapping=aes(x=dist_geo_log, y=pop$pop_specific_fst))+
  geom_point(aes(color=eco_distance), size= 3.5) +
  stat_smooth(method = "lm", color="grey") +
  xlab("Log10 of Distance (km)") +
  ylab("Population Specific Fst") +
  theme_bw() +
  theme(axis.title = element_text(size=10,face="bold"), axis.text = element_text(size=8))+
  scale_color_continuous(name="Environmental Distance", breaks = c(250, 750), 
                         labels = c("Low", "High"))

