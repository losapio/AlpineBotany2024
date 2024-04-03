#packages ------------
install.packages('vegan')
install.packages('FD')
install.packages('glmmTMB', type = 'source')
install.packages('car')
install.packages('ggplot2')
install.packages('DHARMa')
install.packages('effects')
install.packages('effectsize')
install.packages('emmeans')
install.packages('parameters', type = 'source')
install.packages('readxl')
install.packages('AER')
install.packages('Matrix')
install.packages("tidyr")
install.packages('arm')
install.packages('brms')
install.packages('ggpubr')
install.packages('lmerTest')
install.packages('lme4')
install.packages('TMB', type = 'source')
install.packages('rmarkdown')

install.packages('viridis')
install.packages("r2d3")
install.packages("purrr")
install.packages("stringr")
install.packages("tibble")
install.packages("bipartiteD3")
install.packages('RColorBrewer')
install.packages('carData')

library(vegan)
library(FD)
library(glmmTMB)
library(car)
library(ggplot2)
library(bipartite)
library(DHARMa)
library(effects)
library(effectsize)
library(emmeans)
library(parameters)
library(readxl)
library(AER)
library(lme4)
library(lmerTest)
library(tidyr)
library(arm)
library(brms)
library(ggpubr)
library(dplyr)
library(lmerTest)
library(lme4)
library(TMB)
library(rmarkdown)

library(viridis)
library(r2d3)
library(purrr) 
library(stringr)
library(tibble)
library(bipartiteD3)
library(RColorBrewer)
library(carData)

load('ngan_231114.RData')

#######section 1: import data --------------------------
library(data.table)
getwd()
setwd("/Users/tungan/Desktop/NGAN/AlpineJournal")
flora_db <- as.data.frame(read.csv("Flora_corrected.csv", header = TRUE, sep = ";"))
quadratdb <- read.csv("quadrat.csv", header = TRUE, sep = ";")
transectdb <- read.csv("transect.csv", header = TRUE, sep = ";")
quadratdb
str(quadratdb)
colnames(quadratdb)
flora_db #to see the data on r
str(flora_db) #summary of dataset
colnames(flora_db) #to see only column

# create latin binomial for plant species of data to see if something wrong on name
##flora_db
sort(unique(flora_db$Genus)) #sort=put in alphabet to check easier
flora_db$plantsp = paste(flora_db$Genus,flora_db$Species,sep='.') #to create the name
flora_db$pl = paste(flora_db$Stage, flora_db$Plot, sep = '.')
sort(unique(flora_db$plantsp))
sort(unique(flora_db$pl))

flora_db$years = 0 
for(i in 1:nrow(flora_db)){
    if(flora_db$Stage[i]==1) flora_db$years[i] <- (2022-1989)/2
    if(flora_db$Stage[i]==2) flora_db$years[i] <- 2022-((1989+1925)/2)
    if(flora_db$Stage[i]==3) flora_db$years[i] <- 2022-((1925+1900)/2)
    if(flora_db$Stage[i]==4) flora_db$years[i] <- 2022-((1900+1864)/2)
}
flora_db
colnames(flora_db)  

##quadratdb
sort(unique(quadratdb$Genus_F))  ###eg for how to take a look at column 
quadratdb$plantsp = paste(quadratdb$Genus_F,quadratdb$Species_F,sep='.') 
quadratdb$insectfa = paste(quadratdb$Order_I, quadratdb$Family_I, sep = '.') 
quadratdb$insectsp = paste(quadratdb$Genus_I, quadratdb$Species_I, sep = '.')
quadratdb$pl = paste(quadratdb$Stage, quadratdb$Plot, sep = '.')
sort(unique(quadratdb$plantsp))
sort(unique(quadratdb$insectfa))
sort(unique(quadratdb$insectsp))
sort(unique(quadratdb$pl))

quadratdb$years = 0
for(i in 1:nrow(quadratdb)){
    if(quadratdb$Stage[i]==1) quadratdb$years[i] <- (2022-1989)/2
    if(quadratdb$Stage[i]==2) quadratdb$years[i] <- 2022-((1989+1925)/2)
    if(quadratdb$Stage[i]==3) quadratdb$years[i] <- 2022-((1925+1900)/2)
    if(quadratdb$Stage[i]==4) quadratdb$years[i] <- 2022-((1900+1864)/2)
}
quadratdb
colnames(quadratdb)

##transectdb
transectdb$plantsp = paste(transectdb$Genus_F,transectdb$Species_F,sep='.')
transectdb$insectfa = paste(transectdb$Order_I, transectdb$Family_I, sep = '.') 
transectdb$insectsp = paste(transectdb$Genus_I, transectdb$Species_I, sep = '.')
transectdb$pl = paste(transectdb$Stage, transectdb$Plot, sep = '.')
sort(unique(transectdb$plantsp))
sort(unique(transectdb$insectfa))
sort(unique(transectdb$insectsp))
sort(unique(transectdb$pl))

transectdb$years = 0
for(i in 1:nrow(transectdb)){
    if(transectdb$Stage[i]==1) transectdb$years[i] <- (2022-1989)/2
    if(transectdb$Stage[i]==2) transectdb$years[i] <- 2022-((1989+1925)/2)
    if(transectdb$Stage[i]==3) transectdb$years[i] <- 2022-((1925+1900)/2)
    if(transectdb$Stage[i]==4) transectdb$years[i] <- 2022-((1900+1864)/2)
}
transectdb
colnames(transectdb) 

#######section 2: create dataframe --------------------

# create summary table with summary data
# key variables of interest are:
# number of plant species occurring in the plot
# number of plant species in flower (with flowers)
# number of plant species visited by insects
# number of insect visits
# number of insect species
# number of insect genus
# number of insect families
# number of insect orders
# number of insect functional groups
# study design is composed by
# sampling method: quadrant or transect
# stage: 1, 2, 3, 4
# plot: A, B, C, D
# time: time of the day (from 9,30 to 17)
# study design for quadrats is = 4 stages * 4 plots * 3 replicates 
# study design for transect is = 4 stages * 4 plots * 3 replicates 
summary(quadratdb$Time)
summary(transectdb$Time)
# date:
unique(quadratdb$Date)
unique(transectdb$Date)
# Replicate
unique(quadratdb$Replicate)
unique(transectdb$Replicate)

## dataframe for flora_db -----------
colnames(flora_db)
myedgelist <- flora_db[,c(12,13)] #choose columns "plantsp" and "plot"
myedgelist$dummy = "dummy"
flora_matb = frame2webs(myedgelist, varnames = colnames(myedgelist),
                        type.out = "list", emptylist = TRUE)$dummy #combine data
flora_matb = t(flora_matb) #transpose
str(flora_matb)
## plant richness, alpha-diversity (shannon)
pool_df = data.frame(stage = gl(4, 4, 16, label = c('1','2','3','4')),
                     plot = gl(4, 1, 16, label =c('A','B',"C","D")),
                     years = gl(4,4,16, label =unique(flora_db$years)),
                     plant_richness = NA,
                     plant_diversity = NA)
pool_df$years = as.numeric(as.character(pool_df$years))
stages = sort(unique(flora_db$Stage))
plots = sort(unique(flora_db$Plot))

unique(flora_db$plantsp) #number of plant species

for(i in 1:4){ # loop for the stage
    print(i)
    for(j in 1:4){ # loop for the plot
        repstplant = subset(flora_db, Stage==stages[i] & Plot==plots[j])
        #calculate plant species richness as total number of species occurring (present) within a plot (and its buffer)
        pool_df$plant_richness[which(pool_df$stage==stages[i]&pool_df$plot==plots[j])] = 
            length(unique(repstplant$plantsp))
        #calculate plant diversity as number of plant species by shannon index
        pool_df$plant_diversity[which(pool_df$stage==stages[i]&pool_df$plot==plots[j])] = 
            diversity(table(repstplant[,12]), index = "shannon")
    }}
pool_df  #original data
write.csv(pool_df, 'pool_df.csv', row.names = FALSE)

## dataframe for quadratdb ------------
summaryt_quad = data.frame(stage = gl(4, 4*3, 4*4*3,labels = sort(unique(quadratdb$Stage))),
                           plot = gl(4, 3, 4*4*3,labels = sort(unique(quadratdb$Plot))),
                           replicate = gl(3, 1, 4*4*3,labels = sort(unique(quadratdb$Replicate))),
                           years = gl(4,4*3,4*3*4, label = sort(unique(quadratdb$years))),
                           n.plsp.flw = NA,  # number of plant species in flower (with flowers)
                           n.plsp.vis = NA,  # number of plant species visited by insects
                           n.inss.vis = NA, # number of individual insect visits
                           n.insp.vis = NA, # number of insect families
                           insect_abundance = NA,
                           insect_richness = NA,
                           insect_diversity = NA,
                           method = 'quad')
stages = sort(unique(quadratdb$Stage))
plots = sort(unique(quadratdb$Plot))
replicates = sort(unique(quadratdb$Replicate))

for(i in 1:4){ # loop for the stage
    print(i)
    for(j in 1:4) { # loop for the plot
        print(j)
        for(k in 1:3) { # loop for the replicates
            repstpl = subset(quadratdb, Stage==stages[i] & Plot==plots[j] & Replicate==replicates[k])
            # number of plant species visited by insects
            n.plsp.vis = length(unique(repstpl$plantsp))
            summaryt_quad$n.plsp.vis[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&
                                               summaryt_quad$replicate==replicates[k])] <- n.plsp.vis
            # number of individual insect visits
            n.inss.vis = nrow(repstpl)
            summaryt_quad$n.inss.vis[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&
                                               summaryt_quad$replicate==replicates[k])] <- n.inss.vis
            # number of insect species
            n.insp.vis = length(unique(repstpl$insectsp))
            summaryt_quad$n.insp.vis[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&
                                               summaryt_quad$replicate==replicates[k])] <- n.insp.vis
            #calculate insect species abundance as total number of insects occurring (present) within a plot (and its buffer)
            summaryt_quad$insect_abundance[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]
                                                 &summaryt_quad$replicate==replicates[k])] = nrow(repstpl)
            #calculate insect richness as total number of insect species occurring (present) within a plot (and its buffer)
            summaryt_quad$insect_richness[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]
                                                &summaryt_quad$replicate==replicates[k])] = length(unique(repstpl$insectsp))
            #calculate insect diversity as number of insect species by shannon index
            summaryt_quad$insect_diversity[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]
                                                 &summaryt_quad$replicate==replicates[k])] = diversity(table(repstpl[,19]), 
                                                                                                       index = "shannon")
        }
    }
} 

summaryt_quad$years = NA
for(i in 1:nrow(summaryt_quad)){
    if(summaryt_quad$stage[i]==1) summaryt_quad$years[i] <- (2022-1989)/2
    if(summaryt_quad$stage[i]==2) summaryt_quad$years[i] <- 2022-((1989+1925)/2)
    if(summaryt_quad$stage[i]==3) summaryt_quad$years[i] <- 2022-((1925+1900)/2)
    if(summaryt_quad$stage[i]==4) summaryt_quad$years[i] <- 2022-((1900+1864)/2)
}
summaryt_quad

##caculate number of plant flowering for quadrat by combine dataframe
head(summaryt_quad)
head(pool_df)
summaryt_quad = merge(summaryt_quad, pool_df, by.x=c('stage','plot'), by.y=c('stage','plot'))
head(summaryt_quad)
summaryt_quad$plant_richness
summaryt_quad$n.plsp.occ  = NULL
summaryt_quad$years = summaryt_quad$years.y
summaryt_quad$years.y  = NULL
summaryt_quad$years.x  = NULL
head(summaryt_quad)

###get n plant species flowering during sampling => obtain this info from nora's data blooming in flora_db
head(flora_db)
flora_db$Blooming_1[is.na(flora_db$Blooming_1)] = 0 #transfer NA to 0 on data
flora_db$Blooming_2[flora_db$Blooming_2 == 2] = 0
flora_db$Blooming_3[flora_db$Blooming_3 == 2] = 0
flora_db$Blooming_4[flora_db$Blooming_4 == 2] = 0
str(flora_db)

flora_db$Blooming_2 = as.numeric(flora_db$Blooming_2) #to confirm the values are numeric
flora_db$Blooming_3 = as.numeric(flora_db$Blooming_3)
flora_db$Blooming_4 = as.numeric(flora_db$Blooming_4)

stages
plots

for(i in 1:4){
    for(j in 1:4){
        # first blooming
        summaryt_quad$n.plsp.flw[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&summaryt_quad$replicate==1)] = 
            sum(flora_db$Blooming_1[which(flora_db$Stage==stages[i]&flora_db$Plot==plots[j])])
        # second blooming
        summaryt_quad$n.plsp.flw[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&summaryt_quad$replicate==2)] = 
            sum(flora_db$Blooming_2[which(flora_db$Stage==stages[i]&flora_db$Plot==plots[j])])
        # third blooming
        summaryt_quad$n.plsp.flw[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&summaryt_quad$replicate==3)] = 
            sum(flora_db$Blooming_3[which(flora_db$Stage==stages[i]&flora_db$Plot==plots[j])])
        # forth blooming
        summaryt_quad$n.plsp.flw[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&summaryt_quad$replicate==4)] = 
            sum(flora_db$Blooming_4[which(flora_db$Stage==stages[i]&flora_db$Plot==plots[j])])
    }
}

summaryt_quad ###should be the original dataframe
colnames(summaryt_quad)

## dataframe for transectdb --------------
summaryt_tran = data.frame(stage = gl(4, 4*3, 4*4*3,labels = sort(unique(transectdb$Stage))),
                           plot = gl(4, 3, 4*4*3,labels = sort(unique(transectdb$Plot))),
                           replicate = gl(3, 1, 4*4*3,labels = sort(unique(transectdb$Replicate))),
                           years = gl(4,4*3,4*3*4, label = sort(unique(transectdb$years))),
                           n.plsp.flw = NA,  # number of plant species in flowering
                           n.plsp.vis = NA,  # number of plant species visited by insects
                           n.inss.vis = NA, # number of insect visits
                           n.insp.vis = NA, # number of insect species visit
                           insect_abundance = NA,
                           insect_richness = NA,
                           insect_diversity = NA,
                           method = 'tran')
summaryt_tran

stages = sort(unique(transectdb$Stage))
plots = sort(unique(transectdb$Plot))
replicates = sort(unique(transectdb$Replicate))

for(i in 1:4){ # loop for the stage
    print(i)
    for(j in 1:4){ # loop for the plot
        print(j)
        for(k in 1:3){ # loop for the replicates
            repstpl1 = subset(transectdb, Stage==stages[i] & Plot==plots[j] & Replicate==replicates[k])
            # number of plant species visited by insects
            n.plsp.vis = length(unique(repstpl1$plantsp))
            summaryt_tran$n.plsp.vis[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]&
                                               summaryt_tran$replicate==replicates[k])] <- n.plsp.vis
            # number of individual insect visits
            n.inss.vis = nrow(repstpl1)
            summaryt_tran$n.inss.vis[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]&
                                               summaryt_tran$replicate==replicates[k])] <- n.inss.vis
            # number of insect species
            n.insp.vis = length(unique(repstpl1$insectsp))
            summaryt_tran$n.insp.vis[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]&
                                               summaryt_tran$replicate==replicates[k])] <- n.insp.vis
            #calculate insect species abundance as total number of insects occurring (present) within a plot (and its buffer)
            summaryt_tran$insect_abundance[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]
                                                 &summaryt_tran$replicate==replicates[k])] = nrow(repstpl1)
            #calculate insect species richness as total number of insect morphospecies occurring (present) within a plot (and its buffer)
            summaryt_tran$insect_richness[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]
                                                &summaryt_tran$replicate==replicates[k])] = length(unique(repstpl1$insectsp))
            #calculate insect diversity as number of insect species by shannon index
            summaryt_tran$insect_diversity[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]
                                                 &summaryt_tran$replicate==replicates[k])] = diversity(table(repstpl1[,19]), 
                                                                                                       index = "shannon") 
        }
    }
} 

summaryt_tran$years = NA
for(i in 1:nrow(summaryt_tran)){
    if(summaryt_tran$stage[i]==1) summaryt_tran$years[i] <- (2022-1989)/2
    if(summaryt_tran$stage[i]==2) summaryt_tran$years[i] <- 2022-((1989+1925)/2)
    if(summaryt_tran$stage[i]==3) summaryt_tran$years[i] <- 2022-((1925+1900)/2)
    if(summaryt_tran$stage[i]==4) summaryt_tran$years[i] <- 2022-((1900+1864)/2)
}
summaryt_tran
head(summaryt_tran)
head(pool_df)
summaryt_tran = merge(summaryt_tran, pool_df, by.x=c('stage','plot'), by.y=c('stage','plot'))
summaryt_tran
summaryt_tran$years = summaryt_tran$years.y
summaryt_tran$years.y  = NULL
summaryt_tran$years.x  = NULL
colnames(summaryt_tran)
summaryt_tran ##orginal dataframe
#summaryt_tran = summaryt_tran[,c(1:3,5:15,4)] #relocate the column years


## dataframe for both quadrat and transect -----------
colnames(summaryt_quad) == colnames(summaryt_tran)
summaryt = rbind(summaryt_quad,summaryt_tran)
summaryt

write.csv(summaryt, 'summaryt.csv', row.names = FALSE)

#######section 3: statistics+modeling: diversity of flora -------------
#########plant richness + plant diversity (pool_df) ---------------
### model for plant species richness (try with poisson and nbinom1 and nbinom2 but result is similar)
hist(pool_df$plant_richness, 10, col = 'green', las = 1) ### look at distribution (bieu do cot)
mod.plrich_2 = glmmTMB(plant_richness ~ poly(years,2),
                       data = pool_df,
                       family = 'genpois')
mod.plrich_1 = glmmTMB(plant_richness ~ years,
                       data = pool_df,
                       family = 'genpois')
mod.plrich_0 = glmmTMB(plant_richness ~ 1,
                       data = pool_df,
                       family = 'genpois')

summary(mod.plrich_0) # check dispersion
summary(mod.plrich_1) 
summary(mod.plrich_2) 
testDispersion(mod.plrich_0)
testDispersion(mod.plrich_1)
testDispersion(mod.plrich_2)

#log(mean(pool_df$plant_richness)) to compare with Edviance

anova(mod.plrich_0,mod.plrich_1,mod.plrich_2)
Anova(mod.plrich_2)
cohens_f(aov(mod.plrich_2))

##creat table data of Anova and so on
write.csv(round(Anova(mod.plrich_2),3), 'mod.plrich_2.csv')

pdf("1.mod.plrichness.pdf", 5, 5)
ggplot(pool_df, aes(x = years,
                    y = plant_richness)) +
    geom_jitter(width = 2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Plant diversity [n. species]')
dev.off()

###model for plant diversity (shannon index)
hist(pool_df$plant_diversity, 10, col = 'green', las = 1)
mod.pldiversity_2 = lm(plant_diversity ~ poly(years,2),
                             data = pool_df)
mod.pldiversity_1 = lm(plant_diversity ~ years,
                             data = pool_df)
mod.pldiversity_0 = lm(plant_diversity ~ 1,
                             data = pool_df)
summary(mod.pldiversity_0)
summary(mod.pldiversity_1)
summary(mod.pldiversity_2)
testDispersion(mod.pldiversity_0)
testDispersion(mod.pldiversity_1)
testDispersion(mod.pldiversity_2)

anova(mod.pldiversity_0,mod.pldiversity_1,mod.pldiversity_2)
Anova(mod.pldiversity_2)
cohens_f(aov(mod.pldiversity_2))
##creat table data of Anova and so on
write.csv(round(Anova(mod.pldiversity_2),3), 'mod.pldiversity_2.csv')

pdf("2.mod.pldiversity.pdf", 5, 5)
ggplot(pool_df, aes(x = years,
                    y = plant_diversity)) +
    geom_jitter(width = 2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Plant diversity [Shannon index]')
dev.off()

#########plant visited by insect ------------
summaryt$n.plsp.vis
##modelling number of plants visited by insects in summaryt(quadrant+transect)
hist(summaryt$n.plsp.vis, 20, col = 'green', las = 1)
####add replicate to dataframe of summaryt
summaryt$replicate2 = as.numeric(summaryt$replicate)
summaryt$replicate2[which(summaryt$method=='tran')] = summaryt$replicate2[which(summaryt$method=='tran')]+3
summaryt$replicate2 = as.factor(summaryt$replicate2)
#### family poisson is better than negative binomial:(Time+plrich=fixed effects)+(method,rep,plot=random effect) => best results => choose for thesis
mod.plvis.tot_2 = glmer(summaryt$n.plsp.vis ~ poly(years,2)*plant_richness+
                            (1|replicate2)+(1|method)+(1|plot), #add replicate as random
                        data = summaryt,
                        family = "poisson") #quadratic model+plric and year
mod.plvis.tot_2Y = glmer(summaryt$n.plsp.vis ~ poly(years,2)+
                             (1|replicate2)+(1|method)+(1|plot),
                        data = summaryt,
                        family = "poisson") #quadratic model + effect on year
mod.plvis.tot_1 = glmer(summaryt$n.plsp.vis ~ years*plant_richness+(1|replicate2)+(1|method)+(1|plot),
                       data = summaryt,
                       family = "poisson") #linear model + effect on year&plrich
mod.plvis.tot_1R = glmer(summaryt$n.plsp.vis ~ plant_richness+(1|replicate2)+(1|method)+(1|plot),
                         data = summaryt,
                         family = "poisson") #linear model + effect on plrich
#mod.plvis.tot_1Y = glmer(summaryt$n.plsp.vis ~ years+(1|replicate2)+(1|method)+(1|plot),
#                        data = summaryt,
#                        family = "poisson") #linear model + effect on year
mod.plvis.tot_00 = glmer(summaryt$n.plsp.vis ~ 1+(1|replicate2)+(1|method)+(1|plot),
                      data = summaryt,
                      family = "poisson") #null model

summary(mod.plvis.tot_00) # check dispersion
summary(mod.plvis.tot_1) 
summary(mod.plvis.tot_1R)
summary(mod.plvis.tot_2Y) 
summary(mod.plvis.tot_2) 
testDispersion(mod.plvis.tot_00, plot = T)
testDispersion(mod.plvis.tot_1, plot = T)
testDispersion(mod.plvis.tot_1R, plot = T)
testDispersion(mod.plvis.tot_2Y, plot = T)
testDispersion(mod.plvis.tot_2, plot = T) #quadratic is better than linear=>chose

# including year first and then plant richness
anova(mod.plvis.tot_00,mod.plvis.tot_2Y,mod.plvis.tot_2)
# including plant richness first and then years
anova(mod.plvis.tot_00,mod.plvis.tot_1R,mod.plvis.tot_2)

Anova(mod.plvis.tot_2)

pdf("3.mod.plvis.tot.pdf", 5, 5)
ggplot(summaryt, aes(x = years,
                     y = n.plsp.vis,
                     col = method)) +
    geom_jitter(width = 5)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Visited plant [n. species]')
dev.off()

pdf("4.mod.plvis.rich.pdf", 5, 5)  #model with plant diversity
ggplot(summaryt, aes(x = plant_richness,
                     y = n.plsp.vis,
                     col = stage)) +
    geom_jitter(width = 5)+
    geom_smooth(method = lm, formula = y ~ x) +
    theme_classic()+
    labs(x = 'Plant diversity [n. species]',
         y = 'Visited plant [n. species]')
dev.off()

#########plant flowering ------------
#model for plant flowering in quadrat
summaryt_quad$n.plsp.flw 
hist(summaryt_quad$n.plsp.flw, 10, col = 'green', las = 1)
mod.plflw.quad_2 = glmmTMB(n.plsp.flw ~ poly(years,2)+(1|replicate)+(1|plot), # quadratic model
                       data = summaryt_quad,
                       family = "nbinom2")
mod.plflw.quad_1 = glmmTMB(n.plsp.flw ~ years+(1|replicate)+(1|plot), # linear model
                       data = summaryt_quad,
                       family = "nbinom2")
mod.plflw.quad_0 = glmmTMB(n.plsp.flw ~ 1+(1|replicate)+(1|plot),  # null model
                       data = summaryt_quad,
                       family = "nbinom2")
summary(mod.plflw.quad_0) # check dispersion
summary(mod.plflw.quad_1)
summary(mod.plflw.quad_2)  # model coefficients (parameter estimates)
testDispersion(mod.plflw.quad_0, plot = T)
testDispersion(mod.plflw.quad_1, plot = T)
testDispersion(mod.plflw.quad_2, plot = T)

anova(mod.plflw.quad_0,mod.plflw.quad_1,mod.plflw.quad_2) # explained variance
Anova(mod.plflw.quad_2)
cohens_f(aov(mod.plflw.quad_2))

pdf("5.mod.plflw.quad.pdf", 5, 5)
ggplot(summaryt_quad, aes(x = years,
                          y = n.plsp.flw)) +
    geom_jitter(width = 2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Plant flowering [n. species]')
dev.off()

#build combined model
###model for ratio pl.vis/pl.flw
summaryt_quad$std.plvis.plflw = summaryt_quad$n.plsp.vis/(summaryt_quad$n.plsp.flw)
summaryt_quad$std.plvis.plflw[is.nan(summaryt_quad$std.plvis.plflw)] = 0 #transfer NaN to 0 on data
hist(summaryt_quad$std.plvis.plflw , 10, col = 'green', las = 1)

mod.stdplvis.plflw.2 = lmer(std.plvis.plflw ~ poly(years,2)+(1|replicate)+(1|plot),
                       data = summaryt_quad) 
mod.stdplvis.plflw.1 = lmer(std.plvis.plflw ~ years+(1|replicate)+(1|plot),
                       data = summaryt_quad) 
mod.stdplvis.plflw.0 = lmer(std.plvis.plflw ~ 1+(1|replicate),
                       data = summaryt_quad) #boundary singular????

summary(mod.stdplvis.plflw.0) # check dispersion
summary(mod.stdplvis.plflw.1)
summary(mod.stdplvis.plflw.2) 
testDispersion(mod.stdplvis.plflw.0, plot = T)
testDispersion(mod.stdplvis.plflw.1, plot = T)
testDispersion(mod.stdplvis.plflw.2, plot = T)

anova(mod.stdplvis.plflw.0,mod.stdplvis.plflw.1,mod.stdplvis.plflw.2)
Anova(mod.stdplvis.plflw.2)

ranef(mod.stdplvis.plflw.2) #to check intercept of replicate
ranova(mod.stdplvis.plflw.2) #to test the affect of random effect

pdf("6.mod.stdplvis.plflw.quad.pdf", 5, 6) 
ggplot(summaryt_quad, aes(x = years,
                          y = std.plvis.plflw )) +
    geom_point()+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Frequency visited plant [n. sp/n. sp]')
dev.off()

#######section 4: statistic: diversity of insect -------------
#########insect richness -------------
# model for insect species richness in total (summaryt)
summaryt$replicate2 = as.numeric(summaryt$replicate)
summaryt$replicate2[which(summaryt$method=='tran')] = summaryt$replicate2[which(summaryt$method=='tran')]+3
summaryt$replicate2 = as.factor(summaryt$replicate2)
mod.inrich.sp.tot_2 = glmmTMB(insect_richness ~ poly(years,2)*plant_richness+(1|plot/replicate2)+(1|method), #with plot as randomeffect, model convergence problem
                              data = summaryt,
                              family = 'nbinom2')
mod.inrich.sp.tot_2Y = glmmTMB(insect_richness ~ poly(years,2)+(1|plot/replicate2)+(1|method),
                              data = summaryt,
                              family = 'nbinom2')
mod.inrich.sp.tot_2R = glmmTMB(insect_richness ~ plant_richness+(1|plot/replicate2)+(1|method),
                              data = summaryt,
                              family = 'nbinom2')
mod.inrich.sp.tot_1 = glmmTMB(insect_richness ~ years*plant_richness+(1|plot/replicate2)+(1|method),
                              data = summaryt,
                              family = 'nbinom2')
mod.inrich.sp.tot_1Y = glmmTMB(insect_richness ~ years+(1|plot/replicate2)+(1|method),
                              data = summaryt,
                              family = 'nbinom2')
mod.inrich.sp.tot_0 = glmmTMB(insect_richness ~ 1+(1|plot/replicate2)+(1|method),
                               data = summaryt,
                               family = 'nbinom2')

summary(mod.inrich.sp.tot_0)
summary(mod.inrich.sp.tot_1)
summary(mod.inrich.sp.tot_1Y)
summary(mod.inrich.sp.tot_2R)
summary(mod.inrich.sp.tot_2Y)
summary(mod.inrich.sp.tot_2)
testDispersion(mod.inrich.sp.tot_0)
testDispersion(mod.inrich.sp.tot_1)
testDispersion(mod.inrich.sp.tot_1Y)
testDispersion(mod.inrich.sp.tot_2R)
testDispersion(mod.inrich.sp.tot_2Y)
testDispersion(mod.inrich.sp.tot_2)

plot(Effect("years", mod.inrich.sp.tot_2, residuals = T))

# including year first and then plant richness
anova(mod.inrich.sp.tot_0,mod.inrich.sp.tot_2Y,mod.inrich.sp.tot_2)
# including plant richness first and then years
anova(mod.inrich.sp.tot_0,mod.inrich.sp.tot_2R,mod.inrich.sp.tot_2)

Anova(mod.inrich.sp.tot_2)
cohens_f(aov(mod.inrich.sp.tot_2))

pdf("8.mod.inrich.sp.tot.pdf", 5, 5)
ggplot(summaryt, aes(x = years,
                     y = insect_richness, 
                     color = method)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Pollinator diversity [n. species]')
dev.off()

######### insect abundance ------------
summaryt$insect_abundance
# look at distribution
hist(summaryt$insect_abundance, 20, col = 'green', las = 1)
# model for insect species abundance (family) in total, plrich as predictor
mod.inabundance.tot_2 = glmmTMB(insect_abundance ~ poly(years,2)*plant_richness+
                                    (1|replicate2)+(1|method)+(1|plot),
                                data = summaryt,
                                family = 'nbinom2')
mod.inabundance.tot_2Y = glmmTMB(insect_abundance ~ poly(years,2)+(1|replicate2)+(1|method)+(1|plot),
                                data = summaryt,
                                family = 'nbinom2')
mod.inabundance.tot_2R = glmmTMB(insect_abundance ~ plant_richness+(1|replicate2)+(1|method)+(1|plot),
                                data = summaryt,
                                family = 'nbinom2')
mod.inabundance.tot_1 = glmmTMB(insect_abundance ~ years*plant_richness+(1|replicate2)+(1|method)+(1|plot),
                                data = summaryt,
                                family = 'nbinom2')
mod.inabundance.tot_1Y = glmmTMB(insect_abundance ~ years+(1|replicate2)+(1|method)+(1|plot),
                                data = summaryt,
                                family = 'nbinom2')
mod.inabundance.tot_0 = glmmTMB(insect_abundance ~ 1+(1|replicate2)+(1|method)+(1|plot),
                                 data = summaryt,
                                 family = 'nbinom2')

summary(mod.inabundance.tot_0)
summary(mod.inabundance.tot_1Y)
summary(mod.inabundance.tot_1)
summary(mod.inabundance.tot_2R)
summary(mod.inabundance.tot_2Y)
summary(mod.inabundance.tot_2)
testDispersion(mod.inabundance.tot_0)
testDispersion(mod.inabundance.tot_1Y)
testDispersion(mod.inabundance.tot_1)
testDispersion(mod.inabundance.tot_2R)
testDispersion(mod.inabundance.tot_2Y)
testDispersion(mod.inabundance.tot_2)

# including year first and then plant richness
anova(mod.inabundance.tot_0,mod.inabundance.tot_2Y,mod.inabundance.tot_2)
# including plant richness first and then years
anova(mod.inabundance.tot_0,mod.inabundance.tot_2R,mod.inabundance.tot_2)

Anova(mod.inabundance.tot_2)
cohens_f(aov(mod.inabundance.tot_2))

pdf("9.mod.inabundance.tot.pdf", 5, 5)
ggplot(summaryt, aes(x = years,
                     y = insect_abundance, 
                     color = method)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic()+
    labs(x = 'Glacier retreat [years]',
         y = 'Pollinator abundance [n. individual]')
dev.off()

pdf("10.mod.inabund.plrich.pdf", 5, 5) ##add n of plant richness as predictor in model 
ggplot(summaryt, aes(x = plant_richness,
                     y = insect_abundance, 
                     color = stage)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ x) +
    theme_classic() +
    labs(x = 'Plant diversity [n. species]',
         y = 'Insect abundance [n. individual]')
dev.off()

pdf("10.1.mod.indiv.pldiv.pdf", 5, 5) ##add n of plant richness as predictor in model 
ggplot(summaryt, aes(x = plant_richness,
                     y = insect_richness, 
                     color = stage)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ x) +
    theme_classic() +
    labs(x = 'Plant diversity [n. species]',
         y = 'Pollinator diversity [n. species]')
dev.off()

#########insect diversity ------------
colnames(summaryt)
summaryt$insect_diversity = as.numeric(summaryt$insect_diversity)
# look at distribution
hist(summaryt$insect_diversity, 20, col = 'green', las = 1)
# model for insect species diversity in total 
mod.indiversity.tot_2 = lmer(insect_diversity ~ poly(years,2)*plant_richness+
                                 (1|replicate2)+(1|method)+(1|plot),
                                data = summaryt) #problem with fitted mixed model is near singular
mod.indiversity.tot_2Y = lmer(insect_diversity ~ poly(years,2)+(1|replicate2)+(1|method)+(1|plot),
                             data = summaryt)
mod.indiversity.tot_2R = lmer(insect_diversity ~ plant_richness+(1|replicate2)+(1|method)+(1|plot),
                             data = summaryt)
#mod.indiversity.tot_1 = lmer(insect_diversity ~ years*plant_richness+(1|replicate2)(1|method)+(1|plot),
#                                data = summaryt)
mod.indiversity.tot_1Y = lmer(insect_diversity ~ years+(1|replicate2)+(1|method)+(1|plot),
                             data = summaryt)
mod.indiversity.tot_0 = lmer(insect_diversity ~ 1+(1|replicate2)+(1|method)+(1|plot),
                              data = summaryt) 

summary(mod.indiversity.tot_0)
summary(mod.indiversity.tot_1Y)
summary(mod.indiversity.tot_2Y)
summary(mod.indiversity.tot_2R)
summary(mod.indiversity.tot_2)
testDispersion(mod.indiversity.tot_0)
testDispersion(mod.indiversity.tot_1Y)
testDispersion(mod.indiversity.tot_2R)
testDispersion(mod.indiversity.tot_2Y)
testDispersion(mod.indiversity.tot_2)

# including year first and then plant richness
anova(mod.indiversity.tot_0,mod.indiversity.tot_2Y,mod.indiversity.tot_2)
# including plant richness first and then years
anova(mod.indiversity.tot_0,mod.indiversity.tot_2R,mod.indiversity.tot_2)

Anova(mod.indiversity.tot_2)

pdf("15.mod.indiversity.tot.pdf", 5, 5)
ggplot(summaryt, aes(x = years,
                     y = insect_diversity, 
                     color = method)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic() +
    labs(x = 'Glacier retreat [years]',
         y = 'Pollinator diversity [Shannon index]')
dev.off()

######## plant-insect network ----------
#### calculate network diversity for each stage, each plot, and each replicate
# for quadrant
summaryt_quad$replicate = as.numeric(as.character(summaryt_quad$replicate))
replicates = c(1,2,3)

summaryt_quad$intdiv = NA # interaction diversity

for(i in 1:4){ # loop for the stage
    print(i)
    for(j in 1:4){ # loop for the plot
        print(j)
        for(k in 1:3){ # loop for the replicates
            mel = subset(quadratdb, Stage==i & Plot==plots[j] & Replicate==k)
            if(nrow(mel)>0){
                mel <- mel[,c(17,19)] # columns "plantsp","insectsp"
                mel$spint = paste(mel$plantsp,mel$insectsp,sep="-")
                summaryt_quad$intdiv[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&summaryt_quad$replicate==replicates[k])] = 
                    diversity(table(mel$spint), index = "shannon") #careful of replicate
                ### network complexity = connectance
            }
            else summaryt_quad$intdiv[which(summaryt_quad$stage==stages[i]&summaryt_quad$plot==plots[j]&summaryt_quad$replicate==replicates[k])] = 0
        }
    }
} 

summaryt_quad$intdiv
summaryt_quad$insect_richness
#connectance = n.links/n.sp^2 = n.links/n.sp of plants and insects = n.visits/(n. plant sp * n. insect sp)
head(summaryt_quad)
summaryt_quad$connectance = summaryt_quad$insect_abundance/(summaryt_quad$insect_richness*summaryt_quad$plant_richness)
summaryt_quad$connectance[is.nan(summaryt_quad$connectance)] = 0 #transfer NaN to 0 on data

# for transect
summaryt_tran$intdiv = NA # interaction diversity

for(i in 1:4){ # loop for the stage
    print(i)
    for(j in 1:4){ # loop for the plot
        print(j)
        for(k in 1:3){ # loop for the replicates
            meltran = subset(transectdb, Stage==i & Plot==plots[j] & Replicate==k)
            if(nrow(meltran)>0){
                meltran <- meltran[,c(17,19)] # columns "plantsp" and "insectsp"
                meltran$spint = paste(meltran$plantsp,meltran$insectsp,sep="-")
                summaryt_tran$intdiv[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]&summaryt_tran$replicate==replicates[k])] = 
                    diversity(table(meltran$spint), index = "shannon")
                ### network complexity = connectance
            }
            else summaryt_tran$intdiv[which(summaryt_tran$stage==stages[i]&summaryt_tran$plot==plots[j]&summaryt_tran$replicate==replicates[k])] = 0
        }
    }
} 

summaryt_tran$intdiv
#connectance = n.links/n.sp^2 = n.links/n.sp of plants and insects = n.visits/(n. plant sp * n. insect sp)
summaryt_tran$connectance = summaryt_tran$insect_abundance/(summaryt_tran$insect_richness*summaryt_tran$plant_richness)
summaryt_tran$connectance[is.nan(summaryt_tran$connectance)] = 0 #transfer NaN to 0 on data
summaryt_tran

# model for combination quadrant and transect
summaryt_quad_int = summaryt_quad
summaryt_tran_int = summaryt_tran
summaryt_quad_int$std.plvis.plflw = NULL
summaryt_quad_int$std.inab.flw = NULL
summaryt_quad_int$std.inrich.flw = NULL
summaryt_quad_int$std.inabu.plrich = NULL
summaryt_quad_int$std.inrich.plrich = NULL
colnames(summaryt_quad_int) == colnames(summaryt_tran_int)
summaryt_int = rbind(summaryt_quad_int,summaryt_tran_int)
colnames(summaryt_int)
write.csv(summaryt_int, 'summaryt_int.csv', row.names = FALSE)

#model for interaction diversity in total and add plant richness as predictor
summaryt_int$intdiv = as.numeric(summaryt_int$intdiv)
hist(summaryt_int$intdiv, 20, col = 'green', las = 1)

summaryt_int$replicate2 = as.numeric(summaryt_int$replicate)
summaryt_int$replicate2[which(summaryt_int$method=='tran')] = summaryt_int$replicate2[which(summaryt_int$method=='tran')]+3
summaryt_int$replicate2 = as.factor(summaryt_int$replicate2)

mod.intdiv.tot_2 = lmer(intdiv ~ poly(years,2)*plant_richness +(1|replicate2)+(1|method)+(1|plot),
                             data = summaryt_int) ##because of shannon index
mod.intdiv.tot_2Y = lmer(intdiv ~ poly(years,2) + (1|replicate2)+(1|method)+(1|plot),
                             data = summaryt_int)
#mod.intdiv.tot_1 = lmer(intdiv ~ years*plant_richness +(1|replicate2)+(1|method)+(1|plot),
#                         data = summaryt_int)
mod.intdiv.tot_1R = lmer(intdiv ~ plant_richness +(1|replicate2)+(1|method)+(1|plot),
                        data = summaryt_int)
mod.intdiv.tot_1Y = lmer(intdiv ~ years +(1|replicate2)+(1|method)+(1|plot),
                        data = summaryt_int)
mod.intdiv.tot_0 = lmer(intdiv ~ 1+(1|replicate2)+(1|method)+(1|plot),
                              data = summaryt_int) 

summary(mod.intdiv.tot_0)
summary(mod.intdiv.tot_1Y)
summary(mod.intdiv.tot_1R)
summary(mod.intdiv.tot_2Y)
summary(mod.intdiv.tot_2)
testDispersion(mod.intdiv.tot_0)
testDispersion(mod.intdiv.tot_1Y)
testDispersion(mod.intdiv.tot_1R)
testDispersion(mod.intdiv.tot_2Y)
testDispersion(mod.intdiv.tot_2)

# including year first and then plant richness
anova(mod.intdiv.tot_0,mod.intdiv.tot_2Y,mod.intdiv.tot_2)
# including plant richness first and then years
anova(mod.intdiv.tot_0,mod.intdiv.tot_1R,mod.intdiv.tot_2)

Anova(mod.intdiv.tot_2)

pdf("20.mod.intdiv.tot.pdf", 5, 5)
ggplot(summaryt_int, aes(x = years,
                     y = intdiv, 
                     color = method)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic() +
    labs(x = 'Glacier retreat [years]',
         y = 'Interaction diversity [Shannon index]')
dev.off()

pdf("21.mod.intdiv.rich.pdf", 5, 5)  #model with plant richenss
ggplot(summaryt_int, aes(x = plant_richness,
                         y = intdiv,
                         col = stage)) +
    geom_jitter(width = 5)+
    geom_smooth(method = lm, formula = y ~ x) +
    theme_classic()+
    labs(x = 'Plant diversity [n. species]',
         y = 'Interaction diversity [Shannon index]')
dev.off()

#model for interaction connectance in total (dont need to add plant richness as predictor because already have this on ratio)
summaryt_int$connectance = as.numeric(summaryt_int$connectance)
hist(summaryt_int$connectance, 20, col = 'green', las = 1)

summaryt_int$replicate2 = as.numeric(summaryt_int$replicate)
summaryt_int$replicate2[which(summaryt_int$method=='tran')] = summaryt_int$replicate2[which(summaryt_int$method=='tran')]+3
summaryt_int$replicate2 = as.factor(summaryt_int$replicate2)

mod.connectance.tot_2 = lmer(connectance ~ poly(years,2) + (1|replicate)+(1|method)+(1|plot),
                              data = summaryt_int) 
mod.connectance.tot_1 = lmer(connectance ~ years + (1|replicate)+(1|method)+(1|plot),
                              data = summaryt_int)
mod.connectance.tot_0 = lmer(connectance ~ 1 + (1|replicate)+(1|method)+(1|plot),
                               data = summaryt_int)

summary(mod.connectance.tot_0)
summary(mod.connectance.tot_1)
summary(mod.connectance.tot_2)
testDispersion(mod.connectance.tot_0)
testDispersion(mod.connectance.tot_1)
testDispersion(mod.connectance.tot_2)

anova(mod.connectance.tot_0,mod.connectance.tot_1,mod.connectance.tot_2)
Anova(mod.connectance.tot_2)

pdf("22.mod.connectance.tot.pdf", 5, 5)
ggplot(summaryt_int, aes(x = years,
                         y = connectance, 
                         color = method)) +
    geom_jitter(width = 0.2)+
    geom_smooth(method = lm, formula = y ~ poly(x, 2)) +
    theme_classic() +
    labs(x = 'Glacier retreat [years]',
         y = 'Network complexity [Connectance]')
dev.off()


#model of interation and biodiversity (n. plant sp * n. insect sp)
colnames(summaryt_int)
summaryt_int$biodiv = summaryt_int$plant_richness*summaryt_int$insect_richness
summaryt_int$biodiv
pdf("20.1.mod.intdiv.biodiv.pdf", 5, 5)  
ggplot(summaryt_int, aes(x = biodiv,
                         y = intdiv,
                         col = stage)) +
    geom_jitter(width = 5)+
    geom_smooth(method = lm, formula = y ~ x) +
    theme_classic()+
    labs(x = 'Number of species [n.pl sp *n.in sp]',
         y = 'Interaction diversity [Shannon index]')
dev.off()

pdf("20.2.mod.interation.biodiv.pdf", 5, 5)  
ggplot(summaryt_int, aes(x = connectance,
                         y = intdiv,
                         col = stage)) +
    geom_jitter(width = 5)+
    geom_smooth(method = lm, formula = y ~ x) +
    theme_classic()+
    labs(x = 'Number of species [n.pl sp *n.in sp]',
         y = 'Interaction diversity [n. interaction]')
dev.off()


######### build network for insect species ----------
#### network for insect species total
colnames(quadratdb) ==  colnames(transectdb)
web_tot = rbind(quadratdb,transectdb)
colnames(web_tot)
#web_tot$Genus_F = substr(web_tot$Genus_F,1,1)
web_tot$plantsp = paste(web_tot$Genus_F,web_tot$Species_F,sep=' ')
web_tot$insectsp = paste(web_tot$Genus_I,web_tot$Species_I,sep=' ')
#### network for insect species total
myedgelist7 <- web_tot[,c(17,19)] # columns "plantsp"=lower and "insectsp"=higher
myedgelist7$dummy = web_tot$Stage
flora_matb3 = frame2webs(myedgelist7, varnames = colnames(myedgelist7),
                         type.out = "list", emptylist = TRUE)
flora_matb3

pdf("plotweb.stage1.sp.pdf", 9, 5)
plotweb(flora_matb3$'1',text.rot=90, labsize = 0.5, y.lim=c(-1,3), arrow="down.center", col.high = 'lightblue', 
        col.interaction=t(ifelse(flora_matb2$'1'[,] <= 1, adjustcolor('lightblue', alpha.f = 1.5), 
                                 adjustcolor('black', alpha.f = 1))), #interaction<=1:blue
        bor.col.interaction = NA, bor.col.high = NA)
dev.off()

pdf("plotweb.stage2.sp.pdf", 9, 5)
plotweb(flora_matb3$'2',text.rot=90, labsize = 0.5, y.lim=c(-1,3), arrow="down.center", col.high = 'lightblue', 
        col.interaction=t(ifelse(flora_matb2$'2'[,] <= 1, adjustcolor('lightblue', alpha.f = 1.5), 
                                 adjustcolor('black', alpha.f = 1))), 
        bor.col.interaction = NA, bor.col.high = NA)
dev.off()

pdf("plotweb.stage3.sp.pdf", 9, 5)
plotweb(flora_matb3$'3',text.rot=90, labsize = 0.5, y.lim=c(-1,3), arrow="down.center", col.high = 'lightblue', 
        col.interaction=t(ifelse(flora_matb2$'3'[,] <= 1, adjustcolor('lightblue', alpha.f = 1.5), 
                                 adjustcolor('black', alpha.f = 1))), 
        bor.col.interaction = NA, bor.col.high = NA)
dev.off()

pdf("plotweb.stage4.sp.pdf", 9, 5)
plotweb(flora_matb3$'4',text.rot=90, labsize = 0.5, y.lim=c(-1,3), arrow="down.center", col.high = 'lightblue', 
        col.interaction=t(ifelse(flora_matb2$'4'[,] <= 1, adjustcolor('lightblue', alpha.f = 1.5), 
                                 adjustcolor('black', alpha.f = 1))), #interaction<=1:black
        bor.col.interaction = NA, bor.col.high = NA)
dev.off()



###### ----------------
nrow(quadratdb[quadratdb$Order_I == "Diptera",]) #to count number of individual insects
unique(transectdb$Family_I) #to count insect family
nrow(quadrantdb[quadrantdb$Stage == "1",]) #count no of interaction = nrow
write.csv(flora_matb, 'flora_matb.csv', row.names = FALSE)
length(unique(flora_db$Family))
nrow(quadratdb[quadratdb$pl == "4.D",]) #to count insect abundance for each plot
length(unique(quadratdb$insectsp))
nrow(transectdb[transectdb$pl == "1.A",])
summaryt_tran$insect_richness
nrow(quadratdb)
nrow(transectdb)
colnames(quadratdb)


######## save your work #####
save.image('ngan_231114.RData')
