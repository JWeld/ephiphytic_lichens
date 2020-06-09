library(tidyverse)
library(vegan)
library(VIM)

# #load data
# Sp_all <- read_csv("Raw_data/Sp_all.csv")
# Env_all <- read_csv("Raw_data/Env_all.csv")
# 
# Env_all <- as.tibble(Env_all)
# Sp_all <- as.tibble(Sp_all)
# Env_all$X1 <- as.factor(Env_all$X1)
# Sp_all$X1 <- as.factor(Sp_all$X1)
# 
# #visualise data. Requires VIM
# matrixplot(Sp_all, interactive = F, sortby = "X1")
# 
# #looks like missing data for one species `Sphaerophorus globosus`
# apply(is.na(Sp_all), 2, which)
# 
# #remove it
# Sp_all$`Sphaerophorus globosus` <- NULL
# 
# #looks better now 
# matrixplot(Sp_all, interactive = F, sortby = "X1")
# 
# #species richness per sample
# Rich <- rowSums(Sp_all[, 2:56]>0)
# Sp_all <- cbind(Sp_all, Rich)
# 
# #shannon div per sample
# Shan <- diversity(Sp_all[,-1],"shannon")
# Sp_all <- cbind(Sp_all, Shan)


#load data
# spe <- read_csv("Raw_data/Sp_all.csv")
# env <- read_csv("Raw_data/Env_all.csv")
# 
# env <- as_tibble(env)
# spe <- as_tibble(spe)
# #env$X1 <- as.factor(env$X1)
# #spe$X1 <- as.factor(spe$X1)
# 
# #visualise data. Requires VIM
# matrixplot(spe, interactive = F, sortby = "X1")
# 
# #looks like missing data for one species `Sphaerophorus globosus`
# apply(is.na(spe), 2, which)
# 
# #remove it
# spe$`Sphaerophorus globosus` <- NULL
# 
# spe <- rename(spe, ID = X1)
# 
# ##### Lägg till data om individuella träd
# spe$Site<-substr(spe$ID,1,2)
# spe$Year<-substr(spe$ID,3,4)
# spe$Plot<-substr(spe$ID,5,7)
# spe$Tree<-substr(spe$ID,9,nchar(as.character(spe$ID)))
# spe$sp<-substr(spe$Tree,1,2)
# spe$id2<-substr(spe$ID,5,nchar(spe$ID))
# head(spe$id2)
# spe <- dplyr::select(spe, Site, Year, Plot, Tree, sp, id2, ID, everything())
# 
# # Gör om yy till yyyy
# spe$Year <- ifelse(as.numeric(as.character(spe$Year)) > 17,paste("19",spe$Year, sep=""), paste("20",spe$Year, sep="")) 
# spe$Year <- as.factor(spe$Year)
# spe$Site <- as.factor(spe$Site)
# spe$Plot <- as.factor(spe$Plot)
# spe$Tree <- as.factor(spe$Tree)
# spe$id2 <- as.factor(spe$id2)
# spe$ID <- as.factor(spe$ID)
# 
# #species richness per sample
# n_sp <- rowSums(spe[, 8:62]>0)
# spe <- cbind(spe, n_sp)
# 
# #shannon div per sample
# Shan <- diversity(spe[,-c(1:7)],"shannon")
# spe <- cbind(spe, Shan)
# rm(Shan)
# # 
# # #tidy format
# # spe_tidy <- spe %>%
# #   pivot_longer(-ID, names_to = "spe", values_to = "cover")
# # spe_tidy <- spe_tidy %>% dplyr::filter(cover>0)
# # spe_tidy <- spe_tidy %>% group_by(ID) %>% mutate(n_sp = length(unique(spe)))
# # spe_tidy <- spe_tidy %>% mutate(tot_cov = sum(cover))
# # 
# # 
# 
# 
# ############# Ta bort asp, sälg och ek
# keep <- c("PS","PA","BZ","PINU SYL","PICE ABI","BETULA Z")
# spe <- dplyr::filter(spe, sp %in% keep)
# env <- dplyr::filter(env, TreeSp %in% keep)
# 
# 
# #ADD ELLENBERG/HULT VALUES XXXXX
# 
# #common species?
# colSums(spe[,8:64])
# 
# #Make combined responses and env variables dataframe
# Sp2 <- select(spe, c(ID, Shan, n_sp))
# env <- env %>% rename(ID = X1)
# combi <- left_join(env, Sp2, by = "ID")
# #combi$Year <- as.factor(combi$Year)
# combi$Site <- as.factor(combi$Site)
# combi$TreeSp <- as.factor(combi$TreeSp)
# 
# #change year to full (i.e. 1995 not 95)
# library(lubridate)
# #Add 4-digit years. R assumes that 2-digits years are 00-68 = 2000-2068 and 69-99 = 1969-1999
# combi$Year <- parse_datetime(combi$Year, "%y")
# combi$Year <- year(combi$Year)
# 
# #index of year starting 1996
# combi$Year_i <- I(combi$Year - 1996)
# 
# #make plotID
# combi$Tree_ID <- str_sub(combi$ID,-4,-1) %>% str_remove(.,"-")
# combi$Plot <- str_sub(combi$ID,5,14) %>% word(.,1,sep = "\\-") %>% str_remove_all(.," ")
# 
# combi$Plot <- as.factor(combi$Plot)
# combi$Tree_ID <- as.factor(combi$Tree_ID)
# #combi <- combi %>% rename(ID = X1)

#Ulfs data combo####
library(vegdata)#!
library(ggplot2)#!
library(gridExtra)#!
library(nlme)
library(tidyverse)
# Read and transform data
Lav_spe <- read.delim("~/Documents/R/paper_3/Raw_data/Lav_spe.txt", row.names=1)
spe<-Lav_spe
spe<-spe/400
head(spe)
# Nya artdata med överföring av Bryoria capillaris/fuscescens till Bryoria fuscescens
#spe<-read.table("Lav_spe2.txt", header = T, row.names = 1, sep = "\t")
#spe<-spe/400

#Ellebergs from Wirth
wirth <- read_delim("Raw_data/wirth.csv", ";", escape_double = FALSE,
                    col_types = cols(SUB = col_skip()), trim_ws = TRUE)
wirth <- rename(wirth, species = X1)
wirth$species <- gsub(" ", ".", wirth$species) # Lägg till kolumn med artnamn med punkt
wirth$L <- as.numeric(wirth$L)
wirth$T <- as.numeric(wirth$T)
wirth$K <- as.numeric(wirth$K)
wirth$F <- as.numeric(wirth$F)
wirth$R <- as.numeric(wirth$R)
wirth$N <- as.numeric(wirth$N)


Lav_env <- read.csv("~/Documents/R/paper_3/Raw_data/Lav_env.csv")
env <- Lav_env
env <- rename(env, Site = Område)
env$Site<-paste(substring(env$Site,1,2),env$Inv_nr, sep="_")
env$Site<-factor(env$Site)

känslighet <- read.csv("~/Documents/R/paper_3/Raw_data/känslighet.csv")
känsl <- känslighet
känsl <- column_to_rownames(känsl, "Arter")
#känsl$species<-gsub(" ", ".", rownames(känsl)) # Lägg till kolumn med artnamn med punkt
känsl$species<-rownames(känsl) # Lägg till kolumn med artnamn med punkt
känsl$species <- as.character(känsl$species)
känsl$species<-gsub(" ", ".", känsl$species) # Lägg till kolumn med artnamn med punkt

känsl$X.1 <- NULL
känsl$X.2 <- NULL
känsl$X <- NULL

# Ta bort asp, sälg och ek
spe<-subset(spe, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
############# Ta bort asp, sälg och ek
# keep <- c("PS","PA","BZ","PINU SYL","PICE ABI","BETULA Z")
# spe <- dplyr::filter(spe, sp %in% keep)
# env <- dplyr::filter(env, TreeSp %in% keep)
env$TreeSp <- as.factor(env$TreeSp)
env <- as_tibble(env)
env <- column_to_rownames(env, "Art")


#### Hultengrens index and Ellenbergs####
känsl2 <- left_join(känsl,select(wirth, -K))

# Beräkna viktat medel = Känslighetsindex
KInd<-isc(veg=spe, trait.db=känsl2, ivname="K", keyname = "species" ,method ='mean')
pHInd<-isc(veg=spe, trait.db=känsl2, ivname="pH.tal", keyname = "species" ,method ='mean')
NInd<-isc(veg=spe, trait.db=känsl2, ivname="N", keyname = "species" ,method ='mean') # För många NA för att var meningsfull. 0 = NA i ISC

#species richness per sample
n_sp <- rowSums(spe>0)
spe <- cbind(spe, n_sp)

#species div per sample
spe <- rownames_to_column(spe, var = "ID")

Shan <- diversity(spe[, 2:58],"shannon")
div <- cbind(spe, Shan)
rm(Shan)
div <- as_tibble(spe) %>% select(ID,n_sp,Shan)


## Lägg till indexen i miljödatafilen
tmp<-merge(env,KInd, by = "row.names")
row.names(tmp)<-tmp$Row.names
tmp$Row.names<-NULL
head(tmp)
tmp1<-merge(tmp,NInd, by = "row.names")
row.names(tmp1)<-tmp1$Row.names
tmp1$Row.names<-NULL
env<-merge(tmp1,pHInd, by = "row.names")
row.names(env)<-env$Row.names
env$Row.names<-NULL
colnames(env)[5:7] <- c("Sens","N","pH_Ind")

head(env)
tmp<-NULL
tmp1<-NULL

# Ta bort nollor
is.na(env$Sens) <- !env$Sens
is.na(env$N) <- !env$N
is.na(env$pH_Ind) <- !env$pH_Ind

env2 <- rownames_to_column(env, var = "ID")


##### Utveckling för individuella träd####

# Skapa datafil med Siteåde och träd-id i separata kolumner
temp <- env2 %>% column_to_rownames(var="ID")
data <-temp[,c(5:7)]
rm(temp)

names(data)
colnames(data)<-c("Sensitivity","Nitrogen", "pH")
str(data)

data$ID<-row.names(data)
data$Site<-substr(row.names(data),1,2)
data$Year<-as.factor(substr(row.names(data),3,4))
data$Plot<-substr(row.names(data),5,7)
data$Tree<-substr(row.names(data),9,nchar(row.names(data)))
data$sp<-substr(data$Tree,1,2)
data$id2<-substr(row.names(data),5,nchar(row.names(data)))

# Gör om yy till yyyy
data$Year <- ifelse(as.numeric(as.character(data$Year)) > 17,paste("19",data$Year, sep=""), paste("20",data$Year, sep="")) 
data$Year <- as.factor(data$Year)

head(data)
str(data)

# Kontroll av data
names(data)


#combine and rename####
temp <- data
temp <- as_tibble(temp)
#temp <- rownames_to_column(temp, var="ID")%>% as_tibble()
#temp$id <- NULL
temp <- temp %>% rename(Site = Site, Year = Year)
temp$Year <- as.numeric(as.character(temp$Year))
data <- left_join(select(spe, ID, n_sp, Shan), temp,  by= "ID")
#data <- left_join(data, env2, by ="ID")

data$ID <- as_factor(data$ID)
data$Site <- as_factor(data$Site)
data$Plot <- as_factor(data$Plot)
data$Tree <- as_factor(data$Tree)
data$sp <- as_factor(data$sp)
data$id2 <- as_factor(data$id2)

#index of year starting 1996
data$Year_i <- I(data$Year - 1996)
data <- drop_na(data)
data <- select(data, ID, Site, Plot, id2, Tree, sp, Year, Year_i, everything())
str(data)

#plot data
ggplot(data,aes(x=Year,y=Sensitivity,colour=Site))+geom_point()+
  geom_smooth(method="loess",alpha=0.3)

ggplot(data,aes(x=Year,y=Shan,colour=Site))+geom_point()+
  geom_smooth(method="loess",alpha=0.3)

ggplot(data,aes(x=Year,y=pH,colour=Site))+geom_point()+
  geom_smooth(method="loess",alpha=0.3)

ggplot(data,aes(x=Year,y=n_sp,colour=Site))+geom_point()+
  geom_smooth(method="loess",alpha=0.3)

ggplot(data,aes(x=Year,y=Nitrogen,colour=Site))+geom_point()+
  geom_smooth(method="loess",alpha=0.3)#
#
#
#saveRDS(data,"data.RDS")  
# 
# 
#library(DataExplorer)
#create_report(data)


#subset by site
data.An <- dplyr::filter(data, Site == "An")
data.Ga <- dplyr::filter(data, Site == "Ga")
data.Gd <- dplyr::filter(data, Site == "Gd")
data.Ki <- dplyr::filter(data, Site == "Ki")

#Regressions####
######## Sensitivity #################

######## Subset med bara Gårdsjön#####
data_Gd<-subset(data, Site=="Gd")
data_Gd$Year<-factor(data_Gd$Year)
names(data_Gd)
#str(data_Gd)

min(data_Gd$Sensitivity)
(a<-which(is.na(data_Gd$Sensitivity)))
which.min(data_Gd$Sensitivity)
row.names(data_Gd[which(is.na(data_Gd$Sensitivity)),])

# Ett värde har NA i "Sensitivity", men det är ologiskt
# Ersätt med värdet på samma träd inventeringen innan
data_Gd[a,1]
data_Gd[a,1]<-data_Gd[match("Gd06E 3-PS5",rownames(data_Gd)),1]
data_Gd[a,1]


# Boxplot
boxplot(Sensitivity~Year,data=data_Gd, main="Gårdsjön", xlab="Year", ylab="Sensitvity") 

# Claudias regressioner

data_Gd$time <- 1
data_Gd$time[data_Gd$Year=="2001"]<-6
data_Gd$time[data_Gd$Year=="2006"]<-11
data_Gd$time[data_Gd$Year=="2011"]<-16

summary(Gd.mod1<-lme(Sensitivity~time,
                     random=~1|Plot/Tree,
                     correlation=corCAR1(form=~time|Plot/Tree),
                     data=data_Gd))
anova(Gd.mod1)
fixed.effects(Gd.mod1)

summary(Gd.mod1<-lme(Sensitivity~time,
                     random=~1|Plot/Tree,
                     correlation=corCAR1(form=~time|Plot/Tree),
                     data=data_Gd))

Gd.mod1
plot(Gd.mod1, col = c(1:nlevels(as.factor(data_Gd$time))), pch = 16, main="Gårdsjön, lichen sensitvity")

library(broom)
library(broom.mixed)
library(dotwhisker)

summary(Gd.mod3<-lme(Sensitivity~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
anova(Gd.mod3)
fixed.effects(Gd.mod3)

plot(Gd.mod3, col = c(1:nlevels(as.factor(data_Gd$time))), pch = 16, main="Gårdsjön, lichen sensitvity; mod3")
tidy(Gd.mod3) %>% dwplot(.,show_intercept = TRUE)
tidy(Gd.mod3) %>% dwplot(.)

######## Subset med bara Aneboda####

data_An<-subset(data, Site=="An")
data_An$Year<-factor(data_An$Year)
levels(data_An$Year)

min(data_An$Sensitivity)
#which.min(data_An$Sensitivity)

# Boxplot
boxplot(Sensitivity~Year,data=data_An, main="Aneboda",xlab="Year", ylab="Sensitvity") 

## Claudias regressioner

data_An$time<-1
data_An$time[data_An$Year=="2002"]<-6
data_An$time[data_An$Year=="2007"]<-11
data_An$time[data_An$Year=="2012"]<-16

summary(An.mod1<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
anova(An.mod1)
fixed.effects(An.mod1)

plot(An.mod1, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16, main="Aneboda, lichen sensitvity")

tidy(An.mod1) %>% dwplot(.,show_intercept = TRUE)
tidy(An.mod1) %>% dwplot(.)
#summary(An.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
#plot(An.mod2, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)

######## Subset med bara Kindla ####

data_Ki<-subset(data, Site=="Ki")
data_Ki$Year<-factor(data_Ki$Year)
levels(data_Ki$Year)

data_Ki$Plot<-as.factor(data_Ki$Plot)
levels(data_Ki$Plot)

min(data_Ki$Sensitivity)
#which.min(data_Ki$Sensitivity)

# Boxplot
boxplot(Sensitivity~Year,data=data_Ki, main="Kindla", xlab="Year", ylab="Sensitvity") 

## Claudias regressioner

data_Ki$time<-1
data_Ki$time[data_Ki$Year=="2004"]<-6
data_Ki$time[data_Ki$Year=="2008"]<-11
data_Ki$time[data_Ki$Year=="2013"]<-16

summary(Ki.mod1<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
anova(Ki.mod1)
fixed.effects(Ki.mod1)

plot(Ki.mod1, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla, lichen sensitvity")
ACF(Ki.mod1, maxLag = 11)
plot(ACF(Ki.mod1, alpha= .05), main="Kindla, sensitivity")

#summary(Ki.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
#plot(Ki.mod2, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")

######## Subset med bara Gammtratten ####

data_Ga<-subset(data, Site=="Ga")
data_Ga$Year<-factor(data_Ga$Year)
levels(data_Ga$Year)

min(data_Ga$Sensitivity)
which.min(data_Ga$Sensitivity)

# Boxplot
boxplot(Sensitivity~Year,data=data_Ga, main="Gammtraten", xlab="Year", ylab="Sensitvity") 

## Claudias regressioner 
data_Ga$time<-1
data_Ga$time[data_Ga$Year=="2005"]<-6
data_Ga$time[data_Ga$Year=="2010"]<-11
data_Ga$time[data_Ga$Year=="2015"]<-16


summary(Ga.mod1<-lme(Sensitivity~time,
                     random=~1|Plot/Tree,
                     correlation=corCAR1(form=~time|Plot/Tree),
                     data=data_Ga))

lmer(Sensitivity~time + (1|Plot/Tree), data = data_Ga)


anova(Ga.mod1)
fixed.effects(Ga.mod1)

plot(Ga.mod1, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten, lichen sensitvity")
ACF(Ga.mod1, maxLag = 11)
plot(ACF(Ga.mod1, alpha= .1), main="title2")

#summary(Ga.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
#plot(Ga.mod2, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten")

# Alla Siteåden i samma boxplot####
#png(file = "Boxplot_Sensitivity.png", bg = "white", width = 1920, height = 1920, res = 400, pointsize = 12)

par(mfrow = c(2,2))
boxplot(Sensitivity~Year,data=data_Gd, main="Gårdsjön", xlab="Year", ylab="Sensitvity") 
boxplot(Sensitivity~Year,data=data_An, main="Aneboda", xlab="Year", ylab="Sensitvity") 
boxplot(Sensitivity~Year,data=data_Ki, main="Kindla", xlab="Year", ylab="Sensitvity") 
boxplot(Sensitivity~Year,data=data_Ga, main="Gammtraten", xlab="Year", ylab="Sensitvity") 

######## pH #################

# Subset med bara Gårdsjön#####
names(data_Gd)

(a<-which(is.na(data_Gd$pH)))
which(data_Gd$pH==0)
row.names(data_Gd[which(data_Gd$pH==0),])

# Ett vÅrde har 0 i "pH", men det År ologiskt
# ErsÅtt med vÅrdet pÅ samma trÅd inventeringen innan
data_Gd[a,3]
data_Gd[a,3]<-data_Gd[match("Gd06E 3-PS5",rownames(data_Gd)),3]
data_Gd[a,3]

# Boxplot
boxplot(pH~Year,data=data_Gd, main="Gårdsjön", xlab="Year", ylab="pH") 

plot(data_Gd$time, data_Gd$pH)

## Claudias regressioner 

data_Gd$time<-1
data_Gd$time[data_Gd$Year=="2001"]<-6
data_Gd$time[data_Gd$Year=="2006"]<-11
data_Gd$time[data_Gd$Year=="2011"]<-16

summary(Gd_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd[-55,]))
anova(Gd_pH.mod1)
fixed.effects(Gd_pH.mod1)

plot(Gd_pH.mod1, col = c(1:nlevels(as.factor(data_Gd$time))), pch = 16, main="Gårdsjön, pH")

#summary(Gd_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
#plot(Gd_pH.mod2, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)


### Test utan extremvÅden; Gårdsjön
# GÅller endast dÅ alla trÅdslag År med

#which(data_Gd$pH > 3) # Hitta rader med pH > 3
#data_Gd_1<-data_Gd[-c(20,40,60),] # Ta bort rader med pH >3

# Ta bort eventuellt tomma rader och kolumner
#dim(data_Gd_1)
#data_Gd_1<-data_Gd_1[ , which(!apply(data_Gd_1==0,2,all))]# Ta bort tomma kolumner
#data_Gd_1<-data_Gd_1[!rowSums(data_Gd_1, na.rm=TRUE) == 0,] # Ta bort tomma rader
#dim(data_Gd_1)


#data_Gd_1$time<-1
#data_Gd_1$time[data_Gd_1$Year=="2001"]<-6
#data_Gd_1$time[data_Gd_1$Year=="2006"]<-11
#data_Gd_1$time[data_Gd_1$Year=="2011"]<-16

#summary(Gd_pH.mod1_1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd_1[-55,]))
#anova(Gd_pH.mod1_1)
#fixed.effects(Gd_pH.mod1_1)

#plot(Gd_pH.mod1_1, col = c(1:nlevels(as.factor(data_Gd$time))), pch = 16, main="Gårdsjön, utan avvikare")

######## Subset med bara Aneboda ####

# Boxplot
boxplot(pH~Year,data=data_An, main="Aneboda", xlab="Year", ylab="pH") 

## Claudias regressioner

summary(An_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
anova(An_pH.mod1)
fixed.effects(An_pH.mod1)

plot(An_pH.mod1, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16, main="Aneboda, pH")

#summary(An_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
#plot(An_pH.mod2, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)


### Test utan extremvÅden; Aneboda
# Gäller endast då alla trädslag är med

which(data_An$pH > 3) # Hitta rader med pH > 3
data_An_1<-data_An[-50,] # Ta bort rader med pH >3
names(data_An_1)

# Ta bort eventuellt tomma rader och kolumner
dim(data_An_1)
data_An_1<-data_An_1[ , which(!apply(data_An_1==0,2,all))]# Ta bort tomma kolumner
data_An_1<-data_An_1[!rowSums(data_An_1, na.rm=TRUE) == 0,] # Ta bort tomma rader
dim(data_An_1)

summary(An_pH.mod1_1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An_1))
anova(An_pH.mod1_1)
fixed.effects(An_pH.mod1_1)

plot(An_pH.mod1_1, col = c(1:nlevels(as.factor(data_An$time))), pch = 16, main="Aneboda, pH, utan avvikare")

summary(An_pH.mod2_1<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An_1))
plot(An_pH.mod2_1, col = as.numeric(factor(data_An_1$time, levels = unique(data_An_1$time))), pch = 16)

######## Subset med bara Kindla ####

# Boxplot
boxplot(pH~Year,data=data_Ki, main="Kindla", xlab="Year", ylab="pH") 

## Claudias regressioner

summary(Ki_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
anova(Ki_pH.mod1)
fixed.effects(Ki_pH.mod1)

plot(Ki_pH.mod1, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla, pH")

#summary(Ki_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
#plot(Ki_pH.mod2, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")



######## Subset med bara Gammtratten ####

# Boxplot
boxplot(pH~Year,data=data_Ga, main="Gammtratten", xlab="Year", ylab="pH") 

## Claudias regressioner

summary(Ga_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
anova(Ga_pH.mod1)
fixed.effects(Ga_pH.mod1)

plot(Ga_pH.mod1, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten, pH")
ACF(Ga_pH.mod1, maxLag = 11)
plot(ACF(Ga_pH.mod1, alpha= .0001), main="title2")

acf(residuals(Ga_pH.mod1,type="normalized"))


#summary(Ga_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
#plot(Ga_pH.mod2, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")


# Alla Siteåden i samma boxplot####
par(mfrow = c(2, 2))
boxplot(pH~Year,data=data_Gd, main="Gårdsjön, pH", xlab="Year", ylab="pH") #ylim = c(1, 4)
boxplot(pH~Year,data=data_An, main="Aneboda, pH", xlab="Year", ylab="pH") 
boxplot(pH~Year,data=data_Ki, main="Kindla, pH", xlab="Year", ylab="pH") 
boxplot(pH~Year,data=data_Ga, main="Gammtratten, pH", xlab="Year", ylab="pH") 
par(mfrow = c(1, 1))

######## Kväve #################

######## Subset med bara Gårdsjön ####

summary(Gd.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod1, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)

summary(Gd.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree, data=data_Gd))
summary(Gd.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree, data=data_Gd))

summary(Gd.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod2, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)

######## Subset med bara Aneboda ####
## Claudias regressioner

summary(An.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod1, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)

summary(An.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod2, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)


######## Subset med bara Kindla ####

## Claudias regressioner

summary(Ki.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod1, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")

summary(Ki.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod2, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")


######## Subset med bara Gammtratten ####

## Claudias regressioner

summary(Ga.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod1, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")

summary(Ga.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod2, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")


# Alla Siteåden i samma boxplot####
par(mfrow = c(2, 2))
boxplot(Nitrogen ~Year,data=data_Gd, main="Gårdsjön, Nitrogen", xlab="Year", ylab="Nitrogen") #ylim = c(1, 4)
boxplot(Nitrogen ~Year,data=data_An, main="Aneboda, Nitrogen", xlab="Year", ylab="Nitrogen") 
boxplot(Nitrogen ~Year,data=data_Ki, main="Kindla, Nitrogen", xlab="Year", ylab="Nitrogen") 
boxplot(Nitrogen ~Year,data=data_Ga, main="Gammtratten, Nitrogen", xlab="Year", ylab="Nitrogen") 
par(mfrow = c(1, 1))



