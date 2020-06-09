library(tidyverse)
library(vegan)
library(VIM)
library(vegdata)
library(nlme)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(sjPlot)

# Read and transform data
Lav_spe <- read.delim("~/Documents/R/paper_3/Raw_data/Lav_spe.txt", row.names=1)
spe<-Lav_spe
spe<-spe/400
head(spe)
# Nya artdata med överföring av Bryoria capillaris/fuscescens till Bryoria fuscescens
#spe<-read.table("Lav_spe2.txt", header = T, row.names = 1, sep = "\t")
#spe<-spe/400

Lav_env <- read.csv("~/Documents/R/paper_3/Raw_data/Lav_env.csv")
env <- Lav_env
env$Område<-paste(substring(env$Område,1,2),env$Inv_nr, sep="_")
env$Område<-factor(env$Område)

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

# Ta bort asp, sälg och ek???
#spe<-subset(spe, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
#env<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
############# Ta bort asp, sälg och ek
# keep <- c("PS","PA","BZ","PINU SYL","PICE ABI","BETULA Z")
# spe <- filter(spe, sp %in% keep)
# env <- filter(env, TreeSp %in% keep)
env$TreeSp <- as.factor(env$TreeSp)
env <- as_tibble(env)
env <- column_to_rownames(env, "Art")

#### Hultengrens index####
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

#### combine H index and Ellenbergs####
känsl2 <- left_join(känsl,select(wirth, -K))
känsl3 <- full_join(känsl,select(wirth, -K))#more matches with ellenbergs than in hultengren, use for CWM N?
# Beräkna viktat medel = Känslighetsindex
KInd<-isc(veg=spe, trait.db=känsl2, ivname="K", keyname = "species" ,method ='mean')
pHInd<-isc(veg=spe, trait.db=känsl2, ivname="pH.tal", keyname = "species" ,method ='mean')
NInd<-isc(veg=spe, trait.db=känsl2, ivname="N", keyname = "species" ,method ='mean') 
RInd<-isc(veg=spe, trait.db=känsl2, ivname="R", keyname = "species" ,method ='mean') 
TInd<-isc(veg=spe, trait.db=känsl2, ivname="T", keyname = "species" ,method ='mean') 

#species richness per sample
spe$n_sp<-apply(spe>0,1,sum)
#n_spp <- specnumber(spe)


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
colnames(env)[5:7] <- c("Sens","N_Ind","pH_Ind")
head(env)
tmp<-NULL
tmp1<-NULL

# Ta bort nollor
is.na(env$Sens) <- !env$Sens
is.na(env$N_Ind) <- !env$N_Ind
is.na(env$pH_Ind) <- !env$pH_Ind

##### Utveckling för individuella träd####

# Skapa datafil med Område och träd-id i separata kolumner
IM_ind<-env[,c(5:7)]

names(IM_ind)

colnames(IM_ind)<-c("Sensitivity","Nitrogen", "pH")
str(IM_ind)

IM_ind$id<-row.names(IM_ind)
IM_ind$Omr<-substr(row.names(IM_ind),1,2)
IM_ind$År<-as.factor(substr(row.names(IM_ind),3,4))
IM_ind$Plot<-substr(row.names(IM_ind),5,7)
IM_ind$Tree<-substr(row.names(IM_ind),9,nchar(row.names(IM_ind)))
IM_ind$sp<-substr(IM_ind$Tree,1,2)
IM_ind$id2<-substr(row.names(IM_ind),5,nchar(row.names(IM_ind)))

# Gör om yy till yyyy
IM_ind$År <- ifelse(as.numeric(as.character(IM_ind$År)) > 17,paste("19",IM_ind$År, sep=""), paste("20",IM_ind$År, sep="")) 
IM_ind$År <- as.factor(IM_ind$År)

#translate
IM_ind <- rename(IM_ind, ID=id, Site=Omr, Year=År)
#add shannon diversity
#species div per sample
Shan <- diversity(select(spe, -n_sp),"shannon")
div <- cbind(spe, Shan)
div <- rownames_to_column(div,var="ID")
rm(Shan)
div <- as_tibble(div) %>% select(ID,Shan)
data <- left_join(IM_ind, div, by="ID")
#add richness
n_spp <- spe %>% rownames_to_column(var="ID") %>%  select(ID, n_sp)
data <- left_join(data, n_spp, by = "ID")

str(data)
data <- mutate(data, ID2=ID)
data <- column_to_rownames(data, var="ID2")
databak <- data
#Add richness of sensitive species only####
#all species we have values for
summary(känsl2$K) #mean 5.173, median 5, 1stQ 4, 3rdQ 6
summary(känsl2$N) #mean 3.859, median 3, 1stQ 2, 3rdQ 5
summary(känsl3$N) #mean 3.696, median 3, 1stQ 2, 3rdQ 5
spe_list <- colnames(spe)
#for species present in surveys
filter(känsl3, species %in% spe_list) %>% summary #31 spp S: median 4 1stQ 3, 3rdQ 6
filter(känsl3, species %in% spe_list) %>% summary #43 spp N: median 2 1stQ 2. 3rdQ 3

#use 1st and 3rd quartiles to define sensitive/insensitive species,
#as the lower quartile is central to the lower half of the data and
#the upper quartile is central to the upper half of the data.
#Use values for all species we have data for to define "sensitive" or
#only those that are in surveys? 

S_sens_list <- filter(känsl2, K >= 6) %>% select(species)
S_sens_list <- S_sens_list$species
# S_insens_list <- filter(känsl2, K < 3) %>% select(species)
# S_insens_list <- S_insens_list$species

N_sens_list <- filter(känsl3, N <=2) %>% select(species)
N_sens_list <- N_sens_list$species
# N_insens_list <- filter(känsl3, N >3) %>% select(species)
# N_insens_list <- N_insens_list$species
cyanolist <- c("Lobaria.pulmonaria","Peltigera.collina","Leptogium.lichenoides",
               "Parmeliella.triptophylla","Nephroma.parile") 
#only on site one year for one species. Not really worth it.

spe_S_sens <- select(spe, one_of(S_sens_list))
spe_S_sens$n_ss_sp <-apply(spe_S_sens>0,1,sum)
spe_S_sens <- rownames_to_column(spe_S_sens, var = "ID")
S_sens_count <- select(spe_S_sens, ID, n_ss_sp)

# spe_S_insens <- select(spe, one_of(S_insens_list))
# spe_S_insens$n_si_sp <-apply(spe_S_insens>0,1,sum)
# spe_S_insens <- rownames_to_column(spe_S_insens, var = "ID")
# S_insens_count <- select(spe_S_insens, ID, n_si_sp)

spe_N_sens <- select(spe, one_of(N_sens_list))
spe_N_sens$n_ns_sp <-apply(spe_N_sens>0,1,sum)
spe_N_sens <- rownames_to_column(spe_N_sens, var = "ID")
N_sens_count <- select(spe_N_sens, ID, n_ns_sp)

# spe_N_insens <- select(spe, one_of(N_insens_list))
# spe_N_insens$prop_ns <-apply(spe_N_insens>0,1,sum)
# spe_N_insens <- rownames_to_column(spe_N_insens, var = "ID")
# N_insens_count <- select(spe_N_insens, ID, prop_ns)


sens_counts <- full_join(S_sens_count, N_sens_count)
data <- left_join(data, sens_counts)
data$ID <- as.factor(data$ID)

#E 3-PS5 has NA values for Sensitivity,ph and N in 2011
#use values for same tree from the previous survey 2006? or remove 2011 record from dataset?
(a<-which(is.na(data$pH)))
rownames(data)[a]
(a<-which(is.na(data$Nitrogen)))
rownames(data)[a]
data %>% filter(data$id2 == "E 3-PS5")
data[215,1] <- 3.042857
data[215,2] <- 2.757143
data[215,3] <- 5.085714

data$ID <- as.factor(data$ID)
data$id2 <- as.factor(data$id2)
data$Site <- as.factor(data$Site)
data$Plot <- as.factor(data$Plot)
data$Tree <- as.factor(data$Tree)
data$sp <- as.factor(data$sp)

data2 <- cbind(data, TInd)
data2 <- cbind(data2, RInd)
data3 <- filter(data2, ID != "Gd11E 3-PS5") #remove?
data <- data3

#proportion sensitive columns####
data3 <- mutate(data3, prop_ss = n_ss_sp/n_sp)
data3 <- mutate(data3, prop_ns = n_ns_sp/n_sp)
data <- data3

#plot data####
ggplot(data,aes(x=as.numeric(as.character(Year)),y=Sensitivity,colour=Site))+geom_point()+
  geom_smooth(method="loess") 

ggplot(data,aes(x=as.numeric(as.character(Year)),y=Shan,colour=Site))+geom_point()+
  geom_smooth(method="loess")

ggplot(data,aes(x=as.numeric(as.character(Year)),y=n_sp,colour=Site))+geom_point()+
  geom_smooth(method="loess")

ggplot(data,aes(x=as.numeric(as.character(Year)),y=Nitrogen,colour=Site))+geom_point()+
  geom_smooth(method="loess")

ggplot(data,aes(x=as.numeric(as.character(Year)),y=pH,colour=Site))+geom_point()+
  geom_smooth(method="loess")

ggplot(data,aes(x=as.numeric(as.character(Year)),y=RInd,colour=Site))+geom_point()+
  geom_smooth(method="loess")

ggplot(filter(data, TInd > 2),aes(x=as.numeric(as.character(Year)),y=TInd,colour=Site))+geom_point()+
  geom_smooth(method="lm")  #+ facet_wrap(~Site)


ggplot(data,aes(x=as.numeric(as.character(Year)),y=prop_ns,colour=Site))+geom_point()+
  geom_smooth(method="lm")  #+ facet_wrap(~Site)

ggplot(data,aes(x=as.numeric(as.character(Year)),y=prop_ns,colour=Site))+geom_point()+
  geom_smooth(method="lm")  #+ facet_wrap(~Site)
 #+ facet_wrap(~Site)

#
#saveRDS(data,"data.RDS")  
# 
# 
#library(DataExplorer)
#create_report(data)

dep_means <- readRDS("dep_means.RDS")
library(stringr)
deposition <- dep_means %>%  filter(str_detect(ID, "SE"))
#saveRDS(data,"data.RDS")  


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

# Claudias regressioner

data_Gd$time <- 1
data_Gd$time[data_Gd$Year=="2001"]<-6
data_Gd$time[data_Gd$Year=="2006"]<-11
data_Gd$time[data_Gd$Year=="2011"]<-16

summary(Gd.mod1.s<-lme(Sensitivity~time,
                     random=~1|Plot/Tree,
                     correlation=corCAR1(form=~time|Plot/Tree),
                     data=data_Gd))
anova(Gd.mod1.s)
fixed.effects(Gd.mod1.s)

Gd.mod1.s
plot(Gd.mod1.s, col = c(1:nlevels(as.factor(data_Gd$time))), pch = 16, main="Gårdsjön, lichen sensitvity")

# library(broom)
# library(broom.mixed)
# library(dotwhisker)

summary(Gd.mod3.s<-lme(Sensitivity~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
anova(Gd.mod3.s)
fixed.effects(Gd.mod3.s)


######## Subset med bara Aneboda####

data_An<-subset(data, Site=="An")
data_An$Year<-factor(data_An$Year)
levels(data_An$Year)

min(data_An$Sensitivity)
#which.min(data_An$Sensitivity)

## Claudias regressioner

data_An$time<-1
data_An$time[data_An$Year=="2002"]<-6
data_An$time[data_An$Year=="2007"]<-11
data_An$time[data_An$Year=="2012"]<-16

summary(An.mod1.s<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
anova(An.mod1.s)
fixed.effects(An.mod1.s)

plot(An.mod1.s, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16, main="Aneboda, lichen sensitvity")

tidy(An.mod1.s) %>% dwplot(.,show_intercept = TRUE)
tidy(An.mod1.s) %>% dwplot(.)
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

## Claudias regressioner

data_Ki$time<-1
data_Ki$time[data_Ki$Year=="2004"]<-6
data_Ki$time[data_Ki$Year=="2008"]<-11
data_Ki$time[data_Ki$Year=="2013"]<-16

summary(Ki.mod1.s<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
anova(Ki.mod1.s)
fixed.effects(Ki.mod1.s)

plot(Ki.mod1.s, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla, lichen sensitvity")
ACF(Ki.mod1.s, maxLag = 11)
plot(ACF(Ki.mod1.s, alpha= .05), main="Kindla, sensitivity")

#summary(Ki.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
#plot(Ki.mod2, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")

######## Subset med bara Gammtratten ####

data_Ga<-subset(data, Site=="Ga")
data_Ga$Year<-factor(data_Ga$Year)
levels(data_Ga$Year)

min(data_Ga$Sensitivity)
which.min(data_Ga$Sensitivity)

## Claudias regressioner 
data_Ga$time<-1
data_Ga$time[data_Ga$Year=="2005"]<-6
data_Ga$time[data_Ga$Year=="2010"]<-11
data_Ga$time[data_Ga$Year=="2015"]<-16


summary(Ga.mod1.s<-lme(Sensitivity~time,
                     random=~1|Plot/Tree,
                     correlation=corCAR1(form=~time|Plot/Tree),
                     data=data_Ga))

anova(Ga.mod1.s)
fixed.effects(Ga.mod1.s)

plot(Ga.mod1.s, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten, lichen sensitvity")
ACF(Ga.mod1.s, maxLag = 11)
plot(ACF(Ga.mod1.s, alpha= .1), main="title2")

#summary(Ga.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
#plot(Ga.mod2, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten")

#function to plot all sites in ggplot####
customplot <- function(data, x, y){
  name <- deparse(substitute(data)) #for title
  gg <- ggplot(data, aes_q(x = substitute(x),
                         y = substitute(y)
  ))
       gg <- gg +        
                stat_boxplot(geom = 'errorbar') +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter(alpha = 0.15, size = 2, width = 0.07, height = 0) +
                geom_smooth(method="loess",se=FALSE, aes(group=1))+
                ggtitle(name) +
                scale_x_discrete() + 
                scale_y_continuous(limits = c(0, 9))+
                theme(axis.text = element_text(size = 12),
                      axis.title.x = element_text(size=13),
                      axis.title.y = element_text(size=13),
                      plot.title = element_text(size=14, face="bold"),
                      panel.background = element_rect(fill = 'white', colour = "black"))
       return(gg)
}
# Remove outliers when overlaying boxplot with original data points
#p + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2)

#change data names on copies to make labelling plots easier
Gårdsjön <- data_Gd
Aneboda <- data_An
Kindla <- data_Ki
Gammtratten <- data_Ga

a <- customplot(Gårdsjön, x= Year, y=Sensitivity)
b <- customplot(Aneboda, x= Year, y=Sensitivity)
c <- customplot(Kindla,x= Year, y=Sensitivity)
d <- customplot(Gammtratten,x= Year, y=Sensitivity)

library(patchwork)
a+b+c+d

# ggsave(
#   file="sens.tiff",
#   plot = last_plot(),
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm"),
#   dpi = 320,
#   limitsize = TRUE
# )
# Alla Siteåden i samma boxplot####
#png(file = "Boxplot_Sensitivity.png", bg = "white", width = 1920, height = 1920, res = 400, pointsize = 12)

# par(mfrow = c(2,2))
# boxplot(Sensitivity~Year,data=data_Gd, main="Gårdsjön", xlab="Year", ylab="Sensitvity") 
# boxplot(Sensitivity~Year,data=data_An, main="Aneboda", xlab="Year", ylab="Sensitvity") 
# boxplot(Sensitivity~Year,data=data_Ki, main="Kindla", xlab="Year", ylab="Sensitvity") 
# boxplot(Sensitivity~Year,data=data_Ga, main="Gammtraten", xlab="Year", ylab="Sensitvity") 
# par(mfrow = c(1,1))

plot_model(Ga.mod1.s,type = "diag") #sjPlot

# tab_model(Gd.mod1.s,An.mod1.s,Ki.mod1.s,Ga.mod1.s)
# tab_model(Gd.mod1.s,show.icc = FALSE,show.re.var = FALSE,p.val = "kr", show.df = TRUE)

######## pH #################

# Subset med bara Gårdsjön#####
names(data_Gd)

(a<-which(is.na(data_Gd$pH)))
rownames(data_Gd)[a]

## Claudias regressioner 

data_Gd$time<-1
data_Gd$time[data_Gd$Year=="2001"]<-6
data_Gd$time[data_Gd$Year=="2006"]<-11
data_Gd$time[data_Gd$Year=="2011"]<-16

summary(Gd.mod1.ph<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd[-55,]))
anova(Gd.mod1.ph)
fixed.effects(Gd.mod1.ph)

plot(Gd.mod1.ph, col = c(1:nlevels(as.factor(data_Gd$time))), pch = 16, main="Gårdsjön, pH")

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

## Claudias regressioner

summary(An.mod1.ph<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
anova(An.mod1.ph)
fixed.effects(An.mod1.ph)

plot(An.mod1.ph, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16, main="Aneboda, pH")

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

summary(An.mod1.1.ph<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An_1))
anova(An.mod1.1.ph)
fixed.effects(An.mod1.1.ph)

plot(An.mod1.1.ph, col = c(1:nlevels(as.factor(data_An$time))), pch = 16, main="Aneboda, pH, utan avvikare")

summary(An.mod2.1.ph<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An_1))
plot(An.mod2.1.ph, col = as.numeric(factor(data_An_1$time, levels = unique(data_An_1$time))), pch = 16)

######## Subset med bara Kindla ####

## Claudias regressioner

summary(Ki_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
anova(Ki_pH.mod1)
fixed.effects(Ki_pH.mod1)

plot(Ki_pH.mod1, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla, pH")

#summary(Ki_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
#plot(Ki_pH.mod2, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")



######## Subset med bara Gammtratten ####

## Claudias regressioner

summary(Ga.mod1.ph<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
anova(Ga.mod1.ph)
fixed.effects(Ga.mod1.ph)

plot(Ga.mod1.ph, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten, pH")
ACF(Ga.mod1.ph, maxLag = 11)
plot(ACF(Ga.mod1.ph, alpha= .0001), main="title2")

acf(residuals(Ga.mod1.ph,type="normalized"))


#summary(Ga_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
#plot(Ga_pH.mod2, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")

a <- customplot(Gårdsjön, x= Year, y=pH)
b <- customplot(Aneboda, x= Year, y=pH)
c <- customplot(Kindla,x= Year, y=pH)
d <- customplot(Gammtratten,x= Year, y=pH)

#library(patchwork)
a+b+c+d

# Alla Siteåden i samma boxplot####
# par(mfrow = c(2, 2))
# boxplot(pH~Year,data=data_Gd, main="Gårdsjön, pH", xlab="Year", ylab="pH") #ylim = c(1, 4)
# boxplot(pH~Year,data=data_An, main="Aneboda, pH", xlab="Year", ylab="pH") 
# boxplot(pH~Year,data=data_Ki, main="Kindla, pH", xlab="Year", ylab="pH") 
# boxplot(pH~Year,data=data_Ga, main="Gammtratten, pH", xlab="Year", ylab="pH") 
# par(mfrow = c(1, 1))

######## Kväve #################

######## Subset med bara Gårdsjön ####
summary(Gd.mod1.n<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod1.n, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)
anova(Gd.mod1.n)


summary(Gd.mod2.n<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod2.n, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)

######## Subset med bara Aneboda ####
## Claudias regressioner

summary(An.mod1.n <-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod1.n, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)
anova(An.mod1.n)

summary(An.mod2.n <-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod2.n, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)


######## Subset med bara Kindla ####

## Claudias regressioner

summary(Ki.mod1.n<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod1.n, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")
anova(Ki.mod1.n)

summary(Ki.mod2.n<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod2.n, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")


######## Subset med bara Gammtratten ####

## Claudias regressioner

summary(Ga.mod1.n<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod1.n, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")
anova(Ga.mod1.n)

summary(Ga.mod2.n<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod2.n, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")

a <- customplot(Gårdsjön, x= Year, y=Nitrogen)
b <- customplot(Aneboda, x= Year, y=Nitrogen)
c <- customplot(Kindla,x= Year, y=Nitrogen)
d <- customplot(Gammtratten,x= Year, y=Nitrogen)

#library(patchwork)
a+b+c+d



# Alla Siteåden i samma boxplot####
# par(mfrow = c(2, 2))
# boxplot(Nitrogen ~Year,data=data_Gd, main="Gårdsjön, Nitrogen", xlab="Year", ylab="Nitrogen") #ylim = c(1, 4)
# boxplot(Nitrogen ~Year,data=data_An, main="Aneboda, Nitrogen", xlab="Year", ylab="Nitrogen") 
# boxplot(Nitrogen ~Year,data=data_Ki, main="Kindla, Nitrogen", xlab="Year", ylab="Nitrogen") 
# boxplot(Nitrogen ~Year,data=data_Ga, main="Gammtratten, Nitrogen", xlab="Year", ylab="Nitrogen") 
# par(mfrow = c(1, 1))


#Diversity#####


######## Subset med bara Gårdsjön ####
summary(Gd.mod1.d <-lme(Shan~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod1.d, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)
anova(Gd.mod1.d)

summary(Gd.mod1.d <-lme(Shan~time,random=~1|Plot/Tree, data=data_Gd))

summary(Gd.mod2.d <-lme(Shan~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod2.d, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)

######## Subset med bara Aneboda ####
## Claudias regressioner

summary(An.mod1.d <-lme(Shan~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod1.d, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)
anova(An.mod1.d)

summary(An.mod2.d <-lme(Shan~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod2.d, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)


######## Subset med bara Kindla ####

## Claudias regressioner

summary(Ki.mod1.d <-lme(Shan~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod1.d , col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")
anova(Ki.mod1.d)

summary(Ki.mod2.d <-lme(Shan~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod2.d , col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")


######## Subset med bara Gammtratten ####

## Claudias regressioner

summary(Ga.mod1.d <-lme(Shan~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod1.d, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")
anova(Ga.mod1.d)

summary(Ga.mod2.d <-lme(Shan~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod2.d, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten pH")



a <- customplot(Gårdsjön, x= Year, y=Shan)
b <- customplot(Aneboda, x= Year, y=Shan)
c <- customplot(Kindla,x= Year, y=Shan)
d <- customplot(Gammtratten,x= Year, y=Shan)

a+b+c+d

# Alla Siteåden i samma boxplot####
# par(mfrow = c(2, 2))
# boxplot(Shan ~Year,data=data_Gd, main="Gårdsjön, Diversity", xlab="Year", ylab="Shannon index") #ylim = c(1, 4)
# boxplot(Shan ~Year,data=data_An, main="Aneboda, Diversity", xlab="Year", ylab="Shannon index") 
# boxplot(Shan ~Year,data=data_Ki, main="Kindla, Diversity", xlab="Year", ylab="Shannon index") 
# boxplot(Shan ~Year,data=data_Ga, main="Gammtratten, Diversity", xlab="Year", ylab="Shannon index") 
# par(mfrow = c(1, 1))


#Richness#####


######## Subset med bara Gårdsjön ####
summary(Gd.mod1.r <-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod1.r, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)
anova(Gd.mod1.r)

summary(Gd.mod2.r <-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod2.r, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)

######## Subset med bara Aneboda ####
## Claudias regressioner

summary(An.mod1.r <-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod1.r, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)
anova(An.mod1.r)

summary(An.mod2.r <-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod2.r, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)


######## Subset med bara Kindla ####

## Claudias regressioner

summary(Ki.mod1.r <-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod1.r , col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")
anova(Ki.mod1.r)

summary(Ki.mod2.r <-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod2.r , col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")


######## Subset med bara Gammtratten ####

## Claudias regressioner

summary(Ga.mod1.r <-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod1.r, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten")
anova(Ga.mod1.r)

summary(Ga.mod2.r <-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod2.r, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten")

a <- customplot(Gårdsjön, x= Year, y=n_sp)
b <- customplot(Aneboda, x= Year, y=n_sp)
c <- customplot(Kindla,x= Year, y=n_sp)
d <- customplot(Gammtratten,x= Year, y=n_sp)

a+b+c+d
# Alla Siteåden i samma boxplot####
# par(mfrow = c(2, 2))
#boxplot(n_sp ~Year,data=data_Gd, main="Gårdsjön, richness", xlab="Year", ylab="richness") #ylim = c(1, 4)
# boxplot(n_sp ~Year,data=data_An, main="Aneboda, richness", xlab="Year", ylab="richness") 
# boxplot(n_sp ~Year,data=data_Ki, main="Kindla, richness", xlab="Year", ylab="richness") 
# boxplot(n_sp ~Year,data=data_Ga, main="Gammtratten, richness", xlab="Year", ylab="richness") 
# par(mfrow = c(1, 1))


#Richness of S sensitive species####

######## Subset med bara Gårdsjön ####
summary(Gd.mod1.rss <-lme(n_ss_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod1.rss, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)
anova(Gd.mod1.rss) # - ns

######## Subset med bara Aneboda ####
summary(An.mod1.rss <-lme(n_ss_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod1.rss, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)
anova(An.mod1.rss) # - sig

######## Subset med bara Kindla ####
summary(Ki.mod1.rss <-lme(n_ss_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod1.rss , col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")
anova(Ki.mod1.rss) # - ns

######## Subset med bara Gammtratten ####
summary(Ga.mod1.rss <-lme(n_ss_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod1.rss, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten")
anova(Ga.mod1.rss) # - sig

#Boxplot all####
a <- customplot(Gårdsjön, x= Year, y=n_ss_sp)
b <- customplot(Aneboda, x= Year, y=n_ss_sp)
c <- customplot(Kindla,x= Year, y=n_ss_sp)
d <- customplot(Gammtratten,x= Year, y=n_ss_sp)

a+b+c+d

#anova
data.Gd$Year <- as.factor(data.Gd$Year)
anova_S_Gd <- aov(n_ss_sp ~ Year, data = data.Gd)
summary(anova_S_Gd) #ns

data.An$Year <- as.factor(data.An$Year)
anova_S_An <- aov(n_ss_sp ~ Year, data = data.An)
summary(anova_S_An) # sig ***
TukeyHSD(anova_S_An) #97-07, 97-12, 02-12

data.Ki$Year <- as.factor(data.Ki$Year)
anova_S_Ki <- aov(n_ss_sp ~ Year, data = data.Ki)
summary(anova_S_Ki) #ns

data.Ga$Year <- as.factor(data.Ga$Year)
anova_S_Ga <- aov(n_ss_sp ~ Year, data = data.Ga)
summary(anova_S_Ga) # sig *
TukeyHSD(anova_S_Ga) #2000 -2005


#Richness of N sensitive species####

######## Subset med bara Gårdsjön ####
summary(Gd.mod1.rns <-lme(n_ns_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
plot(Gd.mod1.rns, col = as.numeric(factor(data_Gd$time, levels = unique(data_Gd$time))), pch = 16)
anova(Gd.mod1.rns) # - ns

######## Subset med bara Aneboda ####
summary(An.mod1.rns <-lme(n_ns_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_An))
plot(An.mod1.rns, col = as.numeric(factor(data_An$time, levels = unique(data_An$time))), pch = 16)
anova(An.mod1.rns) # - sig

######## Subset med bara Kindla ####
summary(Ki.mod1.rns <-lme(n_ns_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ki))
plot(Ki.mod1.rns, col = as.numeric(factor(data_Ki$time, levels = unique(data_Ki$time))), pch = 16, main = "Kindla")
anova(Ki.mod1.rns) # + ns

######## Subset med bara Gammtratten ####
summary(Ga.mod1.rns <-lme(n_ns_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=data_Ga))
plot(Ga.mod1.rns, col = as.numeric(factor(data_Ga$time, levels = unique(data_Ga$time))), pch = 16, main = "Gammtratten")
anova(Ga.mod1.rns) # + ns

#Boxplot all####
a <- customplot(Gårdsjön, x= Year, y=n_ns_sp)
b <- customplot(Aneboda, x= Year, y=n_ns_sp)
c <- customplot(Kindla,x= Year, y=n_ns_sp)
d <- customplot(Gammtratten,x= Year, y=n_ns_sp)

a+b+c+d

#anova
data.Gd$Year <- as.factor(data.Gd$Year)
anova_N_Gd <- aov(n_ns_sp ~ Year, data = data.Gd)
summary(anova_N_Gd) #ns

data.An$Year <- as.factor(data.An$Year)
anova_N_An <- aov(n_ns_sp ~ Year, data = data.An)
summary(anova_N_An) # sig ***
TukeyHSD(anova_N_An) #97-07,97-12, 02-12

data.Ki$Year <- as.factor(data.Ki$Year)
anova_N_Ki <- aov(n_ns_sp ~ Year, data = data.Ki)
summary(anova_N_Ki) #ns

data.Ga$Year <- as.factor(data.Ga$Year)
anova_N_Ga <- aov(n_ns_sp ~ Year, data = data.Ga)
summary(anova_N_Ga) #ns


#deposition data from ICP IM#####
IM_dep_long <- read_csv(
  "Raw_data/ICPIM_TF_data.csv",
  col_types = cols(
    DataProcessingDate = col_skip(),
    DataQualityFlag = col_skip(),
    Day = col_skip(),
    DeterInfo = col_skip(),
    InstCode = col_skip(),
    PretreCode = col_skip(),
    PretreInfo = col_skip(),
    SamplingLevel = col_skip(),
    SpatialPool = col_skip(),
    ParameterList = col_skip(),
    DeterCode = col_skip(),
    ParameterInfo = col_skip(),
    NeedleAgeCode = col_skip(),
    StatusFlag = col_skip(),
    Unit = col_skip(),
    ListMedium = col_skip(),
    Subprog = col_skip(),
    YearMonth = col_date(format = "%Y%m")
  )
)


#add year column
IM_dep_long$survey_year <-  year(IM_dep_long$YearMonth)
#make line ID
IM_dep_long <-
  transform(IM_dep_long,
            ID = paste0(IM_dep_long$survey_year, sep = "_", IM_dep_long$AreaCode))
IM_dep_long$ID <- as.character(IM_dep_long$ID)

#ws_check function shows some columns need white space trimming out
IM_dep_long$StationCode <- trimws(IM_dep_long$StationCode, "both")
IM_dep_long$Value <- trimws(IM_dep_long$Value, "both")
IM_dep_long$Sample_ID <- trimws(IM_dep_long$Sample_ID, "both")
IM_dep_long$Value <- as.numeric(IM_dep_long$Value)

#add country (same as above but neater code)
IM_dep_long <-
  mutate(IM_dep_long, country = str_sub(IM_dep_long$AreaCode, 1, 2)) #keep first two letters of country code and use to match
IM_dep_long$country <-
  IM_dep_long$country %>% str_replace_all(
    c(
      "AT" = "Austria",
      "DE" = "Germany",
      "DK" = "Denmark",
      "EE" = "Estonia",
      "FI" = "Finland",
      "IT" = "Italy",
      "LT" = "Lithuania",
      "LV" = "Latvia",
      "NO" = "Norway",
      "PL" = "Poland",
      "RU" = "Russia",
      "SE" = "Sweden",
      "CZ" = "Czech",
      "ES" = "Spain",
      "IE" = "Ireland"
    )
  )

#IM data is long, Forests is wide, spread IM to match
parameters <- c("NH4N", "NO3N", "PH", "SO4S")
IM_dep <- filter(IM_dep_long, ParameterCode %in% parameters) %>%
  spread(., ParameterCode, Value)

# rename columns to match FOR_dep
IM_dep <- IM_dep %>% rename(n_nh4 = NH4N, n_no3 = NO3N, ph = PH, s_so4 = SO4S)
IM_dep_SE <-  IM_dep %>%  filter(str_detect(ID, "SE"))
IM_dep_SE_summaries <- IM_dep_SE %>% drop_na() %>%
  group_by(AreaName, survey_year) %>%
  summarise_at(vars(ph, n_nh4, n_no3, s_so4), funs(sum, mean))
IM_dep_SE_summaries$AreaName <-  str_replace(IM_dep_SE_summaries$AreaName, "GÅRDSJÖN F1", "GÅRDSJÖN")
IM_dep_SE_summaries <- IM_dep_SE_summaries %>% rename(Site = AreaName, Year = survey_year)
IM_dep_SE_summaries$Site <- as.factor(IM_dep_SE_summaries$Site)

write.csv(IM_dep_SE, "IM_dep_SE.csv")
# average per site/year combo
IM_dep_means <- IM_dep %>% group_by(ID) %>%
  summarise_at(vars(ph, n_nh4, n_no3, s_so4), funs(mean), na.rm = TRUE) #mean and sd?

deposition <- IM_dep_means %>%  filter(str_detect(ID, "SE"))
deposition$Site <- word(deposition$ID, 2, sep = "_")
deposition$Year <- word(deposition$ID, 1, sep = "_")
deposition$ID <- as.factor(deposition$ID)
deposition$Site <- as.factor(deposition$Site)
deposition$Site <- str_replace(deposition$Site, "SE04", "Gd")
deposition$Site <- str_replace(deposition$Site, "SE14", "An")
deposition$Site <- str_replace(deposition$Site, "SE15", "Ki")
deposition$Site <- str_replace(deposition$Site, "SE16", "Ga")
deposition$ID <- as.factor(deposition$ID)

#make site year combo ID
deposition <-transform(deposition,
                       ID3 = paste0(deposition$Site, sep = "_", deposition$Year))
#deposition$ID3 <- as.character(deposition$ID3)
data <- transform(data,
                  ID3 = paste0(data$Site, sep = "_", data$Year))
#no dep data for An 1997 and Ki 1998, match with nearest available (An 1999 and Ki 1999)
data$ID3 <-  str_replace(data$ID3, "An_1997", "An_1999")
data$ID3 <-  str_replace(data$ID3, "Ki_1998", "Ki_1999")


#Add deposition data to the rest
data4 <- left_join(data, select(deposition, ID3, n_nh4, n_no3, s_so4), by = "ID3")



#plots####
#Positive correlation between all 3 deposiiton types
#deposition over time
fitted_models <-  IM_dep_SE_summaries %>% group_by(Site) %>% do(model = lm(s_so4_mean ~ Year, data = .))
fitted_models$model
library(broom)
fitted_models %>% tidy(model)
fitted_models %>% glance(model)
#fitted_models %>% augment(model)

fitted_models <-  IM_dep_SE_summaries %>% group_by(Site) %>% do(model = lm(n_nh4_mean ~ Year, data = .))
fitted_models %>% tidy(model)
fitted_models %>% glance(model)
#fitted_models %>% augment(model)

fitted_models <-  IM_dep_SE_summaries %>% group_by(Site) %>% do(model = lm(n_no3_mean ~ Year, data = .))
fitted_models %>% tidy(model)
fitted_models %>% glance(model)

IM_dep_SE_summaries <- ungroup(IM_dep_SE_summaries)

a <- ggplot(IM_dep_SE_summaries,aes(x=as.numeric(as.character(Year)),y=s_so4_mean, colour=Site))+geom_point()+
  geom_smooth(method="lm", se= FALSE) + labs(x = "Year", y = "SO4, mg/l") + 
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))

b <- ggplot(IM_dep_SE_summaries,aes(x=as.numeric(as.character(Year)),y=n_nh4_mean,colour=Site))+
  geom_point()+
  labs(x = "Year", y = "NH4, mg/l") + 
  stat_smooth(method="lm", se= FALSE, data = IM_dep_SE_summaries[which(IM_dep_SE_summaries$Site != "GAMMTRATTEN"),] ) +
  stat_smooth(method='loess', se= FALSE, data = IM_dep_SE_summaries[which(IM_dep_SE_summaries$Site == "GAMMTRATTEN"),]) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))


c <- ggplot(IM_dep_SE_summaries,aes(x=as.numeric(as.character(Year)),y=n_no3_mean,colour=Site))+
  geom_point()+
  labs(x = "Year", y = "NO3, mg/l")+
  stat_smooth(method="lm", se= FALSE, data = IM_dep_SE_summaries[which(IM_dep_SE_summaries$Site != "ANEBODA"),] ) +
  stat_smooth(method='loess', se= FALSE, data = IM_dep_SE_summaries[which(IM_dep_SE_summaries$Site == "ANEBODA"),]) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_rect(fill = 'white', colour = "black"))

(a|b)/(c|plot_spacer())

#full monthly data

fitted_models <-  IM_dep_SE %>% group_by(AreaName) %>% do(model = lm(s_so4 ~ survey_year, data = .))
fitted_models$model
library(broom)
fitted_models %>% tidy(model)
fitted_models %>% glance(model)
#fitted_models %>% augment(model)

fitted_models <-  IM_dep_SE %>% group_by(AreaName) %>% do(model = lm(n_nh4 ~ survey_year, data = .))
fitted_models %>% tidy(model)
fitted_models %>% glance(model)
#fitted_models %>% augment(model)

fitted_models <-  IM_dep_SE %>% group_by(AreaName) %>% do(model = lm(n_no3 ~ survey_year, data = .))
fitted_models %>% tidy(model)
fitted_models %>% glance(model)

IM_dep_SE <- ungroup(IM_dep_SE)

a <- ggplot(IM_dep_SE,aes(x=as.numeric(as.character(survey_year)),y=s_so4, colour=AreaName))+
  geom_point()+
  labs(x = "Year", y = "SO4, mg/l") + 
  geom_smooth(method="lm", se= FALSE) + 
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))

b <- ggplot(IM_dep_SE,aes(x=as.numeric(as.character(survey_year)),y=n_nh4,colour=AreaName))+
  geom_point()+
  labs(x = "Year", y = "NH4, mg/l") + 
  stat_smooth(method="lm", se= FALSE, data = IM_dep_SE[which(IM_dep_SE$AreaName != "GAMMTRATTEN"),] ) +
  stat_smooth(method='loess', se= FALSE, data = IM_dep_SE[which(IM_dep_SE$AreaName == "GAMMTRATTEN"),]) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))


c <- ggplot(IM_dep_SE,aes(x=as.numeric(as.character(survey_year)),y=n_no3,colour=AreaName))+
  geom_point()+
  labs(x = "Year", y = "NO3, mg/l")+
  stat_smooth(method="lm", se= FALSE, data = IM_dep_SE[which(IM_dep_SE$AreaName != "ANEBODA"),] ) +
  stat_smooth(method='loess', se= FALSE, data = IM_dep_SE[which(IM_dep_SE$AreaName == "ANEBODA"),]) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_rect(fill = 'white', colour = "black"))

(a|b)/(c|plot_spacer())



fitted_models <-  IM_dep_SE_summaries %>% group_by(Site) %>% do(model = lm(s_so4_mean ~ Year, data = .))
fitted_models$model
library(broom)
fitted_models %>% tidy(model)
fitted_models %>% glance(model)
#fitted_models %>% augment(model)

fitted_models <-  IM_dep_SE_summaries %>% group_by(Site) %>% do(model = lm(n_nh4_mean ~ Year, data = .))
fitted_models %>% tidy(model)
fitted_models %>% glance(model)
#fitted_models %>% augment(model)

fitted_models <-  IM_dep_SE_summaries %>% group_by(Site) %>% do(model = lm(n_no3_mean ~ Year, data = .))
fitted_models %>% tidy(model)
fitted_models %>% glance(model)

IM_dep_SE_summaries <- ungroup(IM_dep_SE_summaries)
library(ggpmisc)

#plot with R2 and p values
group <- IM_dep_SE_summaries$Site
my.formula <- y ~ x 

a <- ggplot(IM_dep_SE_summaries, aes(Year, s_so4_mean, colour = group)) +
  geom_point() +
  labs(x = "Year", y = "SO4, mg/l") +
  geom_smooth(method = "lm",
              se = FALSE,
              formula = my.formula) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))+
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = formula),
    label.x = "middle",
    label.y = "top",
    aes(label = paste(
      "italic(P)*\"-value = \"*",
      signif(..p.value.., digits = 2), sep = ""
    )),
    parse = TRUE
  )  +
  stat_poly_eq(
    formula = my.formula,
    label.y = "top",
    label.x = "right",
    parse = TRUE
  ) 

b <- ggplot(IM_dep_SE_summaries, aes(Year, n_nh4_mean, colour = group)) +
  geom_point() +
  labs(x = "Year", y = "NH4, mg/l")+
  geom_smooth(method = "lm",
              se = FALSE,
              formula = my.formula)  +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black")) +
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = formula),
    label.x = "middle",
    label.y = "top",
    aes(label = paste(
      "italic(P)*\"-value = \"*",
      signif(..p.value.., digits = 2), sep = ""
    )),
    parse = TRUE
  )  +
  stat_poly_eq(
    formula = my.formula,
    label.y = "top",
    label.x = "right",
    parse = TRUE
  ) 

c <- ggplot(IM_dep_SE_summaries, aes(Year, n_no3_mean, colour = group)) +
  geom_point() +
  labs(x = "Year", y = "NO3, mg/l")+
  geom_smooth(method = "lm",
              se = FALSE,
              formula = my.formula) +
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_rect(fill = 'white', colour = "black")) +
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = formula),
    label.x = "middle",
    label.y = "top",
    aes(label = paste(
      "italic(P)*\"-value = \"*",
      signif(..p.value.., digits = 2), sep = ""
    )),
    parse = TRUE
  )  +
  stat_poly_eq(
    formula = my.formula,
    label.y = "top",
    label.x = "right",
    parse = TRUE
  ) 

(a|b)/(c|plot_spacer())


#deposition over time, using all dep data, and annual totals
a <- ggplot(data4,aes(x=as.numeric(as.character(Year)),y=s_so4,colour=Site))+geom_point()+
  geom_smooth(method="lm", se= FALSE) + labs(x = "Year", y = "SO4, mg/l") + 
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))

b <- ggplot(data4,aes(x=as.numeric(as.character(Year)),y=n_nh4,colour=Site))+geom_point()+
  geom_smooth(method="lm", se= FALSE) + labs(x = "Year", y = "NH4, mg/l") + 
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = "black"))

c <- ggplot(data4,aes(x=as.numeric(as.character(Year)),y=n_no3,colour=Site))+geom_point()+
  geom_smooth(method="lm", se= FALSE) + labs(x = "Year", y = "NO3, mg/l")+
  theme(axis.text = element_text(size = 11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_rect(fill = 'white', colour = "black"))

(a|b)/(c|plot_spacer())

#relation deposition and sensitivity
ggplot(data4,aes(x=s_so4,y=Sensitivity,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_no3,y=Sensitivity,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_nh4,y=Sensitivity,colour=Site))+geom_point()+
  geom_smooth(method="lm") 
#sensitive species proportion
ggplot(data4,aes(x=prop_ss,y=prop_ns,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_no3,y=prop_ns,colour=Site))+geom_point()+
  geom_smooth(method="lm") 



#relation deposition and N
ggplot(data4,aes(x=s_so4,y=Nitrogen,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_no3,y=Nitrogen,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_nh4,y=Nitrogen,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

#relation deposition and div
ggplot(data4,aes(x=s_so4,y=Shan,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_no3,y=Shan,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

ggplot(data4,aes(x=n_nh4,y=Shan,colour=Site))+geom_point()+
  geom_smooth(method="lm") 

#boxplot other variables
a <- customplot(Gårdsjön, x= Year, y=RInd)
b <- customplot(Aneboda, x= Year, y=RInd)
c <- customplot(Kindla,x= Year, y=RInd)
d <- customplot(Gammtratten,x= Year, y=RInd)

a+b+c+d
