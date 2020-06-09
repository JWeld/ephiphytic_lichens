# Analysera lavdata från IM

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

# Ta bort asp, sälg och ek
spe<-subset(spe, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
############# Ta bort asp, sälg och ek
# keep <- c("PS","PA","BZ","PINU SYL","PICE ABI","BETULA Z")
# spe <- filter(spe, sp %in% keep)
# env <- filter(env, TreeSp %in% keep)
env$TreeSp <- as.factor(env$TreeSp)
env <- as_tibble(env)
env <- column_to_rownames(env, "Art")

#### Hultengrens index####

# Beräkna viktat medel = Känslighetsindex
KInd<-isc(veg=spe, trait.db=känsl, ivname="K", keyname = "species" ,method ='mean')
pHInd<-isc(veg=spe, trait.db=känsl, ivname="pH.tal", keyname = "species" ,method ='mean')
NInd<-isc(veg=spe, trait.db=känsl, ivname="N.tal", keyname = "species" ,method ='mean') # För många NA för att var meningsfull. 0 = NA i ISC

#species richness per sample
n_sp <- rowSums(spe>0)
spe <- cbind(spe, n_sp)

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

head(IM_ind)
str(IM_ind)

# Kontroll av data
names(IM_ind)

#write.table(IM_ind, file ="IM_ind.txt", quote = F, sep = "\t", na="NA",dec=",", row.names=TRUE, col.names = TRUE)


######## Sensitivity #################

######## Subset med bara Gårdsjön#####
IM_ind_Gd<-subset(IM_ind, Omr=="Gd")
IM_ind_Gd$År<-factor(IM_ind_Gd$År)
names(IM_ind_Gd)
#str(IM_ind_Gd)

min(IM_ind_Gd$Sensitivity)
(a<-which(is.na(IM_ind_Gd$Sensitivity)))
which.min(IM_ind_Gd$Sensitivity)
row.names(IM_ind_Gd[which(is.na(IM_ind_Gd$Sensitivity)),])

# Ett värde har NA i "Sensitivity", men det är ologiskt
# Ersätt med värdet på samma träd inventeringen innan
IM_ind_Gd[a,1]
IM_ind_Gd[a,1]<-IM_ind_Gd[match("Gd06E 3-PS5",rownames(IM_ind_Gd)),1]
IM_ind_Gd[a,1]


# Everything on the same plot
ggplot(IM_ind_Gd, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gårdsjön") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Gd, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gårdsjön") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(Sensitivity~År,data=IM_ind_Gd, main="Gårdsjön", xlab="År", ylab="Sensitvity") 





## Claudias regressioner

IM_ind_Gd$time<-1
IM_ind_Gd$time[IM_ind_Gd$År=="2001"]<-6
IM_ind_Gd$time[IM_ind_Gd$År=="2006"]<-11
IM_ind_Gd$time[IM_ind_Gd$År=="2011"]<-16

summary(Gd.mod1<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Gd))
anova(Gd.mod1)
fixed.effects(Gd.mod1)

Gd.mod1
plot(Gd.mod1, col = c(1:nlevels(as.factor(IM_ind_Gd$time))), pch = 16, main="Gårdsjön, lichen sensitvity")

# library(coefplot)
# coefplot(Gd.mod1)

library(broom)
library(broom.mixed)
library(dotwhisker)
tidy(Gd.mod1) %>% dwplot(.,show_intercept = TRUE)
tidy(Gd.mod1) %>% dwplot(.)


#fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1|Subject)

#newdat <- expand.grid(Sex=unique(Orthodont$Sex),
#                      age=c(min(Orthodont$age),
#                            max(Orthodont$age)))

ggplot(data=IM_ind_Gd, aes(x=time, y=Sensitivity, colour = sp)) +
  geom_point(size=3) +
  geom_line(aes(y=predict(Gd.mod1), group=sp)) #+
  geom_line(data=newdat, aes(y=predict(fm2, level=0, newdata=newdat), size="Population")) +
  scale_size_manual(name="Predictions", values=c("Subjects"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 




#fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1|Subject)

#newdat <- expand.grid(Sex=unique(Orthodont$Sex),
#                      age=c(min(Orthodont$age),
#                            max(Orthodont$age)))

#library(ggplot2)
#p <- ggplot(Orthodont, aes(x=age, y=distance, colour=Sex)) +
#  geom_point(size=3) +
#  geom_line(aes(y=predict(fm2), group=Subject, size="Subjects")) +
#  geom_line(data=newdat, aes(y=predict(fm2, level=0, newdata=newdat), size="Population")) +
#  scale_size_manual(name="Predictions", values=c("Subjects"=0.5, "Population"=3)) +
#  theme_bw(base_size=22) 
#print(p)


#summary(Gd.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Gd))
#plot(Gd.mod2, col = as.numeric(factor(IM_ind_Gd$time, levels = unique(IM_ind_Gd$time))), pch = 16)

summary(Gd.mod3<-lme(Sensitivity~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Gd))
anova(Gd.mod3)
fixed.effects(Gd.mod3)

plot(Gd.mod3, col = c(1:nlevels(as.factor(IM_ind_Gd$time))), pch = 16, main="Gårdsjön, lichen sensitvity; mod3")
tidy(Gd.mod3) %>% dwplot(.,show_intercept = TRUE)
tidy(Gd.mod3) %>% dwplot(.)


## Test
#t<-IM_ind_Gd
#t$Sensitivity[t$År=="2001"]<-t$Sensitivity+2
#t$Sensitivity[t$År=="2006"]<-t$Sensitivity+5
#t$Sensitivity[t$År=="2011"]<-t$Sensitivity+10
#summary(t_GM.mod1<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=t))
#anova(t_GM.mod1)
#fixed.effects(t_GM.mod1)
#plot(t_GM.mod1, col = as.numeric(factor(t$time, levels = unique(t$time))), pch=16)
#ACF(t_GM.mod1, maxLag = 3)
#plot(ACF(t_GM.mod1, alpha= .05), main="Testdata")
#acf(residuals(Ga_pH.mod1,type="normalized"))


######## Subset med bara Aneboda####

IM_ind_An<-subset(IM_ind, Omr=="An")
IM_ind_An$År<-factor(IM_ind_An$År)
levels(IM_ind_An$År)

min(IM_ind_An$Sensitivity)
#which.min(IM_ind_An$Sensitivity)


# Everything on the same plot
ggplot(IM_ind_An, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_An, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(Sensitivity~År,data=IM_ind_An, main="Aneboda",xlab="År", ylab="Sensitvity") 



## Claudias regressioner

IM_ind_An$time<-1
IM_ind_An$time[IM_ind_An$År=="2002"]<-6
IM_ind_An$time[IM_ind_An$År=="2007"]<-11
IM_ind_An$time[IM_ind_An$År=="2012"]<-16

summary(An.mod1<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An))
anova(An.mod1)
fixed.effects(An.mod1)

plot(An.mod1, col = as.numeric(factor(IM_ind_An$time, levels = unique(IM_ind_An$time))), pch = 16, main="Aneboda, lichen sensitvity")

tidy(An.mod1) %>% dwplot(.,show_intercept = TRUE)
tidy(An.mod1) %>% dwplot(.)
#summary(An.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An))
#plot(An.mod2, col = as.numeric(factor(IM_ind_An$time, levels = unique(IM_ind_An$time))), pch = 16)


######## Subset med bara Kindla ####

IM_ind_Ki<-subset(IM_ind, Omr=="Ki")
IM_ind_Ki$År<-factor(IM_ind_Ki$År)
levels(IM_ind_Ki$År)

IM_ind_Ki$Plot<-as.factor(IM_ind_Ki$Plot)
levels(IM_ind_Ki$Plot)

min(IM_ind_Ki$Sensitivity)
#which.min(IM_ind_Ki$Sensitivity)


# Everything on the same plot
ggplot(IM_ind_Ki, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Ki, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(Sensitivity~År,data=IM_ind_Ki, main="Kindla", xlab="År", ylab="Sensitvity") 



## Claudias regressioner

IM_ind_Ki$time<-1
IM_ind_Ki$time[IM_ind_Ki$År=="2004"]<-6
IM_ind_Ki$time[IM_ind_Ki$År=="2008"]<-11
IM_ind_Ki$time[IM_ind_Ki$År=="2013"]<-16

summary(Ki.mod1<-lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ki))
anova(Ki.mod1)
fixed.effects(Ki.mod1)

plot(Ki.mod1, col = as.numeric(factor(IM_ind_Ki$time, levels = unique(IM_ind_Ki$time))), pch = 16, main = "Kindla, lichen sensitvity")
ACF(Ki.mod1, maxLag = 11)
plot(ACF(Ki.mod1, alpha= .05), main="Kindla, sensitivity")

#summary(Ki.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ki))
#plot(Ki.mod2, col = as.numeric(factor(IM_ind_Ki$time, levels = unique(IM_ind_Ki$time))), pch = 16, main = "Kindla")



######## Subset med bara Gammtratten ####

IM_ind_Ga<-subset(IM_ind, Omr=="Ga")
IM_ind_Ga$År<-factor(IM_ind_Ga$År)
levels(IM_ind_Ga$År)

min(IM_ind_Ga$Sensitivity)
which.min(IM_ind_Ga$Sensitivity)


# Everything on the same plot
ggplot(IM_ind_Ga, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Ga, aes(År,Sensitivity, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(Sensitivity~År,data=IM_ind_Ga, main="Gammtraten", xlab="År", ylab="Sensitvity") 



## Claudias regressioner 

IM_ind_Ga$time<-1
IM_ind_Ga$time[IM_ind_Ga$År=="2005"]<-6
IM_ind_Ga$time[IM_ind_Ga$År=="2010"]<-11
IM_ind_Ga$time[IM_ind_Ga$År=="2015"]<-16

summary(Ga.mod1<-lme(Sensitivity~time,
                     random=~1|Plot/Tree,
                     correlation=corCAR1(form=~time|Plot/Tree),
                     data=IM_ind_Ga))
anova(Ga.mod1)
fixed.effects(Ga.mod1)

plot(Ga.mod1, col = as.numeric(factor(IM_ind_Ga$time, levels = unique(IM_ind_Ga$time))), pch = 16, main = "Gammtratten, lichen sensitvity")
ACF(Ga.mod1, maxLag = 11)
plot(ACF(Ga.mod1, alpha= .1), main="title2")

#summary(Ga.mod2<-lme(Sensitivity~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ga))
#plot(Ga.mod2, col = as.numeric(factor(IM_ind_Ga$time, levels = unique(IM_ind_Ga$time))), pch = 16, main = "Gammtratten")



# Alla områden i samma boxplot####
#png(file = "Boxplot_Sensitivity.png", bg = "white", width = 1920, height = 1920, res = 400, pointsize = 12, type = "windows")

par(mfrow = c(2,2))
boxplot(Sensitivity~År,data=IM_ind_Gd, main="Gårdsjön", xlab="År", ylab="Sensitvity") 
boxplot(Sensitivity~År,data=IM_ind_An, main="Aneboda", xlab="År", ylab="Sensitvity") 
boxplot(Sensitivity~År,data=IM_ind_Ki, main="Kindla", xlab="År", ylab="Sensitvity") 
boxplot(Sensitivity~År,data=IM_ind_Ga, main="Gammtraten", xlab="År", ylab="Sensitvity") 





######## pH #################

# Subset med bara Gårdsjön#####
names(IM_ind_Gd)

(a<-which(is.na(IM_ind_Gd$pH)))
which(IM_ind_Gd$pH==0)
row.names(IM_ind_Gd[which(IM_ind_Gd$pH==0),])
rownames(IM_ind_Gd)[a]
# Ett vÅrde har 0 i "pH", men det År ologiskt
# ErsÅtt med vÅrdet pÅ samma trÅd inventeringen innan
IM_ind_Gd[a,3]
IM_ind_Gd[a,3]<-IM_ind_Gd[match("Gd06E 3-PS5",rownames(IM_ind_Gd)),3]
IM_ind_Gd[a,3]



# Everything on the same plot
ggplot(IM_ind_Gd, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gårdsjön") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Gd, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gårdsjön") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(pH~År,data=IM_ind_Gd, main="Gårdsjön", xlab="År", ylab="pH") 

plot(IM_ind_Gd$time, IM_ind_Gd$pH)


## Claudias regressioner 

IM_ind_Gd$time<-1
IM_ind_Gd$time[IM_ind_Gd$År=="2001"]<-6
IM_ind_Gd$time[IM_ind_Gd$År=="2006"]<-11
IM_ind_Gd$time[IM_ind_Gd$År=="2011"]<-16

summary(Gd_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Gd[-55,]))
anova(Gd_pH.mod1)
fixed.effects(Gd_pH.mod1)

plot(Gd_pH.mod1, col = c(1:nlevels(as.factor(IM_ind_Gd$time))), pch = 16, main="Gårdsjön, pH")

#summary(Gd_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Gd))
#plot(Gd_pH.mod2, col = as.numeric(factor(IM_ind_Gd$time, levels = unique(IM_ind_Gd$time))), pch = 16)





### Test utan extremvÅden; Gårdsjön
# GÅller endast dÅ alla trÅdslag År med

#which(IM_ind_Gd$pH > 3) # Hitta rader med pH > 3
#IM_ind_Gd_1<-IM_ind_Gd[-c(20,40,60),] # Ta bort rader med pH >3

# Ta bort eventuellt tomma rader och kolumner
#dim(IM_ind_Gd_1)
#IM_ind_Gd_1<-IM_ind_Gd_1[ , which(!apply(IM_ind_Gd_1==0,2,all))]# Ta bort tomma kolumner
#IM_ind_Gd_1<-IM_ind_Gd_1[!rowSums(IM_ind_Gd_1, na.rm=TRUE) == 0,] # Ta bort tomma rader
#dim(IM_ind_Gd_1)


#IM_ind_Gd_1$time<-1
#IM_ind_Gd_1$time[IM_ind_Gd_1$År=="2001"]<-6
#IM_ind_Gd_1$time[IM_ind_Gd_1$År=="2006"]<-11
#IM_ind_Gd_1$time[IM_ind_Gd_1$År=="2011"]<-16

#summary(Gd_pH.mod1_1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Gd_1[-55,]))
#anova(Gd_pH.mod1_1)
#fixed.effects(Gd_pH.mod1_1)

#plot(Gd_pH.mod1_1, col = c(1:nlevels(as.factor(IM_ind_Gd$time))), pch = 16, main="Gårdsjön, utan avvikare")



######## Subset med bara Aneboda ####

# Everything on the same plot
ggplot(IM_ind_An, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_An, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(pH~År,data=IM_ind_An, main="Aneboda", xlab="År", ylab="pH") 



## Claudias regressioner

IM_ind_An$time<-1
IM_ind_An$time[IM_ind_An$År=="2002"]<-6
IM_ind_An$time[IM_ind_An$År=="2007"]<-11
IM_ind_An$time[IM_ind_An$År=="2012"]<-16

summary(An_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An))
anova(An_pH.mod1)
fixed.effects(An_pH.mod1)

plot(An_pH.mod1, col = as.numeric(factor(IM_ind_An$time, levels = unique(IM_ind_An$time))), pch = 16, main="Aneboda, pH")

#summary(An_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An))
#plot(An_pH.mod2, col = as.numeric(factor(IM_ind_An$time, levels = unique(IM_ind_An$time))), pch = 16)


### Test utan extremvÅden; Aneboda
# Gäller endast då alla trädslag är med

which(IM_ind_An$pH > 3) # Hitta rader med pH > 3
IM_ind_An_1<-IM_ind_An[-50,] # Ta bort rader med pH >3
names(IM_ind_An_1)

# Ta bort eventuellt tomma rader och kolumner
dim(IM_ind_An_1)
IM_ind_An_1<-IM_ind_An_1[ , which(!apply(IM_ind_An_1==0,2,all))]# Ta bort tomma kolumner
IM_ind_An_1<-IM_ind_An_1[!rowSums(IM_ind_An_1, na.rm=TRUE) == 0,] # Ta bort tomma rader
dim(IM_ind_An_1)

summary(An_pH.mod1_1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An_1))
anova(An_pH.mod1_1)
fixed.effects(An_pH.mod1_1)

plot(An_pH.mod1_1, col = c(1:nlevels(as.factor(IM_ind_An$time))), pch = 16, main="Aneboda, pH, utan avvikare")

summary(An_pH.mod2_1<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An_1))
plot(An_pH.mod2_1, col = as.numeric(factor(IM_ind_An_1$time, levels = unique(IM_ind_An_1$time))), pch = 16)




######## Subset med bara Kindla ####

# Everything on the same plot
ggplot(IM_ind_Ki, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Ki, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(pH~År,data=IM_ind_Ki, main="Kindla", xlab="År", ylab="pH") 



## Claudias regressioner

IM_ind_Ki$time<-1
IM_ind_Ki$time[IM_ind_Ki$År=="2004"]<-6
IM_ind_Ki$time[IM_ind_Ki$År=="2008"]<-11
IM_ind_Ki$time[IM_ind_Ki$År=="2013"]<-16

summary(Ki_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ki))
anova(Ki_pH.mod1)
fixed.effects(Ki_pH.mod1)

plot(Ki_pH.mod1, col = as.numeric(factor(IM_ind_Ki$time, levels = unique(IM_ind_Ki$time))), pch = 16, main = "Kindla, pH")

#summary(Ki_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ki))
#plot(Ki_pH.mod2, col = as.numeric(factor(IM_ind_Ki$time, levels = unique(IM_ind_Ki$time))), pch = 16, main = "Kindla")



######## Subset med bara Gammtratten ####


# Everything on the same plot
ggplot(IM_ind_Ga, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Ga, aes(År,pH, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(pH~År,data=IM_ind_Ga, main="Gammtratten", xlab="År", ylab="pH") 



## Claudias regressioner

IM_ind_Ga$time<-1
IM_ind_Ga$time[IM_ind_Ga$År=="2005"]<-6
IM_ind_Ga$time[IM_ind_Ga$År=="2010"]<-11
IM_ind_Ga$time[IM_ind_Ga$År=="2015"]<-16

summary(Ga_pH.mod1<-lme(pH~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ga))
anova(Ga_pH.mod1)
fixed.effects(Ga_pH.mod1)

plot(Ga_pH.mod1, col = as.numeric(factor(IM_ind_Ga$time, levels = unique(IM_ind_Ga$time))), pch = 16, main = "Gammtratten, pH")
ACF(Ga_pH.mod1, maxLag = 11)
plot(ACF(Ga_pH.mod1, alpha= .0001), main="title2")

acf(residuals(Ga_pH.mod1,type="normalized"))


#summary(Ga_pH.mod2<-lme(pH~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ga))
#plot(Ga_pH.mod2, col = as.numeric(factor(IM_ind_Ga$time, levels = unique(IM_ind_Ga$time))), pch = 16, main = "Gammtratten pH")


# Alla områden i samma boxplot####
par(mfrow = c(2, 2))
boxplot(pH~År,data=IM_ind_Gd, main="Gårdsjön, pH", xlab="År", ylab="pH") #ylim = c(1, 4)
boxplot(pH~År,data=IM_ind_An, main="Aneboda, pH", xlab="År", ylab="pH") 
boxplot(pH~År,data=IM_ind_Ki, main="Kindla, pH", xlab="År", ylab="pH") 
boxplot(pH~År,data=IM_ind_Ga, main="Gammtratten, pH", xlab="År", ylab="pH") 
par(mfrow = c(1, 1))




######## Kväve #################





######## Subset med bara Gårdsjön ####
IM_ind_Gd<-subset(IM_ind, Omr=="Gd")
IM_ind_Gd$År <- as.numeric(as.character(IM_ind_Gd$År))
IM_ind_Gd$År<-factor(IM_ind_Gd$År, labels=c("96", "01", "06", "11"))
#names(IM_ind_Gd)
#str(IM_ind_Gd)

#min(IM_ind_Gd$Nitrogen)
#which.min(IM_ind_Gd$Nitrogen)



# Everything on the same plot
ggplot(IM_ind_Gd, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gårdsjön") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Gd, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gårdsjön") +
  labs(col = "Tree sp") +
  xlab("Year")



## Claudias regressioner

IM_ind_Gd$time<-1
IM_ind_Gd$time[IM_ind_Gd$År=="01"]<-6
IM_ind_Gd$time[IM_ind_Gd$År=="06"]<-11
IM_ind_Gd$time[IM_ind_Gd$År=="11"]<-16

summary(Gd.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=drop_na(IM_ind_Gd)))
plot(Gd.mod1, col = c(1:nlevels(as.factor(IM_ind_Gd$time))), pch = 16)
legend(x=2, y=2, legend=c(levels(IM_ind_Gd$År)), col = c(1:4),cex= 0.8, pch = 16)

summary(Gd.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=drop_na(IM_ind_Gd)))
plot(Gd.mod2, col = as.numeric(factor(IM_ind_Gd$time, levels = unique(IM_ind_Gd$time))), pch = 16)

######## Subset med bara Aneboda ####

IM_ind_An<-subset(IM_ind, Omr=="An")
IM_ind_An$År<-factor(IM_ind_An$År)
levels(IM_ind_An$År)

min(IM_ind_An$Nitrogen)
#which.min(IM_ind_An$Nitrogen)


# Everything on the same plot
ggplot(IM_ind_An, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_An, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")



## Claudias regressioner

IM_ind_An$time<-1
IM_ind_An$time[IM_ind_An$År=="02"]<-6
IM_ind_An$time[IM_ind_An$År=="07"]<-11
IM_ind_An$time[IM_ind_An$År=="12"]<-16

summary(An.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An))
plot(An.mod1, col = as.numeric(factor(IM_ind_An$time, levels = unique(IM_ind_An$time))), pch = 16)

summary(An.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_An))
plot(An.mod2, col = as.numeric(factor(IM_ind_An$time, levels = unique(IM_ind_An$time))), pch = 16)


######## Subset med bara Kindla ####

IM_ind_Ki<-subset(IM_ind, Omr=="Ki")
IM_ind_Ki$År<-factor(IM_ind_Ki$År)
levels(IM_ind_Ki$År)

IM_ind_Ki$Plot<-as.factor(IM_ind_Ki$Plot)

min(IM_ind_Ki$Nitrogen)
#which.min(IM_ind_Ki$Nitrogen)


# Everything on the same plot
ggplot(IM_ind_Ki, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Ki, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")


## Claudias regressioner

IM_ind_Ki$time<-1
IM_ind_Ki$time[IM_ind_Ki$År=="04"]<-6
IM_ind_Ki$time[IM_ind_Ki$År=="08"]<-11
IM_ind_Ki$time[IM_ind_Ki$År=="13"]<-16

summary(Ki.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ki))
plot(Ki.mod1, col = as.numeric(factor(IM_ind_Ki$time, levels = unique(IM_ind_Ki$time))), pch = 16, main = "Kindla")

summary(Ki.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ki))
plot(Ki.mod2, col = as.numeric(factor(IM_ind_Ki$time, levels = unique(IM_ind_Ki$time))), pch = 16, main = "Kindla")



######## Subset med bara Gammtratten ####

IM_ind_Ga<-subset(IM_ind, Omr=="Ga")
IM_ind_Ga$År<-factor(IM_ind_Ga$År)
#levels(IM_ind_Ga$År)
#head(IM_ind_Ga)

#min(IM_ind_Ga$pH)


# Everything on the same plot
ggplot(IM_ind_Ga, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(IM_ind_Ga, aes(År,Nitrogen, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")



## Claudias regressioner

IM_ind_Ga$time<-1
IM_ind_Ga$time[IM_ind_Ga$År=="05"]<-6
IM_ind_Ga$time[IM_ind_Ga$År=="10"]<-11
IM_ind_Ga$time[IM_ind_Ga$År=="15"]<-16

summary(Ga.mod1<-lme(Nitrogen~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ga))
plot(Ga.mod1, col = as.numeric(factor(IM_ind_Ga$time, levels = unique(IM_ind_Ga$time))), pch = 16, main = "Gammtratten pH")

summary(Ga.mod2<-lme(Nitrogen~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=IM_ind_Ga))
plot(Ga.mod2, col = as.numeric(factor(IM_ind_Ga$time, levels = unique(IM_ind_Ga$time))), pch = 16, main = "Gammtratten pH")


# Alla områden i samma boxplot####
par(mfrow = c(2, 2))
boxplot(n_sp ~År,data=IM_ind_Gd, main="Gårdsjön, Nitrogen", xlab="År", ylab="Nitrogen") #ylim = c(1, 4)
boxplot(Nitrogen ~År,data=IM_ind_An, main="Aneboda, Nitrogen", xlab="År", ylab="Nitrogen") 
boxplot(Nitrogen ~År,data=IM_ind_Ki, main="Kindla, Nitrogen", xlab="År", ylab="Nitrogen") 
boxplot(Nitrogen ~År,data=IM_ind_Ga, main="Gammtratten, Nitrogen", xlab="År", ylab="Nitrogen") 
par(mfrow = c(1, 1))

