# Scritp fÅr att analysera lavdata fr?n IM

library(ggplot2)
library(nlme)

# Read and transform data
#spe<-read.table("Lav_spe.txt", header = T, row.names = 1, sep = "\t")
#spe<-spe/400

# Nya artdata med ?verfÅring av Bryoria capillaris/fuscescens till Bryoria fuscescens
Lav_spe <- read.delim("~/Documents/R/paper_3/Raw_data/Lav_spe.txt", row.names=1)
spe<-Lav_spe
spe<-spe/400
head(spe)
##### Ber?kna antal arter per tr?d
spe$n_sp<-apply(spe>0,1,sum)


##### Ber?kna summa t?ckning per tr?d
spe$täck<-apply(spe[,c(1:55)],1,sum)
head(spe)



##### L?gg till data om individuella tr?d

spe$id<-row.names(spe)
spe$Omr<-substr(row.names(spe),1,2)
spe$År<-as.factor(substr(row.names(spe),3,4))
spe$Plot<-substr(row.names(spe),5,7)
spe$Tree<-substr(row.names(spe),9,nchar(row.names(spe)))
spe$sp<-substr(spe$Tree,1,2)
spe$id2<-substr(row.names(spe),5,nchar(row.names(spe)))
head(spe$id2)

# GÅr om yy till yyyy
spe$År <- ifelse(as.numeric(as.character(spe$År)) > 17,paste("19",spe$År, sep=""), paste("20",spe$År, sep="")) 
spe$År <- as.factor(spe$År)

head(spe[,57:ncol(spe)])

#write.table(spe, file ="spe_individ.txt", quote = F, sep = "\t", na="NA",dec=",", row.names=TRUE, col.names = TRUE)



##########################################
############# Ta bort asp, s?lg och ek
spe<-subset(spe, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env$TreeSp<-factor(env$TreeSp)
###########################################

# Plotta alla År och omr?den
ggplot(spe, aes(x=År, y=n_sp, color=Omr, group= id2)) + 
  geom_point() +
  geom_smooth(method=lm, se=F)

# Frekvenstabell fÅr alla vÅrdtr?d
ant_träd<-as.data.frame(table(spe$sp, spe$År, by = spe$Omr))
#write.table(ant_tr?d, file = "antal_tr?d.txt", sep = "\t", quote = FALSE)


##################################################

######################################
######## Regressioner per omr?de #################
######################################

######## Subset med bara GÅrdsj?n
spe_Gd<-subset(spe, Omr=="Gd")
spe_Gd$År<-factor(spe_Gd$År)
names(spe_Gd)
head(spe_Gd)

#min(spe_Gd$n_sp)
#(a<-which(is.na(spe_Gd$n_sp)))
#which.min(spe_Gd$n_sp)
#row.names(spe_Gd[which(is.na(spe_Gd$n_sp)),])


# Everything on the same plot
ggplot(spe_Gd, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("GÅrdsj?n") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_Gd, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("GÅrdsj?n") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~År,data=spe_Gd, main="GÅrdsj?n", xlab="År", ylab="Number of species") 




##########################
## Claudias regressioner

spe_Gd$time<-1
spe_Gd$time[spe_Gd$År=="2001"]<-6
spe_Gd$time[spe_Gd$År=="2006"]<-11
spe_Gd$time[spe_Gd$År=="2011"]<-16

summary(Gd.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
anova(Gd.mod1)
fixed.effects(Gd.mod1)

plot(Gd.mod1, col = c(1:nlevels(as.factor(spe_Gd$time))), pch = 16, main="GÅrdsj?n, number of species")

#summary(Gd.mod2<-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
#plot(Gd.mod2, col = as.numeric(factor(spe_Gd$time, levels = unique(spe_Gd$time))), pch = 16)

summary(Gd.mod3<-lme(n_sp~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
anova(Gd.mod3)
fixed.effects(Gd.mod3)

plot(Gd.mod3, col = c(1:nlevels(as.factor(spe_Gd$time))), pch = 16, main="GÅrdsj?n, number of species; mod3")


## Summat?ckning
summary(Gd.mod3<-lme(t?ck~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
anova(Gd.mod3)
fixed.effects(Gd.mod3)

plot(Gd.mod3, col = c(1:nlevels(as.factor(spe_Gd$time))), pch = 16, main="GÅrdsj?n, sum cover")

# Boxplot
boxplot(täck~År,data=spe_Gd, main="GÅrdsj?n", xlab="År", ylab="Sum cover") 



## Test
#t<-spe_Gd
#t$n_sp[t$År=="2001"]<-t$n_sp+2
#t$n_sp[t$År=="2006"]<-t$n_sp+5
#t$n_sp[t$År=="2011"]<-t$n_sp+10
#summary(t_GM.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=t))
#anova(t_GM.mod1)
#fixed.effects(t_GM.mod1)
#plot(t_GM.mod1, col = as.numeric(factor(t$time, levels = unique(t$time))), pch=16)
#ACF(t_GM.mod1, maxLag = 3)
#plot(ACF(t_GM.mod1, alpha= .05), main="Testdata")
#acf(residuals(t_GM.mod1,type="normalized"))


######## Subset med bara Aneboda

spe_An<-subset(spe, Omr=="An")
spe_An$År<-factor(spe_An$År)
levels(spe_An$År)

min(spe_An$n_sp)
#which.min(spe_An$n_sp)


# Everything on the same plot
ggplot(spe_An, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_An, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~År,data=spe_An, main="Aneboda",xlab="År", ylab="Number of species") 


##########################
## Claudias regressioner

spe_An$time<-1
spe_An$time[spe_An$År=="2002"]<-6
spe_An$time[spe_An$År=="2007"]<-11
spe_An$time[spe_An$År=="2012"]<-16

summary(An.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_An))
anova(An.mod1)
fixed.effects(An.mod1)

plot(An.mod1, col = as.numeric(factor(spe_An$time, levels = unique(spe_An$time))), pch = 16, main="Aneboda, number of species")

#summary(An.mod2<-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_An))
#plot(An.mod2, col = as.numeric(factor(spe_An$time, levels = unique(spe_An$time))), pch = 16)

summary(An.mod3<-lme(n_sp~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_An))
anova(An.mod3)
fixed.effects(An.mod3)
plot(An.mod3, col = c(1:nlevels(as.factor(spe_An$time))), pch = 16, main="Aneboda, number of species; mod3")



######## Subseet med bara Kindla

spe_Ki<-subset(spe, Omr=="Ki")
spe_Ki$År<-factor(spe_Ki$År)
levels(spe_Ki$År)

spe_Ki$Plot<-as.factor(spe_Ki$Plot)
levels(spe_Ki$Plot)

min(spe_Ki$n_sp)
#which.min(spe_Ki$n_sp)


# Everything on the same plot
ggplot(spe_Ki, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_Ki, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~År,data=spe_Ki, main="Kindla", xlab="År", ylab="Kindla, number of species") 


##########################
## Claudias regressioner

spe_Ki$time<-1
spe_Ki$time[spe_Ki$År=="2004"]<-6
spe_Ki$time[spe_Ki$År=="2008"]<-11
spe_Ki$time[spe_Ki$År=="2013"]<-16

summary(Ki.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Ki))
anova(Ki.mod1)
fixed.effects(Ki.mod1)

plot(Ki.mod1, col = as.numeric(factor(spe_Ki$time, levels = unique(spe_Ki$time))), pch = 16, main = "Kindla, number of species")
ACF(Ki.mod1, maxLag = 3)
plot(ACF(Ki.mod1, alpha= .05), main="Kindla, n_sp")

#summary(Ki.mod2<-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Ki))
#plot(Ki.mod2, col = as.numeric(factor(spe_Ki$time, levels = unique(spe_Ki$time))), pch = 16, main = "Kindla")

summary(Ki.mod3<-lme(n_sp~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Ki))
anova(Ki.mod3)
fixed.effects(Ki.mod3)
plot(Ki.mod3, col = c(1:nlevels(as.factor(spe_Ki$time))), pch = 16, main="Kindla, number of species; mod3")


######## Subset med bara Gammtratten

spe_Ga<-subset(spe, Omr=="Ga")
spe_Ga$År<-factor(spe_Ga$År)
levels(spe_Ga$År)

min(spe_Ga$n_sp)
#which.min(spe_Ga$n_sp)


# Everything on the same plot
ggplot(spe_Ga, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_Ga, aes(År,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~År,data=spe_Ga, main="Gammtraten", xlab="Year", ylab="Number of species") 


##########################
## Claudias regressioner

spe_Ga$time<-1
spe_Ga$time[spe_Ga$År=="2005"]<-6
spe_Ga$time[spe_Ga$År=="2010"]<-11
spe_Ga$time[spe_Ga$År=="2015"]<-16

summary(Ga.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Ga))
anova(Ga.mod1)
fixed.effects(Ga.mod1)

plot(Ga.mod1, col = as.numeric(factor(spe_Ga$time, levels = unique(spe_Ga$time))), pch = 16, main = "Gammtratten, number of species")
ACF(Ga.mod1, maxLag = 3)
plot(ACF(Ga.mod1, alpha= .1), main="title2")

#summary(Ga.mod2<-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Ga))
#plot(Ga.mod2, col = as.numeric(factor(spe_Ga$time, levels = unique(spe_Ga$time))), pch = 16, main = "Gammtratten")

summary(Ga.mod3<-lme(n_sp~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Ga))
anova(Ga.mod3)
fixed.effects(Ga.mod3)
plot(Ga.mod3, col = c(1:nlevels(as.factor(spe_Ga$time))), pch = 16, main="Gammtraten, number of species; mod3")



# Alla omr?den i samma boxplot
#png(file = "Boxplot_n_sp.png", bg = "white", width = 1920, height = 1920, res = 400, pointsize = 12, type = "windows")

par(mfrow = c(2,2))
boxplot(n_sp~År,data=spe_Gd, main="GÅrdsj?n", xlab="Year", ylab="Number of species") 
boxplot(n_sp~År,data=spe_An, main="Aneboda", xlab="Year", ylab="Number of species") 
boxplot(n_sp~År,data=spe_Ki, main="Kindla", xlab="Year", ylab="Number of species") 
boxplot(n_sp~År,data=spe_Ga, main="Gammtraten", xlab="Year", ylab="Number of species") 

#dev.off(which = dev.cur())
#dev.off(which = dev.cur())
#dev.off(which = dev.cur())
par(mfrow = c(1,1))


