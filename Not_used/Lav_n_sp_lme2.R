# Scritp f?r att analysera lavdata fr?n IM

library(ggplot2)
library(nlme)


#load data
spe <- read_csv("Raw_data/Sp_all.csv")
env <- read_csv("Raw_data/Env_all.csv")

env <- as_tibble(env)
spe <- as_tibble(spe)
#env$X1 <- as.factor(env$X1)
#spe$X1 <- as.factor(spe$X1)

#visualise data. Requires VIM
matrixplot(spe, interactive = F, sortby = "X1")

#looks like missing data for one species `Sphaerophorus globosus`
apply(is.na(spe), 2, which)

#remove it
spe$`Sphaerophorus globosus` <- NULL

spe <- rename(spe, ID = X1)

##### Lägg till data om individuella träd
spe$Site<-substr(spe$ID,1,2)
spe$Year<-substr(spe$ID,3,4)
spe$Plot<-substr(spe$ID,5,7)
spe$Tree<-substr(spe$ID,9,nchar(as.character(spe$ID)))
spe$sp<-substr(spe$Tree,1,2)
spe$id2<-substr(spe$ID,5,nchar(spe$ID))
head(spe$id2)
spe <- select(spe, Site, Year, Plot, Tree, sp, id2, everything())

# Gör om yy till yyyy
spe$Year <- ifelse(as.numeric(as.character(spe$Year)) > 17,paste("19",spe$Year, sep=""), paste("20",spe$Year, sep="")) 
spe$Year <- as.factor(spe$Year)
spe$Site <- as.factor(spe$Site)
spe$Plot <- as.factor(spe$Plot)
spe$Tree <- as.factor(spe$Tree)
spe$id2 <- as.factor(spe$id2)
spe$ID <- as.factor(spe$ID)

#species richness per sample
n_sp <- rowSums(spe[, 8:62]>0)
spe <- cbind(spe, n_sp)

#shannon div per sample
Shan <- diversity(spe[,-c(1:7)],"shannon")
spe <- cbind(spe, Shan)

# 
# #tidy format
# spe_tidy <- spe %>%
#   pivot_longer(-ID, names_to = "spe", values_to = "cover")
# spe_tidy <- spe_tidy %>% filter(cover>0)
# spe_tidy <- spe_tidy %>% group_by(ID) %>% mutate(n_sp = length(unique(spe)))
# spe_tidy <- spe_tidy %>% mutate(tot_cov = sum(cover))
# 
# 



##########################################
############# Ta bort asp, s?lg och ek
#spe<-subset(spe, spe$TreeSp == "PINU SYL" | spe$TreeSp == "PICE ABI" | spe$TreeSp == "BETULA Z")
#env<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
#env$TreeSp<-factor(env$TreeSp)
keep <- c("PS","PA","BZ","PINU SYL","PICE ABI","BETULA Z")
spe <- filter(spe, sp %in% keep)
env <- filter(env, TreeSp %in% keep)

###########################################

# Plotta alla Year och Site?den
ggplot(spe, aes(x=Year, y=n_sp, colour = Site)) + 
  geom_point() +
  geom_smooth(method=lm, se=F)

# Frekvenstabell f?r alla v?rdtr?d
ant_träd<-as.data.frame(table(spe$sp, spe$Year, by = spe$Site))
write.table(ant_tr?d, file = "antal_tr?d.txt", sep = "\t", quote = FALSE)


##################################################

######################################
######## Regressioner per Site?de #################
######################################

######## Subset med bara G?rdsj?n
spe_Gd<-subset(spe, Site=="Gd")
spe_Gd$Year<-factor(spe_Gd$Year)
names(spe_Gd)
head(spe_Gd)

#min(spe_Gd$n_sp)
#(a<-which(is.na(spe_Gd$n_sp)))
#which.min(spe_Gd$n_sp)
#row.names(spe_Gd[which(is.na(spe_Gd$n_sp)),])


# Everything on the same plot
ggplot(spe_Gd, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("G?rdsj?n") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_Gd, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("G?rdsj?n") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~Year,data=spe_Gd, main="G?rdsj?n", xlab="Year", ylab="Number of species") 




##########################
## Claudias regressioner

spe_Gd$time<-1
spe_Gd$time[spe_Gd$Year=="2001"]<-6
spe_Gd$time[spe_Gd$Year=="2006"]<-11
spe_Gd$time[spe_Gd$Year=="2011"]<-16

summary(Gd.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
anova(Gd.mod1)
fixed.effects(Gd.mod1)

plot(Gd.mod1, col = c(1:nlevels(as.factor(spe_Gd$time))), pch = 16, main="G?rdsj?n, number of species")

#summary(Gd.mod2<-lme(n_sp~as.factor(time),random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
#plot(Gd.mod2, col = as.numeric(factor(spe_Gd$time, levels = unique(spe_Gd$time))), pch = 16)

summary(Gd.mod3<-lme(n_sp~time+sp,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
anova(Gd.mod3)
fixed.effects(Gd.mod3)

plot(Gd.mod3, col = c(1:nlevels(as.factor(spe_Gd$time))), pch = 16, main="G?rdsj?n, number of species; mod3")


## Summat?ckning
summary(Gd.mod3<-lme(t?ck~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=spe_Gd))
anova(Gd.mod3)
fixed.effects(Gd.mod3)

plot(Gd.mod3, col = c(1:nlevels(as.factor(spe_Gd$time))), pch = 16, main="G?rdsj?n, sum cover")

# Boxplot
boxplot(t?ck~Year,data=spe_Gd, main="G?rdsj?n", xlab="Year", ylab="Sum cover") 



## Test
#t<-spe_Gd
#t$n_sp[t$Year=="2001"]<-t$n_sp+2
#t$n_sp[t$Year=="2006"]<-t$n_sp+5
#t$n_sp[t$Year=="2011"]<-t$n_sp+10
#summary(t_GM.mod1<-lme(n_sp~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=t))
#anova(t_GM.mod1)
#fixed.effects(t_GM.mod1)
#plot(t_GM.mod1, col = as.numeric(factor(t$time, levels = unique(t$time))), pch=16)
#ACF(t_GM.mod1, maxLag = 3)
#plot(ACF(t_GM.mod1, alpha= .05), main="Testdata")
#acf(residuals(t_GM.mod1,type="normalized"))


######## Subset med bara Aneboda

spe_An<-subset(spe, Site=="An")
spe_An$Year<-factor(spe_An$Year)
levels(spe_An$Year)

min(spe_An$n_sp)
#which.min(spe_An$n_sp)


# Everything on the same plot
ggplot(spe_An, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_An, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Aneboda") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~Year,data=spe_An, main="Aneboda",xlab="Year", ylab="Number of species") 


##########################
## Claudias regressioner

spe_An$time<-1
spe_An$time[spe_An$Year=="2002"]<-6
spe_An$time[spe_An$Year=="2007"]<-11
spe_An$time[spe_An$Year=="2012"]<-16

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

spe_Ki<-subset(spe, Site=="Ki")
spe_Ki$Year<-factor(spe_Ki$Year)
levels(spe_Ki$Year)

spe_Ki$Plot<-as.factor(spe_Ki$Plot)
levels(spe_Ki$Plot)

min(spe_Ki$n_sp)
#which.min(spe_Ki$n_sp)


# Everything on the same plot
ggplot(spe_Ki, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_Ki, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Kindla") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~Year,data=spe_Ki, main="Kindla", xlab="Year", ylab="Kindla, number of species") 


##########################
## Claudias regressioner

spe_Ki$time<-1
spe_Ki$time[spe_Ki$Year=="2004"]<-6
spe_Ki$time[spe_Ki$Year=="2008"]<-11
spe_Ki$time[spe_Ki$Year=="2013"]<-16

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

spe_Ga<-subset(spe, Site=="Ga")
spe_Ga$Year<-factor(spe_Ga$Year)
levels(spe_Ga$Year)

min(spe_Ga$n_sp)
#which.min(spe_Ga$n_sp)


# Everything on the same plot
ggplot(spe_Ga, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_line(aes(col=substr(Tree,1,2))) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

ggplot(spe_Ga, aes(Year,n_sp, col=substr(Tree,1,2), group=id2)) + 
  geom_point() + 
  geom_smooth(method=lm, se=F) +
  ggtitle("Gammtratten") +
  labs(col = "Tree sp") +
  xlab("Year")

# Boxplot
boxplot(n_sp~Year,data=spe_Ga, main="Gammtraten", xlab="Year", ylab="Number of species") 


##########################
## Claudias regressioner

spe_Ga$time<-1
spe_Ga$time[spe_Ga$Year=="2005"]<-6
spe_Ga$time[spe_Ga$Year=="2010"]<-11
spe_Ga$time[spe_Ga$Year=="2015"]<-16

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



# Alla Site?den i samma boxplot
#png(file = "Boxplot_n_sp.png", bg = "white", width = 1920, height = 1920, res = 400, pointsize = 12, type = "windows")

par(mfrow = c(2,2))
boxplot(n_sp~Year,data=spe_Gd, main="Gårdsjön", xlab="Year", ylab="Number of species") 
boxplot(n_sp~Year,data=spe_An, main="Aneboda", xlab="Year", ylab="Number of species") 
boxplot(n_sp~Year,data=spe_Ki, main="Kindla", xlab="Year", ylab="Number of species") 
boxplot(n_sp~Year,data=spe_Ga, main="Gammtraten", xlab="Year", ylab="Number of species") 

#dev.off(which = dev.cur())
#dev.off(which = dev.cur())
#dev.off(which = dev.cur())
par(mfrow = c(1,1))


