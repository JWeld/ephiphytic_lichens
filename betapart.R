#BETAPART STUFF###############################################################
#Beta diversity (variation of the species composition of assemblages) may reflect two 
#different phenomena. There are two potential ways in which two species assemblages 
#can be ‘different’. One is species replacement (i.e. turnover), which consists in the 
#substitution of species in one site by different species in the other site. 
#The second way is species loss (or gain), which implies the elimination (or addition) 
#of species in only one of the sites, and leads to the poorest assemblage being a strict 
#subset of the richest one (a pattern called nestedness).
library(betapart)
#Use sign function to convert abundance data 
#into prescence abscence matrix that can be fed to Betapart. Can't include factors (non numeric)
dat <- select(spe, -n_sp)
dat.pa=sign(dat) #exclude them. works but lose first factor columns
dat.pa=cbind(Site=0, dat.pa)
dat.pa=cbind(Tree=0, dat.pa) 
dat.pa=cbind(Plot=0, dat.pa) #add back some empty columns in that case...
dat.pa=cbind(Year=0, dat.pa)
dat.pa=cbind(ID=0, dat.pa) 
dat.pa$Year <- data2$Year #copy the values from object previously produced in epiphytes.R
dat.pa$Site=data2$Site
dat.pa$Plot=data2$Plot
dat.pa$Tree=data2$Tree
dat.pa$ID=data2$ID#to put them back in.
str(dat.pa)

#Betapart doesn't like year and plot columns. Subset by year...
#fieldpresabs_1 created above
period1 <- c(1996,1997,1998,2000)
period2 <- c(2001,2002,2004,2005)
period3 <- c(2006,2007,2008,2010)
period4 <- c(2011,2012,2013,2015)

dat.pa.1 <- dat.pa %>% filter(Year %in% period1)  
dat.pa.2 <- dat.pa %>% filter(Year %in% period2)  
dat.pa.3 <- dat.pa %>% filter(Year %in% period3)   
dat.pa.4 <- dat.pa %>% filter(Year %in% period4)  

#break down by site
Gd1 <- filter(dat.pa.1, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Gd2 <- filter(dat.pa.2, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID) 
Gd3 <- filter(dat.pa.3, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Gd4 <- filter(dat.pa.4, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID)

An1 <- filter(dat.pa.1, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)
An2 <- filter(dat.pa.2, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)
An3 <- filter(dat.pa.3, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)
An4 <- filter(dat.pa.4, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)

Ki1 <- filter(dat.pa.1, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ki2 <- filter(dat.pa.2, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ki3 <- filter(dat.pa.3, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ki4 <- filter(dat.pa.4, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)

Ga1 <- filter(dat.pa.1, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ga2 <- filter(dat.pa.2, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ga3 <- filter(dat.pa.3, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ga4 <- filter(dat.pa.4, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)

#check things look right 
str(Gd1)

# make betapart objects 
Gd1c <- betapart.core(Gd1) 
Gd2c <- betapart.core(Gd2) 
Gd3c <- betapart.core(Gd3) 
Gd4c <- betapart.core(Gd4) 

An1c <- betapart.core(An1) 
An2c <- betapart.core(An2) 
An3c <- betapart.core(An3) 
An4c <- betapart.core(An4) 

Ki1c <- betapart.core(Ki1) 
Ki2c <- betapart.core(Ki2) 
Ki3c <- betapart.core(Ki3) 
Ki4c <- betapart.core(Ki4) 

Ga1c <- betapart.core(Ga1) 
Ga2c <- betapart.core(Ga2) 
Ga3c <- betapart.core(Ga3) 
Ga4c <- betapart.core(Ga4) 

#Returns three values-turnover and nestedness components, + total multi-plot dissimilarity across the site
(Ga1m <- beta.multi(Gd1c))#turnover=0.7406015, nestedness=0.09151528, betadiv=0.8321168
(Ga2m <- beta.multi(Gd2c))#turnover=0.7419355 nestedness=0.07624633, betadiv=0.8181818
(Ga3m <- beta.multi(Gd3c))#turnover=0.7472527 nestedness=0.06972338, betadiv=0.8169761
(Ga4m <- beta.multi(Gd4c))#turnover=0.6455696 nestedness=0.1439041, betadiv=0.7894737

Ga1p = beta.pair(Ga1c)
Ga2p = beta.pair(Ga2c)
Ga3p = beta.pair(Ga3c)
Ga4p = beta.pair(Ga4c)

#boxplot turnover (pass 1 to function), nestedness (2), or total (3)
fbeta.pair= function(x){
  y=c(Ga1p[x],Ga2p[x], Ga3p[x], Ga4p[x])
  boxplot(y)
  title(main="Beta.pair,all 4 sample years")
}

par(mfrow=c(1,3))
fbeta.pair(1) 

#boxplot turnover (pass 1 to function), nestedness (2), or total (3)
fbeta.multi= function(x){
  y=c(Ga1p[x], Ga2p[x], Ga3p[x], Ga4p[x])
  boxplot(y)
  title(main="Beta.multi, all 4 sample years")
}

fbeta.multi(1)

#temporal change between years
#apply betapart temp (must be same sized matrices)
betapart_t_Gd <-  beta.temp(Gd1c, Gd4c)
betapart_t_An <-  beta.temp(An1c, An4c)
betapart_t_Ki <-  beta.temp(Ki1c, Ki4c)
betapart_t_Ga <-  beta.temp(Ga1c, Ga4c)

par(mfrow=c(1,3))
boxplot(betapart_t_Gd)
title("betapart_t_Gd")
boxplot(betapart_t_An)
title("betapart_t_An")
boxplot(betapart_t_Ki)
title("betapart_t_Ki")
boxplot(betapart_t_Ga)
title("betapart_t_Ga")

#plot probability distributions####
# sampling across equal sites (CHANGE to the beta.core object for the site of interest!)
samp_1 <- beta.sample(Ki1c, sites=10, samples=100, index.family="sor")
samp_2 <- beta.sample(Ki2c, sites=10, samples=100, index.family="sor") 
samp_3 <- beta.sample(Ki3c, sites=10, samples=100, index.family="sor") 
samp_4 <- beta.sample(Ki4c, sites=10, samples=100, index.family="sor") 

# plotting the distributions of components 
dist_1 <- samp_1$sampled.values
dist_2 <- samp_2$sampled.values 
dist_3 <- samp_3$sampled.values 
dist_4 <- samp_4$sampled.values 

#Aneboda
par(mfrow=c(1,1))
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0, 17), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Aneboda" )
# Add a legend
legend(0.6, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.6, 12, legend=c("1997", "2002", "2007", "2012"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

#Gårdsjön
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0, 17), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Gårdsjön" )
# Add a legend
legend(0.7, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.7, 12, legend=c("1996", "2001", "2006", "2011"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

#Gammtratten
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0, 18), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Gammtratten" )
# Add a legend
legend(0.5, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.5, 12, legend=c("2000", "2005", "2010", "2015"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

#Kindla
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.8), ylim=c(0, 17), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Kindla" )
# Add a legend
legend(0.5, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.5, 12, legend=c("1998", "2004", "2008", "2013"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

#Not Used####

# #alternative Jaccard
# par(mfrow=c(1,1))
# plot(density(
#   dist_1$beta.JNE), xlim=c(0,0.9), ylim=c(0, 22), xlab='Beta diversity', main='', col='grey60',lwd=2)
# #lines(density(dist_1$beta.JAC),col='grey60', lwd=2) 
# lines(density(dist_1$beta.JNE),col='grey60', lty=1, lwd=2) 
# lines(density(dist_1$beta.JTU),col='grey60', lty=2, lwd=2)
# #lines(density(dist_2$beta.JAC),col='red', lwd=2) 
# lines(density(dist_2$beta.JNE),col='red', lty=1,lwd=2) 
# lines(density(dist_2$beta.JTU),col='red', lty=2,lwd=2) 
# #lines(density(dist_3$beta.JAC),col='green', lwd=2) 
# lines(density(dist_3$beta.JNE),col='green', lty=1, lwd=2) 
# lines(density(dist_3$beta.JTU),col='green', lty=2, lwd=2)
# #lines(density(dist_4$beta.JAC),col='blue', lwd=2) 
# lines(density(dist_4$beta.JNE),col='blue', lty=1, lwd=2) 
# lines(density(dist_4$beta.JTU),col='blue', lty=2, lwd=2)
# #title("gray=1,r=2, g=3, b=4, right=total, left=nest, dash=turn" )
# title("Gårdsjön")
# # Add a legend
# legend(0.7, 18, legend=c("1997", "2002", "2007", "2012"),
#        col=c("grey60","red","green", "blue"), lty=1, cex=0.8)
# text(0.15, 18, labels = "Turnover")
# text(0.65, 10, labels = "Nestedness")

# # Now, we'll calculate the Jaccard index and its partitions of turnover and nestedness. We can calculate Sorensen index instead by using the argument     index.family="sorensen"    .
# dist<-beta.pair(select(dat.pa.4,-Year, -Plot, -Site, -Tree, -ID), index.family="sorensen")
# dist<-bray.part(dat)
# # To get the pairwise Jaccard index turnover partition between communities, type: dist[[1]]. To get nestedness partition, type: dist[[2]]. To get all beta diversity: dist[[3]].
# groups <- data$Site
# bd<-betadisper(dist[[3]],groups)
# plot(bd)
# boxplot(bd)
# anova(bd)



#exploratory PCA####
# group <- data4$Site
# dat <- select(data4, -c(ID, Site, Year, Plot, Tree, TInd, RInd, pH, id2, ID3, sp))
# 
# fit <- rda(dat, scale = TRUE)
# pl <- ordiplot(fit, type = "none")
# points(pl, "sites", pch=21, col="red", bg="yellow")
# text(pl, "species", col="blue", cex=0.9)
# ordihull(fit, groups = data4$Site, label = TRUE)






