# Scritp för att analysera lavdata från IM

library(vegan)
library(labdsv)
library(plotrix)
library(betapart)
library(vegdata)
library(ggplot2)
library(gridExtra)
library(reshape2)


setwd("C:/Users/grandin/Documents/Ulf IVM/IM/Lavar Utvärd 2016/R")

spe<-read.table("Lav_spe.txt", header = T, row.names = 1, sep = "\t")
spe<-spe/400

env<-read.table("Lav_env.txt", header = T, row.names = 1, sep = "\t")
as.factor(env$År)


# Kombinera Område och år
env$OmrÅr<-paste(substring(env$Område,1,2),env$Inv_nr, sep="_")
env$OmrÅr<-factor(env$OmrÅr)
names(env)
#levels(env$OmrÅr)

# Ta bort asp, sälg och ek
spe1<-subset(spe, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env1<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env1$TreeSp<-factor(env1$TreeSp)
names(env1)

# Skapa subset med bara tall, gran och björk
#oädel<-env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z" 


#############################################################
# Slå samman spe och env, och dela up på de olika områdena
#######################################################################################

spenv<-merge(spe1,env1, all = TRUE, by = 0)
row.names(spenv)<-spenv[[1]] # Fixa radnamn
spenv<-spenv[,-1] # ta bort kolumn med radnamn
levels(spenv$Område)
names(spenv)

############# Gårdsjön #####
Gård<-subset(spenv, subset = spenv$Område == "Gd")
# Check number of rows with rowSums = 0
sum(rowSums(Gård[,c(1:(ncol(Gård)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-Gård[!rowSums(Gård[,c(1:(ncol(Gård)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(Gård[,c(1:(ncol(Gård)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
Gård<-Gård[ , which(!apply(Gård==0,2,all))]
sum(colSums(Gård[,c(1:(ncol(Gård)-5))], na.rm=TRUE) == 0)

# Ändra alla Cladonia till Cladonia.sp. för OmrÅr = Gd_4
# 2011 noterades Cladonia som C.sp, alla tidigare år som Cladonia.coniocraea
# Detta gav ett väldigt utslag för 2011 i NMDS och adonis

#Cladonia.coniocraea	Cladonia.digitata	Cladonia.pyxidata	Cladonia.sp.	Cladonia.squamosa

Gård$Cladonia.sp<-Gård$Cladonia.coniocraea+Gård$Cladonia.digitata+Gård$Cladonia.pyxidata+Gård$Cladonia.sp.+Gård$Cladonia.squamosa
Gård$Cladonia.coniocraea<-NULL
Gård$Cladonia.digitata<-NULL
Gård$Cladonia.pyxidata<-NULL
Gård$Cladonia.sp.<-NULL
Gård$Cladonia.squamosa<-NULL

names(Gård)
Gård<-Gård[,c(1:5,27,6:26)] # Så att artena komme i bokstavsorning (och inte Clad sp sist)
names(Gård)

#write.table(Gård, "C:/Users/Ulf/Ulf IVM/IM/Lavar Utvärd 2016/Gårdsjön/Gård0.txt", quote=FALSE, sep = "\t", dec = ",")#, row.names = TRUE, col.names = TRUE)


############# Aneboda #####
Aneb<-subset(spenv, subset = spenv$Område == "An")
# Check number of rows with rowSums = 0
sum(rowSums(Aneb[,c(1:(ncol(Aneb)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-Aneb[!rowSums(Aneb[,c(1:(ncol(Aneb)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(Aneb[,c(1:(ncol(Aneb)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
Aneb<-Aneb[ , which(!apply(Aneb==0,2,all))]
sum(colSums(Aneb[,c(1:(ncol(Aneb)-5))], na.rm=TRUE) == 0)


# Aneboda utan avvikare

# 2 avvikande ytor: An07I 7-BZ3 och An12J 7-PS4
# Båda ytorna var nya resp. år.
Ane_avv<-c("An07I 7-BZ3","An12J 7-PS4")
match(Ane_avv, rownames(Aneb))
Aneb1<-Aneb[-c(26,59),]

# Check number of rows with rowSums = 0
sum(rowSums(Aneb1[,c(1:(ncol(Aneb1)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-Aneb1[!rowSums(Aneb1[,c(1:(ncol(Aneb1)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(Aneb1[,c(1:(ncol(Aneb1)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
Aneb1<-Aneb1[ , which(!apply(Aneb1==0,2,all))]
sum(colSums(Aneb1[,c(1:(ncol(Aneb1)-5))], na.rm=TRUE) == 0)




############# Kindla #####
Kind<-subset(spenv, subset = spenv$Område == "Ki")
# Check number of rows with rowSums = 0
sum(rowSums(Kind[,c(1:(ncol(Kind)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-Kind[!rowSums(Kind[,c(1:(ncol(Kind)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(Kind[,c(1:(ncol(Kind)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
Kind<-Kind[ , which(!apply(Kind==0,2,all))]
sum(colSums(Kind[,c(1:(ncol(Kind)-5))], na.rm=TRUE) == 0)

############# Gammratten #####
Gam<-subset(spenv, subset = spenv$Område == "Ga")
# Check number of rows with rowSums = 0
sum(rowSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-Gam[!rowSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
Gam<-Gam[ , which(!apply(Gam==0,2,all))]
sum(colSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0)


#Skapa data för beräknigar


########### Aneboda #########################

# Skapa dataset där 0 kodas om till NA
Aneb0<-Aneb
is.na(Aneb0[,c(1:(ncol(Aneb0)-5))])<-!Aneb0[,c(1:(ncol(Aneb0)-5))]
names(Aneb0)

Aneb0m<-melt(Aneb0[,c(1:20,22)], id.vars="År", variable.name="Species", value.name="Cover",na.rm=T)
head(Aneb0m)

# Alternativt Aneboda med avvikare borttagna
Aneb1_0<-Aneb1
is.na(Aneb1_0[,c(1:(ncol(Aneb1_0)-5))])<-!Aneb1_0[,c(1:(ncol(Aneb1_0)-5))]
names(Aneb1_0)

Aneb1_0m<-melt(Aneb1_0[,c(1:18,20)], id.vars="År", variable.name="Species", value.name="Cover",na.rm=T)
head(Aneb1_0m)


########### Gammtratten #########################

# Skapa dataset där 0 kodas om till NA
Gam0<-Gam
is.na(Gam0[,c(1:(ncol(Gam0)-5))])<-!Gam0[,c(1:(ncol(Gam0)-5))]
names(Gam0)

Gam0m<-melt(Gam0[,c(1:25,27)], id.vars="År", variable.name="Species", value.name="Cover",na.rm=T)
head(Gam0m)

########### Gårdsjön #########################

# Skapa dataset där 0 kodas om till NA
Gård0<-Gård
is.na(Gård0[,c(1:(ncol(Gård0)-5))])<-!Gård0[,c(1:(ncol(Gård0)-5))]
names(Gård0)

Gård0m<-melt(Gård0[,c(1:22,24)], id.vars="År", variable.name="Species", value.name="Cover",na.rm=T)
head(Gård0m)

########### Kindla #########################

# Skapa dataset där 0 kodas om till NA
Kind0<-Kind
is.na(Kind0[,c(1:(ncol(Kind0)-5))])<-!Kind0[,c(1:(ncol(Kind0)-5))]
names(Kind0)

Kind0m<-melt(Kind0[,c(1:20,22)], id.vars="År", variable.name="Species", value.name="Cover",na.rm=T)
head(Kind0m)



##### Plottar

par(mfrow=c(2,2),oma=c(0,0,2,0))

# 1. Aneboda
Tä_An<-lm(formula=Cover~År, data=Aneb0m)
summary(Tä_An)

plot(Cover~År, data=Aneb0m, main = "Aneboda")
# Add fit lines
abline(lm(Cover~År, data=Aneb0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(Tä_An)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')


# 2. Gammtratten

Tä_Gm<-lm(formula=Cover~År, data=Gam0m)
summary(Tä_Gm)

plot(Cover~År, data=Gam0m, main = "Gammtratten")
# Add fit lines
abline(lm(Cover~År, data=Gam0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(Tä_Gm)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')



#3. Gårdsjön
Tä_Gd<-lm(formula=Cover~År, data=Gård0m)
summary(Tä_Gd)

plot(Cover~År, data=Gård0m, main = "Gårdsjön")
# Add fit lines
abline(lm(Cover~År, data=Gård0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(Tä_Gd)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

# 4. Kindla
Tä_Ki<-lm(formula=Cover~År, data=Kind0m)
summary(Tä_Ki)

plot(Cover~År, data=Kind0m, main = "Kindla")
# Add fit lines
abline(lm(Cover~År, data=Kind0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(Tä_Ki)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n')



title("Lichen cover", outer=TRUE)


par(mfrow = c(1, 1))


# Jämförelse mellan de båa Aneboda-varianterna

par(mfrow=c(2,1),oma=c(0,0,2,0))

# 1. Aneboda, med avvikare
Tä_An<-lm(formula=Cover~År, data=Aneb0m)
summary(Tä_An)

plot(Cover~År, data=Aneb0m, main = "All trees")
# Add fit lines
abline(lm(Cover~År, data=Aneb0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(Tä_An)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

# 2. Aneboda, utan avvikare
Tä_An1<-lm(formula=Cover~År, data=Aneb1_0m)
summary(Tä_An1)

plot(Cover~År, data=Aneb1_0m, main = "Without outliers")
# Add fit lines
abline(lm(Cover~År, data=Aneb1_0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(Tä_An1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 19, y = 2.5, labels = mylabel)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

title("Lichen cover, Aneboda", outer=TRUE)


par(mfrow = c(1, 1))
