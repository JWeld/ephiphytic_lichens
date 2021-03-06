# Scritp f�r att analysera lavdata fr�n IM

library(vegan)
library(labdsv)
library(plotrix)
library(betapart)
library(vegdata)
library(ggplot2)
library(gridExtra)
library(reshape2)


setwd("C:/Users/grandin/Documents/Ulf IVM/IM/Lavar Utv�rd 2016/R")

spe<-read.table("Lav_spe.txt", header = T, row.names = 1, sep = "\t")
spe<-spe/400

env<-read.table("Lav_env.txt", header = T, row.names = 1, sep = "\t")
as.factor(env$�r)


# Kombinera Omr�de och �r
env$Omr�r<-paste(substring(env$Omr�de,1,2),env$Inv_nr, sep="_")
env$Omr�r<-factor(env$Omr�r)
names(env)
#levels(env$Omr�r)

# Ta bort asp, s�lg och ek
spe1<-subset(spe, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env1<-subset(env, env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z")
env1$TreeSp<-factor(env1$TreeSp)
names(env1)

# Skapa subset med bara tall, gran och bj�rk
#o�del<-env$TreeSp == "PINU SYL" | env$TreeSp == "PICE ABI" | env$TreeSp == "BETULA Z" 


#############################################################
# Sl� samman spe och env, och dela up p� de olika omr�dena
#######################################################################################

spenv<-merge(spe1,env1, all = TRUE, by = 0)
row.names(spenv)<-spenv[[1]] # Fixa radnamn
spenv<-spenv[,-1] # ta bort kolumn med radnamn
levels(spenv$Omr�de)
names(spenv)

############# G�rdsj�n #####
G�rd<-subset(spenv, subset = spenv$Omr�de == "Gd")
# Check number of rows with rowSums = 0
sum(rowSums(G�rd[,c(1:(ncol(G�rd)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-G�rd[!rowSums(G�rd[,c(1:(ncol(G�rd)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(G�rd[,c(1:(ncol(G�rd)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
G�rd<-G�rd[ , which(!apply(G�rd==0,2,all))]
sum(colSums(G�rd[,c(1:(ncol(G�rd)-5))], na.rm=TRUE) == 0)

# �ndra alla Cladonia till Cladonia.sp. f�r Omr�r = Gd_4
# 2011 noterades Cladonia som C.sp, alla tidigare �r som Cladonia.coniocraea
# Detta gav ett v�ldigt utslag f�r 2011 i NMDS och adonis

#Cladonia.coniocraea	Cladonia.digitata	Cladonia.pyxidata	Cladonia.sp.	Cladonia.squamosa

G�rd$Cladonia.sp<-G�rd$Cladonia.coniocraea+G�rd$Cladonia.digitata+G�rd$Cladonia.pyxidata+G�rd$Cladonia.sp.+G�rd$Cladonia.squamosa
G�rd$Cladonia.coniocraea<-NULL
G�rd$Cladonia.digitata<-NULL
G�rd$Cladonia.pyxidata<-NULL
G�rd$Cladonia.sp.<-NULL
G�rd$Cladonia.squamosa<-NULL

names(G�rd)
G�rd<-G�rd[,c(1:5,27,6:26)] # S� att artena komme i bokstavsorning (och inte Clad sp sist)
names(G�rd)

#write.table(G�rd, "C:/Users/Ulf/Ulf IVM/IM/Lavar Utv�rd 2016/G�rdsj�n/G�rd0.txt", quote=FALSE, sep = "\t", dec = ",")#, row.names = TRUE, col.names = TRUE)


############# Aneboda #####
Aneb<-subset(spenv, subset = spenv$Omr�de == "An")
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
# B�da ytorna var nya resp. �r.
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
Kind<-subset(spenv, subset = spenv$Omr�de == "Ki")
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
Gam<-subset(spenv, subset = spenv$Omr�de == "Ga")
# Check number of rows with rowSums = 0
sum(rowSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0)
# Remove rows with only zeros
#gras<-Gam[!rowSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0,]
# Check number of columns with colSums = 0
sum(colSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0)
# Remove columns with only zeros
Gam<-Gam[ , which(!apply(Gam==0,2,all))]
sum(colSums(Gam[,c(1:(ncol(Gam)-5))], na.rm=TRUE) == 0)


#Skapa data f�r ber�knigar


########### Aneboda #########################

# Skapa dataset d�r 0 kodas om till NA
Aneb0<-Aneb
is.na(Aneb0[,c(1:(ncol(Aneb0)-5))])<-!Aneb0[,c(1:(ncol(Aneb0)-5))]
names(Aneb0)

Aneb0m<-melt(Aneb0[,c(1:20,22)], id.vars="�r", variable.name="Species", value.name="Cover",na.rm=T)
head(Aneb0m)

# Alternativt Aneboda med avvikare borttagna
Aneb1_0<-Aneb1
is.na(Aneb1_0[,c(1:(ncol(Aneb1_0)-5))])<-!Aneb1_0[,c(1:(ncol(Aneb1_0)-5))]
names(Aneb1_0)

Aneb1_0m<-melt(Aneb1_0[,c(1:18,20)], id.vars="�r", variable.name="Species", value.name="Cover",na.rm=T)
head(Aneb1_0m)


########### Gammtratten #########################

# Skapa dataset d�r 0 kodas om till NA
Gam0<-Gam
is.na(Gam0[,c(1:(ncol(Gam0)-5))])<-!Gam0[,c(1:(ncol(Gam0)-5))]
names(Gam0)

Gam0m<-melt(Gam0[,c(1:25,27)], id.vars="�r", variable.name="Species", value.name="Cover",na.rm=T)
head(Gam0m)

########### G�rdsj�n #########################

# Skapa dataset d�r 0 kodas om till NA
G�rd0<-G�rd
is.na(G�rd0[,c(1:(ncol(G�rd0)-5))])<-!G�rd0[,c(1:(ncol(G�rd0)-5))]
names(G�rd0)

G�rd0m<-melt(G�rd0[,c(1:22,24)], id.vars="�r", variable.name="Species", value.name="Cover",na.rm=T)
head(G�rd0m)

########### Kindla #########################

# Skapa dataset d�r 0 kodas om till NA
Kind0<-Kind
is.na(Kind0[,c(1:(ncol(Kind0)-5))])<-!Kind0[,c(1:(ncol(Kind0)-5))]
names(Kind0)

Kind0m<-melt(Kind0[,c(1:20,22)], id.vars="�r", variable.name="Species", value.name="Cover",na.rm=T)
head(Kind0m)



##### Plottar

par(mfrow=c(2,2),oma=c(0,0,2,0))

# 1. Aneboda
T�_An<-lm(formula=Cover~�r, data=Aneb0m)
summary(T�_An)

plot(Cover~�r, data=Aneb0m, main = "Aneboda")
# Add fit lines
abline(lm(Cover~�r, data=Aneb0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(T�_An)
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

T�_Gm<-lm(formula=Cover~�r, data=Gam0m)
summary(T�_Gm)

plot(Cover~�r, data=Gam0m, main = "Gammtratten")
# Add fit lines
abline(lm(Cover~�r, data=Gam0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(T�_Gm)
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



#3. G�rdsj�n
T�_Gd<-lm(formula=Cover~�r, data=G�rd0m)
summary(T�_Gd)

plot(Cover~�r, data=G�rd0m, main = "G�rdsj�n")
# Add fit lines
abline(lm(Cover~�r, data=G�rd0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(T�_Gd)
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
T�_Ki<-lm(formula=Cover~�r, data=Kind0m)
summary(T�_Ki)

plot(Cover~�r, data=Kind0m, main = "Kindla")
# Add fit lines
abline(lm(Cover~�r, data=Kind0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(T�_Ki)
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


# J�mf�relse mellan de b�a Aneboda-varianterna

par(mfrow=c(2,1),oma=c(0,0,2,0))

# 1. Aneboda, med avvikare
T�_An<-lm(formula=Cover~�r, data=Aneb0m)
summary(T�_An)

plot(Cover~�r, data=Aneb0m, main = "All trees")
# Add fit lines
abline(lm(Cover~�r, data=Aneb0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(T�_An)
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
T�_An1<-lm(formula=Cover~�r, data=Aneb1_0m)
summary(T�_An1)

plot(Cover~�r, data=Aneb1_0m, main = "Without outliers")
# Add fit lines
abline(lm(Cover~�r, data=Aneb1_0m), col="red", lwd =2)

# Add r2 and p-value to plot
modsum = summary(T�_An1)
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
