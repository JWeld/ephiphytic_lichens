install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(devtools)
install_github("julianfaraway/brinla")
install_github("gfalbery/ggregplot")

library(INLA)
library(faraway)
library(gridExtra)
library(brinla)
library(ggregplot)
library(lme4)
library(nlme)

library(interactions)


dat <- data
dat$Year_i <- as.numeric(dat$Year_i)

#explicit nesting
dat$siteplot <- factor(paste0(dat$Site,dat$Plot))
dat$siteplottree <- factor(paste0(dat$Site,dat$Plot,dat$Tree))
dat$Sensitivity.p <- dat$Sensitivity/9

dat.An <- filter(dat, Site == "An")
dat.Ga <- filter(dat, Site == "Ga")
dat.Gd <- filter(dat, Site == "Gd")
dat.Ki <- filter(dat, Site == "Ki")

#linear model for ref
lmod <- glmer(Nitrogen ~  Year_i +(1|Site/Plot/Tree),family="poisson", data)

lme(rich ~ Year,
    random = ~ 1 + Site, 
    #correlation=corCAR1(form=~Year|Site|plot|Tree_ID),
    data=data, method="REML")

lmer1 <- lmer(Nitrogen ~ Year_s + (1|Site/plot/Tree_ID), data = data)      
lmer2 <- lmer(Nitrogen ~ Year_s + (1 + Year_s|Site/plot), data = data) 

glmer1 <- glmer(Nitrogen ~  Year_s + (1|Site/plot/Tree_ID),family="poisson", data)

fit <- glmer1

(mm_plot <- ggplot(combi, aes(x = Year, y = Nitrogen, colour = plot)) +
        facet_wrap(~Site, nrow=2) +   # a panel for each site
        geom_point(alpha = 0.5) +
        geom_line(data = cbind(combi, pred = predict(fit)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
        theme(legend.position = "none",
              panel.spacing = unit(2, "lines"))  # adding space between panels
)

summary(fit)
ranef(fit) #ranef returns the estimated deviation
#to get the fitted average response per site
response <- fixef(fit) + ranef(fit)$Site
response$Site <-rownames(response)
names(response)[1]<-"Intercept"
response <- response[,c(2,1)]
#plot
ggplot(response,aes(x=Site,y=Intercept))+geom_point()

sim_slopes(fit)

fit <- lmod6
acf(resid(fit), main = "ACF")
pacf(resid(fit), main = "pACF")




lme.S.An <- -lme(Sensitivity~time,random=~1|Plot/Tree,correlation=corCAR1(form=~time|Plot/Tree), data=dat.An)
#longitudanal
psid$slperson <- psid$person
formula <- log(income) ~ cyear*sex+age+educ + f(person, model="iid") + f(slperson, cyear , model="iid")


#INLA
formula <- Sensitivity.p ~ 1 + Year_i  + f(Site, model ="iid") + f(siteplot, model ="iid") + f(siteplottree, model ="iid")
imod <- inla(formula, family="beta", data=dat, control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(imod)

cbind(imod$summary.fixed[,1:2],summary(An.mod1)$coef[,1:2])

cbind(imod$summary.fixed[,1:2],imod2$summary.fixed[,1:2])

bri.hyperpar.summary(imod)
summary(lmod)$sigma

#plot SD
bri.hyperpar.plot(imod)
#plot fixed parameters
bri.fixed.plot(imod)

plot(imod, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)

Efxplot(imod)

resp <- "rich"
covar <- c("n_nh4","n_no3", "latitude", "longitude", "survey_year")
INLAModelSel(resp, covar, "ID_siteplot", "iid", "poisson", data)


# Set prior on precision
prec.prior <- list(prec = list(param = c(0.001, 0.001)))
#nested effects, OBS! MEMORY SUCK!
Zlt <- as(model.matrix( ~ 0 + ID_site:plot.x, data = data), "Matrix")

# Index for technician
data$IDs <- 1:nrow(data)
# Index for technician:sample
data$IDsp <- 1:nrow(data)

imod2 <- inla(rich ~ 1 + n_nh4 + n_no3 + f(survey_year, model="ar1") +
                f(IDs, model = "z", Z = Zlt, hyper = prec.prior) +
                f(IDsp, model = "z", Z = Zlt, hyper = prec.prior),
              data = data, control.predictor = list(compute = TRUE))



#Aneboda only####
formula <- Sensitivity.p ~ 1 + Year_i  + f(siteplot, model ="iid") + f(siteplottree, model ="iid")
imod <- inla(formula, family="beta", data=dat.An, control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(imod)

bri.hyperpar.summary(imod)

#plot SD
bri.hyperpar.plot(imod)
#plot fixed parameters
bri.fixed.plot(imod)

plot(imod, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)

Efxplot(imod) #decrease with time

#Gårdsjön only####
formula <- Sensitivity.p ~ 1 + Year_i  + f(siteplot, model ="iid") + f(siteplottree, model ="iid")
imod <- inla(formula, family="beta", data=dat.Gd, control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(imod)

bri.hyperpar.summary(imod)

#plot SD
bri.hyperpar.plot(imod)
#plot fixed parameters
bri.fixed.plot(imod)

plot(imod, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)

Efxplot(imod) #small increase


#Gammtratten only####
formula <- Sensitivity.p ~ 1 + Year_i  + f(siteplot, model ="iid") + f(siteplottree, model ="iid")
imod <- inla(formula, family="beta", data=dat.Ga, control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(imod)

bri.hyperpar.summary(imod)

#plot SD
bri.hyperpar.plot(imod)
#plot fixed parameters
bri.fixed.plot(imod)

plot(imod, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)

Efxplot(imod) #decrease


#Kindla only####
formula <- Sensitivity.p ~ 1 + Year_i  + f(siteplot, model ="iid") + f(siteplottree, model ="iid")
imod <- inla(formula, family="beta", data=dat.Ki, control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
summary(imod)

bri.hyperpar.summary(imod)

#plot SD
bri.hyperpar.plot(imod)
#plot fixed parameters
bri.fixed.plot(imod)

plot(imod, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)

Efxplot(imod) #small decrease
