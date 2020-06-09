library(mgcv)
library(mgcViz)
dat <- data

dat$Year_i <- as.ordered(dat$Year_i)
dat$Site <- as.ordered(dat$Site)
dat$Plot <- as.ordered(dat$Plot)
dat$Tree <- as.ordered(dat$Tree)
dat$id2 <- as.ordered(dat$id2)
dat$ID <- as.ordered(dat$ID)

dat$Year <- as.factor(dat$Year_i)

dat$Year <- as.numeric(as.character(dat$Year))
dat$Tree <- as.character(dat$Tree)

dat$Site <- as.factor(dat$Site)
dat$Plot <- as.factor(dat$Plot)
dat$Tree <- as.factor(dat$Tree)
dat$id2 <- as.factor(dat$id2)
dat$sp <- as.factor(dat$sp)
dat$ID <- as.factor(dat$ID)


gamm.s <-
  gamm(
    Sensitivity/9 ~ s(Year,Site,bs="fs",k=10),
      #s(Site/Plot/Tree, bs = "re"),
    correlation = corCAR1(form=~Year|Site/Plot/Tree),
    family = betar,
    data = data
  )

gamm.s <-
  gamm(
    Sensitivity/9 ~ s(Year,k=10, by = Site),
    #s(Site/Plot/Tree, bs = "re"),
    correlation = corCAR1(form=~Year|Site/Plot/Tree),
    family = betar,
    data = dat
  )

#s(id2, bs="re")+ s(Site, bs="re")

gamm.sp <-
  gamm(
    n_sp ~ s(Year,k=10, by = Site) +
      s(Site/Plot/Tree, bs = "re"),
    correlation = corCAR1(form=~Year|Site/Plot/Tree),
    family = poisson,
    data = dat
  )

gamm.shan <-
  gamm(
    Shan ~s(Year,k=10, by = Site) +
      s(Site/Plot/Tree, bs = "re"),
    correlation = corCAR1(form=~Year|Site/Plot/Tree),
    family = gaussian,
    data = dat
  )

gamm.n <-
  gamm(
    Nitrogen/9 ~s(Year,k=10, by = Site) +
      s(Site/Plot/Tree, bs = "re"),
    correlation = corCAR1(form=~Year|Site/Plot/Tree),
    niterPQL=50,
    family = betar,
    data = dat
  )

gamm.ph <-
  gamm(
    pH ~ s(Year,k=10, by = Site) +
      s(Site/id2, bs = "re"),
    correlation = corCAR1(form=~Year|Site/id2),
    family = gaussian,
    data = dat
  )


gam.s0 <- gam(Sensitivity ~ s(Year, k=5, m=2, bs="tp") +
      s(Year, by= Site, k=5, m=1, bs="tp") +
      s(Site, bs="re", k=12),
      data=dat, method="REML")


gam.s <- gam(Sensitivity ~ 
                   s(Year, by=Site, m=1, bs="tp") +
                   s(id2, bs="re")+ s(Site, bs="re"), data=dat, method="REML")

gam.sp <- gam(n_sp ~ 
               s(Year, by=Site, m=1, bs="tp") +
               s(Site/Plot/Tree, bs="re"), data=dat, method="REML")

gam.shan <- gam(Shan ~ 
               s(Year, by=Site, m=1, bs="tp") +
               s(Site/Plot/Tree, bs="re", k=10), data=dat, method="REML")

gam.n <- gam(Nitrogen ~ 
               s(Year, by=Site, m=1, bs="tp") +
               s(Site/Plot/Tree, bs="re", k=20), data=dat, method="REML")

gam.ph <- gam(pH ~ 
               s(Year, by=Site, m=1, bs="tp") +
               s(Site/Plot/Tree, bs="re", k=20), data=dat, method="REML")


fit <- gamm.n$gam
fit <- gam.s
summary(fit)
plot(fit, pages = 1, all.terms = TRUE, shade = TRUE, shift = coef(fit)[1])
gam.check(fit, pages = 1)
acf(resid(fit), main = "ACF")
pacf(resid(fit), main = "pACF")

b <- getViz(fit)
print(plot(b, allTerms = T), pages = 1) # Calls print.plotGam()


#qgam####
library(qgam)

qgam.s <-
  qgam(
    Nitrogen ~ s(Year,Site,bs="fs",k=5)+
      s(Site/Plot/Tree, bs = "re"),
    data = dat, qu = 0.5
  )

qgam.s <-
  qgam(
    Nitrogen ~ s(Year,Site,bs="fs",k=5),
    data = dat, qu = 0.5
  )
fit <- qgam.s
summary(fit)
plot(fit, pages = 1, all.terms = TRUE, shade = TRUE, shift = coef(fit)[1])
check.qgam(fit)
gam.check(fit)

b <- getViz(fit)
print(plot(b, allTerms = T), pages = 1) # Calls print.plotGam()

