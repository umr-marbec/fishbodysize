require(metafor) #metafor
require(phytools) #physignal
require(geiger) #name.check
require(picante) #phylosignal
require(MuMIn) #AICc
require(car) #qqPlot
require(ggplot2)
library(MCMCglmm)
rm(list=ls())

dir=dirname(rstudioapi::getActiveDocumentContext()$path)

# import data
#---------------------------------------

dataset=read.csv2(paste(dir, "/Appendix_1.csv", sep=""))

# subsetting data
#---------------------------------------
dataset <- subset(dataset, dataset$meta.analysis=="body_size") #choose body size data



# editing data for analyses
#---------------------------------------
#Converting r into Z
datasetok <- escalc(measure="ZCOR", ri=r.effect.size, ni=total.N, data=dataset, append=TRUE)
dataset = as.data.frame(datasetok)



# phylogeny
#----------------------------------------
#import tree
tree<-read.newick("FISH_tree.txt")
tree2=collapse.singles(tree, root.edge = FALSE)
tree2


#excluding spp not included in our data set
tree2 = drop.tip(tree2, "Culaea_inconstans")
tree2 = drop.tip(tree2, "Pimephales_promelas")
tree2 = drop.tip(tree2, "Stegastes_planifrons")
tree2 = drop.tip(tree2, "Clupea_harengus")
tree2 = drop.tip(tree2, "Stegastes_partitus")
tree2 = drop.tip(tree2, "Micropterus_salmoides")
tree2 = drop.tip(tree2, "Pastinachus_sephen")




#checking if tree is OK
is.rooted(tree2)
is.binary.tree(tree2) #it's FALSE because there are polytomies
is.ultrametric(tree2)
tree2
plot(tree2, cex=1, lwd=3, edge.width=3)
axisPhylo()



#Converting tree into a variance-covariance matrix:
#---------------------------------------
varcor = vcv(corBrownian(1, tree2), corr=TRUE)
varcor



# creates random factor "id"
#---------------------------------------
id = as.factor(1:(length(dataset$r.effect.size)))
length(id) #131





#=====================================================
# Testing the random-effects
#=====================================================

# 1 random
#-------------
# id only |--> AICc = 331.7047
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id), data=dataset, method="ML")
summary(meta)


# ref |--> AICc = 324.5497 ***BEST MODEL!
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), data=dataset, method="ML")
summary(meta)


# species |--> AICc = 327.9936
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|species), data=dataset, method="ML")
summary(meta)


# phylo |--> AICc = 331.2894
dataset$phylo <- dataset$species
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|phylo), R = list(phylo = varcor), data=dataset, method="ML")
summary(meta)



# 2 random
#-------------
# ref + phylo |--> AICc =  326.6782 
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref, ~1|phylo), R = list(phylo = varcor), data=dataset, method="ML")
summary(meta)


# ref + species |--> AICc =  326.2811 
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref, ~1|species), data=dataset, method="ML")
summary(meta)


# phylo + species |--> AICc =  330.1221 
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|phylo, ~1|species), R = list(phylo = varcor), data=dataset, method="ML")
summary(meta)




# 3 random
#----------------
# ref + phylo + species |--> AICc =  328.4436 
meta <- rma.mv(yi=yi, V=vi, random = list(~1|id, ~1|ref, ~1|phylo, ~1|species), R = list(phylo = varcor), data=dataset, method="ML")
summary(meta)






# The best model includes only "ref" as random factors:
#-------------
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), data=dataset)
summary(meta)



# forest plot
par(mar=c(4, 4.5, 1, 2))

forest.rma(meta, cex=.35, cex.axis=1, cex.lab=1.1, col="gray", xlab="Fisher's Z",  
           xlim=c(-4.5, 4.5), annotate=F, showweights=T, at=c(-4,-3,-2,-1,0,1,2,3,4),
           border="black", order="obs", slab=NA, mlab="Overall effect size", psize=2, 
           pch=16, efac=c(0,0), lty="solid")
par(mar=c(5, 4, 4, 2))







#=====================================================
# I^2 multi-level and H^2
#=====================================================

meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML")
summary(meta)

#variation within study
s2m <- sum(1/meta$vi) * (meta$k-1) / (sum(1/meta$vi)^2 - sum((1/meta$vi)^2))
s2m #s^2_m: variation within study = random error= 0.04941828

s2t <- sum(meta$sigma2) + s2m
s2t #s^2_t= 0.7052767

#the "total" variation
total_I2 = sum(meta$sigma2)/s2t
total_I2 # "total" variation = 92.99%


#Partitioning the "total" variation
meta$sigma2[1] / s2t # I^2 between effect sizes (id) = 0.8089954
meta$sigma2[2] / s2t # I^2 ref = 0.1209352





# Measuring the phylogenetic signal
#--------------------------------------
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref, ~1|phylo), 
               data=dataset, method="REML")
summary(meta)

#variation within study
s2m <- sum(1/meta$vi) * (meta$k-1) / (sum(1/meta$vi)^2 - sum((1/meta$vi)^2))
s2m #s^2_m: variation within study = random error= 0.04941828

s2t <- sum(meta$sigma2) + s2m
s2t #s^2_t= 0.7029109

#the "total" variation
total_I2 = sum(meta$sigma2)/s2t
total_I2 # "total" variation = 0.9296948


#Partitioning the "total" variation
meta$sigma2[1] / s2t # I^2 between effect sizes (id) = 0.7888426
meta$sigma2[2] / s2t # I^2 ref = 0.113723
meta$sigma2[3] / s2t # I^2 phylo = 0.02712931

H2<-meta$sigma2[3] / sum(meta$sigma2)
H2 # 2.92%







#=====================================================
# Egger's regressions and Funnel Plots using "metafor"
#=====================================================

# Egger's regression using the raw effect sizes
#------------------------------
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML")
summary(meta)


s.e. = sqrt(dataset$vi)
precision  = 1/s.e.
Egger1 = lm(dataset$yi*precision ~ precision)
summary(Egger1) #Intercept=-1.2559 P=0.25 (No evidence of publication bias)




# Egger's regression using the Meta-Analytic Residuals
# We should use the best model with moderators
#------------------------------
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ schoaling)
summary(meta)


s.e. = sqrt(dataset$vi)
precision  = 1/s.e.
resid = rstandard(meta)$resid
head(resid)
Egger.resid = lm(resid*precision ~ precision)
summary(Egger.resid) #Intercept = -0.9094; P = 0.393 (No evidence of publication bias)





#Plotting the Funnel Plots
#------------------------------
par(mfrow=c(1,2))
par(mar=c(5, 4.5, 4, 2))


#funnel plot of raw effect sizes
plot(dataset$yi, precision, cex=1.8, pch = 21, bg="hotpink4", cex.axis=1.3, 
     cex.lab=1.3, xlab="Fisher's z", ylab="Precision (1/S.E.)")
abline(v=0.7769,lwd=3)


#funnel plot of Meta-Analytic Residuals
plot(resid, precision, cex=1.8, pch = 21, bg="lightblue4", 
     cex.axis=1.3, cex.lab=1.3, xlab="Meta-analytic residuals",
     ylab="Precision (1/S.E.)")
abline(v=0,lwd=3)


par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 2))
#dev.off()



#=====================================================
# Egger's regression and Funnel Plot using MCMCglmm
#=====================================================

# Running the best model with MCMCglmm to obtain the 
#residuals free from random factors (metafor do not do this!)
#-----------------------------------------
library(MCMCglmm)

#prior

prior2 <- list(G = list(G1 = list(V = 1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list(V = 1, nu = 0.02))
dataset$family<-NULL
dataset$REF= as.factor(dataset$ref)
dataset$IDI= as.factor(dataset$ID)

# one should do this with the "best model"

metaMC=MCMCglmm(yi ~ schoaling, random=~REF +IDI,
                data=dataset,family="gaussian", prior=prior2,
                nitt=130000,burnin=30000,thin=100, verbose=F, pr=T)

# see mixing graphically 
plot(metaMC$Sol[,1]) # fixed effects
plot(metaMC$VCV) # random effects


Prediction1<-predict(metaMC, marginal=~REF+IDI)
s.e. = sqrt(dataset$vi)
precision  = 1/s.e.
MR1<-dataset$yi-Prediction1
zMR1<-MR1*precision

# Egger's regression
Egger1<-lm(zMR1~precision)

# No significant intercept estimate 
summary(Egger1)


# Egger's regression
#--------------------------------------------------
# the intercept is not signficant - no evidnece for funnel asymetry

# lm(formula = zMR1 ~ precision)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7.1138 -2.3562 -0.6235  1.6604 18.3514 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -0.9025     1.0609  -0.851    0.397
# precision     0.3310     0.2353   1.407    0.162


#Intercept = -0.9025; P = 0.397 (No evidence of publication bias)

# funnel plot - plotting orignal and meta-analytic residuals
#--------------------------------------------------


plot(MR1, precision, cex=1.8, pch = 21, bg="lightblue4", 
     cex.axis=1.3, cex.lab=1.3, xlab="Meta-analytic residuals",
     ylab="Precision (1/S.E.)")
abline(v=0,lwd=3)


#dev.off()






#=====================================================
# Testing for the effect of approach speed and SD on effect sizes
#=====================================================

# approach speed
#---------------------------
mean(dataset$approach.speed.cm.sec, na.rm=T)
standDev = sd(dataset$approach.speed.cm.sec, na.rm=T)
StandErr = standDev/sqrt(120)
StandErr


#formating data
data2 <- subset(dataset, dataset$approach.speed.cm.sec >=1)
id = as.factor(1:(length(data2$r.effect.size)))
length(id) #120



# model
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=data2, method="REML",
               mod= ~ approach.speed.cm.sec)
summary(meta)





# SD
#---------------------------
mean(dataset$starting.distance.m, na.rm=T)
standDev = sd(dataset$starting.distance.m, na.rm=T)
StandErr = standDev/sqrt(67)
StandErr


#formating data
data2 <- subset(dataset, dataset$starting.distance.m >=1)
id = as.factor(1:(length(data2$r.effect.size)))
length(id) #67



#model
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=data2, method="REML",
               mod= ~ starting.distance.m)
summary(meta)








#=====================================================
# Collinearity
#=====================================================
source("HighstatLib.R") #Load code of Zuur


dataset$schoaling_coded = as.numeric(dataset$schoaling)
dataset$environment_coded = as.numeric(dataset$environment)
dataset$protected.area_coded = as.numeric(dataset$protected.area)


# creates a data frame only with the covariates
cov <- dataset[,c('mean.bodysize', 
                'schoaling_coded',
                'environment_coded',
                'longevity.years',
                'trophic.level', 
                'protected.area_coded')] 

cov$mean.bodysize <- sqrt(cov$mean.bodysize)
cov$longevity.years <- log10(cov$longevity.years)
cov$trophic.level <- log10(cov$trophic.level)
head(cov)


# VIF 
corvif(cov) # VIF each moderator should have a value < 3. 
is.data.frame(cov)
names = names(cov)


#displays a multiplot figure showing the correlation coefficients and their associated scatter plot
pairs(cov[,],lower.panel = panel.cor, labels=names,  cex.labels=1)








#=====================================================
# Multi-model inference of the covariates
#=====================================================


# creates random factor "id"
#---------------------------------------
id = as.factor(1:(length(dataset$r.effect.size)))
length(id) #131


# changing order of level in schoaling factor
#---------------------------------------
dataset$schoaling = factor(dataset$schoaling, level = c('solitary', 'groups'))



# creates a new function to run in MuMIn
updated.rma.mv <- updateable(rma.mv)
updated.rma.mv


full.model <- updated.rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
                             data=dataset, method="ML",
                             mod= ~ 
                               sqrt(mean.bodysize) + 
                               schoaling +
                               environment +
                               log10(longevity.years) +
                               log10(trophic.level) + 
                               protected.area)
summary(full.model)


# testing call of new function
getCall(full.model)

# update the new function to run in MuMIn
update(full.model)


#=============================
### additional methods for "rma.mv" class (made by Kamil Barton)
### we need this to run model selection with rma.mv in MuMIn
#=============================
formula.rma.mv <- function (x, ...) return(eval(getCall(x)$mods))

makeArgs.rma.mv <-
  function (obj, termNames, comb, opt, ...) {
    ret <- MuMIn:::makeArgs.default(obj, termNames, comb, opt)
    names(ret)[1L] <- "mods"
    ret
  }

nobs.rma.mv <-
  function (object, ...)
    attr(logLik(object), "nall")

coefTable.rma.mv <- function (model, ...)
  MuMIn:::.makeCoefTable(model$b, model$se, coefNames = rownames(model$b))
#=============================


# testing dredge
#dredge(full.model, evaluate=F) # show all candidate models
candidate.models = dredge(full.model, REML = FALSE)
write.csv(candidate.models,'candidate.models_ma_bodysize.csv')


# displays delta AICc <2
subset(candidate.models, delta < 2)


# model averaging
summary(model.avg(candidate.models))


# importance of each predictor
importance(candidate.models)




#================================================================
# coefficients estimated by REML, but models estimated by ML
#================================================================

#---
Model-averaged coefficients:  
  (full average) 
Estimate Std. Error z value Pr(>|z|)
intrcpt                 0.506856   0.496436   1.021    0.307
schoalinggroups         0.241187   0.226803   1.063    0.288
sqrt(mean.bodysize)     0.027716   0.071843   0.386    0.700
log10(longevity.years) -0.063415   0.197253   0.321    0.748
log10(trophic.level)    0.134354   0.436746   0.308    0.758
environmentpelagic     -0.024436   0.112992   0.216    0.829
protected.areayes      -0.005061   0.075789   0.067    0.947

(conditional average) 
Estimate Std. Error z value Pr(>|z|)  
intrcpt                 0.50686    0.49644   1.021   0.3073  
schoalinggroups         0.36076    0.18386   1.962   0.0497 *
  sqrt(mean.bodysize)     0.08657    0.10501   0.824   0.4097  
log10(longevity.years) -0.21323    0.31446   0.678   0.4977  
log10(trophic.level)    0.45827    0.70865   0.647   0.5178  
environmentpelagic     -0.09139    0.20403   0.448   0.6542  
protected.areayes      -0.01993    0.14942   0.133   0.8939  
---
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Relative variable importance: 
                    schoaling sqrt(mean.bodysize) log10(longevity.years) log10(trophic.level) environment protected.area
Importance:          0.67      0.32                0.30                   0.29                 0.27        0.25          
N containing models:   32        32                  32                     32                   32          32          
#---




  
  




#==============================================================
#################### Plotting Moderators #########################
#==============================================================

#Full model
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 sqrt(mean.bodysize) + 
                 schoaling +
                 environment +
                 log10(longevity.years) +
                 log10(trophic.level) + 
                 protected.area)
summary(meta)





#=====================================================
#  schoaling
#=====================================================
dataset$schoaling = factor(dataset$schoaling, level = c('solitary', 'groups'))

data2 = dataset

data2$schoaling = as.factor(data2$schoaling)



#minimal model
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 schoaling-1)
summary(meta)

summary(data2$schoaling) #shows sample size of each level


# organizing a data frame
cols = cbind(meta[[1]], meta[[6]], meta[[7]]); cols
level.names = c("Solitary", "Grouped"); level.names
dat.sum = data.frame(level.names, cols); dat.sum
colnames(dat.sum) = c("levels", "r", "ci.lb", "ci.ub"); dat.sum
is.numeric(dat.sum$r)
factor(dat.sum$levels)
dat.sum$levels = factor(dat.sum$levels, level = c("Solitary", "Grouped"))
factor(dat.sum$levels)
dat.sum



#plot
ggplot(dat.sum, aes(x=levels, y=r)) + 
  geom_abline(intercept=0, slope=0, colour="grey", linetype=2, size=1)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), colour="black", width=.0001, size=1.5) +
  geom_point(size=22, colour="black", show_guide = FALSE )+
  scale_size( guide = "none" )+
  geom_point(size=20, fill = "cornflowerblue", pch=21) +
  xlab("Schoaling behavior") +
  ylab(expression("Fisher's z")) +
  theme_bw((base_size = 26)) +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank ())+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0.005))+
  annotate("text", size = 8, label="35", x = "Solitary", y = -0.5)+
  annotate("text", size = 8, label="96", x = "Grouped", y = -0.5) #+






#=====================================================
#  body size
#=====================================================
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 sqrt(mean.bodysize))
summary(meta)


#all points
ggplot(dataset, aes(x=sqrt(mean.bodysize), y=yi)) + 
  geom_point(size=9, colour="black", show_guide = FALSE )+
  scale_size(guide = "none" )+
  geom_point(size=8, fill = "#D55E00", pch=21) + 
  xlab("Squared root of body size (cm)") +
  ylab(expression("Fisher's z")) +
  theme_classic(base_size = 26)+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank ())+
  geom_abline(intercept=0.4034, slope=0.0796, colour="black", linetype=1, size=1)+
  theme(axis.title.x=element_text(vjust=-0.74))+
  theme(axis.title.y=element_text(vjust=1.2))+
  scale_y_continuous(breaks=c(-1, 0, 1, 2, 3))







#=====================================================
#  longevity
#=====================================================
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 log10(longevity.years))
summary(meta)


#all points
ggplot(dataset, aes(x=log10(longevity.years), y=yi)) + 
  geom_point(size=9, colour="black", show_guide = FALSE )+
  scale_size( guide = "none" )+
  geom_point(size=8, fill = "darkolivegreen", pch=21) + 
  xlab(expression("Log"[10]*" of longevity (years)")) +
  ylab(expression("Fisher's z")) +
  theme_classic(base_size = 26)+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank ())+
  geom_abline(intercept=0.9617, slope=-0.1690, colour="black", linetype=1, size=1)+
  theme(axis.title.x=element_text(vjust=-0.74))+
  theme(axis.title.y=element_text(vjust=1.2))+
  scale_y_continuous(breaks=c(-1, 0, 1, 2, 3))







#=====================================================
#  trophic level
#=====================================================
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 log10(trophic.level))
summary(meta)


#all points
ggplot(dataset, aes(x=log10(trophic.level), y=yi)) + 
  geom_point(size=9, colour="black", show_guide = FALSE )+
  scale_size( guide = "none" )+
  geom_point(size=8, fill = "coral4", pch=21) + 
  xlab(expression("Log"[10]*" of trophic level")) +
  ylab(expression("Fisher's z")) +
  theme_classic(base_size = 26)+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank ())+
  geom_abline(intercept=0.6800, slope=0.2387, colour="black", linetype=1, size=1)+
  theme(axis.title.x=element_text(vjust=-0.74))+
  theme(axis.title.y=element_text(vjust=1.2))+
  scale_y_continuous(breaks=c(-1, 0, 1, 2, 3))







#=====================================================
#  Environment
#=====================================================
data2 = dataset

data2$environment = as.factor(data2$environment)


#minimal model
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 environment-1)
summary(meta)

summary(data2$environment) #shows sample size of each level


# organizing a data frame
cols = cbind(meta[[1]], meta[[6]], meta[[7]]); cols
level.names = c("Demersal", "Pelagic"); level.names
dat.sum = data.frame(level.names, cols); dat.sum
colnames(dat.sum) = c("levels", "r", "ci.lb", "ci.ub"); dat.sum
is.numeric(dat.sum$r)
factor(dat.sum$levels)
dat.sum$levels = factor(dat.sum$levels, level = c("Demersal", "Pelagic"))
factor(dat.sum$levels)
dat.sum



#plot
ggplot(dat.sum, aes(x=levels, y=r)) + 
  geom_abline(intercept=0, slope=0, colour="grey", linetype=2, size=1)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), colour="black", width=.0001, size=1.5) +
  geom_point(size=22, colour="black", show_guide = FALSE )+
  scale_size( guide = "none" )+
  geom_point(size=20, fill = "burlywood4", pch=21) +
  xlab("Environment") +
  ylab(expression("Fisher's z")) +
  theme_bw((base_size = 27)) +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank ())+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0.005))+
  annotate("text", size = 8, label="104", x = "Demersal", y = -0.5)+
  annotate("text", size = 8, label="27", x = "Pelagic", y = -0.5) #+






#=====================================================
#  Protected area
#=====================================================
data2 = dataset

data2$environment = as.factor(data2$protected.area)


#minimal model
meta <- rma.mv(yi=yi, V=vi, random = list(~ 1|id, ~1|ref), 
               data=dataset, method="REML",
               mod= ~ 
                 protected.area-1)
summary(meta)

summary(data2$protected.area) #shows sample size of each level


# organizing a data frame
cols = cbind(meta[[1]], meta[[6]], meta[[7]]); cols
level.names = c("Unprotected", "Protected"); level.names
dat.sum = data.frame(level.names, cols); dat.sum
colnames(dat.sum) = c("levels", "r", "ci.lb", "ci.ub"); dat.sum
is.numeric(dat.sum$r)
factor(dat.sum$levels)
dat.sum$levels = factor(dat.sum$levels, level = c("Unprotected", "Protected"))
factor(dat.sum$levels)
dat.sum



#plot
ggplot(dat.sum, aes(x=levels, y=r)) + 
  geom_abline(intercept=0, slope=0, colour="grey", linetype=2, size=1)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), colour="black", width=.0001, size=1.5) +
  geom_point(size=22, colour="black", show_guide = FALSE )+
  scale_size( guide = "none" )+
  geom_point(size=20, fill = "darkgoldenrod2", pch=21) +
  xlab("Area protection status") +
  ylab(expression("Fisher's z")) +
  theme_bw((base_size = 27)) +
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank ())+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0.005))+
  annotate("text", size = 8, label="64", x = "Unprotected", y = -0.5)+
  annotate("text", size = 8, label="67", x = "Protected", y = -0.5) #+


#-------- END -----------