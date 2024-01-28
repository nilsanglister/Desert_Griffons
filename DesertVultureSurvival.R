######Desert Survival analysis########
setwd('C:/Users/admin/Desktop/R codes')
library(readxl)

#Path1="C:/Users/admin/Desktop/R codes/Hever_survival.xlsx"
Path1="C:/Users/admin/Desktop/R codes/Hever_survival_unique.xlsx"
#insert path to first excel file "Hever_survival_unique.xlsx"
####Survival for desert site####
heverSurv <- read_excel(Path1)
DF=as.data.frame(heverSurv)
names(DF)
un<-unique(DF$Nili_id)
length(un)

getwd()
dim(DF)
head(DF)
#To see the first five rows, type: 

DF[1:5,]

####packages####

library(survival) # this is the cornerstone command for survival analysis in R
library(ggplot2) # newer package that does nice plots
library(survminer)#for graphs
library(dplyr)
install.packages("survminer")
library(wesanderson)

####defining variables####
DF$sex <- as.factor(DF[,"sex"]) # R calls categorical variables factors
DF$fu_time <- DF[,"fu_time"] # continuous variable (numeric) 
#fu_time (follow-up time, i.e. time in days since admission to hospital) 
DF$death <- DF[,"dead"] # binary variable (numeric) 
DF$status_det<-as.factor(DF[,"status_det"])
hist(DF$fu_time,breaks=52)
DF$season <- as.factor(DF[,"season"]) # R calls categorical variables factors
DF$age_release <- DF[,"age_release"] # continuous variable (numeric) 
DF$origin<- as.factor(DF[,"origin"]) # R calls categorical variables factors
DF$days_acclim <- DF[,"days_acclim"] # continuous variable (numeric) 
DF$pre_2020<- as.factor(DF[,"pre_2020"]) # R calls categorical variables factors
DF$releaseNum<- as.factor(DF[,"repeat"]) # R calls categorical variables factors


#for 90 days
km_fit <- survfit(Surv(DF$fu_ninety, DF$dead_ninety) ~ 1)

plot(km_fit)




summary(km_fit, times = c(7*(1:15))) #the values used for curve weekly
summary(km_fit, times = c(7*(1:8),30*(3:12)))
summary(km_fit, times = c(1:7,14, 21,  30*(1:12)))

################################33
#The above code asks for output every day for the first week, then at 30, 60 and 90 days, and then every 90 days thereafter. Here's the output:

#plot kaplan maier for different predictors

km_season2<- survfit(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$season) 
plot(km_season2)

km_group2<- survfit(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$pre_2020) 
plot(km_group2)#can't use this group 1 is only those that died quickly

km_age<- survfit(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$age_release) 
plot(km_age)#can't use this group too many at age 2

km_releaseNo<- survfit(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$releaseNum) 
plot(km_releaseNo)#can't use this group, group un equal

km_origin<- survfit(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$origin) 
plot(km_origin)#can't use this group, group un equal


survdiff(Surv(time=DF$fu_ninety, event=DF$dead_ninety) ~ DF$season, rho=0) 
survdiff(Surv(time=DF$fu_ninety, event=DF$dead_ninety) ~ DF$age_release, rho=0) 
survdiff(Surv(time=DF$fu_ninety, event=DF$dead_ninety) ~ DF$releaseNum, rho=0) 
survdiff(Surv(time=DF$fu_ninety, event=DF$dead_ninety) ~ DF$origin, rho=0) 

#for releaseNum
DF$releaseNum[DF$releaseNum=='3']='2'

DF$releaseNum=droplevels(DF$releaseNum)
levels(DF$releaseNum)
survdiff(Surv(time=DF$fu_ninety, event=DF$dead_ninety) ~ DF$releaseNum, rho=0) 
names(DF)
table(DF$day_one)
SurvObj<-Surv(time=DF$fu_time, event=DF$death)
fit1 <- survfit(SurvObj ~ season, data = DF)
summary(fit1)
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90))

####graph km for season 90 days USE THIS ONE!!!!!!####
SurvObj<-Surv(time=DF$fu_ninety, event=DF$dead_ninety)
fit1 <- survfit(SurvObj ~ season, data = DF)
summary(fit1)
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90),  break.x.by=10,censor=F)

#with confidence intervals# do not use!
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90),break.x.by=10, conf.int=T,censor=F)

#graph km for age no difference
SurvObj<-Surv(time=DF$fu_ninety, event=DF$dead_ninety)
fit1 <- survfit(SurvObj ~ age_release, data = DF)
summary(fit1)
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90),conf.int=T,censor=F)

table(DF$season)

#graph km for no. of release no difference
SurvObj<-Surv(time=DF$fu_ninety, event=DF$dead_ninety)
fit1 <- survfit(SurvObj ~ releaseNum, data = DF)
summary(fit1)
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90),conf.int=T,censor=F)

#graph km for origin no difference
SurvObj<-Surv(time=DF$fu_ninety, event=DF$dead_ninety)
fit1 <- survfit(SurvObj ~ origin, data = DF)
summary(fit1)
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90),conf.int=T,censor=F)
fit1 <- survfit(SurvObj ~ origin, data = DF)
summary(fit1)
ggsurvplot(fit1, data = DF, pval = TRUE, xlim=c(0,90),  break.x.by=10,censor=F)
#yields the log-rank or Mantel-Haenszel test

####km plot for age####
names(DF)
str(DF$age)
#adding dichtomic age class

                                         
DF$age_release <- as.factor(DF[,"age_release"]) 


####check these variables for models
#Age (assumed to be continuous)
summary(DF$age)
hist(DF$age)
shapiro.test(DF$age)
DF$age[DF$age=='3']='2'
DF$age[DF$age=='4']='2'
DF$age=as.factor(DF$age)
DF$age=droplevels(DF$age)
levels(DF$age)

table(DF$age)
#Gender
names(DF)
t<-table(DF$sex, exclude=NULL)
addmargins(t)
round(100*prop.table(t),digits=1) 

#days in acclimation
str(DF$days_acclim)
summary(DF$days_acclim)
hist(DF$days_acclim)
shapiro.test(DF$days_acclim)

#season
t <- table(DF$season, exclude=NULL) 
addmargins(t) # adds the total (a "sum" column) 
round(100*prop.table(t),digits=1)


#origin
str(DF$origin)
DF$origin<-as.factor(DF$origin)
table(DF$origin, exclude=NULL)




#################################
####Multiple cox model####
cox <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$age_release  + DF$days_acclim + DF$season + DF$origin)

summary(cox)

#only season significant
cox <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~  DF$season )

summary(cox)

cox <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~  DF$season+ DF$days_acclim )

summary(cox)#0.186 coefficient for season and days of acclimation only season 0.27


cox <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~  DF$season+ DF$days_acclim+ DF$age )

summary(cox)#exp of coefiicient not different from without age

#
#for forest graph season + days acclim

################look at 90 days only######
# for 90 days
surv_object <- Surv(time = DF$fu_ninety, event = DF$dead_ninety)
fit.coxph <- coxph(surv_object ~ season, 
                   data = DF)

ggforest(fit.coxph, data = DF)


####Checking forthe proportionality assumption####
fitSex <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$sex) # fit the desired model
fitSeason <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$season) 

fitAge <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$age) 

fitAcclim <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$days_acclim) 

fitOrigin<-coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$origin)

fitSum <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$age_release  + DF$days_acclim + DF$season + DF$origin+DF$reg)

fitSeasonAge<-coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$age_release  +  DF$season )


####checking for other residuals in cox regression#####
res.cox <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$season) 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

               

#####for all sites #####
####################################################3333
######Survival analysis for all release sites########

#for all releases

#insert path for second excel file "released_survival_unique.xlsx"

ReleasedSurv <- read_excel(Path1)
DF=as.data.frame(ReleasedSurv)
names(DF)
un<-unique(DF$Nili_id)
length(un)

getwd()
dim(DF)
head(DF)

#To see the first five rows, type: 

DF[1:5,]


####defining variables####
names(DF)
dim(DF)#107

length(unique(DF$Nili_id))#107
DF$fu_time <- DF[,"fu_time"] # continuous variable (numeric) 
DF$fu_ninety <- DF[,"fu_ninety"] # continuous variable (numeric) 
#fu_time (follow-up time, i.e. time in days since admission to hospital) 
DF$death <- DF[,"dead"] # binary variable (numeric) 
DF$dead_ninety <- DF[,"dead_ninety"] # binary variable (numeric) 

hist(DF$fu_time,breaks=52)
hist(DF$fu_ninety,breaks=13)
DF$region_release<-as.factor(DF[,"region_release"])
DF$season <- as.factor(DF[,"season"]) # R calls categorical variables factors
DF$age_release <- DF[,"age_release"] # continuous variable (numeric) 
DF$origin<- as.factor(DF[,"origin"]) # R calls categorical variables factors
DF$days_acclim <- DF[,"days_acclim"] # continuous variable (numeric) 
DF$rep<- as.factor(DF[,"repeat"]) # R calls categorical variables factors
DF$status_det<- as.factor(DF[,"status_det"]) # R calls categorical variables factors
DF$status<- as.factor(DF[,"status"]) # R calls categorical variables factors
DF$sex=as.factor(DF$sex)
DF$age=as.factor(DF$age)


####check these variables for models
#Age (assumed to be continuous)
DF$age=as.factor(DF$age)
summary(DF$age)

shapiro.test(DF$age)


DF$age[DF$age=='6']='3'
DF$age[DF$age=='5']='3'
DF$age[DF$age=='4']='3'
DF$age=droplevels(DF$age)
levels(DF$age)
t<-table(DF$age, exclude=NULL)
addmargins(t)



#######cox models for 90 days######

###########################################################3
#####Now without 2020-2022 without negev these years comparitive for different sites#####
levels(DF$region_release)
table(DF$region_release)
#without the wild and releases in negev
DFRelease<-DF%>% 
  filter(region_release!='sde boker')%>%
  filter(region_release!="negev")%>%
  filter(release_year<"2020")

#without released in negev  
DF2019<-DFRelease%>% 
  
  filter(release_year<"2020")


table(DF2019$region_release)
 DF2019$region_release=droplevels(DF2019$region_release)
levels(DF2019$region_release)
DF2019$release_year=as.factor(DF2019$release_year)
DF2019$release_year=droplevels(DF2019$release_year)
levels(DF2019$release_year)




####defining variables####
names(DF2019)
un<-unique(DF2019$Nili_id)
length(un)
DF2019$fu_ninety <- DF2019[,"fu_ninety"] # continuous variable (numeric) 
#fu_time (follow-up time, i.e. time in days since admission to hospital) 
DF2019$death <- DF2019[,"dead"] # binary variable (numeric) 
DF2019$dead_ninety <- DF2019[,"dead_ninety"] # binary variable (numeric) 


hist(DF2019$fu_ninety,breaks=13)
DF2019$region_release<-as.factor(DF2019[,"region_release"])
DF2019$season <- as.factor(DF2019[,"season"]) # R calls categorical variables factors
DF2019$age_release <- DF2019[,"age_release"] # continuous variable (numeric) 
DF2019$origin<- as.factor(DF2019[,"origin"]) # R calls categorical variables factors
DF2019$days_acclim <- DF2019[,"days_acclim"] # continuous variable (numeric) 
DF2019$rep<- as.factor(DF2019[,"repeat"]) # R calls categorical variables factors
DF2019$status_det<- as.factor(DF2019[,"status_det"]) # R calls categorical variables factors

DF2019$sex=as.factor(DF2019$sex)
DF2019$age=as.factor(DF2019$age)

#####for 90 days####


#####Cox Models########
cox<-coxph(Surv(fu_ninety,dead_ninety)~season, data=DF2019)
summary(cox)

cox <- coxph(Surv(fu_ninety, dead_ninety) ~ 1, data = DF2019) # take variables straight from DF
summary(cox)

cox <- coxph(Surv(fu_ninety, dead_ninety) ~ region_release, data = DF2019) # take variables straight from DF

summary(cox)

cox <- coxph(Surv(fu_ninety, dead_ninety) ~ age, data = DF2019) # take variables straight from DF
summary(cox)

####check these variables for models
#Age (assumed to be continuous)
DF2019$age=as.factor(DF$age)
summary(DF2019$age)



DF2019$age[DF2019$age=='0']='1'
DF$age[DF$age=='6']='3'
DF$age[DF$age=='5']='3'
DF$age[DF$age=='4']='3'
DF2019$age=droplevels(DF2019$age)
levels(DF2019$age)
t<-table(DF2019$age, exclude=NULL)
addmargins(t)
#Gender
names(DF)
t<-table(DF2019$sex, exclude=NULL)
addmargins(t)
round(100*prop.table(t),digits=1) 

#days in acclimation
str(DF2019$days_acclim)
summary(DF2019$days_acclim)
hist(DF2019$days_acclim)
shapiro.test(DF$days_acclim)

#season
t <- table(DF2019$season, exclude=NULL) 
addmargins(t) # adds the total (a "sum" column) 
round(100*prop.table(t),digits=1)


#origin
str(DF2019$origin)
DF2019$origin<-as.factor(DF2019$origin)
table(DF2019$origin, exclude=NULL)

#region
str(DF2019$region_release)
DF2019$region_release<-as.factor(DF2019$region_release)
table(DF2019$region_release, exclude=NULL)

droplevels(DF2019$region_release)
#reorder to check the relations between golan and hever (can change back)

#DF2019$region_release<-factor(DF2019$region_release,levels=c('Golan','Carmel','Hever'))

DF2019$region_release<-factor(DF2019$region_release,levels=c('Carmel','Golan','Hever','sde boker'))
levels(DF2019$region_release)
table(DF2019$region_release)
#repeat
str(DF2019$rep)
DF2019$rep<-as.factor(DF2019$rep)
table(DF2019$rep, exclude=NULL)
DF2019<- DF2019%>%
  filter(rep!='2')
DF2019$rep=droplevels(DF2019$rep)
table(DF2019$rep)

####survival models#####

cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$origin) 
summary(cox) 

cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age) 
summary(cox) 

cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$rep) 
summary(cox) 

cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$season) 
summary(cox) 

cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$region_release) 
summary(cox) 


####Multiple cox model####
FitSum <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age  + DF2019$days_acclim + DF2019$season + DF2019$origin+  DF2019$region_release+DF2019$age)

summary(FitSum)

FitAge <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age) 

FitSeason <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$season)

FitRegion <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$region_release)

FitOrigin<- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$origin)

FitAcclim<- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$days_acclim)

FitAgeSeason <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age  +  DF2019$season )

FitAgeRegion <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age  +  DF2019$region_release )

FitAgeOrigin <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age  +  DF2019$origin)

FitAgeAcclim <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$age  +  DF2019$days_acclim)

FitSeasonRegion <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$season  +  DF2019$region_release )


FitSeasonOrigin <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$season  +  DF2019$origin)

FitSeasonAcclim <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$season  +  DF2019$days_acclim)

FitRegionOrigin <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$region_release  +  DF2019$origin)

FitRegionAcclim <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$region_release  +  DF2019$days_acclim)

FitOriginAcclim <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$origin  +  DF2019$days_acclim)

ModelNames=c('FitSum','FitAge','FitSeason','FitOrigin','FitRegion', 'FitAcclim','FitAgeSeason', 'FitAgeRegion', 'FitAgeOrigin','FitAgeAcclim','FitSeasonRegion', 'FitSeasonOrigin','FitSeasonAcclim','FitRegionOrigin','FitRegionAcclim', 'FitOriginAcclim')
ModelsList=list(FitSum,FitAge,FitSeason,FitOrigin,FitRegion, FitAcclim,FitAgeSeason, FitAgeRegion, FitAgeOrigin,FitAgeAcclim,FitSeasonRegion, FitSeasonOrigin,FitSeasonAcclim,FitRegionOrigin,FitRegionAcclim, FitOriginAcclim)


#comparing aic
TableAicMix=aictab(cand.set=ModelsList,modnames=ModelNames,weights = TRUE)

print(TableAicMix)



# only region and releaseseason are significant
cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~  DF2019$region_release+ DF2019$season )

summary(cox)#0.186 coefficient for season and days of acclimation only season 0.27




####Checking forthe proportionality assumption####
fit <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$region_release) # fit the desired model
fit <- coxph(Surv(DF$fu_ninety, DF$dead_ninety) ~ DF$season) 
temp <- cox.zph(fit)# apply the cox.zph function to the desired model


#
#for forest graph season + days acclim
surv_object <- Surv(time = DF2019$fu_ninety, event = DF2019$dead_ninety)
fit.coxph <- coxph(surv_object ~ season+ region_release, 
                   data = DF2019)

ggforest(fit.coxph, data = DF2019)
print(temp) # display the results

plot(temp) # plot the curves

#nicer plot
ggcoxzph(temp) 

#kaplan myer plot
km_fit <- survfit(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$region_release) 

plot(km_fit)
plot(km_fit, xlab = "time", ylab = "Survival probability") # label the axes 
#install.packages("xfun")

####checking for other residuals in cox regression#####
res.cox <- coxph(Surv(DF2019$fu_ninety, DF2019$dead_ninety) ~ DF2019$season) 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 


#pretty kaplan meier plots
#graph km for origin no difference
SurvObj<-Surv(time=DF2019$fu_ninety, event=DF2019$dead_ninety)
fit1 <- survfit(SurvObj ~ season, data = DF2019)
summary(fit1)
#with out CI
ggsurvplot(fit1, data = DF2019, pval = TRUE, xlim=c(0,90),censor=F, break.x.by=10)
#with CI
ggsurvplot(fit1, data = DF2019, pval = TRUE, xlim=c(0,90),conf.int=T,censor=F, break.x.by=10)

#this one for 90 days
SurvObj<-Surv(time=DF2019$fu_ninety, event=DF2019$dead_ninety)
fit1 <- survfit(SurvObj ~ region_release, data = DF2019)
summary(fit1)
ggsurvplot(fit1, data = DF2019, pval = TRUE, xlim=c(0,90),conf.int=T,censor=F)

ReleaseTable<-table(DF2019$region_release, DF2019$season)
addmargins(ReleaseTable)


ggsurvplot(fit1, data = DF2019, pval = TRUE, palette = c("#990099", "#339900","#FFCC00",'midnightblue'),xlim=c(0,90), break.x.by=7, censor=F)

#use this one with CI
ggsurvplot(fit1, data = DF2019, pval = TRUE, palette = c("#990099", "#339900","#FFCC00"),xlim=c(0,90), conf.int=T,censor=F, break.x.by=10)






#############For cause of death graph#############

#######Morbidity and Mortality#######
#load packages
library(ggplot2)
library(readxl)
library(tidyverse)
library(hrbrthemes)

library(viridis)
library(dplyr)

Path1="C:/Users/admin/Desktop/R codes/mmForChe.xlsx"
#insert path to upload morbidity and mortality causes from free-roaming vultures
mmVultures <- read_excel(paste(Path1,sep=''))
mmVultures=as.data.frame(mmVultures)


#For all data
#See how many rows and columns you have in your dataset: 

dim(mmVultures)

#looking at how many individuals are sampled
length(unique(mmVultures$ID))#how many different wing tags are there?
#how many individuals died
table (mmVultures$STATUS)

#how many in each year?
table (mmVultures$Year)

#how many from every region?
table (mmVultures$Region)

#how many from Known how many unknown?
table (mmVultures$COM_KNOWN)

#how many in different age_groups??
table (mmVultures$Age_group)

#how many combined poisoning?
table (mmVultures$Cause_o_mm)

#how many from different causes?
table (mmVultures$cause_combined)

#defining data
mmVultures$Age_group=as.factor(mmVultures$Age_group)
mmVultures$Year=as.factor(mmVultures$Year)
mmVultures$Region=as.factor(mmVultures$Region)
mmVultures$COM_KNOWN=as.factor(mmVultures$COM_KNOWN)
mmVultures$cause_combined=as.factor(mmVultures$cause_combined)
mmVultures$cause_combined=reorder(mmVultures$cause_combined, +count(mmVultures$cause_combined))
mmVultures$poison_combined=as.factor(mmVultures$poison_combined)
mmVultures$pestisideVsOther=as.factor(mmVultures$pestisideVsOther)
mmVultures$released_wild=as.factor(mmVultures$released_wild)

###########################
###for all wild vultures from 2015-2019
#if we do use 2020 we have 7 more vultures
#observations=40 individuals=40

mmWild<-mmVultures%>% 
   filter(released_wild!="released")%>%
  filter(Year!="2022")%>%
  filter(Year!="2021")%>%
  filter(Year!="2020")%>%
  filter(Year!="2010")%>%
  filter(Year!="2011")%>%
  filter(Year!="2012")%>%
  filter(Year!="2013")%>%
  filter(Year!="2014")%>%
  
  filter(cause_combined!="NA")

mmWild$Year<-droplevels(mmWild$Year)
mmWild$released_wild=droplevels(mmWild$released_wild)
levels(mmWild$released_wild)
dim(mmWild)
levels(mmWild$Year)
#looking at how many individuals are sampled
length(unique(mmWild$ID))#how many different wing tags are there?

#how many from different causes?
summary_wild <- mmWild %>% group_by(cause_combined) %>%
  summarise(number.birds = n())
t<-table (mmWild$cause_combined)
addmargins(t)

mmWild$cause_combined<-dplyr::recode(mmWild$cause_combined, "illness" = "illness", "infrastructure" = "infrastructure", "lead" = "poisoning", "pesticides" = "poisoning", "other" = "other" , "persecution"="persecution", "uk"="uk")
mmWild$cause_combined<-droplevels(mmWild$cause_combined)
levels(mmWild$cause_combined)

table(mmWild$cause_combined)
####### adding released vultures####
#filter released data
#add path to 3rd excel file released_survival_unique.xlsx"


ReleasedMM <- read_excel(paste(Path2,sep=''))
ReleasedMM=as.data.frame(ReleasedMM)
ReleasedMM$release_year<-as.factor(ReleasedMM$release_year)
ReleasedMM$cause_combined<-as.factor(ReleasedMM$cause_combined)
levels(ReleasedMM$cause_combined)
table(ReleasedMM$cause_combined)

#combine poisoning groups
ReleasedMM$cause_combined<-dplyr::recode(ReleasedMM$cause_combined, "illness" = "illness", "infrastructure" = "infrastructure", "lead" = "poisoning", "pesticides" = "poisoning", "other" = "other" , "persecution"="persecution", "uk"="uk")
ReleasedMM$cause_combined<-droplevels(ReleasedMM$cause_combined)
levels(ReleasedMM$cause_combined)

table(ReleasedMM$cause_combined)
#for 2016-2020 
#if we want to add 2020 
mmRelease<-ReleasedMM%>% 
  filter(release_year!="2022")%>%
  filter(release_year!="2021")%>%
  filter(release_year!="2020")%>%
  filter(cause_combined!="NA")
mmRelease$release_year<-droplevels(mmRelease$release_year)
mmRelease$cause_combined<-droplevels(mmRelease$cause_combined)
levels(mmRelease$release_year)
levels(mmRelease$cause_combined)
#For released
#See how many rows and columns you have in your dataset: 

dim(mmRelease)
names(mmRelease)
#looking at how many individuals are sampled
length(unique(mmRelease$Nili_id))#how many different wing tags are there?
#how many individuals died
table (mmRelease$region_release)

#how many in each year?
table (mmRelease$release_year)


#how many in different age_groups??
table (mmRelease$age_release)



#how many from different causes?
table (mmRelease$cause_combined)


#####################3
#For hever
mmHever$cause_combined=as.factor(mmHever$cause_combined)

mmHever<-mmRelease%>% 
  filter(region_release=="Hever")%>%
  filter(cause_combined!="NA")

mmHever$cause_combined<-droplevels(mmHever$cause_combined)

levels(mmHever$cause_combined)
levels(mmRelease$region_release)

mmHever$region_release=droplevels(mmHever$region_release)
levels(mmHever$region_release)

dim(mmHever)
length(unique(mmHever$Nili_id))
table (mmHever$cause_combined)

mmCarmelGolan<-mmRelease%>% 
  filter(region_release!="Hever")%>%
  filter(region_release!="negev")%>%
  filter(cause_combined!="NA")

#mmCarmel<-mmRelease%>% 
  filter(region_release=="Carmel")%>%
    filter(cause_combined!="NA")
#mmGolan<-mmRelease%>% 
  filter(region_release=="Golan")%>%
  filter(cause_combined!="NA")

dim(mmCarmelGolan)
length(unique(mmCarmelGolan$Nili_id))
cg<-table (mmCarmelGolan$cause_combined)

length(unique(mmCarmel$Nili_id))
table (mmCarmel$cause_combined)


length(unique(mmGolan$Nili_id))
table (mmGolan$cause_combined)
####combining data bases#####

wild_mortality <- mmWild %>% group_by(cause_combined) %>%
  summarise(mortality = n())
wild_mortality<-as.data.frame(wild_mortality)
x<-c("wild")
wild_mortality$origin<-as.factor(x)
sum(wild_mortality$mortality)
wild_mortality$prop<-(wild_mortality$mortality)/49

####ALL released, better  to use carmel golan together######
release_mortality <- mmRelease %>% group_by(cause_combined) %>%
  summarise(mortality = n())
release_mortality<-as.data.frame(release_mortality)
y<-c("released")
release_mortality$origin<-as.factor(y)
sum(release_mortality$mortality)
release_mortality$prop<-(release_mortality$mortality)/55
##########################

hever_mortality<-mmHever %>% group_by(cause_combined) %>%
  summarise(mortality = n())
hever_mortality<-as.data.frame(hever_mortality)
z<-c("hever")
hever_mortality$origin<-as.factor(z)
sum(hever_mortality$mortality)
hever_mortality$prop<-(hever_mortality$mortality)/13


CarmelGolan_mortality<-mmCarmelGolan %>% group_by(cause_combined) %>%
  summarise(mortality = n())
CarmelGolan_mortality<-as.data.frame(CarmelGolan_mortality)
y<-c("carmel_golan")
CarmelGolan_mortality$origin<-as.factor(y)
sum(CarmelGolan_mortality$mortality)
CarmelGolan_mortality$prop<-(CarmelGolan_mortality$mortality)/38


summary.table <- rbind(wild_mortality, CarmelGolan_mortality,hever_mortality) 

summary.table$cause_combined<-factor(summary.table$cause_combined, levels=c('uk','illness','poisoning','infrastructure','persecution','other'))
summary.table<-as.data.frame(summary.table)

####chi.square####
# Create the contingency table
death_proportions <- matrix(c(12, 0,28,9,0,4,0,9,7,2,5,24), nrow = 3, byrow = TRUE)
colnames(death_proportions) <- c("trauma", "illness", "poisoning","uk") # Add the column names
rownames(death_proportions) <- c("Wild", "Judea","Carmel" ) # Add the row names
#Next, you can use the chisq.test() function to perform the chi-square test:
  
 
# Perform chi-square test of independence
result <- chisq.test(death_proportions)

# Print the test result
print(result)
#does not work

#try with only unknown known
death_proportions <- matrix(c(40,9,4,9,14,24), nrow = 3, byrow = TRUE)
colnames(death_proportions) <- c("known","uk") # Add the column names
rownames(death_proportions) <- c("Wild", "Judea","Carmel" ) # Add the row names
#Next, you can use the chisq.test() function to perform the chi-square test:


# Perform chi-square test of independence
result <- chisq.test(death_proportions)

# Print the test result
print(result)
death_proportions<-as.data.frame(death_proportions)
kruskal.test(death_proportions)


# Create the contingency table
death_proportions <- matrix(c(x1, x2, ..., y1, y2, ..., z1, z2, ...), nrow = 3, byrow = TRUE)
colnames(death_proportions) <- c("Cause1", "Cause2", ...) # Add the column names
rownames(death_proportions) <- c("Region1", "Region2", "Region3") # Add the row names

#The chisq.test() function will return the test statistic, degrees of freedom, and the p-value associated with the test. The p-value represents the probability of observing the data given that the null hypothesis is true. If the p-value is below your chosen significance level (e.g., 0.05), you can reject the null hypothesis and conclude that there is a significant association between the death causes and the regions.

#You can also extract specific values from the test result, such as the chi-square statistic, degrees of freedom, and p-value, using the following code:
  
 
# Extract specific values from the test result
chi_square <- result$statistic
degrees_of_freedom <- result$parameter
p_value <- result$p.value

# Print the extracted values
print(chi_square)
print(degrees_of_freedom)
print(p_value)





####graph for causes of mortality####

summary.table$cause_combined<-factor(summary.table$cause_combined, levels=c('uk','illness','pesticides','infrastructure','persecution','lead','nsaids','other'))

#my color scale from wes anderson darjeeling
ggplot(summary.table, aes(x = origin , y = prop, fill = cause_combined)) +   geom_bar(stat = 'identity') + 
  scale_fill_manual(values=c("#5BBCD6","#F2AD00","#FF4444","#F98400", "#00A080", "#046C9A"))+
  ggtitle('Causes of Morbidity and mortality') + 
  
  xlab('Group') + 
  ylab('Proportion') + 
  coord_flip() + 
  theme_minimal()

