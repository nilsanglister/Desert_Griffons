########For Hever article combining several codes together########
######questions and and average temp. for days of death for released vultures compared to control wild vultures judea, negev
#2) Height above ground
#3) Total distance/Max. displacement
#4) Flight yes or no



###loading packages###

library(glmmTMB)
library(caret)
library(emmeans)
library(FSA)
library(readxl)
library(ggplot2)
library(nlme)
library(AICcmodavg)
library('suncalc')
library("raster")
library('lme4')
library('viridis')
library(prediction)
library(DHARMa)
library(lmerTest)
###loading packages shorthand###
pcks=list('tidyverse','sp','dplyr',"move","mapproj","ggmap","lubridate",
          "htmltools","leaflet.extras","geosphere") 
sapply(pcks, require, char=TRUE)
#does not download adehabitat

#setting directory
#set directory


#####for temps at time of death#####
#This compares the maximum temperature per day of GPS transmitter before hour of death of a released vulture compared to the maximum temp per day for a wild born vulture in the same time at the same area

#this one is for paired design (6 released, 6 wild)
#create path 1 to upload excel file of temperatures

#checking predictors
hevertemps <- read_excel(Path1)
hevertemps=as.data.frame(hevertemps)
names(hevertemps)
un<-unique(hevertemps$vult_tag)
length(un)
HeverWild<- hevertemps%>%
  filter(captive=="Wild") 
un<-unique(HeverWild$vult_tag)
length(un)


hevertemps$vult_tag=as.factor(hevertemps$vult_tag)
levels(hevertemps$vult_tag)
hevertemps$group_release=as.factor(hevertemps$group_release)
levels(hevertemps$group_release)
hevertemps$captive=as.factor(hevertemps$captive)
levels(hevertemps$captive)
str(hevertemps$maxtemp)
table(hevertemps$captive)


#plotting
hist(hevertemps$maxtemp)
ggplot(data=hevertemps, aes(x=captive, y=maxtemp, fill=captive))+ geom_boxplot() + geom_smooth(se=FALSE, method="lm")+theme_bw(base_size = 10)

#regular statistics
summary(hevertemps$maxtemp)

group_by(hevertemps, captive) %>%
  dplyr::summarise(    count = n(),
                       mean = mean(maxtemp, na.rm = TRUE),
                       sd = sd(maxtemp, na.rm = TRUE),
                       se=sd/sqrt(count),
                       median = median(maxtemp, na.rm = TRUE),
                       IQR = IQR(maxtemp, na.rm = TRUE),
                      
                       Minimum=min(maxtemp, na.rm = TRUE),
                       Maximum=max(maxtemp,na.rm = TRUE))


shapiro.test(hevertemps$maxtemp)#is normaly distributed p=0.05


Modeltemps=glmmTMB(maxtemp  ~ captive + (1|group_release) ,family=Gamma(link="log"), data=hevertemps, verbose = F)
summary(Modeltemps)

#checking model assumptions
x=simulateResiduals(Modeltemps)
plot(x)




####getting data from movebank and xcel####
#for the released vs. all the wild vultures in the region
#only for 10 days after release excluding day one for released animals and wild-born in the same region, downloaded movement data
#data has been shifted in order to avoide critical data for conservation and we only used data in the region of Israel and Jordan to avoid long-range forrays

#upload csv of movement data for griffons:"mergeHever.csv"

#create path 2
mergeHever <- read.csv(Path2)

mergeHever<-as.data.frame(mergeHever)

unique(mergeHever$local_identifier)
mergeHever$released_wild<-mergeHever$released.wild

####entering functions####


####height calculations function####
getDEM <- function(type="file",box,filename=DEM_name,resoluton=NULL) # a subfunction used within GroudHeight
{
  # producing a DEM file in the "box" limits based on the provided data and type
  # other DEM files are available at http://sendimage.whu.edu.cn/res/DEM_share/SRTM1/N30,E30/
  
  # concatenate two raster files: 
  # imported_raster1 <- raster(str_name)
  # str_name<-"DEM_Files/N33E035/n33_e035_1arc_v3.tif" # DEM on geographic grid
  # imported_raster2 <- raster(str_name)
  # str_name<-"DEM_Files/from_michal/n32_e035_1arc_v3.tif.ovr" # DEM on 0:1800 grid 
  # imported_raster <- merge(imported_raster1,imported_raster2)
  # save a merged raster file locally for further use: 
  # writeRaster(imported_raster, filename="DEM_Harod.tif", format="GTiff", overwrite=TRUE)
  wgs84 <<- "+proj=longlat +ellps=WGS84 +datum=WGS84"
  coordinates(box) <- ~lon+lat
  proj4string(box) <- CRS(wgs84)
  if (type=="file")
  {imported_raster <- raster(filename)
  DEM <- crop(imported_raster, extent(box), snap="out")
  print(sprintf("the source resolution is (%3.1fm lat, %3.1fm lon)", res(DEM)[1]*111e3,res(DEM)[2]*94e3))
  if (!is.null(resoluton)) { 
    requiredRes=resoluton/111e3  #resolution in degrees
    DEM <- projectRaster(DEM,crs=CRS(wgs84),res=requiredRes)
    print(sprintf("the final resolution is (%3.1fm lat, %3.1fm lon)", res(DEM)[1]*111e3,res(DEM)[2]*94e3))
  }
  }
  if (type=="web")
  {
    # getData('ISO3')
    # getData("alt", country = 'ISR')
    DEM1 <- raster::getData('SRTM', lon=box$lon[1], lat=box$lat[1])
    DEM2 <- raster::getData('SRTM', lon=box$lon[2], lat=box$lat[2])
    DEM3 <- raster::getData('SRTM', lon=box$lon[1], lat=box$lat[2])
    DEM4 <- raster::getData('SRTM', lon=box$lon[2], lat=box$lat[1])
    DEM <- raster::merge(DEM1,DEM2,DEM3,DEM4)
    DEM <- raster::crop(DEM, extent(box), snap="out")
    print(sprintf("the source resolution is (%3.1fm lat, %3.1fm lon)", res(DEM)[1]*111e3,res(DEM)[2]*94e3))
    if (!is.null(resoluton)) { 
      requiredRes=resoluton/111e3  #resolution in degrees
      DEM <- projectRaster(DEM,crs=CRS(wgs84),res=requiredRes)
      print(sprintf("the final resolution is (%3.1fm lat, %3.1fm lon)", res(DEM)[1]*111e3,res(DEM)[2]*94e3))
    }
  }
  return(DEM)
}
#Use data frame with location_lat & long
#####function for adding ground height####
addGroundHeight <- function(df)
{
  wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
  box <- data.frame(lat=c(min(df$location_lat),max(df$location_lat)),lon=c(min(df$location_long),max(df$location_long)))
  DEM <- getDEM(type="web",box)
  # ANTS.df <- ANTS.df[which(!is.na(ANTS.df$LAT)&!is.na(ANTS.df$LON)),]
  coordinates(df) <- ~location_long+location_lat
  proj4string(df) <- CRS(wgs84)
  ANTS <- SpatialPoints(df, proj4string=CRS(wgs84))
  df$GroundASL <- raster::extract(DEM,ANTS)
  # return(as.data.frame(ANTS.df))
  return(as.data.frame(df))
}

#### function for visualising data in leaflet####
atl_mapleaf <- function(llpd1,Year_Month=NULL)
{
  # llpd1 <- llpd1 %>% dplyr::filter(tag_local_identifier==tagIdentifier)
  if (!is.null(Year_Month))
    llpd1 <- llpd1 %>% dplyr::filter(YearMonth==Year_Month)
  itm<-"+init=epsg:2039 "
  wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  coordinates(llpd1)<-~location_long+location_lat
  proj4string(llpd1)<-CRS(wgs84)
  # proj4string(dd1)<-CRS(itm)
  # llpd1 <- spTransform(dd1, wgs84)
  
  require("RColorBrewer")
  require("leaflet.extras")
  require("htmltools")
  # display.brewer.all()
  # display.brewer.pal(n = 4, name = 'RdYlBu')
  # col=brewer.pal(n = 4, name = 'RdYlBu')
  col=brewer.pal(n = 6, name = 'Dark2')
  
  ll<-leaflet() %>% 
    addProviderTiles('Esri.WorldImagery') %>% # 'CartoDB.Positron' 'OpenStreetMap.Mapnik' 'Stadia.AlidadeSmooth','CartoDB.Positron'
    addCircles(data=llpd1, weight = 5, fillOpacity = 1,color = col[4],group="original",
               popup = ~htmlEscape(paste0("time=",as.character((llpd1$timestamp)),
                                          ", identifier=",as.character((llpd1$local_identifier)),
                                          ", sat_count=",as.character((llpd1$gps_satellite_count)),
                                          ", ground_speed=",as.character((llpd1$ground_speed)),
                                          ", height=",as.character(round(llpd1$height_above_msl)),
                                          ", groundElev=",as.character(round(llpd1$GroundASL)),
                                          ", temperature=",as.character(round(llpd1$external_temperature)),
                                          ", voltage=",llpd1$tag_voltage,
                                          ", charge=",llpd1$battery_charge_percent  
               ))) %>%
    addPolylines(data=llpd1@coords, weight = 1, opacity = 1,col=col[4],group="original") %>% 
    addScaleBar(position = c("bottomleft"), 
                options = scaleBarOptions(imperial=FALSE,maxWidth=200)) 
  # %>% 
  #   addLayersControl(
  #     overlayGroups = c("original","GPS"),
  #     options = layersControlOptions(collapsed = FALSE, autoZIndex=TRUE)) 
  ll
}
####Movement daily parameters####
##eitam I also need temperature and height not during flight, see if ok!!!#
MovementVar <- function(data,flight_speed_threshold=4,Rename=F)
{
  if(Rename)#needed if taken from csv data
    data <- data %>%  rename(local_identifier='individual-local-identifier',
                             location_long='location-long', location_lat='location-lat',
                             ground_speed='ground-speed',height_above_msl='height-above-msl')
  
  data <- data %>%dplyr::mutate(date=date(timestamp)) %>% 
    dplyr::select('local_identifier','location_long','location_lat','ground_speed','external_temperature',
                  'height_above_msl','date','timestamp',"released_wild", "repeat."  ,'GroundASL','origin','age_release',"days_acclim" ,
                  "status",  'pre_2020','group_release','day_one', 'day_ten','timefromday1','season','status', 'death_date', 'days_to_death')
  
  #if we are using other data and other predictors check in the selected columns above
  
  
  #get times of sunset sunsrise for every bout###
  sun =  getSunlightTimes(data=data.frame(date=date(data$timestamp),lat=data$`location_lat`,lon=data$`location_long`),keep =c("sunset","sunrise"))
  data$sunrise <- sun$sunrise
  data$sunset <- sun$sunset
  #get one hour before sunrise, one hour after sunset####
  dailymov<-data %>% 
    arrange(timestamp) %>% 
    filter(timestamp>(sunrise-duration(1, "hours"))) %>% 
    filter(timestamp<(sunset+duration(1, "hours")))  %>%
    group_by(local_identifier,date,group_release) %>%
    dplyr::mutate(pointsPerDay=n(),#how many gps points I have for every day for this tag
           step_time=as.vector(as.numeric(timestamp-lag(timestamp),units ='mins')),#time difference between two consecutive points
           step_time25=quantile(step_time,probs=0.25,na.rm=TRUE),
           flightHeight=ifelse(height_above_msl-GroundASL<0,0,height_above_msl-GroundASL+13)) %>% #calculating the time difference quartiles 
    # filter(step_time25<25) %>% # here we cut off data that has a time difference larger than 25 minutes between points
    #change to 32 if too much data is lost
    #filter(pointsPerDay>30) %>% # here we cut off data that has less than 30 bouts per day
    
    dplyr::mutate(displacement=distHaversine(cbind(location_long[1],location_lat[1]),
                                      cbind(location_long,location_lat)),
           deltaHeight0=abs(height_above_msl-height_above_msl[1]),
           step_dist=distHaversine(cbind(location_long,location_lat),
                                   cbind(lead(location_long),lead(location_lat))),
           flight_stepTime=ifelse(ground_speed>flight_speed_threshold,step_time,NA),#time between bouts of flight
           flightHeight=ifelse(ground_speed>flight_speed_threshold,flightHeight,NA),#height above ground if flying with correction negative valuse changed to zero. 13 meters added
           flightspeed=ifelse(ground_speed>flight_speed_threshold,ground_speed,NA),
           flightTime=(ifelse(ground_speed>flight_speed_threshold,timestamp,NA)),#time at flight (defined in previous line above 2 meter per second)
           segmentType=(ifelse(ground_speed>flight_speed_threshold,1,0)),#flying=1, notflying segment=0
           Diff=segmentType-lag(segmentType),#here we want to identify transition between flight and non-flight state
           set=cumsum(ifelse(is.na(Diff),0,as.numeric(Diff!=0)))# gives numbers to flight segments in a day
           
    ) %>%
    group_by(local_identifier,date,set) %>% 
    dplyr::mutate(
      segTotal_distance=sum(step_dist, na.rm=TRUE),
      segFlight_time=sum(flight_stepTime,na.rm=TRUE),
    )  %>% 
    group_by(local_identifier,date) %>%
    dplyr::mutate( numberOfFlighetSegments=max(set)/2,
            avgSegmentDistance=mean(segTotal_distance,na.rm=TRUE),
            sdSegmentDistance=sd(segTotal_distance,na.rm=TRUE),
            avgSegmentTime=mean(segFlight_time,na.rm=TRUE),
            sdSegmentTime=sd(segFlight_time,na.rm=TRUE),
            locationsPerDay=n(),
            net_displacement=tail(displacement,1),
            max_displecement=max(displacement, na.rm=TRUE),
            max_deltaHeight0=max(deltaHeight0, na.rm=TRUE),
            mean_height=mean(height_above_msl),
            deltaHeight0_75=quantile(deltaHeight0,probs=0.75,na.rm=TRUE),
            deltaHeight0_50=quantile(deltaHeight0,probs=0.5,na.rm=TRUE),
               max_temp=max(external_temperature,na.rm=TRUE),
            median_temp=median(external_temperature,na.rm=TRUE),
            total_distance=sum(step_dist, na.rm=TRUE),
            maxFlight_heightAGL=max(flightHeight, na.rm=TRUE),
            avgFlight_heightAGL=mean(flightHeight,na.rm=TRUE),
            sdFlight_heightAGL=sd(flightHeight,na.rm=TRUE),
            flight_time=sum(flight_stepTime,na.rm=TRUE),
            maxFlight_speed=max(flightspeed, na.rm=TRUE),
            avgFlight_speed=mean(flightspeed,na.rm=TRUE),
            sdFlight_speed= sd(flightspeed,na.rm=TRUE),
            setOffTimeMinute=(min(flightTime,na.rm=TRUE)-as.numeric(sunrise))/60,
            returnBeforesunsetMinute=(as.numeric(sunset)-max(flightTime,na.rm=TRUE))/60,
            
    )  %>%
    slice(1) %>% 
    dplyr:: select('local_identifier','date','timestamp','location_long','location_lat','locationsPerDay','net_displacement','max_displecement',
                   'total_distance','avgFlight_heightAGL','avgFlight_heightAGL','sdFlight_heightAGL','flight_time',
                   'maxFlight_speed','avgFlight_speed','sdFlight_speed','setOffTimeMinute','returnBeforesunsetMinute','mean_height',
                   'numberOfFlighetSegments','avgSegmentDistance','sdSegmentDistance','avgSegmentTime','sdSegmentTime',"released_wild"  ,'GroundASL','origin','age_release',"days_acclim" ,
                   "status",  'pre_2020','group_release','day_one', 'day_ten','timefromday1','season','status', 'death_date', 'days_to_death','max_deltaHeight0','released_wild','repeat.','origin',
                   'max_temp','median_temp','deltaHeight0_75', 'deltaHeight0_50') %>% 
    # #for the mycoplasma data set I enter the predictors as pos_neg, origin, age_class
    #if we are using other data and other predictors check in the selected columns above
    
    ungroup()
  return(dailymov)}



###move function for flight
### constant parameters ###
flight_speed_threshold <- 4



####working with functions####

mergeHever_ten_d<-as.data.frame(mergeHever)
mergeHever_ten_d$timestamp <- as.POSIXct(mergeHever_ten_d$timestamp, format = c("%Y-%m-%d %H:%M:%S"), tz = "UTC")
mergeHever_ten_d <- addGroundHeight(mergeHever_ten_d)
mergeHever_ten_d <- mergeHever_ten_d %>% dplyr::mutate(flightHeight=height_above_msl-GroundASL,
                      flightHeight_corrected=ifelse(flightHeight<0,0,flightHeight+13))

names(mergeHever_ten_d)
mergeHeverMov<-MovementVar(mergeHever_ten_d,flight_speed_threshold = 4)#can change flightspeed threshold here

un<-unique(mergeHeverMov$local_identifier)
length(un)



####for statistics, check and clean up data####
#first YEAR
str(mergeHeverMov)
mergeHeverMov[complete.cases(mergeHeverMov), ]


mergeHeverMov$released_wild=as.factor(mergeHeverMov$released_wild)
levels(mergeHeverMov$released_wild)

mergeHeverMov$repeat.=as.factor(mergeHeverMov$repeat.)
levels(mergeHeverMov$repeat.)

mergeHeverMov$origin=as.factor(mergeHeverMov$origin)
levels(mergeHeverMov$origin)

##categorizing age into groups###
mergeHeverMov$age_class <- ifelse(mergeHeverMov$age_release < 1, "juvenile", 
                                 ifelse(mergeHeverMov$age_release>=1 & mergeHeverMov$age_release   <=4, "subadult", "adult"))
                                       
mergeHeverMov$age_class <-as.factor(mergeHeverMov$age_class)  
levels(mergeHeverMov$age_class)                                             

mergeHeverMov$status=as.factor(mergeHeverMov$status)
levels(mergeHeverMov$status)

mergeHeverMov$pre_2020=as.factor(mergeHeverMov$pre_2020)
levels(mergeHeverMov$pre_2020)

mergeHeverMov$group_release=as.factor(mergeHeverMov$group_release)
levels(mergeHeverMov$group_release)

mergeHeverMov$season=as.factor(mergeHeverMov$season)
levels(mergeHeverMov$season)

mergeHeverMov$timefromday1=as.integer(mergeHeverMov$timefromday1)
levels(mergeHeverMov$timefromday1)


####adding columns with distance in km, spped km/hr, time in hr####
HeverMov<-mergeHeverMov%>%
  mutate(net_km=net_displacement/1000, max_km=max_displecement/1000, 
         totdiskm=total_distance/1000, maxspeedkm=maxFlight_speed*3.6, 
         avgspeedkm=avgFlight_speed*3.6, avgsegkm=avgSegmentDistance/1000, avgseghr=avgSegmentTime/60, setOffhr=setOffTimeMinute/60,
         returnhr=returnBeforesunsetMinute/60,flightHr=flight_time/60,
         dayFromDayOne=as.Date(date)-as.Date(day_one))
HeverMov<-as.data.frame(HeverMov)
#########################################################



####filtering and defining factors####
#only more than 12 locations perday
HeverCompare<- HeverMov%>%
  filter(locationsPerDay>12) 


HeverCompare$pre_2020=droplevels(HeverCompare$pre_2020)
levels(HeverCompare$pre_2020)
un<-unique(HeverCompare$local_identifier)
length(un)
table(HeverCompare$local_identifier,HeverCompare$pre_2020 )
HeverCompare$local_identifier=as.factor(HeverCompare$local_identifier)
HeverCompare$local_identifier=droplevels(HeverCompare$local_identifier)

HeverCompare$group_release=as.factor(HeverCompare$group_release)
HeverCompare$date <- as.Date(HeverCompare$date)

#finding how many days for individual
n_days_df <- HeverCompare %>% group_by(local_identifier, group_release) %>% mutate(n_days = as.numeric(difftime(max(date), min(date), units="days")))%>%
  ungroup()
#filtering out individual releases with less than 3 days of data
HeverCompare<-n_days_df%>% 
  filter(n_days>=2)


HeverCompare$released_wild=droplevels(HeverCompare$released_wild)
HeverCompare$date=as.factor(HeverCompare$date)



#########################
#can start here for the statistics#####
#can upload data after calculation of the movement functions
#create path 3 for csv file of movement calculations per day "HeverCompare.csv"
HeverCompare <- read.csv(Path3)

un<-unique(HeverCompare$local_identifier)
length(un)

levels(HeverCompare$age_class)
HeverSubadult<- HeverCompare%>%
  filter(age_class=="subadult") %>%
  filter(released_wild=="released")
un<-unique(HeverSubadult$local_identifier)
length(un)
HeverAdult<- HeverCompare%>%
  filter(age_class=="adult") %>%
  filter(released_wild=="released")
un<-unique(HeverAdult$local_identifier)
length(un)
HeverJuvenile<- HeverCompare%>%
  filter(age_class=="juvenile") %>%
  filter(released_wild=="released")
un<-unique(HeverAdult$local_identifier)
length(un)

#####comparing only subadults for both released and wild####
HeverCompare1<- HeverCompare%>%
  filter(age_class=="subadult") 
un<-unique(HeverCompare1$local_identifier)
length(un)
HeverCompare<-HeverCompare1
####comparing no. of individuals for release wild in every season####
# Released vs. Wild, Summer vs. Winter Groups
HeverSumCat<- HeverCompare%>%
  filter(season=="summer") %>%
  filter(origin=="cat")
HeverSumCat$season=droplevels(HeverSumCat$season)
HeverSumCat$released_wild=droplevels(HeverSumCat$released_wild)
HeverSumCatlocal_identifier=droplevels(HeverSumCat$local_identifier)
un<-unique(HeverSumCat$local_identifier)
length(un)

HeverSumbp<- HeverCompare%>%
  filter(season=="summer") %>%
  filter(origin=="bp")
HeverSumbp$season=droplevels(HeverSumbp$season)
HeverSumbp$released_wild=droplevels(HeverSumbp$released_wild)
HeverSumbp$local_identifier=droplevels(HeverSumbp$local_identifier)
un<-unique(HeverSumbp$local_identifier)
length(un)
#for Summer released n=5 individuals 29 days for imported
# for summer n=2 individuals 10 days for breeding program do not use

# Released vs. Wild, Summer vs. Winter Groups
HeverSumWild<- HeverCompare%>%
  filter(season=="summer") %>%
  filter(origin=="wild")
un<-unique(HeverSumWild$local_identifier)
length(un)
# for summer wild N=2/n=193 days


HeverSumWild$local_identifier=droplevels(HeverSumWild$local_identifier)
HeverSumWild$season=droplevels(HeverSumWild$season)
HeverSumWild$released_wild=droplevels(HeverSumWild$released_wild)
HeverSumWild$group_release=droplevels(HeverSumWild$group_release)
SumWild<-table(HeverSumWild$local_identifier,HeverSumWild$group_release)
addmargins(SumWild)


# Released vs. Wild, Summer vs. Winter Groups
HeverWinRelease<- HeverCompare%>%
  filter(season=="winter") %>%
  filter(origin=="bp")
un<-unique(HeverWinRelease$local_identifier)
levels(HeverCompare$origin)
length(un)




 #####Movement parameters##### 
 
  ####fly today- probability of flying#####
#adding fly today yes or no
HeverCompare$flytoday3<-cut(HeverCompare$flightHr,breaks=c(0,1,10), labels=c('0','1'))
HeverCompare$flytoday3=as.factor(HeverCompare$flytoday3)
levels(HeverCompare$flytoday3)
Did_not_fly<-HeverCompare[which(HeverCompare$flytoday3=='0'),]
Did_fly<-HeverCompare[which(HeverCompare$flytoday3=='1'),]

#statistics#

#####models for flight today####

##use age if relevant##
NullModel        =glmmTMB( flytoday3 ~ 1 +            (1|local_identifier) +(1|group_release) ,family=binomial(link="logit"), data=HeverCompare)


Release           =glmmTMB( flytoday3 ~ 1 + released_wild +         (1|local_identifier)+(1|group_release) ,family=binomial(link="logit"), data=HeverCompare)



Origin=glmmTMB( flytoday3 ~ 1 + origin+ (1|group_release)+         (1|local_identifier) ,family=binomial(link="logit"), data=HeverCompare)


Season    =glmmTMB( flytoday3 ~ 1 + season+  (1|group_release)+ (1|local_identifier) ,family=binomial(link="logit"), data=HeverCompare)

ReleaseSeason=glmmTMB(flytoday3~  released_wild+season+ (1|local_identifier)+(1|group_release),family=binomial(link="logit"), data=HeverCompare)

ReleaseCorSeason=glmmTMB(flytoday3~  released_wild+season+ (1|local_identifier)+(1|group_release),family=binomial(link="logit"), data=HeverCompare)


OriginSeason=glmmTMB(flytoday3~  origin+season+ (1|local_identifier)+(1|group_release),family=binomial(link="logit"), data=HeverCompare)

#OriginCorSeason=glmmTMB(flytoday3~  origin+season+origin*season+ (1|local_identifier)+(1|group_release),family=binomial(link="logit"), data=HeverCompare)
#cannot use last model due to only two released in summer from breeding program


#comparing aic 

ModelNames=c('NullModel','Release','Season','Origin',
             'ReleaseSeason', 'ReleaseCorSeason',  'OriginSeason')
ModelsList=list(NullModel,Release,Season, Origin,
                ReleaseSeason, ReleaseCorSeason,  OriginSeason)

#comparing aic
TableAicMix=aictab(cand.set=ModelsList,modnames=ModelNames,weights = TRUE)
print(TableAicMix)

summary(Origin)
x=simulateResiduals(Origin)
plot(x)
########emmeans####



##Top model origin

library(effects)
allEffects(Origin)
plot(allEffects(Origin))
##plotting effects##
summaryOrigin<-summary(Origin)

coef<-coefficients(Origin)

efOrigin<-effect('origin',Origin)
summary(efOrigin)
EffOrigin<-as.data.frame(efOrigin)





ggplot(EffOrigin, aes(origin, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw()




  
  

###############################################
#####check parameters for vultures that did fly######
mergMovfly<-HeverCompare%>% 
  filter(flytoday3!="0") 
mergMovfly$flytoday3=droplevels(mergMovfly$flytoday3)


levels(mergMovfly$flytoday3)
shapiro.test(mergMovfly$max_km)#not normal distribution
shapiro.test(log(mergMovfly$max_km))#not normal distribution

ggplot(data=mergMovfly, aes(x=origin, y=max_km))+ geom_boxplot() + geom_smooth(se=FALSE, method="lm")+theme_bw(base_size = 10)


#################################################3
#####model for vultures that did fly####
#can use origin, season,release vs. wild seperately but not together, season with all repeats are not significant




#without interactions without origin (only 1 bp) no group release
#####model for vultures that did fly####
#can use origin, season,release vs. wild seperately but not together, season with all repeats are not significant

NullModel=glmmTMB(max_km~  1+ (1|local_identifier),family=Gamma(link="log"), data=mergMovfly)

Season=glmmTMB(max_km~  season+ (1|local_identifier),family=Gamma(link="log"), data=mergMovfly)

Release=glmmTMB(max_km~  released_wild+ (1|local_identifier),family=Gamma(link="log"), data=mergMovfly)

ModelNames=c('NullModel','Release','Season'
             )
ModelsList=list(NullModel,Release,Season
                )
#comparing aic
TableAicMix=aictab(cand.set=ModelsList,modnames=ModelNames,weights = TRUE)
print(TableAicMix)
###################

#comes out 
summary(Release)
allEffects(Release)
x=simulateResiduals(Release)
plot(x)
summaryRelease<-summary(Release)

efRelease<-effect('released_wild',Release)
summary(efRelease)
EffRelease<-as.data.frame(efRelease)

coef<-coefficients(summaryRelease)



ggplot(EffRelease, aes(released_wild, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw()
#########################################3
####ststistics max km####
group_by(mergMovfly, released_wild) %>%
  dplyr::summarise(    count = n(),
                       mean = mean(max_km, na.rm = TRUE),
                       sd = sd(max_km, na.rm = TRUE),
                       se=sd/sqrt(count),
                       median = median(max_km, na.rm = TRUE),
                       IQR = IQR(max_km, na.rm = TRUE),
                       
                       Minimum=min(max_km, na.rm = TRUE),
                       Maximum=max(max_km,na.rm = TRUE))


###For Ns####
Heverreleased<- mergMovfly%>%
  filter(released_wild=="released") 
un<-unique(Heverreleased$local_identifier)
length(un)

HeverWild<- mergMovfly%>%
  filter(released_wild=="wild") 
un<-unique(HeverWild$local_identifier)
length(un)

ggplot(data=mergMovfly, aes(x=season, y=max_km, fill=released_wild))+ geom_boxplot() + geom_smooth(se=FALSE, method="lm")+theme_bw(base_size = 10)
ggplot(data=mergMovfly, aes(x=released_wild, y=max_km, fill=released_wild))+ geom_boxplot() + geom_smooth(se=FALSE, method="lm")+theme_bw(base_size = 10)


library(ggplot2)

##the general statistics##
group_by(mergMovfly, origin) %>%
  dplyr::summarise(    count = n(),
                       mean = mean(max_km, na.rm = TRUE),
                       sd = sd(max_km, na.rm = TRUE),
                       se=sd/sqrt(count),
                       median = median(max_km, na.rm = TRUE),
                       IQR = IQR(max_km, na.rm = TRUE),
                       
                       Minimum=min(max_km, na.rm = TRUE),
                       Maximum=max(max_km,na.rm = TRUE))

HeverCat<- mergMovfly%>%
  filter(origin=="cat") 
un<-unique(HeverCat$local_identifier)
length(un)

HeverCatSum<- HeverCat%>%
  filter(season=="summer") 
un<-unique(HeverCatSum$local_identifier)
length(un)

HeverCatWin<- HeverCat%>%
  filter(season=="winter") 
un<-unique(HeverCatWin$local_identifier)
length(un)

HeverBP<- mergMovfly%>%
  filter(origin=="bp") 
un<-unique(HeverBP$local_identifier)
length(un)

HeverBPSum<- HeverBP%>%
  filter(season=="summer") 
un<-unique(HeverBPSum$local_identifier)
length(un)

HeverBPWin<- HeverBP%>%
  filter(season=="winter") 
un<-unique(HeverBPWin$local_identifier)
length(un)

HeverWild<- mergMovfly%>%
  filter(origin=="wild") 
un<-unique(HeverWild$local_identifier)
length(un)

HeverWildWin<- HeverWild%>%
  filter(season=="winter") 
un<-unique(HeverWildWin$local_identifier)
length(un)

HeverWildSum<- HeverWild%>%
  filter(season=="summer") 
un<-unique(HeverWildSum$local_identifier)
length(un)

group_by(HeverCat, season) %>%
  dplyr::summarise(    count = n(),
                       mean = mean(max_km, na.rm = TRUE),
                       sd = sd(max_km, na.rm = TRUE),
                       se=sd/sqrt(count),
                       median = median(max_km, na.rm = TRUE),
                       IQR = IQR(max_km, na.rm = TRUE),
                       
                       Minimum=min(max_km, na.rm = TRUE),
                       Maximum=max(max_km,na.rm = TRUE))

table<-table(mergMovfly$season, mergMovfly$origin)

addmargins(table)
#Heverbp=7 individuals all sub adult
HeverBP<- HeverCompare%>%
  filter(origin=="bp") 
un<-unique(HeverBP$local_identifier)
length(un)

HeverBPSum<- HeverBP%>%
  filter(season=="summer") 
un<-unique(HeverBPSum$local_identifier)
length(un)

group_by(HeverBP, season) %>%
  dplyr::summarise(    count = n(),
                       mean = mean(max_km, na.rm = TRUE),
                       sd = sd(max_km, na.rm = TRUE),
                       se=sd/sqrt(count),
                       median = median(max_km, na.rm = TRUE),
                       IQR = IQR(max_km, na.rm = TRUE),
                       
                       Minimum=min(max_km, na.rm = TRUE),
                       Maximum=max(max_km,na.rm = TRUE))
#HeverWild n=15, 1 sub adult, 14 adult
HeverWild<- HeverCompare%>%
  filter(origin=="wild") 
un<-unique(HeverWild$local_identifier)
length(un)

HeverWildSum<- HeverWild%>%
  filter(season=="summer") 
un<-unique(HeverWildSum$local_identifier)
length(un)

group_by(HeverWild, season) %>%
  dplyr::summarise(    count = n(),
                       mean = mean(max_km, na.rm = TRUE),
                       sd = sd(max_km, na.rm = TRUE),
                       se=sd/sqrt(count),
                       median = median(max_km, na.rm = TRUE),
                       IQR = IQR(max_km, na.rm = TRUE),
                       
                       Minimum=min(max_km, na.rm = TRUE),
                       Maximum=max(max_km,na.rm = TRUE))

