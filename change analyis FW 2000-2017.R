# This script is written to utilize the change.analysis section of the spsurvey package.
# Author Jay Silvanima
# 02/11/2019
# Tony Olsen incorporated additional code during review 
# 03/01/ - 03/06/2019
# Load psurvey.analysis library
library("spsurvey")

# import the combined Exclusion and Results file that will be used for the analysis

# This is the Flowing Waters file for 2000-2017

setwd ('Z:/data analysis/Change Analysis/Cycle1 FW 2000-2017')

##Descriptive Stats tmdl basins
##Load combined file containing information from both time periods.
##File created independently from this script.

comb <- read.csv('C1C32017COMB.csv')

##comb <- read.csv('Original_Data/C1C32017COMB.csv')

################################################################
###NNC and DO Calculations
################################################################

names(comb)

comb$TN<-(comb$Kjeldahl.Nitrogen..Total..as.N+comb$Nitrate.Nitrite..Total..as.N)



### Pass=1 AND Fail=0
## Using the NNC maximum values

comb$TNcat<-ifelse((comb$TN_NNC >= comb$TN),1,0) 

comb$TPcat<-ifelse((comb$TP_NNC >= comb$Phosphorus..Total..as.P),1,0)

comb$ChlaCat<- ifelse((comb$CHLA_NNC >= comb$Chlorophyll.A..Monochromatic) ,0,1)

#load lubridate
#need the date in the yyyy-mm-dd format

library(lubridate)

# Generate the files used in change.analysis

# do equal area projection for variance estimation
tmp <- marinus(comb$LAT_DEC, comb$LONG_DEC)
comb$xmarinus <- tmp[,'x']
comb$ymarinus <- tmp[,'y']

nr<-nrow(comb)

## The following is the targetsize in miles for the flowing waters in 
## Florida's 29  drainage basins. Excludes canals within the Northwest Florida
## and Suwannee River Water Management Districts as these canals were never
## included in our statistical surveys due to their limited extent.
## The resource extents for three basins are exluded because no sites were
## selected during cycle one from these basins.

targetsizeGroup <- c("Apalachicola - Chipola" = 1814.74,
                "Caloosahatchee" = 249.40,
                ##"Charlotte Harbor" = 52.45,
                "Choctawhatchee - St. Andrew" = 2681.42,
                "Everglades" = 456.09,
                "Everglades West Coast" = 326.23,
                "Fisheating Creek" = 248.35,
                ## "Florida Keys" =  0.09,
                "Indian River Lagoon" = 117.09,
                "Kissimmee River" = 736.68,
                "Lake Okeechobee" = 157.98,
                "Lake Worth Lagoon - Palm Beach Coast" = 215.41,
                "Lower St. Johns" = 1268.08,
                "Middle St. Johns" = 616.44,
                "Nassau - St. Marys" = 520.24,
                "Ochlockonee - St. Marks" = 1474.90,
                "Ocklawaha" = 380.42,
                "Pensacola" = 2325.68,
                "Perdido" = 272.02,
                "Sarasota Bay - Peace - Myakka" = 1637.37,
                "Southeast Coast - Biscayne Bay" = 317.48,
                "Springs Coast" = 147.83,
                "St. Lucie - Loxahatchee" = 216.34,
                "Suwannee" = 1969.26,
                "Tampa Bay" = 146.52,
                "Tampa Bay Tributaries" = 1011.69,
                ##"Upper East Coast" = 121.93,
                "Upper St. Johns" = 537.52,
                "Withlacoochee" = 436.31)


## The following is the targetsize in miles for the flowing waters in 
## Florida's 6 reporting units per 2017 coverage. Excludes canals within the Northwest Florida
## and Suwannee River Water Management Districts as these canals were never
## included in our statistical surveys due to their limited extent.

##targetsizeZone <- c("ZONE 1"=8629.32 ,"ZONE 2"=1901.23,
##                "ZONE 3"=3615.71,"ZONE 4"=3547.18,
##                "ZONE 5"=1455.62,"ZONE 6"=1384.50)

## The following is the targetsize in miles for the flowing waters in 
## Florida's 6 reporting units. Excludes canals within the Northwest Florida
## and Suwannee River Water Management Districts as these canals were never
## included in our statistical surveys due to their limited extent.
## The resource extents for three TMDL basins (Charlotte Harbor, Florida Keys,
## and Upper East Coast are exluded because no sites were
## selected during cycle one from these basins.

targetsizeZone <- c("ZONE 1"=8629.32 ,"ZONE 2"=1901.23,
                    "ZONE 3"=3497.13,"ZONE 4"=3509.88,
                    "ZONE 5"=1438.75,"ZONE 6"=1384.50)

# Look at Cycle, Zone and Group
addmargins(table(comb$GROUP_NAME, comb$C3_REPORTING_UNIT, comb$CYCLE,
                 useNA = 'ifany'))
# have 768 observations from Cycle 1 and 711 from Cycle 2. No missing values.
# All sites are target and sampled.  Assume this is all sites evaluated.
# Do not know if this is true. Makes assumption that sample frame matches the
# target population and that all sites could be sampled.

comb$ZoneGroup <- factor(paste(comb$C3_REPORTING_UNIT, comb$GROUP_NAME, sep = "_"))
levels(comb$ZoneGroup)
# note that Suwannee Group is split between Zone 1 and Zone 2. Is this an issue?

##### Weight adjustment

# Cycle 1: Adjust weights to 2017 framesize with 29 basins used in the design
comb$wgt_1 <- adjwgt(comb$CYCLE == 1, comb$WGT, comb$GROUP_NAME,
                     framesize = targetsizeGroup)

# Cycle 3: Adjust weights to 2017 framesize with 6 zones used in the design
comb$wgt_3 <- adjwgt(comb$CYCLE == 3, comb$WGT, comb$C3_REPORTING_UNIT,
                     framesize = targetsizeZone)

addmargins(table(comb$wgt_1 == 0, comb$wgt_3 == 0))

comb$wgt <- comb$wgt_1 + comb$wgt_3

# look at weight sums for Cycle 1
tmp <- tapply(comb$wgt_1, list(comb$GROUP_NAME, comb$C3_REPORTING_UNIT), sum)
tmp[is.na(tmp)] <- 0
round(addmargins(tmp), 1)
# Note that Charlotte Harbor and Upper East Coast do not have sites.
# Also Group length totals 20281.5

# look at weight sums for Cycle 3
tmp <- tapply(comb$wgt_3, list(comb$GROUP_NAME, comb$C3_REPORTING_UNIT), sum)
tmp[is.na(tmp)] <- 0
round(addmargins(tmp), 1)
# Zone Length totals 20360.8 which is close to length for Cycle 1.
# Note differences in GIS layer scale and how extracted account for differences.
# Not unusual for this to happen.

## write.csv(comb,file='comb.csv')

#Now to get the dataframe for the sites

mysites <- data.frame(siteID=comb$PK_RANDOM, 
                      Survey1=ifelse(comb$CYCLE=="1",TRUE,FALSE),
                      Survey2=ifelse(comb$CYCLE=="3", TRUE,FALSE))

##################################################################
### Define subpopulations (reporting units) and create design for
### Florida's TMDL basins
##################################################################

mysubpop <- data.frame(siteID=comb$PK_RANDOM,
                       Combined=rep("Basins Combined", nr), 
                       Zone=comb$C3_REPORTING_UNIT, # do Zone estimates as well
                       Basin=comb$GROUP_NAME)

## mydsgn <- data.frame(siteID=comb$PK_RANDOM, 
##                     wgt=comb$wgt,
##                     xcoord=comb$xmarinus ,
##                     ycoord=comb$ymarinus ,
##                     stratum=comb$GROUP_NAME)
# since use local nhbd variance estimator (in most cases) likely
# okay to ignore stratum so that nbhd uses nearby sites regardless if they
# are in neighboring GROUP or ZONE. have very small sample sizes otherwise.
mydsgn <- data.frame(siteID=comb$PK_RANDOM, 
                     wgt=comb$wgt,
                     xcoord=comb$xmarinus ,
                     ycoord=comb$ymarinus)

########################################################################
### Create categorical and continuous distribution of specific indicator
### to run change analysis on
########################################################################

### Total Nitrogen

mydata.cat <- data.frame(siteID=comb$PK_RANDOM, 
                         TNcategory=comb$TNcat)

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          TN=comb$TN)

#### Total Nitrogen change.analysis

CategoryTNCAT <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                   data.cat=mydata.cat,data.cont=mydata.cont, test=c("mean","median"))

########################################################################
##### write out total nitrogen change analysis results
########################################################################

CategoryTNCAT

write.csv(CategoryTNCAT$catsum,file='TNcatsumC1C3.csv')
write.csv(CategoryTNCAT$contsum_mean,file='TNcontsummeanC1C3.csv')
write.csv(CategoryTNCAT$contsum_median,file='TNcontsummedianC1C3.csv')

### Nitrate + Nitrite as N (NOx)


mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          NOx=comb$Nitrate.Nitrite..Total..as.N)

#### NOx change.analysis

CategoryNOx <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                               data.cont=mydata.cont, test=c("mean","median"))

########################################################################
##### write out NOx change analysis results
########################################################################

CategoryNOx

write.csv(CategoryNOx$contsum_mean,file='NOxcontsummeanC1C3.csv')
write.csv(CategoryNOx$contsum_median,file='NOxcontsummedianC1C3.csv')


### TKN


mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          TKN=comb$Kjeldahl.Nitrogen..Total..as.N)

#### TKN change.analysis

CategoryTKN <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                               data.cont=mydata.cont, test=c("mean","median"))

########################################################################
##### write out TKN change analysis results
########################################################################

CategoryTKN

write.csv(CategoryTKN$contsum_mean,file='TKNcontsummeanC1C3.csv')
write.csv(CategoryTKN$contsum_median,file='TKNcontsummedianC1C3.csv')


#### Total Phosphorus

mydata.cat <- data.frame(siteID=comb$PK_RANDOM, 
                         TPcategory=comb$TPcat)

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          TP=comb$Phosphorus..Total..as.P)

#### Total Phosphorus change.analysis

CategoryTPCAT <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                 data.cat=mydata.cat,data.cont=mydata.cont, test=c("mean","median"))


########################################################################
##### write out total phosphorus change analysis results
########################################################################

CategoryTPCAT

write.csv(CategoryTPCAT$catsum,file='TPcatsumC1C3.csv')
write.csv(CategoryTPCAT$contsum_mean,file='TPcontsummeanC1C3.csv')
write.csv(CategoryTPCAT$contsum_median,file='TPcontsummedianC1C3.csv')

#### Chlorophyll A

mydata.cat <- data.frame(siteID=comb$PK_RANDOM, 
                         Chlacategory=comb$ChlaCat)

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          ChlA=comb$Chlorophyll.A..Monochromatic)

#### Chlorophyll A corrected change.analysis
 

CategoryChlA <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                data.cat=mydata.cat, data.cont=mydata.cont, test=c("mean","median"))

########################################################################
##### write out chlorophyll a corrected change analysis results
########################################################################

CategoryChlA

write.csv(CategoryChlA$catsum,file='ChlAcatsumC1C3.csv')
write.csv(CategoryChlA$contsum_mean,file='ChlAcontsummeanC1C3.csv')
write.csv(CategoryChlA$contsum_median,file='ChlAcontsummedianC1C3.csv')

#### Total Chloride

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          CL=comb$Chloride..Total)

#### Total Chloride change.analysis

CategoryCL <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                data.cont=mydata.cont, test=c("mean","median"))

########################################################################
##### write out Chloride corrected change analysis results
########################################################################

CategoryCL

write.csv(CategoryCL$contsum_mean,file='CLcontsummeanC1C3.csv')
write.csv(CategoryCL$contsum_median,file='CLcontsummedianC1C3.csv')

# Specific conductance field
# Create Specific Conductance categories of 0, 1000, 10000 microsiemens 

SCCat <- cut(comb$Specific.Conductance..Field, breaks=c(0,1000,10000), include.lowest=TRUE)
comb$SCCat <- SCCat
comb$SCCat <- as.factor(comb$SCCat)

#### Specific conductance change.analysis

mydata.cat <- data.frame(siteID=comb$PK_RANDOM, 
                         speccond=comb$SCCat)

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          SCField=comb$Specific.Conductance..Field)

CategorySpCond <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                  data.cat=mydata.cat,data.cont=mydata.cont, test=c("mean","median"))


########################################################################
##### write out Specific Conductance change analysis results
########################################################################

CategorySpCond

write.csv(CategorySpCond$catsum,file='SpCondcatsumC1C3.csv')
write.csv(CategorySpCond$contsum_mean,file='SpCondMeanC1C3.csv')
write.csv(CategorySpCond$contsum_median,file='SpCondMedianC1C3.csv')



# pH field

#### pH change.analysis

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          pHField=comb$pH..Field)

CategorypH <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                  data.cont=mydata.cont, test=c("mean","median"))


########################################################################
##### write out pH change analysis results
########################################################################

CategorypH

write.csv(CategorypH$contsum_mean,file='pHCondMeanC1C3.csv')
write.csv(CategorypH$contsum_median,file='pHCondMedianC1C3.csv')


# Temperature

#### Temperature change.analysis

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          TempField=comb$Water.Temperature)

CategoryTemp <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                  data.cont=mydata.cont, test=c("mean","median"))


########################################################################
##### write out temperature change analysis results
########################################################################

CategoryTemp

write.csv(CategoryTemp$contsum_mean,file='TempCondMeanC1C3.csv')
write.csv(CategoryTemp$contsum_median,file='TempCondMedianC1C3.csv')

# Dissolved Oxygen

#### Dissolved Oxygen change.analysis

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          DO=comb$Oxygen..Dissolved..Field)

CategoryDO <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                                data.cont=mydata.cont, test=c("mean","median"))


########################################################################
##### write out Dissolved Oxygen analysis results
########################################################################

CategoryDO

write.csv(CategoryDO$contsum_mean,file='DOCondMeanC1C3.csv')
write.csv(CategoryDO$contsum_median,file='DOCondMedianC1C3.csv')



# Total Organic Carbon

#### Total Organic Carbon change.analysis

mydata.cont <- data.frame(siteID=comb$PK_RANDOM, 
                          TOC=comb$Organic.Carbon..Total)

CategoryTOC <- change.analysis(sites=mysites, subpop=mysubpop, design=mydsgn,
                              data.cont=mydata.cont, test=c("mean","median"))


########################################################################
##### write out Dissolved Oxygen analysis results
########################################################################

CategoryTOC

write.csv(CategoryTOC$contsum_mean,file='TOCCondMeanC1C3.csv')
write.csv(CategoryTOC$contsum_median,file='TOCCondMedianC1C3.csv')


##########################################################################
########### exploratory data analyses  Added 06/30/2020
##########################################################################

library(EnvStats)
library(MASS)
library(fitdistrplus)


Fn <- ecdf(comb$Water.Temperature)
Fn(comb$Water.Temperature)

plot(comb$Water.Temperature, , ylab="Fn(comb$Water.Temperature)", verticals = FALSE,
     col.01line = "gray70", pch = 19)

plot(Fn, xlab = 'Sample Quantiles of water temperature', ylab = '', 
     main = 'Empirical Cumluative Distribution Lake Water Temperature FL')

plot(ecdf(comb[comb$CYCLE=="1",]$pH..Field))

plot(ecdf(comb[comb$CYCLE=="1",]$pH..Field),col='green')
lines(ecdf(comb[comb$CYCLE=="3",]$pH..Field),col='red')

plot(ecdf(comb[comb$CYCLE=="1",]$Specific.Conductance..Field),col='green')
lines(ecdf(comb[comb$CYCLE=="3",]$Specific.Conductance..Field),col='red')

plot(ecdf(comb[comb$CYCLE=="1",]$Chlorophyll.A..Monochromatic.),col='green')
lines(ecdf(comb[comb$CYCLE=="3",]$Chlorophyll.A..Monochromatic.),col='red')

plot(ecdf(comb[comb$CYCLE=="1",]$Water.Temperature),col='green')
lines(ecdf(comb[comb$CYCLE=="3",]$Water.Temperature),col='red')



CategoryTemp$contsum_mean

