###################################################################
#BIOMOD_FormatingData Initialise the datasets for usage in biomod2
###################################################################

#open required libraries

library(sp)
library(raster)
library(rgdal)
library(dismo)
library(ecospat)
library(SDMTools)
library(biomod2)
library(earth)
library(gam)
library(mgcv)
library(usdm)
library(ape)
library(ncf)
library(AICcmodavg)
library(foreign)
library(spThin)
library(ENMeval)
library(dplyr) #filter function for spatial thinning
library(randomForest) #for random forest model analysis
library(ggplot2) #better plots
library(viridis) #color pallets for ggplots
library(reshape2) #help reorganize tables for ggplots
library(gridExtra) #help with ploting with ggplot2

## Load required data from GitHub

t <- tempdir()
setwd(t) # set temporal wd

# download data from github
url <- "https://github.com/MJCSRibeiro/Lanius_biomod2/archive/refs/heads/main.zip"
download.file(url,destfile="Lanius_biomod2-main.zip")
unzip(zipfile="Lanius_biomod2-main.zip")

setwd(dir = file.path(t,"Lanius_biomod2-main"))
newfolder <- "Results"
dir.create(file.path(dirname(paste(getwd(),"/Maxent",sep="")),newfolder))


#.................................................................................................
## Data:
#.................................................................................................

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

set.seed(27) # set seed to obtain same results

# OPEN THE DATA


sp<-"Lanius collurio" # set name of target species


# Original Pseudoabsences study area
lcdata<-raster(paste(getwd(),"/Final qgis products/manelFULLdata.asc", sep="")) #call raster file of full data created in manel_data_geographics script
plot(lcdata) #plot data
lcdata<-as.data.frame(rasterToPoints(lcdata)) #save data in data frame
colnames(lcdata)<-c('x','y','pres')
lcdata$pres[lcdata$pres==0]<-NA #transform absences into pseudoabsences
dados<-"FULL"
condicao<-"original"

# Breeding Range Polygon buffer 5 degree (500km)
lcdata<-raster(paste(getwd(),"/Final qgis products/manelpoly5.asc", sep="")) #call raster file of full data created in manel_data_geographics script
plot(lcdata) #plot data
lcdata<-as.data.frame(rasterToPoints(lcdata)) #save data in data frame
colnames(lcdata)<-c('x','y','pres')
lcdata$pres[lcdata$pres==0]<-NA #transform absences into pseudoabsences
dados<-"FULL"
condicao<-"poly5"

# Breeding Range Polygon buffer 10 degree (1000km)
lcdata<-raster(paste(getwd(),"/Final qgis products/manelpoly10.asc", sep="")) #call raster file of full data created in manel_data_geographics script
plot(lcdata) #plot data
lcdata<-as.data.frame(rasterToPoints(lcdata)) #save data in data frame
colnames(lcdata)<-c('x','y','pres')
lcdata$pres[lcdata$pres==0]<-NA #transform absences into pseudoabsences
dados<-"FULL"
condicao<-"poly10"

# Rectangular area max and min long/lat coords of presences
lcdata<-raster(paste(getwd(),"/Final qgis products/manelsquared.asc", sep="")) #call raster file of full data created in manel_data_geographics script
plot(lcdata) #plot data
lcdata<-as.data.frame(rasterToPoints(lcdata)) #save data in data frame
colnames(lcdata)<-c('x','y','pres')
lcdata$pres[lcdata$pres==0]<-NA #transform absences into pseudoabsences
dados<-"FULL"
condicao<-"squared"


## Iberian Peninsula independent data
lcpi.ras<-raster(paste(getwd(),"/Final qgis products/laniusPIindependentcensus.asc", sep="")) #call raster of PI census
plot(crop(lcpi.ras,extent(-10, 2.7, 34.8, 44.5))) #plot PI data

pi.dt<-as.data.frame(rasterToPoints(lcpi.ras)) #extract to data frame to correct extent issue
ext<-extent(-12, 105, 9, 80) #extent of all other rasters
rr<-raster(ext,res=0.5) #base raster
rr<-rasterize(pi.dt[,1:2], rr, field=pi.dt[,3], fun="max", na.rm=T) #create new raster with corrected extent
NAvalue(rr)<--9999
lcpi.ras<-rr; rm(rr) #replace lc pi raster with new raster with corrected extent

lcpi<-as.data.frame(rasterToPoints(lcpi.ras)) #save PI census in data frame
colnames(lcpi)<-c('x','y','pres') #change column names of PI census data frame

# Remove PI from calibration data - Original Pseudoabsences study area
lcdata<-raster(paste(getwd(),"/Final qgis products/manelFULLdata.asc", sep="")) #call raster file of full data created in manel_data_geographics script
plot(lcdata) #plot data
lcdata<-as.data.frame(rasterToPoints(lcdata)) #save data in data frame
colnames(lcdata)<-c('x','y','pres')
lcdata<-merge(lcdata,lcpi, by=c('x','y'), all.x=T) #merge full range and PI into one data frame
lcdata<-lcdata[-which(!is.na(lcdata$pres.y)),] #remove rows with PI data
lcdata<-lcdata[,-4] #remove PI pres info column
colnames(lcdata)<-c('x','y','pres')
lcdata$pres[lcdata$pres==0]<-NA #transform absences into pseudoabsences
dados<-"noPI"
condicao<-"original"

#................................................................................
# Environmental variables to be used
#................................................................................
bias=crop(raster(paste(getwd(),"/current/bias.asc", sep="")),extent(-12, 105, 9, 80)) # density of same niche species
tmin =crop(raster(paste(getwd(),"/current/TMIN.asc", sep="")),extent(-12, 105, 9, 80)) # mean temperature of coldest month
tmax=crop(raster(paste(getwd(),"/current/TMAX.asc", sep="")),extent(-12, 105, 9, 80)) # mean temperature of warmest month
tmean=crop(raster(paste(getwd(),"/current/TMEAN.asc", sep="")),extent(-12, 105, 9, 80)) # annual mean temperature
prec=crop(raster(paste(getwd(),"/current/PREC.asc",sep="")),extent(-12, 105, 9, 80)) # annual precipitation
precsd=crop(raster(paste(getwd(),"/current/PRECSD.asc",sep="")),extent(-12, 105, 9, 80)) # precipitation seasonality (standard deviation)
tsd=crop(raster(paste(getwd(),"/current/TSD.asc",sep="")),extent(-12, 105, 9, 80)) # temperature seasonality (standard deviation)

## Round variables (?)

bias<-round(bias, digits = 0) #bias density is rounded to integer

tmin<-round(tmin, digits = 2) #temperatures are rounded to 2 decimals
tmean<-round(tmean, digits = 2)
tmax<-round(tmax, digits = 2)
tsd<-round(tsd, digits = 2)

prec<-round(prec, digits = 0) #annual precipitation is rounded to integer
precsd<-round(precsd, digits = 2) #coefficient of variation of precipitation as a percentage is rounded to 2 decimals (might change?)

bias@data@names<-'bias'; tmin@data@names<-'tmin'; tmean@data@names<-'tmean'; tmax@data@names<-'tmax'; tsd@data@names<-'tsd'
prec@data@names<-'prec'; precsd@data@names<-'precsd'

# Reclassify landuses (using landcover data with 9 landuse categories as starting point)
#old classification: 1=forest, 2=shrubland, 3=savannah, 4=grassland, 5=wetland, 6=agriculture, 7=urban, 8=ice, 9=barren
#new classification: 1=agriculture, 2=shrub, 3=open habitats (savanna, grassland), 4=forest, 5=other (ice+hot desert+urban+wetland)

land=crop(raster(paste(getwd(),"/current/land.asc", sep="")),extent(-12, 105, 9, 80))
land=as.factor(land)

# use land as a factor variable
land<-reclassify(land, matrix(c(5.5,6,5, #level 5 is agriculture
                                1.5,2,2, #level 2 is shrubland
                                2.5,4,3, #level 3 is open
                                0,1,4, #level 4 is forest
                                4.5,5,1,
                                6.5,9,1), #level 1 is other
                              ncol = 3, byrow = T)) #reclassifies land variable from old values (first 2 columns in the matrix) to new value (third column)

land<-as.factor(land)
rat<-levels(land)[[1]] #gets matrix from land levels
rat$landcover<-c('other','shrub','open','forest','agri') #renames factors to facilitate interpretation
levels(land)<-rat; land #renames levels for land
plot(land)


### Check variables
#myvars<-as.data.frame(rasterToPoints(stack(tmax, tsd, prec, precsd, bias, agri, shrub, open, forest, other))) #land as binary
myvars<-as.data.frame(rasterToPoints(stack(tmax, tsd, prec, bias, land))) #land as factor and no precsd
myebba2<-raster(paste(getwd(),"/Final qgis products/manelEBBA2.asc",sep=""))
myebba2<-as.data.frame(rasterToPoints(myebba2))
mycompletedata<-merge(lcdata,myvars,by=c('x','y'), all.x=T) #merge lc occurrence to variables in each coordinate
mycompletedata<-merge(mycompletedata,myebba2,by=c('x','y'),all.x=T) #merge ebba2 data to data frame
#mycompletedata<-mycompletedata[complete.cases(mycompletedata[,-c(3,14)]),] #remove rows where explanatory variables have NA
mycompletedata<-mycompletedata[complete.cases(mycompletedata[,-c(3,9)]),] #land factor and no precsd

bias_data<-as.data.frame(rasterToPoints(bias)) #create bias data frame
#bias_data<-merge(bias_data,mycompletedata[,c(1,2,14)], by=c('x','y'), all.x=T) #merge with ebba2 data with land binary
bias_data<-merge(bias_data,mycompletedata[,c(1,2,9)], by=c('x','y'), all.x=T) #merge with ebba2 data with land factor
bias_data$bias[which(!is.na(bias_data$manelEBBA2))]<-10 #maximize sampling effort in areas with ebba2 sampling
mycompletedata$bias[which(!is.na(mycompletedata$manelEBBA2))]<-10 #apply same to mycompletedata table

ext<-extent(-12, 105, 9, 80)
rrrr<-raster(ext, res=0.5)
bias<-rasterize(bias_data[,c('x','y')], rrrr, field=bias_data$bias) #create new sampling bias grid where EBBA2 is accounted for
plot(bias)
bias@data@names<-'bias'
rm(rrrr,bias_data) #remove unnecessary objects


#### verify correlation and multicollinearity

#correlation graph
ecotable<-as.data.frame(rasterToPoints(myExpl)); View(ecotable) # transformar os rasters das variaveis num data frame
cord<-data.frame()
h<-1
v<-c()
for(j in 1:length(ecotable[1,])){
  for(i in 1:length(ecotable[,1])){
    if(is.na(ecotable[i,j])==T & i%in%v==F){
      v<-append(v,c(i))
      cord[h,1]<-ecotable[i,1]
      cord[h,2]<-ecotable[i,2]
      h<-h+1
    }
  }
} # removes NA's
colnames(cord)<-c('x','y')
ecotable<-ecotable[-v,]

ecospat.cor.plot(ecotable[,c(3:length(ecotable[1,]))])
ecospat.cor.plot(ecotable[,c(3,4,5,6,7)])
ecospat.cor.plot(ecotable[,c(8,9)])

pairs(stack(tmin,tmax,tmean,tsd,prec,precsd,bias))
pairs(stack(tmax,tsd,precsd,land))

layerStats(stack(tmin,tmax,tmean,tsd,prec,precsd, land), stat = 'pearson', na.rm = T)
layerStats(stack(tmax,tsd,prec,precsd,land), stat = 'pearson', na.rm = T)


#variance inflation factor (VIF)

vifstep(stack(tmin,tmax,tmean,tsd,prec,precsd,bias)) # necessary to remove TMIN and TMEAN from model (and maybe PRECSD if transformed)
vifcor(stack(tmin,tmax,tmean,tsd,prec,precsd,bias), th = 0.7 ) # remove variables with high VIF score and above correlation threshold

vifcor(stack(tmax,tsd,prec,precsd,bias), th=0.5) # maybe run models with only these three variables (DISCUSS)

boxplot(stack(tmin,tmax,tmean,tsd),col='lightblue',ylim=c(-4,4),xlim=c(0,5),main='Temperature Boxplots')
hist(land, col = 'lightblue', main = 'Land cover', xlab = 'Levels')


### Spatial thinning

## by using filter function quantiles and 'spthin' package with "dplyr" to filter

mypresences<-mycompletedata[,c('x','y','pres','manelEBBA2')]
mypresences$ebba<-rep(0) #create ebba column that represents 1 where ebba data exists and 0 where it doesn't
mypresences$ebba[which(!is.na(mypresences$manelEBBA2))]<-1
mypresences$thin<-rep(1) #create thin placeholder column for thinning

quant_thresh<-0.05 #exclude coords that are in the low and top percent, avoiding edges from being thinned

datathin<-as.data.frame(filter(mypresences[which(mypresences$ebba==1),], (mypresences[which(mypresences$ebba==1),2] >= quantile(mypresences[which(mypresences$ebba==1),2], quant_thresh) 
                                                          & mypresences[which(mypresences$ebba==1),2] <= quantile(mypresences[which(mypresences$ebba==1),2], 1 - quant_thresh)
                                                          & mypresences[which(mypresences$ebba==1),1] >= quantile(mypresences[which(mypresences$ebba==1),1], quant_thresh) 
                                                          & mypresences[which(mypresences$ebba==1),1] <= quantile(mypresences[which(mypresences$ebba==1),1], 1 - quant_thresh)))) #filter coords to be thinned
edge_points<-as.data.frame(filter(mypresences[which(mypresences$ebba==1),c('x','y')], (mypresences[which(mypresences$ebba==1),2] <= quantile(mypresences[which(mypresences$ebba==1),2], quant_thresh)
                                                               | mypresences[which(mypresences$ebba==1),2] >= quantile(mypresences[which(mypresences$ebba==1),2], 1 - quant_thresh)
                                                               | mypresences[which(mypresences$ebba==1),1] <= quantile(mypresences[which(mypresences$ebba==1),1], quant_thresh) 
                                                               | mypresences[which(mypresences$ebba==1),1] >= quantile(mypresences[which(mypresences$ebba==1),1], 1 - quant_thresh)))) #filter edge areas
edge_points$edge<-rep(1) #identify edge points

thin_data<-thin(loc.data = datathin[datathin$ebba==1,  #only ebba2 presence data is thinned
                                      c('x','y','thin')],
                lat.col = "y",
                long.col = 'x',
                spec.col = 'thin',
                thin.par = 50, #distance in km to separate each occ
                reps = 1, #number of replicates of thinned data
                locs.thinned.list.return = T,
                write.files = F,
                write.log.file = F)

datathin <- thin_data[[which.max(sapply(thin_data, nrow))]]
datathin$datathin<-rep(1) #create new presence column
colnames(datathin)<-c('x','y','datathin') #rename to be able to merge


mypresences<-merge(mypresences,datathin,by=c('x','y'),all.x=T) #merge columns with new thinned presence data
mypresences<-merge(mypresences,edge_points,by=c('x','y'),all.x=T) #merge columns to identify edge data
mypresences<-mypresences[-which(is.na(mypresences$datathin) #remove NA from thinned data so ebba2 data isn't turned to PA
                                & mypresences$ebba==1 #only data sampled in ebba2
                                & mypresences$manelEBBA2==1 #only presences
                                & is.na(mypresences$edge)) #coords that are not in edge areas
                         ,] #make thinned pres


#####

#stack all variables
# myExpl = stack(tmax, tsd, prec, precsd, bias, agri, shrub, open, forest, other) #stack all layers
# names(myExpl)<-c("tmax",  "tsd",  "prec", "precsd", "bias", "agri", "shrub", "open", "forest", "other") #corresponding layer names
# NAvalue(myExpl)=-9999
# plot(myExpl)

# land as factor
#myExpl = stack(tmax, tsd, prec, precsd, bias, land) #stack all layers
#names(myExpl)<-c("tmax",  "tsd",  "prec", "precsd", "bias", "land") #corresponding layer names
#NAvalue(myExpl)=-9999
#plot(myExpl)

# only bioclimatic variables
#myExpl = stack(tmax, tsd, prec, precsd) #junta todas as layers num só
#names(myExpl)<-c("tmax", "tsd", "prec", "precsd") #muda os nomes das layers
#NAvalue(myExpl)=-9999
#plot(myExpl)

# no precsd variable
myExpl = stack(tmax, tsd, prec, bias, land) #stack all layers
names(myExpl)<-c("tmax",  "tsd",  "prec", "bias", "land") #corresponding layer names
NAvalue(myExpl)=-9999
plot(myExpl)



#.....................................
#### MODELLING ####
#.....................................

# 1. Formatting Data
length(mypresences$pres) #number of data points
sum(mypresences$pres, na.rm = T) #number of presences in occurrence data
sum(is.na(mypresences$pres)) #number of pseudo-absences
myBiomodData <- BIOMOD_FormatingData(resp.var = mypresences$pres, 
                                     expl.var = myExpl, #explanatory variables to be used
                                     resp.xy = mypresences[,c('x','y')], #coordinates of presences
                                     eval.resp.var = lcpi$pres, #PI census occurrence data
                                     eval.resp.xy = lcpi[,c('x','y')], #PI census coordinates
                                     resp.name = sp,
                                     PA.nb.rep=1, #create a unique set of pseudo absences
                                     PA.nb.absences = 10000, #a total of 10000 pseudo absences were selected
                                     PA.strategy = "random" #how to create pseudo absences
                                     )

plot(myBiomodData) #mostra os pontos usados, assim como as pseudo-ausencias
plot(myBiomodData@data.mask$input_data) #pontos usados
plot(myBiomodData@data.mask$PA1) #pseudo-ausencias
plot(crop(myBiomodData@data.mask$PA1,extent(-10, 2.7, 34.8, 44.5))) #area de interesse


# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips =                                    
                                           list(path_to_maxent.jar = paste(getwd(),'/Maxent',sep=""),
                                                memory_allocated = 512,
                                                background_data_dir = 'default',
                                                #maximumbackground = 'default',
                                                #maximumiterations = 200,
                                                visible = FALSE,
                                                linear = TRUE,
                                                quadratic = T,
                                                product = T,
                                                threshold = T,
                                                hinge = T,
                                                lq2lqptthreshold = 80,
                                                l2lqthreshold = 10,
                                                hingethreshold = 15,
                                                beta_threshold = -1,
                                                beta_categorical = -1,
                                                beta_lqp = -1,
                                                beta_hinge = -1,
                                                defaultprevalence = 0.5),
                                         GLM=NULL,
                                         RF=NULL,
                                         GBM=NULL)

# 3. Doing Modelisation

#models using subsample (10 replicates  with 70% of points for calibration)
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                     models = c("MAXENT.Phillips","GLM", "GBM", "RF"),
                                     models.options = myBiomodOption,
                                     NbRunEval=10,#NUMBER OF REPLICATES
                                     DataSplit=70,#PERCENTAGE OF THE DATA TO CALIBRATE THE MODEL
                                    Prevalence=0.5,#WEIGHT FOR PRESENCES
                                     models.eval.meth = c('TSS', 'ROC', 'KAPPA'),
                                     do.full.models = F,
                                    #Yweights = myPresPAdf_full$weight,
                                     rescal.all.models=F,
                                     modeling.id='test',
                                    VarImport = 1,
                                    SaveObj = T)
#help("BIOMOD_Modeling")


#models using full data
myBiomodModelOut_full <- BIOMOD_Modeling(myBiomodData,
                                     models = c('MAXENT.Phillips', 'GLM', 'GBM', 'RF'),
                                     models.options = myBiomodOption,
                                     NbRunEval=1,#NUMBER OF REPLICATES
                                     DataSplit=100,#PERCENTAGE OF THE DATA TO CALIBRATE THE MODEL
                                     Prevalence=0.5,#WEIGHT FOR PRESENCES
                                     models.eval.meth = c('TSS','ROC', 'KAPPA'),
                                     do.full.models = T,
                                     #Yweights = myPresPAdf$weight,
                                     rescal.all.models=F,
                                     modeling.id='test',
                                     VarImport = 1, #CHECK VARIABLE IMPORTANCE ('0' FOR TURN OFF)
                                     SaveObj = T)


##################

#to save model performance
slot=slot(myBiomodModelOut@models.evaluation, "val")
str(slot)
slot.full=slot(myBiomodModelOut_full@models.evaluation, "val")
str(slot.full)

table=as.data.frame(matrix(NA, ncol=8, nrow=132))
colnames(table)=c("model", "run", "metric", "test", "eval", "cutoff", "sensitivity", "specificity")
table[,1]=c(rep(c(rep("Maxent",3),rep("GLM", 3), rep("GBM",3), rep("RF",3)), 11))
table[,2]=c(rep(1,12), rep(2,12), rep(3,12), rep(4,12), rep(5,12), rep(6,12), rep(7,12), rep(8,12), rep(9,12), rep(10,12), rep("full",12))
table[,3]=c(rep(rep(c("TSS", "ROC", "KAPPA"),10),4), rep(c("TSS", "ROC", "KAPPA"),4))

for(i in 0:40){
  table[(1+(i*3)):(3+(i*3)),4:8]=slot[(1+(i*15)):(15+(i*15))]
}
 
for(i in 0:3){
  table[(121+(i*3)):(123+(i*3)),4:8]=slot.full[(1+(i*15)):(15+(i*15))]
}

View(table)

#save in a file

write.table(table, file = paste(getwd(),"/Results/", dados, condicao, ".txt", sep=""), col.names=T,  row.names=F, sep="\t")


#representação gráfica dos modelos em relação às variáveis

myglm<-BIOMOD_LoadModels(myBiomodModelOut, models='GLM')
myRespPlot<-response.plot2(
  models = myglm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

#dim(myRespPlot) ; dimnames(myRespPlot)

mymaxent<-BIOMOD_LoadModels(myBiomodModelOut, models='MAXENT.Phillips')
myRespPlot<-response.plot2(
  models = mymaxent,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

mygbm<-BIOMOD_LoadModels(myBiomodModelOut, models='GBM')
myRespPlot<-response.plot2(
  models = mygbm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

myrf<-BIOMOD_LoadModels(myBiomodModelOut, models='RF')
myRespPlot<-response.plot2(
  models = myrf,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))


###-------------------------------------------

# verificar as importancias das variáveis - NOTA: poderá ser necessário normalizar valores das importancias, para podermos comparar umas com as outras
# NOTA: os valores das importancia são a média de 1-cor(Var.x,Var.y) para todas as variaveis, e são meramente informativos, não sendo valores absolutos da contribuição das variáveis para os modelos
impor.variab<-t(as.data.frame(get_variables_importance(myBiomodModelOut))) #table for variable importance of split models
impor.variab
impor.variab.full<-t(as.data.frame(get_variables_importance(myBiomodModelOut_full))) #table for variable importance of full models
impor.variab.full

write.table(impor.variab, file=paste(getwd(),"/Results/",dados,condicao,"_variable_importance.txt", sep=""), col.names=T, row.names = F, sep="\t")
write.table(impor.variab.full, file=paste(getwd(),"/Results/",dados,condicao,"_variable_importance_FULL.txt", sep=""), col.names=T, row.names = F, sep="\t")


###-------------------------------------------


# 4. Projection on different environmental conditions
scenario=c("current_scn", "2050_A1B_MR", "2050_A1B_CS", "2050_A2_MR", "2050_A2_CS")
scn<-1
m_name<-c('Maxent','GLM','GBM','RF')
for (scn in 1:length(scenario)){
  bias=crop(raster(paste(getwd(),"/",scenario[scn],"/bias.asc", sep="")),extent(-12, 105, 9, 80))
  tmin =crop(raster(paste(getwd(),"/",scenario[scn],"/TMIN.asc", sep="")),extent(-12, 105, 9, 80))
  tmax=crop(raster(paste(getwd(),"/",scenario[scn],"/TMAX.asc", sep="")),extent(-12, 105, 9, 80))
  tmean=crop(raster(paste(getwd(),"/",scenario[scn],"/TMEAN.asc", sep="")),extent(-12, 105, 9, 80))
  prec=crop(raster(paste(getwd(),"/",scenario[scn],"/PREC.asc",sep="")),extent(-12, 105, 9, 80))
  land=crop(raster(paste(getwd(),"/",scenario[scn],"/land.asc", sep="")),extent(-12, 105, 9, 80))
  precsd=crop(raster(paste(getwd(),"/",scenario[scn],"/PRECSD.asc",sep="")),extent(-12, 105, 9, 80))
  tsd=crop(raster(paste(getwd(),"/",scenario[scn],"/TSD.asc",sep="")),extent(-12, 105, 9, 80))
  land=as.factor(land)
# Reclassify landuses for projection
  if (scenario [scn]== "current_scn"){ #the current scenario builds on layers from landcover and then needs a differnet classification to the rest of scenarios which are based on layers from IMAGE project.
    land<-reclassify(land, matrix(c(0,1,4, #level 4 is forest
                                    1.5,2,2, #level 2 is shrubland
                                    2.5,4,3, #level 3 is open
                                    4.5,5,1, #level 1 is other
                                    5.5,6,5, #level 5 is agriculture
                                    6.5,9,1), #level 1 is other
                                  ncol = 3, byrow = T)) #reclassifies land variable from interval of old values (first 2 columns in the matrix) to new value (third column)
    land<-as.factor(land)
    rat<-levels(land)[[1]]
    rat$landcover<-c('other','shrub','open','forest','agri')
    levels(land)<-rat

    # agri=reclassify(land, matrix(c(0,5,0,
    #                              5.5,6, 1,#only category six refer to agricultural habitats
    #                              6.5,9,0), ncol=3, byrow=TRUE))
    # shrub=reclassify(land, matrix(c(0,1,0,
    #                                1.5,2,1,#categories 2 (shrubland)
    #                                2.5,9, 0), ncol=3, byrow=TRUE))
    # 
    # open=reclassify(land, matrix(c(0,2,0,
    #                              2.5,4,1,#3 (savannas) & 4 (grassland) correspond to open habitats
    #                              4.5,9,0), ncol=3, byrow=TRUE))
    # 
    # forest=reclassify(land, matrix(c(0,1,1, #only category one refer to forests
    #                                1.5,9, 0), ncol=3, byrow=TRUE))
    # other=reclassify(land, matrix(c(0,4,0,
    #                               4.5,5,1,#category 5 (wetlands) correspond to other
    #                               5.5,6,0,
    #                               6.5,9,1), ncol=3, byrow=TRUE))# categories 7 (urban), 8(snow&ice) and 9 (barren)

    
    ## Round variables (?)
    
    bias<-round(bias, digits = 0) #bias density is rounded to integer
    
    tmin<-round(tmin, digits = 2) #temperatures are rounded to 2 decimals
    tmean<-round(tmean, digits = 2)
    tmax<-round(tmax, digits = 2)
    tsd<-round(tsd, digits = 2)
    
    prec<-round(prec, digits = 0) #anual precipitation is rounded to integer
    precsd<-round(precsd, digits = 2) #coefficient of variation of precipitation as a percentage is rounded to 2 decimals (might change?)
    
    bias@data@names<-'bias'; tmin@data@names<-'tmin'; tmean@data@names<-'tmean'; tmax@data@names<-'tmax'; tsd@data@names<-'tsd'
    prec@data@names<-'prec'; precsd@data@names<-'precsd'
    
    #standardize variables
    #tmin<-scale(tmin)
    #tmax<-scale(tmax)
    #tmean<-scale(tmean)
    #tsd<-scale(tsd)
    #prec<-scale(prec)
    #precsd<-scale(precsd)
    
    
    #create independent variable stack
    # myExpl.new = stack(tmax, tsd, prec, precsd, bias, agri, shrub, open, forest, other) #stack ensemble variables
    # names(myExpl.new)<-c("tmax", "tsd", "prec", "precsd", "bias", "agri", "shrub", "open", "forest", "other")
    # NAvalue(myExpl.new)=-9999
    # plot(myExpl.new)
    
    #land as factor
    #myExpl.new = stack(tmax, tsd, prec, precsd, bias, land) #stack ensemble variables
    #names(myExpl.new)<-c("tmax", "tsd", "prec", "precsd", "bias", "land")
    #NAvalue(myExpl.new)=-9999
    #plot(myExpl.new)
    
    
    
    #no precsd
    myExpl.new = stack(tmax, tsd, prec, bias, land) #stack ensemble variables
    names(myExpl.new)<-c("tmax", "tsd", "prec", "bias", "land")
    NAvalue(myExpl.new)=-9999
    
    #New ensemble table
    ens.table=as.data.frame(matrix(NA, ncol=10, nrow=44))
    colnames(ens.table)=c("model", "run", "threshold","AUC","omission.rate", "sensitivity", "specificity","prop.correct","Kappa","TSS")
    ens.table[,1]=c(rep(c(rep("Maxent",1),rep("GLM", 1), rep("GBM",1), rep("RF",1)), 11))
    ens.table[,2]=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), rep(10,4), rep("full",4))
    
    # For cycle to extract fitted values from projections and build ensemble table
    modelcomp<-1
    for(modelcomp in 1:length(myBiomodModelOut@models.computed)){
      
      myBiomodProjection_all <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                  new.env = myExpl.new,
                                                  proj.name = scenario[scn],
                                                  selected.models = myBiomodModelOut@models.computed[modelcomp],
                                                  binary.meth = 'TSS',
                                                  do.stack=T,
                                                  compress = FALSE,
                                                  build.clamping.mask = F)
      
      
      
      # 5. Build ensemble-models that will be taken as reference and Project 
      #ensamble
      myBiomodEM_all <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                                chosen.models = myBiomodModelOut@models.computed[modelcomp],
                                                em.by = 'all',
                                                eval.metric = c('TSS'),
                                                eval.metric.quality.threshold = -1,
                                                prob.mean = T,
                                                prob.median = F,
                                                prob.mean.weight=F,
                                                committee.averaging=F)
      #project
      myEnsembleForecasting_all=BIOMOD_EnsembleForecasting( myBiomodEM_all,
                                                            projection.output = myBiomodProjection_all)
      
      
      u_all=slot(myEnsembleForecasting_all@proj, "val")
      writeRaster(u_all, paste(getwd(),"/Results/trial", sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T)
      
      currentpred<-raster(paste(getwd(),"/Results/trial_1.asc",sep=""))
      
      #Testing model fit for the ensembled model
      #------------------------------------------
      
      #Extract model predictions for test data coordinates
      ep=raster::extract(currentpred, cbind(lcpi$x, lcpi$y))
      
      #Merge test data, with model predictions in the same points
      p=cbind(lcpi, ep)
      colnames(p)=c("x", "y", "pres", "p")#change field names
      
      #Find optimum threshold for accuracy analyses
      thresh=optim.thresh(lcpi$pres, ep/1000)
      t=as.numeric(unlist(thresh))#convert to numeric vector. This are differnt threshold. I will select the 5 becasue is often used (maximum sensitivity+specificity)
      #t<-thresh[modelcomp,]
      
      #Calcultae ensembled model accuracy
      ac=accuracy(lcpi$pres, ep/1000, threshold = t[5])#I use the threshhold in position 5
      
      #For single models
      ens.table[modelcomp,3:9]<-round(ac[1:7],digits = 3)
      ens.table[modelcomp,10]<-round(ac[4]+ac[5]-1, digits = 3)
      
    } #single models analysis
    
    View(ens.table)
    
    # NOTA: para tornar varias linhas em 'comentarios' é seleccionar e carregar Ctrl+Shift+C
    
    modelcomp<-1
    for(modelcomp in 1:length(myBiomodModelOut_full@models.computed)){
      myBiomodProjection_full <- BIOMOD_Projection(modeling.output = myBiomodModelOut_full,
                                                   new.env = myExpl.new,
                                                   proj.name = scenario[scn],
                                                   selected.models = myBiomodModelOut_full@models.computed[modelcomp],
                                                   binary.meth = 'TSS',
                                                   do.stack=T,
                                                   compress = FALSE,
                                                   build.clamping.mask = F)
      
      #ensamble
      myBiomodEM_full<- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut_full,
                                                chosen.models = myBiomodModelOut_full@models.computed[modelcomp],
                                                em.by = 'all',
                                                eval.metric = c('TSS'),
                                                eval.metric.quality.threshold = c(-1),
                                                prob.mean = T,
                                                prob.median = F,
                                                prob.mean.weight=F,
                                                committee.averaging=F)
      
      
      # #project
      myEnsembleForecasting_full=BIOMOD_EnsembleForecasting( myBiomodEM_full,
                                                             projection.output = myBiomodProjection_full)
      
      u_full=slot(myEnsembleForecasting_full@proj, "val")
      writeRaster(u_full, paste(getwd(),"/Results/",dados,condicao,m_name[modelcomp],scenario[scn], sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T)
      
      currentpred<-raster(paste(getwd(),"/Results/",dados,condicao,m_name[modelcomp],scenario[scn],"_1.asc",sep=""))
      plot(currentpred)
      #Testing model fit for the ensembled model
      #----------------------------------------------------------------------------------
      
      
      #Extract model predictions for test data coordinates
      ep=raster::extract(currentpred, cbind(lcpi$x,lcpi$y))
      
      #Merge test data, with model predictions in the same points
      p=cbind(lcpi, ep)
      colnames(p)=c("x", "y", "pres", "p")#change field names
      
      #Find optimum threshold for accuracy analyses
      thresh=optim.thresh(lcpi$pres, ep/1000)
      t=as.numeric(unlist(thresh))#convert to numeric vector. This are different threshold. I will select the 5 becasue is often used (maximum sensitivity+specificity)
      #t<-thresh[modelcomp+40,]
      
      
      #Calcultate ensembled model accuracy
      ac=accuracy(lcpi$pres, ep/1000, threshold = t[5])#I use the threshhold in position 5
      
      #For full models
      ens.table[modelcomp+40,3:9]<-round(ac[1:7],digits = 3)
      ens.table[modelcomp+40,10]<-round(ac[4]+ac[5]-1, digits = 3)
      
      #create binary transformation map
      bin_current<-BinaryTransformation(currentpred,ens.table$threshold[modelcomp+40]*1000)
      plot(bin_current)
      
      writeRaster(bin_current, paste(getwd(),"/Results/",dados,condicao,m_name[modelcomp],'bin', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T)
      
    } #full models analysis
    
    View(ens.table)
    
    
    
    ## Run full Ensemble
    # Select models
    
    measure<-'TSS' #select accuracy filter (AUC, Kappa or TSS)
    sel.models<-c()
    limit<-0.7    #select accuracy quality threshold
    for (i in 1:40){
      for(j in c(4,9,10)){
        if(ens.table[i,j]>=limit & colnames(ens.table[j])==measure){
          sel.models<-append(sel.models,i)
        }
      }
    }; sel.models; length(sel.models) #show selected models for ensemble and total number
    
    #Running projections and creating ensemble
    #----------------------------------------------------------------------------------
    # Projection
    myBiomodProjection_all <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                new.env = myExpl.new,
                                                proj.name = scenario[scn],
                                                selected.models = myBiomodModelOut@models.computed[sel.models],
                                                binary.meth = 'TSS',
                                                do.stack=T,
                                                compress = FALSE,
                                                build.clamping.mask = F)
    
    
    # 5. Build ensemble-models that will be taken as reference and Project 
    #ensamble
    myBiomodEM_all <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                              chosen.models = myBiomodModelOut@models.computed[sel.models],
                                              em.by = 'all',
                                              eval.metric = c('TSS'),
                                              eval.metric.quality.threshold = -1,
                                              prob.mean = T,
                                              prob.median = F,
                                              prob.mean.weight=F,
                                              committee.averaging=F,
                                              prob.cv = T,
                                              VarImport = 1)
    #project
    myEnsembleForecasting_all=BIOMOD_EnsembleForecasting(myBiomodEM_all,
                                                        projection.output = myBiomodProjection_all)
    
    #----------------------------------------------------------------------------------
    #Testing model fit for the ensembled model
    #----------------------------------------------------------------------------------
    
    #accuracy measures
    impor.variab.ens<-t(as.data.frame(get_variables_importance(myBiomodEM_all)))
    write.table(impor.variab.ens, file=paste(getwd(),"/Results/",dados,condicao,scenario[scn],"_variable_importance_ensemble.txt", sep=""), col.names=T, row.names = F, sep="\t")
    
    
    u_all=slot(myEnsembleForecasting_all@proj, "val")
    writeRaster(u_all, paste(getwd(),"/Results/",dados,condicao,scenario[scn],'ensemble', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T)
    
    currentpred<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[scn],"ensemble_1.asc",sep=""))
    plot(currentpred)
    #....................................
    
    #Extract model predictions for test data coordinates
    ep=raster::extract(currentpred, cbind(lcpi$x, lcpi$y))
    
    #Merge test data, with model predictions in the same points
    p=cbind(lcpi, ep)
    colnames(p)=c("x", "y", "pres", "p")#change field names
    
    #Find optimum threshold for accuracy analyses
    thresh=optim.thresh(lcpi$pres, ep/1000)
    t=as.numeric(unlist(thresh))#convert to numeric vector. This are differnt threshold. I will select the 5 because is often used (maximum sensitivity+specificity)
    
    
    #Calculate ensembled model accuracy
    ac=accuracy(lcpi$pres, ep/1000, threshold = t[5])#I use the threshhold in position 5
    ens.table[length(ens.table[,1])+1,]<-c('Ensemble',paste(measure,'>=',limit,sep = ""),round(ac,digits = 3),round(ac[4]+ac[5]-1,digits = 3))
    tail(ens.table)
    
    #create binary transformation map
    bin_current<-BinaryTransformation(currentpred,ens.table$threshold[ens.table$model=='Ensemble']*1000)
    plot(bin_current)
    
    writeRaster(bin_current, paste(getwd(),"/Results/",dados,condicao,scenario[scn],'bin', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T)
    
    #----------------------------------------------------------------------------------
    
    #write ensemble table
    write.table(ens.table,file=paste(getwd(),"/Results/",dados,condicao,scenario[scn],"_ensemble_model.txt",sep=""), col.names = T, row.names = F, sep = "\t") #Writes ensemble table
    
    #results table for average reults of models
    results.table<-as.data.frame(matrix(NA, ncol = 13, nrow=4)) #results table for all models
    colnames(results.table)<-c('models','tss.sd','tss.mean','tss.min','tss.max','sens.sd','sens.mean','sens.min','sens.max','spec.sd','spec.mean','spec.min','spec.max')
    m<-unique(ens.table$model)[-5] #model names except ensemble
    for(i in 1:length(m)){
      mod.table<-subset(ens.table,model==m[i] & run!='full')
      results.table[i,]<-c(m[i],
                           sd(mod.table$TSS),
                           mean(mod.table$TSS),
                           min(mod.table$TSS),
                           max(mod.table$TSS),
                           sd(mod.table$sensitivity),
                           mean(mod.table$sensitivity),
                           min(mod.table$sensitivity),
                           max(mod.table$sensitivity),
                           sd(mod.table$specificity),
                           mean(mod.table$specificity),
                           min(mod.table$specificity),
                           max(mod.table$specificity))
    } #calculate values for the metrics
    write.table(results.table, file=paste(getwd(),"/Results/",dados,condicao,scenario[scn],"_results.txt", sep=""), col.names = T, row.names = F, sep = "\t") #Writes results table
    
    
  } #current conditions ensemble modeling
  
  else{ #reclassification for future scenarios with landuses derived from IMAGE project
    
    land<-reclassify(land, matrix(c(0,1,5, #level 5 is agriculture
                                    15.5,16,2, #level 2 is shrubland
                                    1.5,2,3,
                                    6.5,7,3,
                                    13.5,14,3,
                                    16.5,17,3,#level 3 is open
                                    3.5,5,4,
                                    7.5,13,4,
                                    17.5,19,4,#level 4 is forest
                                    5.5,6,1,
                                    14.5,15,1#level 1 is other
                                    ),
                                  ncol = 3, byrow = T)) #reclassifies land variable from interval of old values (first 2 columns in the matrix) to new value (third column)
    #NOTE: Layer 3 of Land raster doesn't seem to have any values in it, so it doesn't have to be reclassified
    land<-as.factor(land)
    rat<-levels(land)[[1]]
    rat$landcover<-c('other','shrub','open','forest','agri')
    levels(land)<-rat

    # agri=reclassify(land, matrix(c(0,1,1,#only category 1 refers to agriculture
    #                          1.5,19, 0), ncol=3, byrow=TRUE))
    # shrub=reclassify(land,matrix(c(0,15,0,
    #                            15.5,16,1, # category 16 is shrubland
    #                            16.5,19,0),ncol = 3,byrow = TRUE))
    # open=reclassify(land, matrix(c(0,1,0,
    #                          1.5,2,1,# category 2 (extensive grassland) correpond to open
    #                          3.5,6, 0,
    #                          6.5,7,1,# category 7 (tundra) correpond to open
    #                          7.5,13, 0,
    #                          13.5, 14, 1, # category 14 (grassland steppe) correpond to open
    #                          14.5,16,0,
    #                          16.5,17, 1,#categories 16 (scrubland) and 17(savanna) also are open habitat
    #                          17.5, 19, 0), ncol=3, byrow=TRUE))
    # forest=reclassify(land, matrix(c(0,2,0,
    #                            3.5,5, 1, #categories 4 (regrowth forest - abandoning) and 5 (regrowth forest - timber) are forests
    #                            5.5,7,0,
    #                            7.5,13, 1, #categories 8 (wooded tundra), 9 (boreal forest), 10 (cool conifer), 11 (temp mixed forest), 12(temp decid forest) and 13 (warm mixed forest) are forests
    #                            13.5, 17, 0,
    #                            17.5, 19, 1), ncol=3, byrow=TRUE)) #categories 18 (tropical woodland) and 19 (tropical forests) are also forests
    # other=reclassify(land, matrix(c(0,5,0,
    #                           5.5,6,1,#category 6 (ice) correspond to other
    #                           6.5,17,0,
    #                           14.5,15,1, #catgory 15 (hot desert) correpond to other
    #                           15.5,19, 0), ncol=3, byrow=TRUE))
    
    
    ## Round variables (?)
    
    bias<-round(bias, digits = 0) #bias density is rounded to integer
    
    tmin<-round(tmin, digits = 2) #temperatures are rounded to 2 decimals
    tmean<-round(tmean, digits = 2)
    tmax<-round(tmax, digits = 2)
    tsd<-round(tsd, digits = 2)
    
    prec<-round(prec, digits = 0) #anual precipitation is rounded to integer
    precsd<-round(precsd, digits = 2) #coefficient of variation of precipitation as a percentage is rounded to 2 decimals (might change?)
    
    bias@data@names<-'bias'; tmin@data@names<-'tmin'; tmean@data@names<-'tmean'; tmax@data@names<-'tmax'; tsd@data@names<-'tsd'
    prec@data@names<-'prec'; precsd@data@names<-'precsd'
    
    #standardize variables
    #tmin<-scale(tmin)
    #tmax<-scale(tmax)
    #tmean<-scale(tmean)
    #tsd<-scale(tsd)
    #prec<-scale(prec)
    #precsd<-scale(precsd)
    
    
    #create independent variable stack
    # myExpl.new = stack(tmax, tsd, prec, precsd, bias, agri, shrub, open, forest, other) #stack ensemble variables
    # names(myExpl.new)<-c("tmax", "tsd", "prec", "precsd", "bias", "agri", "shrub", "open", "forest", "other")
    # NAvalue(myExpl.new)=-9999
    # plot(myExpl.new)
    
    #land as factor
    #myExpl.new = stack(tmax, tsd, prec, precsd, bias, land) #stack ensemble variables
    #names(myExpl.new)<-c("tmax", "tsd", "prec", "precsd", "bias", "land")
    #NAvalue(myExpl.new)=-9999
    #plot(myExpl.new)
    
    
    #no precsd
    myExpl.new = stack(tmax, tsd, prec, bias, land) #stack ensemble variables
    names(myExpl.new)<-c("tmax", "tsd", "prec", "bias", "land")
    NAvalue(myExpl.new)=-9999
    
    
    # Create new selected models projections based on future scenarios
    myFutureProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut, #call calibrated models
                                        new.env = myExpl.new,
                                        proj.name = scenario[scn],
                                        selected.models = myBiomodModelOut@models.computed[sel.models], #select only TSS>0.7
                                        do.stack=T,
                                        compress = F,
                                        build.clamping.mask = F)
    
    
    # 5. Build ensemble-models that will be taken as reference and Project 

    #project
    myFutureEnsembleProj <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_all, #call ensemble model
                                                        projection.output = myFutureProj) #use new model future projections
    
    #create raster
    
    u_future=slot(myFutureEnsembleProj@proj, "val") #store future projections probability values
    writeRaster(u_future, paste(getwd(),"/Results/",dados,condicao,scenario[scn],'ensemble', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T)
    
    currentpred<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[scn],"ensemble_1.asc",sep=""))
    plot(currentpred)
    
    #binary transformation
    bintrans<-BinaryTransformation(currentpred,ens.table$threshold[ens.table$model=='Ensemble']*1000) #transform probability predictions of future projections into binary pres-abs using current condition threshold
    plot(bintrans)
    
    writeRaster(bintrans, paste(getwd(),"/Results/",dados,condicao,scenario[scn],'bin', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T) #write binary future projection raster
    
    #range shift
    rangeshift<-BIOMOD_RangeSize(bin_current, bintrans) #calculate range shift
    write.table(rangeshift$Compt.By.Models, file=paste(getwd(),"/Results/",dados,condicao,scenario[scn],"_rangeshift.txt", sep=""), col.names = T, row.names = F, sep = "\t") #Writes range shift table
    
    plot(rangeshift$Diff.By.Pixel) #1=gain in presence; 0=equal absence; -1=equal presence; -2=loss of presence
    writeRaster(rangeshift$Diff.By.Pixel, paste(getwd(),"/Results/",dados,condicao,scenario[scn],'rangeshift', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T) #write range shift raster
    } #future projections
  
  if(scn==6){
    #create an ensemble from the two global climate change models predictions (maybe?)
    a1b_miroc_2050<-raster(paste(getwd(),"/Results/EBBA2FINAL_noprecsd2050_A1B_MRensemble_1.asc",sep="")) #a1b eco with miroc climate
    a1b_csiro_2050<-raster(paste(getwd(),"/Results/EBBA2FINAL_noprecsd2050_A1B_CSensemble_1.asc",sep="")) #a1b eco with csiro climate
    a2_miroc_2050<-raster(paste(getwd(),"/Results/EBBA2FINAL_noprecsd2050_A2_MRensemble_1.asc",sep="")) #a2 eco with miroc climate
    a2_csiro_2050<-raster(paste(getwd(),"/Results/EBBA2FINAL_noprecsd2050_A2_CSensemble_1.asc",sep="")) #a2 eco with csiro climate
    
    a1b_models<-stack(a1b_miroc_2050,a1b_csiro_2050) #combine both a1b economic projections
    a1b_ensemble<-mean(a1b_models); names(a1b_ensemble@data@names)<-"a1b ensemble" #mean ensemble of a1b model predictions
    plot(a1b_ensemble) #plot ensemble predictions of a1b projections
    a1b_bin<-BinaryTransformation(a1b_ensemble,ens.table$threshold[ens.table$model=='Ensemble']*1000) #binary projection of a1b    
    plot(a1b_bin) #plot pres-abs projections of a1b
    writeRaster(a1b_bin, paste(getwd(),"/Results/",dados,condicao,'a1b_bin', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T) #write a1b binary raster
    a1b_rs<-BIOMOD_RangeSize(bin_current,a1b_bin) #calculate range shift of a1b
    write.table(a1b_rs$Compt.By.Models, file=paste(getwd(),"/Results/",dados,condicao,"a1b_rangeshift.txt", sep=""), col.names = T, row.names = F, sep = "\t") #Writes a1b range shift table
    plot(a1b_rs$Diff.By.Pixel) #1=gain in presence; 0=equal absence; -1=equal presence; -2=loss of presence
    writeRaster(a1b_rs$Diff.By.Pixel, paste(getwd(),"/Results/",dados,condicao,'a1b_rangeshift', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T) #write a1b range shift raster
    
    a2_models<-stack(a2_miroc_2050,a2_csiro_2050) #combine both a2 economic projections
    a2_ensemble<-mean(a2_models); names(a2_ensemble@data@names)<-"a2 ensemble" #mean ensemble of a2 model predictions
    plot(a2_ensemble) #plot ensemble predictions of a2 projections
    a2_bin<-BinaryTransformation(a2_ensemble,ens.table$threshold[ens.table$model=='Ensemble']*1000) #binary projection of a2
    plot(a2_bin) #plot pres-abs projections of a2
    writeRaster(a2_bin, paste(getwd(),"/Results/",dados,condicao,'a2_bin', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T) #write a2 binary raster
    a2_rs<-BIOMOD_RangeSize(bin_current,a2_bin) #calculate range shift of a2
    write.table(a2_rs$Compt.By.Models, file=paste(getwd(),"/Results/",dados,condicao,"a2_rangeshift.txt", sep=""), col.names = T, row.names = F, sep = "\t") #Writes a2 range shift table
    plot(a2_rs$Diff.By.Pixel) #1=gain in presence; 0=equal absence; -1=equal presence; -2=loss of presence
    writeRaster(a2_rs$Diff.By.Pixel, paste(getwd(),"/Results/",dados,condicao,'a2_rangeshift', sep=""), format ='ascii', bylayer=T, suffix='numbers',overwrite=T) #write a2 range shift raster
    
    } #global climate model ensemble
  
} #ALL SCENARIOS MODELING

###-------------------------------------------

### Modeling algorithms analysis

## Load Full models for various evaluations
BIOMOD_LoadModels(myBiomodModelOut_full)
lc_MAXFull<-Lanius.collurio_PA1_Full_MAXENT.Phillips
lc_GLMFull<-Lanius.collurio_PA1_Full_GLM@model
lc_GBMFull<-Lanius.collurio_PA1_Full_GBM@model
lc_RFFull<-Lanius.collurio_PA1_Full_RF@model

## GLM
lc_GLMFull[["call"]] # Check GLM formula
summary(lc_GLMFull) # Check coefficients of GLM
lc_GLMFull[["anova"]] # Check ANOVA table


myRespPlot<-response.plot2(
  models = Lanius.collurio_PA1_Full_GLM@model_name, # name of model to be called
  Data = get_formal_data(myBiomodModelOut_full,'expl.var'), # table of expl variables
  show.variables = c('tmax','tsd','prec','bias','land'), # set of expl variables to plot
  col=c('blue'), #colour of plot
  #fixed.var.metric = 'max', # choose metric to set other variables values
  data_species=get_formal_data(myBiomodModelOut_full,'resp.var'), # only plot response curve for presence points
  #do.bivariate = T, # plot 3 dimension response curve for all pairs of variables
  legend=F # have legend for plot
); par(mfrow=c(1,1))


### Check overdispersion (if c-hat>1 then theres overdispersion, if c-hat>4 theres lack-of-fit)
lc_GLMFull$deviance/lc_GLMFull$df.residual #the residual deviance divided by degrees of freedom from residuals gives an estimate of c-hat

### Check Spatial Autocorrelation
lcenv<-as.data.frame(rasterToPoints(stack(myBiomodData@data.mask$PA1,tmax,tsd,prec,bias,land))) #create table with covariate values in each coord
lcenv<-na.omit(lcenv)
lcenv<-subset(lcenv, PA1!=-1) #remove -1 from data, which represent points without presence-PA data
length(lc_GLMFull$residuals)==length(lcenv[,1]) #must be true, so we are sure to work with same coordinates as used in model
dists<-as.matrix(dist(lcenv[,1:2])) #creates matrix with distance between each coord (!! VERY BIG FILE !!)
dists.inv<-1/dists #invert the matrix
diag(dists.inv)<-0 #put zeros in the diagonal of the matrix
Moran.I(lc_GLMFull$residuals, dists.inv) #run Moran test for spatial autocorrelation (if p.value<0.05 there's autocorrelation)
rm(dists, dists.inv)
glmrsd<-lc_GLMFull$residuals #puts GLM model residuals into object
rnd<-sample(1:length(glmrsd), 500, replace = T) #sample randomly 500 points
spat.cor<-correlog(lcenv[rnd,'x'], lcenv[rnd,'y'], glmrsd[rnd], increment = 2, resamp = 10)
plot(spat.cor); abline(a=0, b=0, col = 'red') #plot correlogram of residuals
plot(lcenv[,1:2][order(glmrsd),], pch = 7, col = rev(heat.colors(5900)), cex = .48, main = 'Residuals GLM', xlab = 'Latitude', ylab = 'Longitude', font.lab=2)

## MAXENT

myRespPlot<-response.plot2(
  models = Lanius.collurio_PA1_Full_MAXENT.Phillips@model_name, # name of model to be called
  Data = get_formal_data(myBiomodModelOut_full,'expl.var'), # table of expl variables
  show.variables = c('tmax','tsd','prec','bias','land'), # set of expl variables to plot
  col=c('blue'), #colour of plot
  #fixed.var.metric = 'median', # choose metric to set other variables
  data_species=get_formal_data(myBiomodModelOut_full,'resp.var'), # only plot response curve for presence points
  #do.bivariate = T, # plot 3 dimension response curve for all pairs of variables
  legend=F # have legend for plot
); par(mfrow=c(1,1))

## GBM
gbm.perf(lc_GBMFull) # Give number of trees for optimal fit (can be total number of trees if not enough were computed)
gbm.perf(lc_GBMFull, method = "cv", plot.it = T) # Graph error of the model in function of total number of trees
summary(lc_GBMFull, method = relative.influence, plotit = T) # Relative influence of variables (similar to AIC)
summary(lc_GBMFull, method = permutation.test.gbm, plotit = F) # Reduction in predictive performance when variable is permuted (similar to RF)


myRespPlot<-response.plot2(
  models = Lanius.collurio_PA1_Full_GBM@model_name, # name of model to be called
  Data = get_formal_data(myBiomodModelOut_full,'expl.var'), # table of expl variables
  show.variables = c('tmax','tsd','prec','bias','land'), # set of expl variables to plot
  col=c('blue'), #colour of plot
  #fixed.var.metric = 'max', # choose metric to set other variables
  data_species=get_formal_data(myBiomodModelOut_full,'resp.var'), # only plot response curve for presence points
  #do.bivariate = T, # plot 3 dimension response curve for all pairs of variables
  legend=F # have legend for plot
); par(mfrow=c(1,1))


## RF
importance(lc_RFFull) # Extract Gini Impurity from variables
plot(margin(lc_RFFull)) # Plot margin function (negtive values are incorrect predictions, positive values are correct predictions; red - absence, blue - presence)

myRespPlot<-response.plot2(
  models = Lanius.collurio_PA1_Full_RF@model_name, # name of model to be called
  Data = get_formal_data(myBiomodModelOut_full,'expl.var'), # table of expl variables
  show.variables = c('tmax','tsd','prec','bias','land'), # set of expl variables to plot
  col=c('blue'), #colour of plot
  #fixed.var.metric = 'median', # choose metric to set other variables
  data_species=get_formal_data(myBiomodModelOut_full,'resp.var'), # only plot response curve for presence points
  #do.bivariate = T, # plot 3 dimension response curve for all pairs of variables
  legend=F # have legend for plot
); par(mfrow=c(1,1))

## Model accuracy plot

x.metric<-"sensitivity" # set which metric x axis corresponds to
y.metric<-"specificity" #set which metric y axis corresponds to

ggplot(ens.table, aes(x = sensitivity, y = specificity, color = model)) +
  geom_point(size = 4) +  # Size of points
  labs(title = "Model Accuracy Values", #plot name
       x = x.metric, #name of x axis
       y = y.metric) + #name of y axis
  scale_color_manual(values = c("GLM" = "blue2", "Maxent" = "chocolate","GBM" = "aquamarine3", "RF" = "darkolivegreen", "Ensemble" = "darkgoldenrod1")) + # Customize colors of points
  theme_minimal() + #set a minimal theme to plot
  theme(legend.position = "right") + #where to put legend
  stat_ellipse(aes(fill = model),type = "t",level = 0.4, alpha = 0.1, geom = "polygon", color = NA) + #create cluster polygon for each model
  scale_fill_manual(values = c("GLM" = "blue2", "Maxent" = "chocolate","GBM" = "aquamarine3", "RF" = "darkolivegreen", "Ensemble" = "darkgoldenrod1")) + #change color of polygons
  xlim(min(ens.table$sensitivity),1) + #x axis value range
  ylim(min(ens.table$specificity),1) + #y axis value range
  geom_abline(slope = 1, intercept = 0, linetype = 1, color = "grey") # add line that follows a balanced Sensitivity+Specificity


# Ranks plot
x_ranked <- ens.table[,c('model','run','sensitivity')] %>%
  arrange(desc(sensitivity)) %>%
  mutate(Overall_Rank = row_number()) # create table with ranked metric for x axis of plot
y_ranked <- ens.table[,c('model','run','specificity')] %>%
  arrange(desc(specificity)) %>%
  mutate(Overall_Rank = row_number()) # create table with ranked metric for y axis of plot
metric_ranks<-merge(x_ranked,y_ranked, by=c('model','run')) #combine tables to plot


ggplot(metric_ranks, aes(x = Overall_Rank.x, y = Overall_Rank.y, color = model)) +
  geom_point(size = 4) +  # Size of points
  labs(title = "Model Accuracy Ranks",
       x = x.metric,
       y = y.metric) +
  scale_color_manual(values = c("GLM" = "blue2", "Maxent" = "chocolate","GBM" = "aquamarine3", "RF" = "darkolivegreen", "Ensemble" = "darkgoldenrod1")) + # Customize colors of points
  theme_minimal() + #set a minimal theme to plot
  theme(legend.position = "right") + #where to put legend
  scale_x_reverse() + #reverse x so lower ranks appear on the right
  scale_y_reverse() + #reverse y so lower ranks appear on the top
  stat_ellipse(aes(fill = model),type = "t",level = 0.45, alpha = 0.1, geom = "polygon", color = NA) + #create cluster polygon for each model
  scale_fill_manual(values = c("GLM" = "blue2", "Maxent" = "chocolate","GBM" = "aquamarine3", "RF" = "darkolivegreen", "Ensemble" = "darkgoldenrod1")) + #change color of polygons
  xlim(length(metric_ranks[,1]),1) +
  ylim(length(metric_ranks[,1]),1) +
  geom_abline(slope = 1, intercept = 0, linetype = 1, color = "grey") # add line that follows a balanced Sensitivity+Specificity

# Single models comparison
sm_table<-ens.table[which(ens.table$model!='Ensemble' & ens.table$run!='full'),] #single models table

ggplot(sm_table, aes(x = model, y = TSS, fill = model)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.2), color = "black", alpha = 0.4) +
  geom_abline(slope = 0, intercept = 0.7, color = 'red', linetype = 1) +
  labs(title = "Distribution of TSS by Algorithm",
       y = "TSS",
       x = "Model") +
  theme_minimal()
  
ggplot(sm_table, aes(x = model, y = TSS, fill = model)) +
  geom_boxplot() +
  labs(title = "Model Performance (TSS) by Algorithm",
        y = "TSS",
        x = "Algorithm") +
  geom_abline(slope = 0, intercept = 0.7, color = 'red', linetype = 1) +
  theme_minimal() +
  theme(legend.position = "none") #remove legend


## Model maps
# Continuous predictions for the whole area
dados<-"FULL" #set type of data the calls model predictions
condicao<-"original" #set occurrence map design


hs_pred<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[1],"ensemble_1.asc",sep="")) #get model predictions for current conditions
plot(hs_pred)

suitability_df <- as.data.frame(hs_pred, xy = TRUE) #put raster into data frame
colnames(suitability_df) <- c("longitude", "latitude", "suitability") #change column names
suitability_df$suitability<-suitability_df$suitability/10
world_map <- map_data("world")

ggplot() +
  # Plot the habitat suitability predictions as a raster
  geom_raster(data = suitability_df, aes(x = longitude, y = latitude, fill = suitability)) +
  # Set a color gradient for habitat suitability
  scale_fill_viridis(name = "Habitat Suitability", option = "turbo", limits = c(0, 100)) +
  # Customize map limits to the Palearctic region
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +
  # Country borders layer, drawn over the suitability map
  geom_path(data = world_map, aes(x = long, y = lat, group = group), color = "black", size = 0.2) +
  # Additional aesthetic elements
  labs(
    title = "Habitat Suitability Map for Red-backed Shrike in the whole range",
    subtitle = "Ensemble Model Predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# Occurrence prediction for the whole area
pa_pred<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[1],"bin.asc",sep="")) #call raster with bin pred for current
occurrence_df <- as.data.frame(pa_pred, xy = TRUE) #put raster into data frame
colnames(occurrence_df) <- c("longitude", "latitude", "occurrence") #change column names
occurrence_df$occurrence<-as.factor(occurrence_df$occurrence) #change occurrences to factor
world_map <- map_data("world")

ggplot() +
  # Plot the occurrence predictions as a raster
  geom_raster(data = occurrence_df, aes(x = longitude, y = latitude, fill = occurrence)) +
  # Set a color gradient for occurrence
  scale_fill_manual(values = c("lightgray", "darkred"), 
                    name = "Occurrence", 
                    labels = c("Absence", "Presence")) +
  # Customize map limits to the breeding range
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +
  # Country borders layer, drawn over the occurrence map
  geom_path(data = world_map, aes(x = long, y = lat, group = group), color = "black", size = 0.7) +
  # Additional aesthetic elements
  labs(
    title = "Occurrence Map for Red-backed Shrike in the whole range",
    subtitle = "Ensemble Model Predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# Habitat Suitability Iberian Peninsula
iberian_map<-world_map[which(world_map$region%in%c('Portugal','Spain') & is.na(world_map$subregion)),] #extract world map only for Iberian Peninsula
masked_habitat <- mask(hs_pred, lcpi.ras) #mask predictions raster for only the iberian peninsula
suitability_iberian<-as.data.frame(masked_habitat, xy=T) #extract values from raster
colnames(suitability_iberian) <- c("longitude", "latitude", "suitability") #change column names
suitability_iberian$suitability<-suitability_iberian$suitability/10

iberian_pres<-lcpi[lcpi$pres==1,] #extract presences from census data
iberian_pres<-merge(suitability_iberian,iberian_pres,by=c(1,2),all.y=T) #merge data to remove presences where suitability NA
iberian_pres<-iberian_pres[-which(is.na(iberian_pres$suitability)),] #remove presences where suitability is NA
iberian_pres<-iberian_pres[,-3] #remove suitability column

ggplot() +
  # Plot the habitat suitability predictions as a raster
  geom_raster(data = suitability_iberian, aes(x = longitude, y = latitude, fill = suitability)) +
  # Set a color gradient for habitat suitability
  scale_fill_viridis(name = "Habitat Suitability", option = "turbo", limits = c(0, 100), alpha=1) +
  # Customize map limits to the Iberian Peninsula
  coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
  # Country borders layer, drawn over the suitability map
  geom_path(data = iberian_map, aes(x = long, y = lat, group = group), color = "black", size = 0.7) +
  # Add presences
  geom_point(data = iberian_pres, aes(x=longitude,y=latitude),color='black',size=1.5) +
  # Additional aesthetic elements
  labs(
    title = "Habitat Suitability Map for Red-backed Shrike in the Iberian Peninsula",
    subtitle = "Ensemble Model Predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Map Binary Presence-Absence predictions
iberian_map<-world_map[which(world_map$region%in%c('Portugal','Spain') & is.na(world_map$subregion)),] #extract world map only for Iberian Peninsula
pa_pred<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[1],"bin.asc",sep="")) #call raster with bin pred
pa_pred_ip<-mask(pa_pred,lcpi.ras) #mask areas outside Iberian Peninsula
occurrence_ip<-as.data.frame(pa_pred_ip,xy=T) #extract raster values
colnames(occurrence_ip) <- c("longitude", "latitude", "occurrence") #change column names
occurrence_ip$occurrence<-as.factor(occurrence_ip$occurrence)

ggplot() +
  # Plot the occurrence predictions as a raster
  geom_raster(data = occurrence_ip, aes(x = longitude, y = latitude, fill = occurrence)) +
  # Set a color gradient for occurrence
  scale_fill_manual(values = c("lightgray", "darkred"), 
                    name = "Occurrence", 
                    labels = c("Absence", "Presence")) +
  # Customize map limits to the Iberian Peninsula
  coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
  # Country borders layer, drawn over the occurrence map
  geom_path(data = iberian_map, aes(x = long, y = lat, group = group), color = "black", size = 0.7) +
  # Add presences
  geom_point(data = iberian_pres, aes(x=longitude,y=latitude), color='black',size=1.5) +
  # Additional aesthetic elements
  labs(
    title = "Occurrence Map for Red-backed Shrike in the Iberian Peninsula",
    subtitle = "Ensemble Model Predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Model algorithm map probability predictions for whole range
dados<-"FULL" #set type of data the calls model predictions
condicao<-"original" #set occurrence map design


hs_glm<-raster(paste(getwd(),"/Results/",dados,condicao,'GLM',scenario[1],"_1.asc",sep="")) #get model predictions
hs_maxent<-raster(paste(getwd(),"/Results/",dados,condicao,'Maxent',scenario[1],"_1.asc",sep="")) #get model predictions
hs_gbm<-raster(paste(getwd(),"/Results/",dados,condicao,'GBM',scenario[1],"_1.asc",sep="")) #get model predictions
hs_rf<-raster(paste(getwd(),"/Results/",dados,condicao,'RF',scenario[1],"_1.asc",sep="")) #get model predictions

suitability_models <- as.data.frame(stack(hs_glm,hs_maxent,hs_gbm,hs_rf), xy = TRUE) #put raster into data frame
colnames(suitability_models) <- c("longitude", "latitude", "GLM",'Maxent','GBM','RF') #change column names
suit_models_long <- melt(suitability_models, id.vars = c("longitude", "latitude"),
                variable.name = "model",
                value.name = "suitability") #reshape table for ggplot
suit_models_long$suitability<-suit_models_long$suitability/10

world_map <- map_data("world")

ggplot(suit_models_long, aes(x = longitude, y = latitude, fill = suitability)) +
  geom_raster() +
  scale_fill_viridis(name = "Habitat Suitability", option = "turbo") +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +  # Adjust for the whole range
  labs(
    title = "Habitat Suitability Predictions by Model",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Model algorithm occurrence for Iberian Peninsula
dados<-"FULL" #set type of data the calls model predictions
condicao<-"original" #set occurrence map design

iberian_map<-world_map[which(world_map$region%in%c('Portugal','Spain') & is.na(world_map$subregion)),] #extract world map only for Iberian Peninsula

pa_glm<-raster(paste(getwd(),"/Results/",dados,condicao,'GLM',"bin.asc",sep="")) #get model predictions
pa_maxent<-raster(paste(getwd(),"/Results/",dados,condicao,'Maxent',"bin.asc",sep="")) #get model predictions
pa_gbm<-raster(paste(getwd(),"/Results/",dados,condicao,'GBM',"bin.asc",sep="")) #get model predictions
pa_rf<-raster(paste(getwd(),"/Results/",dados,condicao,'RF',"bin.asc",sep="")) #get model predictions

pa_glm<-mask(pa_glm,lcpi.ras)
pa_maxent<-mask(pa_maxent,lcpi.ras)
pa_gbm<-mask(pa_gbm,lcpi.ras)
pa_rf<-mask(pa_rf,lcpi.ras)

occurrence_models <- as.data.frame(stack(pa_glm,pa_maxent,pa_gbm,pa_rf), xy = TRUE) #put raster into data frame
colnames(occurrence_models) <- c("longitude", "latitude", "GLM",'Maxent','GBM','RF') #change column names
occ_models_long <- melt(occurrence_models, id.vars = c("longitude", "latitude"),
                         variable.name = "model",
                         value.name = "occurrence") #reshape table for ggplot
occ_models_long$occurrence<-as.factor(occ_models_long$occurrence) #change to factor
world_map <- map_data("world")

ggplot() +
  geom_raster(data = occ_models_long, aes(x = longitude, y = latitude, fill = occurrence)) +
  # Set a color gradient for occurrence
  scale_fill_manual(values = c("lightgray", "darkred"), 
                    name = "Occurrence", 
                    labels = c("Absence", "Presence")) +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  # Customize map limits to the Iberian Peninsula
  coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
  # Country borders layer, drawn over the occurrence map
  geom_path(data = iberian_map, aes(x = long, y = lat, group = group), color = "black", size = 1) +
  # Add presences
  geom_point(data = iberian_pres, aes(x=longitude,y=latitude), color='black',size=1.5) +
  labs(
    title = "Occurrence Predictions by Model",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Map habitat suitability for future scenarios
dados<-"FULL" #set type of data the calls model predictions
condicao<-"original" #set occurrence map design


hs_2050_A1B_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[2],"ensemble_1.asc",sep="")) #get model predictions
hs_2050_A1B_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[3],"ensemble_1.asc",sep="")) #get model predictions
hs_2050_A2_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[4],"ensemble_1.asc",sep="")) #get model predictions
hs_2050_A2_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[5],"ensemble_1.asc",sep="")) #get model predictions

suitability_future <- as.data.frame(stack(hs_2050_A1B_MR,hs_2050_A1B_CS,hs_2050_A2_MR,hs_2050_A2_CS), xy = TRUE) #put raster into data frame
colnames(suitability_future) <- c("longitude", "latitude", "A1B_MR",'A1B_CS','A2_MR','A2_CS') #change column names
suit_future_long <- melt(suitability_future, id.vars = c("longitude", "latitude"),
                         variable.name = "model",
                         value.name = "suitability") #reshape table for ggplot
suit_future_long$suitability<-suit_future_long$suitability/10

world_map <- map_data("world")

ggplot(suit_future_long, aes(x = longitude, y = latitude, fill = suitability)) +
  geom_raster() +
  scale_fill_viridis(name = "Habitat Suitability", option = "turbo", limits = c(0, 100)) +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +  # Adjust for the whole range
  labs(
    title = "Habitat Suitability Predictions by Scenario",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Map occurrence for future scenarios
pa_2050_A1B_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[2],"bin.asc",sep="")) #get model predictions
pa_2050_A1B_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[3],"bin.asc",sep="")) #get model predictions
pa_2050_A2_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[4],"bin.asc",sep="")) #get model predictions
pa_2050_A2_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[5],"bin.asc",sep="")) #get model predictions

pa_2050_A1B_MR<-mask(pa_2050_A1B_MR,lcpi.ras)
pa_2050_A1B_CS<-mask(pa_2050_A1B_CS,lcpi.ras)
pa_2050_A2_MR<-mask(pa_2050_A2_MR,lcpi.ras)
pa_2050_A2_CS<-mask(pa_2050_A2_CS,lcpi.ras)


occurrence_future <- as.data.frame(stack(pa_2050_A1B_MR,pa_2050_A1B_CS,pa_2050_A2_MR,pa_2050_A2_CS), xy = TRUE) #put raster into data frame
colnames(occurrence_future) <- c("longitude", "latitude", "A1B_MR",'A1B_CS','A2_MR','A2_CS') #change column names
occ_future_long <- melt(occurrence_future, id.vars = c("longitude", "latitude"),
                         variable.name = "model",
                         value.name = "occurrence") #reshape table for ggplot
occ_future_long$occurrence<-as.factor(occ_future_long$occurrence)

ggplot() +
  geom_raster(data = occ_future_long, aes(x = longitude, y = latitude, fill = occurrence)) +
  # Set a color gradient for occurrence
  scale_fill_manual(values = c("lightgray", "darkred"), 
                    name = "Occurrence", 
                    labels = c("Absence", "Presence")) +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  # Customize map limits to the Iberian Peninsula
  coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
  # Country borders layer, drawn over the occurrence map
  geom_path(data = iberian_map, aes(x = long, y = lat, group = group), color = "black", size = 1) +
  # Add presences
  geom_point(data = iberian_pres, aes(x=longitude,y=latitude), color='black',size=1.5) +
  labs(
    title = "Occurrence Predictions by Scenario",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Map range shift for future scenarios
rs_2050_A1B_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[2],"rangeshift_1.asc",sep="")) #get model predictions
rs_2050_A1B_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[3],"rangeshift_1.asc",sep="")) #get model predictions
rs_2050_A2_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[4],"rangeshift_1.asc",sep="")) #get model predictions
rs_2050_A2_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[5],"rangeshift_1.asc",sep="")) #get model predictions


rangeshift_future <- as.data.frame(stack(rs_2050_A1B_MR,rs_2050_A1B_CS,rs_2050_A2_MR,rs_2050_A2_CS), xy = TRUE) #put raster into data frame
colnames(rangeshift_future) <- c("longitude", "latitude", "A1B_MR",'A1B_CS','A2_MR','A2_CS') #change column names
rs_future_long <- melt(rangeshift_future, id.vars = c("longitude", "latitude"),
                        variable.name = "model",
                        value.name = "shift") #reshape table for ggplot
rs_future_long$shift<-as.factor(rs_future_long$shift)

ggplot() +
  geom_raster(data = rs_future_long, aes(x = longitude, y = latitude, fill = shift)) +
  # Set a color gradient for shift
  scale_fill_manual(values = c("darkred", "black",'lightgrey','darkgreen'), 
                    name = "Shift", 
                    labels = c("Loss", "Presence (Stable)","Absence (Stable)","Gain")) +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  # Customize map limits to the Iberian Peninsula
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +
  labs(
    title = "Range Shift Predictions by Scenario",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

rs_ip_A1B_MR<-mask(rs_2050_A1B_MR,lcpi.ras)
rs_ip_A1B_CS<-mask(rs_2050_A1B_CS,lcpi.ras)
rs_ip_A2_MR<-mask(rs_2050_A2_MR,lcpi.ras)
rs_ip_A2_CS<-mask(rs_2050_A2_CS,lcpi.ras)

rangeshift_ip <- as.data.frame(stack(rs_ip_A1B_MR,rs_ip_A1B_CS,rs_ip_A2_MR,rs_ip_A2_CS), xy = TRUE) #put raster into data frame
colnames(rangeshift_ip) <- c("longitude", "latitude", "A1B_MR",'A1B_CS','A2_MR','A2_CS') #change column names
rs_ip_long <- melt(rangeshift_ip, id.vars = c("longitude", "latitude"),
                       variable.name = "model",
                       value.name = "shift") #reshape table for ggplot
rs_ip_long$shift<-as.factor(rs_ip_long$shift)

ggplot() +
  geom_raster(data = rs_ip_long, aes(x = longitude, y = latitude, fill = shift)) +
  # Set a color gradient for shift
  scale_fill_manual(values = c("darkred", "black",'lightgrey','darkgreen'), 
                    name = "Shift", 
                    labels = c("Loss", "Presence (Stable)","Absence (Stable)","Gain")) +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  # Customize map limits to the Iberian Peninsula
  coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
  # Country borders layer, drawn over the occurrence map
  geom_path(data = iberian_map, aes(x = long, y = lat, group = group), color = "black", size = 1) +
  labs(
    title = "Range Shift Predictions by Scenario for the Iberian Peninsula",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Map coefficient of variation (uncertainty)
cv_ensemble<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[1],"ensemble_2.asc",sep="")) #get model predictions
cv_2050_A1B_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[2],"ensemble_2.asc",sep="")) #get model predictions
cv_2050_A1B_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[3],"ensemble_2.asc",sep="")) #get model predictions
cv_2050_A2_MR<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[4],"ensemble_2.asc",sep="")) #get model predictions
cv_2050_A2_CS<-raster(paste(getwd(),"/Results/",dados,condicao,scenario[5],"ensemble_2.asc",sep="")) #get model predictions

uncertainty_ensemble<-as.data.frame(cv_ensemble, xy =T) #extract cv of ensemble
colnames(uncertainty_ensemble)<-c("longitude", "latitude","cv")

ggplot() +
  # Plot the uncertainty predictions as a raster
  geom_raster(data = uncertainty_ensemble, aes(x = longitude, y = latitude, fill = cv)) +
  # Set a color gradient for uncertainty
  scale_fill_viridis(name = "Uncertainty", option = "turbo") +
  # Customize map limits to the Palearctic region
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +
  # Country borders layer, drawn over the uncertainty map
  geom_path(data = world_map, aes(x = long, y = lat, group = group), color = "black", size = 0.2) +
  # Additional aesthetic elements
  labs(
    title = "Uncertainty Map for Red-backed Shrike in the whole range",
    subtitle = "Ensemble Model Predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )


uncertainty_future <- as.data.frame(stack(cv_2050_A1B_MR,cv_2050_A1B_CS,cv_2050_A2_MR,cv_2050_A2_CS), xy = TRUE) #put raster into data frame
colnames(uncertainty_future) <- c("longitude", "latitude","A1B_MR",'A1B_CS','A2_MR','A2_CS') #change column names
cv_future_long <- melt(uncertainty_future, id.vars = c("longitude", "latitude"),
                         variable.name = "model",
                         value.name = "uncertainty") #reshape table for ggplot


world_map <- map_data("world")

ggplot(cv_future_long, aes(x = longitude, y = latitude, fill = uncertainty)) +
  geom_raster() +
  scale_fill_viridis(name = "Uncertainty", option = "turbo") +
  facet_wrap(~ model, ncol = 2) +  # Change ncol to adjust the layout as needed
  coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +  # Adjust for the whole range
  labs(
    title = "Uncertainty Predictions by Scenario",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

## Plot occurrence data

occ_data<-myBiomodData@data.mask$input_data
pi_data<-mask(occ_data,lcpi.ras)
plot(stack(pi_data,lcpi.ras), xlim=c(-10,4), ylim=c(35,45))

## Decision tree plot
# Load required libraries
library(raster)
library(rpart) #create decision tree
library(rpart.plot) #plot decision tree

# Load your raster files (example file paths)
tree_data<-stack(occ_data,myExpl)
tree_data <- as.data.frame(na.omit(values(tree_data))) # Remove NA values


# Rename columns for clarity
colnames(tree_data) <- c("occurrence","tmin","tmean","tmax", "tsd", "prec", "precsd","bias","land")

# Convert land use to a factor
tree_data$land <- as.factor(tree_data$land)
tree_data$occurrence[which(tree_data$occurrence==-1)]<-0
tree_data$occurrence <- as.factor(tree_data$occurrence)

# Train a decision tree model
tree_model <- rpart(
  occurrence ~ tmax + tsd + prec + land ,
  data = tree_data,
  method = "class",
  control = rpart.control(maxdepth = 5)
)

# Plot the decision tree
rpart.plot(tree_model, type = 0, extra = 0, fallen.leaves = F, main = "Decision Tree example")


## Ensemble response curves
myglm<-BIOMOD_LoadModels(myBiomodModelOut, models='GLM')
myRespPlot_GLM<-response.plot2(
  models = myglm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

mymaxent<-BIOMOD_LoadModels(myBiomodModelOut, models='MAXENT.Phillips')
myRespPlot_MaxEnt<-response.plot2(
  models = mymaxent,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

mygbm<-BIOMOD_LoadModels(myBiomodModelOut, models='GBM')
myRespPlot_GBM<-response.plot2(
  models = mygbm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

myrf<-BIOMOD_LoadModels(myBiomodModelOut, models='RF')
myRespPlot_RF<-response.plot2(
  models = myrf,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = c('tmax','tsd','prec','bias','land'),
  col=c('blue'),
  legend=F,
  data_species=get_formal_data(myBiomodModelOut,'resp.var')
); par(mfrow=c(1,1))

# Assuming you have separate dataframes for each algorithm and the ensemble:
# rf_data, gbm_data, and ensemble_data with columns: Variable, x, prediction

# Add a column to identify the algorithm
rf_data <- myRespPlot_RF[,c(2,3,5)] %>% mutate(Algorithm = "RF")
gbm_data <- myRespPlot_GBM[,c(2,3,5)] %>% mutate(Algorithm = "GBM")
glm_data <- myRespPlot_GLM[,c(2,3,5)] %>% mutate(Algorithm = "GLM")
max_data <- myRespPlot_MaxEnt[,c(2,3,5)] %>% mutate(Algorithm = "MaxEnt")

# Combine all datasets
combined_data <- bind_rows(rf_data, gbm_data, glm_data, max_data)

# Average predictions by Variable and x for each Algorithm
combined_data<-combined_data[-which(combined_data$expl.name=='bias'),] #remove bias so it doesnt get plotted
averaged_combined_data <- combined_data %>%
  group_by(Algorithm, expl.name, expl.val) %>%
  summarize(prediction = mean(pred.val), .groups = "drop") # do average of predictions (single line)
averaged_combined_data <- combined_data %>%
  group_by(Algorithm, expl.name, expl.val) %>%
  summarize(
    mean = mean(pred.val),
    lower = mean(pred.val) - 1.96 * sd(pred.val) / sqrt(n()),
    upper = mean(pred.val) + 1.96 * sd(pred.val) / sqrt(n()),
    .groups = "drop"
  ) # do confidence intervals 95%
averaged_combined_data$expl.name <- factor(
  averaged_combined_data$expl.name,
  levels = c("tmax","tsd", "prec", 'land') # Specify the desired order
) #change expl vars names to factors for ggplot2
averaged_combined_data$Algorithm <- factor(
  averaged_combined_data$Algorithm,
  levels = c("GLM","MaxEnt", "GBM", 'RF') # Specify the desired order
) #change model names to factors for ggplot2


# Plot all models' response curves
ggplot(averaged_combined_data, aes(x = expl.val, y = prediction, color = Algorithm)) +
  geom_line(data = subset(averaged_combined_data, expl.name != "land"), size = 1.2) +
  geom_point(data = subset(averaged_combined_data, expl.name == "land"), size = 2) +
  facet_wrap(~expl.name, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Environmental Variable",
    y = "Predicted Suitability",
    color = "Model",
    title = "Response Curves for All Models"
  ) +
  theme_minimal()

# Plot confidence intervals
ggplot(averaged_combined_data, aes(x = expl.val, y = mean, color = Algorithm, fill = Algorithm)) +
  geom_line(data = subset(averaged_combined_data, expl.name != "land"), size = 1.2) +
  geom_ribbon(
    data = subset(averaged_combined_data, expl.name != "land"),
    aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA
  ) +
  geom_point(
    data = subset(averaged_combined_data, expl.name == "land"),
    size = 2
  ) +
  facet_wrap(~expl.name, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Environmental Variable",
    y = "Predicted Suitability",
    color = "Model",
    fill = "Model",
    title = "Response Curves with Confidence Intervals"
  ) +
  theme_minimal()

response_data<-merge(myRespPlot_RF,myRespPlot_GBM, by=c('expl.name','expl.val')) #merge both resp plot info from GBM and RF (used in ensemble model) 
response_data$ensemble<-rowMeans(response_data[,c('pred.val.x','pred.val.y')]) # create ensemble predictions out of average
response_data<-response_data[,-c(3:(length(response_data[1,])-1))] #remove unecessary colums
response_long <- melt(
  response_data,
  id.vars = c('expl.name','expl.val'),
  variable.name = "model",
  value.name = "prediction"
) #transform into long format for ggplo2
response_long<-response_long[-which(response_long$expl.name=='bias'),] #remove bias so it doesnt get plotted
response_long$expl.name <- factor(
  response_long$expl.name,
  levels = c("tmax","tsd", "prec", 'land') # Specify the desired order
) #change expl vars names to factors for ggplot2

response_long <- response_long %>%
  group_by(expl.name, expl.val) %>%
  summarize(prediction = mean(prediction), .groups = "drop") #remove duplicate values
response_long$VarType <- ifelse(response_long$expl.name == "land", "categorical", "continuous") #create extra column to separate bioclimatic (continuous) and land (cathegorical)

ggplot(response_long, aes(x = expl.val, y = prediction)) +
  geom_line(data = subset(response_long, VarType == "continuous"), color = "blue", size = 1) +
  geom_point(data = subset(response_long, VarType == "categorical"), color = "blue", size = 2) +
  facet_wrap(~expl.name, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Environmental Variable",
    y = "Predicted Suitability",
    title = "Response Curves for Ensemble Model"
  ) +
  theme_minimal()

## Create new maps for environmental covariates

# Whole breeding range
scenario=c("current_scn", "2050_A1B_MR", "2050_A1B_CS", "2050_A2_MR", "2050_A2_CS")
scn<-1

env_limits<-as.data.frame(matrix(NA, nrow = 3, ncol=2)) #store limit values for all plots to remain consistent colour gradient
rownames(env_limits)<-c('tmax','tsd','prec')
colnames(env_limits)<-c('max','min')
env_limits[1,1]<-50.97;env_limits[2,1]<-18.40;env_limits[3,1]<-5329
env_limits[1,2]<-2.42;env_limits[2,2]<-0.51;env_limits[3,2]<-0

for (scn in 1:length(scenario)){
  tmax=crop(raster(paste(getwd(),"", scenario[scn],"/TMAX.asc", sep="")),extent(-12, 105, 9, 80))
  tsd=crop(raster(paste(getwd(),"", scenario[scn], "/TSD.asc",sep="")),extent(-12, 105, 9, 80))
  prec=crop(raster(paste(getwd(),"",scenario[scn], "/PREC.asc",sep="")),extent(-12, 105, 9, 80))
  land=crop(raster(paste(getwd(),"",scenario[scn],"/land.asc", sep="")),extent(-12, 105, 9, 80))
  land=as.factor(land)
  # Reclassify landuses for projection
  if (scenario [scn]== "current_scn"){
    # reclassify land for current conditions
    land<-reclassify(land, matrix(c(0,1,4, #level 4 is forest
                                    1.5,2,2, #level 2 is shrubland
                                    2.5,4,3, #level 3 is open
                                    4.5,5,1, #level 1 is other
                                    5.5,6,5, #level 5 is agriculture
                                    6.5,9,1), #level 1 is other
                                  ncol = 3, byrow = T)) #reclassifies land variable from interval of old values (first 2 columns in the matrix) to new value (third column)
    land<-as.factor(land)
    rat<-levels(land)[[1]]
    rat$landcover<-c('other','shrub','open','forest','agri')
    levels(land)<-rat
    
  }else{
    # reclassify land for future conditions
    land<-reclassify(land, matrix(c(0,1,5, #level 5 is agriculture
                                    15.5,16,2, #level 2 is shrubland
                                    1.5,2,3,
                                    6.5,7,3,
                                    13.5,14,3,
                                    16.5,17,3,#level 3 is open
                                    3.5,5,4,
                                    7.5,13,4,
                                    17.5,19,4,#level 4 is forest
                                    5.5,6,1,
                                    14.5,15,1),#level 1 is other
                                  ncol = 3, byrow = T)) #reclassifies land variable from interval of old values (first 2 columns in the matrix) to new value (third column)
    #NOTE: Layer 3 of Land raster doesn't seem to have any values in it, so it doesn't have to be reclassified
    land<-as.factor(land)
    rat<-levels(land)[[1]]
    rat$landcover<-c('other','shrub','open','forest','agri')
    levels(land)<-rat
  }
  
  mycovars<-as.data.frame(rasterToPoints(stack(tmax,tsd,prec,land)))
  colnames(mycovars)<-c('x','y','tmax','tsd','prec','land')
  mycovars<-na.omit(mycovars)
  
  mycovars$tmax <- as.numeric(mycovars$tmax)
  mycovars$prec <- as.numeric(mycovars$prec)
  mycovars$tsd <- as.numeric(mycovars$tsd)
  mycovars$land<-factor(mycovars$land)

  # Separate continuous variables
  df_continuous <- mycovars[, c("x", "y", "tmax", "tsd", "prec")]

  # Separate categorical variables
  df_categorical <- mycovars[, c("x", "y", "land")]

  # Melt continuous variables
  df_continuous_long <- melt(df_continuous, id.vars = c("x", "y"))

  # Melt categorical variables (land)
  df_categorical_long <- melt(df_categorical, id.vars = c("x", "y"))

  # Convert 'land' to factor
  df_categorical_long$value <- as.factor(df_categorical_long$value)

  # Plot for continuous variables
  plots <- list()
  for (var in c("tmax", "tsd", "prec")) {
    p <- ggplot(subset(df_continuous_long, variable == var), aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_viridis(option = "turbo", limit=c(env_limits[var,2],env_limits[var,1])) +
      coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +  # Adjust for the whole range
      theme_minimal() +
      ggtitle(paste(var))
    plots[[var]] <- p
  }

  land_colors <- c(
    "1" = "#B0B0B0",  # Gray for 'Other'
    "2" = "#9ACD32",  # Light Green for 'Shrub'
    "3" = "#FFD700",  # Yellow for 'Open'
    "4" = "#228B22",  # Dark Green for 'Forest'
    "5" = "#A0522D"   # Brown for 'Agriculture'
  )

  # Plot for categorical variable (land)
  p_land <- ggplot(subset(df_categorical_long, variable == "land"), aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_manual(values = land_colors) +
    coord_cartesian(xlim = c(-6, 90), ylim = c(32,70)) +  # Adjust for the whole range
    theme_minimal() +
    ggtitle(paste("land"))
  plots[["land"]] <- p_land

  # Arrange all plots in a 2x2 grid
  grid.arrange(
    plots[["tmax"]],
    plots[["tsd"]],
    plots[["prec"]],
    plots[["land"]],
    ncol = 2
  )

}

# Iberian Peninsula
iberian_map<-world_map[which(world_map$region%in%c('Portugal','Spain') & is.na(world_map$subregion)),]

scenario=c("current_scn", "2050_A1B_MR", "2050_A1B_CS", "2050_A2_MR", "2050_A2_CS")
scn<-1

env_limits<-as.data.frame(matrix(NA, nrow = 3, ncol=2)) #store limit values for all plots to remain consistent colour gradient
rownames(env_limits)<-c('tmax','tsd','prec')
colnames(env_limits)<-c('max','min')
env_limits[1,1]<-50.97;env_limits[2,1]<-18.40;env_limits[3,1]<-5329
env_limits[1,2]<-2.42;env_limits[2,2]<-0.51;env_limits[3,2]<-0

for (scn in 1:length(scenario)){
  tmax=crop(raster(paste(getwd(),"", scenario[scn],"/TMAX.asc", sep="")),extent(-12, 105, 9, 80))
  tsd=crop(raster(paste(getwd(),"", scenario[scn], "/TSD.asc",sep="")),extent(-12, 105, 9, 80))
  prec=crop(raster(paste(getwd(),"",scenario[scn], "/PREC.asc",sep="")),extent(-12, 105, 9, 80))
  land=crop(raster(paste(getwd(),"",scenario[scn],"/land.asc", sep="")),extent(-12, 105, 9, 80))
  land=as.factor(land)
  # Reclassify landuses for projection
  if (scenario [scn]== "current_scn"){
    # reclassify land for current conditions
    land<-reclassify(land, matrix(c(0,1,4, #level 4 is forest
                                    1.5,2,2, #level 2 is shrubland
                                    2.5,4,3, #level 3 is open
                                    4.5,5,1, #level 1 is other
                                    5.5,6,5, #level 5 is agriculture
                                    6.5,9,1), #level 1 is other
                                  ncol = 3, byrow = T)) #reclassifies land variable from interval of old values (first 2 columns in the matrix) to new value (third column)
    land<-as.factor(land)
    rat<-levels(land)[[1]]
    rat$landcover<-c('other','shrub','open','forest','agri')
    levels(land)<-rat
    
  }else{
    # reclassify land for future conditions
    land<-reclassify(land, matrix(c(0,1,5, #level 5 is agriculture
                                    15.5,16,2, #level 2 is shrubland
                                    1.5,2,3,
                                    6.5,7,3,
                                    13.5,14,3,
                                    16.5,17,3,#level 3 is open
                                    3.5,5,4,
                                    7.5,13,4,
                                    17.5,19,4,#level 4 is forest
                                    5.5,6,1,
                                    14.5,15,1),#level 1 is other
                                  ncol = 3, byrow = T)) #reclassifies land variable from interval of old values (first 2 columns in the matrix) to new value (third column)
    #NOTE: Layer 3 of Land raster doesn't seem to have any values in it, so it doesn't have to be reclassified
    land<-as.factor(land)
    rat<-levels(land)[[1]]
    rat$landcover<-c('other','shrub','open','forest','agri')
    levels(land)<-rat
  }
  
  mycovars<-as.data.frame(rasterToPoints(stack(tmax,tsd,prec,land)))
  colnames(mycovars)<-c('x','y','tmax','tsd','prec','land')
  mycovars<-na.omit(mycovars)
  
  mycovars$tmax <- as.numeric(mycovars$tmax)
  mycovars$prec <- as.numeric(mycovars$prec)
  mycovars$tsd <- as.numeric(mycovars$tsd)
  mycovars$land<-factor(mycovars$land)
  
  # Separate continuous variables
  df_continuous <- mycovars[, c("x", "y", "tmax", "tsd", "prec")]
  
  # Separate categorical variables
  df_categorical <- mycovars[, c("x", "y", "land")]
  
  # Melt continuous variables
  df_continuous_long <- melt(df_continuous, id.vars = c("x", "y"))
  
  # Melt categorical variables (land)
  df_categorical_long <- melt(df_categorical, id.vars = c("x", "y"))
  
  # Convert 'land' to factor
  df_categorical_long$value <- as.factor(df_categorical_long$value)
  
  # Plot for continuous variables
  plots <- list()
  for (var in c("tmax", "tsd", "prec")) {
    p <- ggplot(subset(df_continuous_long, variable == var), aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_viridis(option = "turbo", limit=c(env_limits[var,2],env_limits[var,1])) +
      coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
      theme_minimal() +
      ggtitle(paste(var))
    plots[[var]] <- p
  }
  
  land_colors <- c(
    "1" = "#B0B0B0",  # Gray for 'Other'
    "2" = "#9ACD32",  # Light Green for 'Shrub'
    "3" = "#FFD700",  # Yellow for 'Open'
    "4" = "#228B22",  # Dark Green for 'Forest'
    "5" = "#A0522D"   # Brown for 'Agriculture'
  )
  
  # Plot for categorical variable (land)
  p_land <- ggplot(subset(df_categorical_long, variable == "land"), aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_manual(values = land_colors) +
    coord_cartesian(xlim = c(min(iberian_map$long), max(iberian_map$long)), ylim = c(min(iberian_map$lat),max(iberian_map$lat))) +
    theme_minimal() +
    ggtitle(paste("land"))
  plots[["land"]] <- p_land
  
  # Arrange all plots in a 2x2 grid
  grid.arrange(
    plots[["tmax"]],
    plots[["tsd"]],
    plots[["prec"]],
    plots[["land"]],
    ncol = 2
  )
  
}