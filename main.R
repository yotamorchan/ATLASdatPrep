# kingfisher tutorial 

library(toolsForAtlas) # useful functions for atlas users
library(dplyr)         # give nice grammar for datamanipulations
options("digits"=14)   # show long numbers till 14 

source("functions\\visualMaps.R")

dataFolderName="data"
birdCode="kf1845"
srceFileName<-sprintf("%s\\%s.sqlite",dataFolderName, birdCode)
raw.df<-loadFromSQLite(dbn=srceFileName,
                       tbl="localizations")
colnames(raw.df)
dataColumns<-which(colnames(raw.df) %in% c("TAG","TIME","X","Y","VARX","VARY","COVXY"))

head(raw.df[, dataColumns])
raw.df<-raw.df[, dataColumns]
## data preparation
raw.df$TIME<-as.double(raw.df$TIME)
# as.POSIXct() return epoch time in seconds. converting to milisec by x1000
tagTime<-as.numeric(as.POSIXct("2017-03-19 10:00:00", tz="UTC"))*1000
# using dplyr filter to remove rows with early time
raw.df<-raw.df%>%filter(TIME>=tagTime)

# filter validation : show the minimume time
min(raw.df$TIME)
#count "localizations" with round X, Y (in hula are X=257000, Y=78000) - they are not "real localizations"
fake.ix<-which(((raw.df$X-trunc(raw.df$X))+(raw.df$Y-trunc(raw.df$Y)))==0)
if (length(fake.ix) >0) {
  print ("in raw data %d localizations, suspeced as fake and removed")
  raw.df<-raw.df[-fake.ix,]
}

## monitor filters yield
recordAnalysis.lst<-list()
recordAnalysis.lst[["raw data"]]<-data.frame("nrow"=nrow(raw.df))

# get number of rows raw data frame
recordAnalysis.lst[["raw data"]]$nrow
df1<-addLocAttribute(raw.df, locAttributs=c("locQuality"))

# show first rows with added attributes
format.data.frame(head(df1), digits=5)
quantile(df1$stdErrXY, c(0.5,0.6,0.7,0.8,0.9, 0.95,0.99))
df1<-df1%>%filter(stdErrXY<100)
recordAnalysis.lst[["filter stdErrXY"]]<-data.frame("nrow"=nrow(df1))
df1<-addLocAttribute(df1, locAttributs=c("speed"))
quantile(df1$spd, c(0.5,0.6,0.7,0.8,0.9,0.95,0.99), na.rm=TRUE)

# define parameters thresholds:
#  v_filter_max_v =  maximum speed (m/sec)
#  min_time_dif_for_stats = minimum number of good localizations - refered as "reliable"
#  v_filter_time_dif = time gap the algorithm is rubust to.
optionsArg=c("v_filter_max_v"=20,
             "min_time_dif_for_stats"=10,
             "v_filter_time_dif"=10)
df2<-filterByVelocity(df1,
                      dataLegend=c("X","Y","TIME"),
                      options=optionsArg)
maps.1<-map.ll(df.lst=list("stdErrXY fltr"=df1,
                           "velocity fltr"=df2),
               timeDuration=2*3600*1000)
names(maps.1)
maps.1[["1845"]]

xyt<-triangularSmooth(df2$X, df2$Y, df2$TIME)
df3<-df2[,c("TAG","X","Y","TIME")]
df3[,c("X","Y","TIME")]<-xyt[,c("X","Y","TIME")]
df3<-addLocAttribute(df3, locAttributs=c("speed","angle"))
maps.2<-map.ll(df.lst=list("velocity fltr"=df2,
                           "triangularSmooth"=df3),
               timeDuration=2*3600*1000)
(maps.2[["1845"]])

# check if exists subfolder names by this day
# if not exists, yet  - create one.
dateStr<-as.character(Sys.Date(), "%Y-%b-%d")
saveFolderName<-sprintf("results\\stg1\\%s", dateStr)
if (!dir.exists(saveFolderName)){
  dir.create(file.path(saveFolderName))
}

#save csv format
saveFileName<-sprintf("%s\\stg1_%s.csv",saveFolderName, birdCode )
write.csv(df3, saveFileName,  row.names = FALSE )

loadFileName<-sprintf("results\\stg1_%s.csv",birdCode )
kf.df<-read.csv(loadFileName)

# for each tag should run AdpFixedPoint
kf.adp<-AdpFixedPoint(
  time.vec=kf.df$TIME,
  x=kf.df$X,
  y=kf.df$Y,
  adp_rng=20,
  smp_rte =4000,
  obs_min=12,
  p_lim=6,
  time_gap =5 * 60 * 1000 )
kf.adp<-kf.adp[which(kf.adp$duration>1),]  # when duration =1 is not real ADP
head(kf.adp)



