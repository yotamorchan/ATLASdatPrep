
SplitNights<-function(Raw,times, obs_start, obs_end){

  # 1. Trim to HOURS of night only (based on obs_start and obs_end) for all tags:
  df<-Raw
  df$dateTime<-as.POSIXct(df$TIME/1000, tz="UTC",origin="1970-01-01")
  df$date<-as.Date(df$dateTime)
  df<-df%>%filter((hour(dateTime)>=obs_start) |(hour(dateTime)<obs_end))
  df$Night<-NA
  
  # 2. Trim EACH tag for it's corresponding capture and off times (from "times" table):  
  TagList<-unique(df$TAG)
  dnight_all<-NULL
  for (j in 1:length(TagList)){
    capture<-as.double(times$capture_unix[times$TAG==TagList[j]]) # time to start recording localizations withouth the car tracks (capture night)
    off<-as.double(times$off_unix[times$TAG==TagList[j]]) # time I recorded the bat no longer moves/the tag has stopped or fallen.
    TAG<-filter(df,(TAG==TagList[j]) & (TIME>=capture) & (TIME<=off) )
    refDate<-as.Date(as.POSIXct(capture/1000, tz="UTC",origin="1970-01-01"))
    TAG$Night<-as.numeric((TAG$date-refDate)+1)
    nextNight.ix<-which(hour(TAG$dateTime)<=obs_end)
    TAG$Night[nextNight.ix]<-TAG$Night[nextNight.ix]-1
    
    dnight_all<-rbind(TAG,dnight_all)
  }
  
  return(dnight_all)
}

