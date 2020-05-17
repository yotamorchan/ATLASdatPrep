library(dplyr)

# merge close ADP locations where their distance is less than adp_rng and time gap less than time_gap
mergeCloseAdp<-function(adp.df,
                        adp_rng=20,
                        smp_rte=4000,
                        time_gap=5*60*1000){
  adp.df<-adp.df%>%arrange(start)
  nn<-nrow(adp.df)
  adp.df$dT<-NA
  adp.df$dT[1:(nn-1)]<-(adp.df$start[2:nn]-adp.df$end[1:(nn-1)])
  adp.df$distance<-NA
  adp.df$distance[1:(nn-1)]<-with(adp.df,
                                    expr={
                                      sqrt((medX[2:nn]-medX[1:(nn-1)])^2+(medY[2:nn]-medY[1:(nn-1)])^2)
                                    })
  
  # gaps big nough between stopping locations
  gap.ix<-which((adp.df$dT[1:(nn-1)]>time_gap) | (adp.df$distance[1:(nn-1)]>adp_rng))
  if ((nn-1) %in% gap.ix){
    gap.ix<-c(gap.ix, nn)
  }
  gap.adp.df<-adp.df[gap.ix,]
  
  # trying to merge close consequent stopping locations
  prev_i<--1
  if (nrow(gap.adp.df)<nn){
    close.adps<-which((adp.df$dT<=time_gap) & (adp.df$distance<=adp_rng))
    merged.adp<-adp.df[0, ]
    isSeq=FALSE
    j<-1
    for (i in close.adps) {
      if (isSeq==TRUE){
        adp.row<-merged.adp[j,]
      }else{
        adp.row<-adp.df[i,]
      }
    
      start<-adp.row$start
      end <- adp.df$end[i+1]
      duration <- end-start
      num_loc <- adp.row$num_loc+adp.df$num_loc[i+1]
      position_qlt <- round((num_loc/(1+duration/smp_rte)),3) # position quality
      medX=mean(adp.row$medX, adp.df$medX[i+1])
      medY=mean(adp.row$medY, adp.df$medY[i+1])
      if (i+2<=nn){
        dT<-adp.df$start[i+2]-end
        distance<- sqrt((adp.df$medX[i+2]-medX)^2+(adp.df$medY[i+2]-medY)^2)
      }else{
        dT<-NA
        distance<-NA
      }
      merged.adp[j,]=c("start"=start,
                       "end"=end,
                       "duration"=duration,
                       "num_loc"=num_loc,
                       "position_qlt"=position_qlt,
                       "medX"=medX,
                       "medy"=medY,
                       "dT"=dT,
                       "distance"=distance)
      
     
      
      if (!is.na(dT)){
        if ((dT>time_gap) | (distance>adp_rng)){
          j<-j+1  
          isSeq=FALSE
        }
        else{
          isSeq=TRUE
        }
      }
      
      prev_i=i
      
    }
    
  }
  ret.adp.df<-rbind(gap.adp.df, merged.adp)
  ret.adp.df<-ret.adp.df%>%arrange(start)
  return(ret.adp.df)
}

# assin ADP information to each localization which fall withing ADP time range
assignADPLoc<-function(df=NULL,
                       adp.df=NULL){
  
  if (is.null(df)){
    print("assignADPLoc(): missing data frames of localizations")
    return(NULL)
  }
  if (is.null(adp.df)){
    print("assignADPLoc(): missing data frames for ADP locations")
    
  }
  df$ADP.ID<-NA
  df$ADP.X<-NA
  df$ADP.Y<-NA
  df$ADP.quality<-NA
  
  for (i in 1:nrow(adp.df)){
    locsubset.ix<-which(between(df$TIME, adp.df$start[i], adp.df$end[i]))
    df$ADP.ID[locsubset.ix]<-i
    df$ADP.X[locsubset.ix]<-adp.df$medX[i]
    df$ADP.Y[locsubset.ix]<-adp.df$medY[i]
    df$ADP.quality[locsubset.ix]<-adp.df$position_qlt[i]
  }
  return(df)
}