library(toolsForAtlas)
library(dplyr)
library(leaflet)
library(RColorBrewer)


options(digits=14)

is.POSIXct <- function(x) inherits(x, "POSIXct")

genColorMap<-function(){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  return(col_vector)
}


# generate leaflet map to view one of more localizations set(s) on same map
# if the data frames list has names - it will be use for the layers descrition
map.ll<-function(df.lst=NULL,
                       tags.list=NULL,
                       startTime=NULL,
                       timeDuration=3*3600*1000){
  
  if (is.null(df.lst)){
    print(sprintf("mapLayers.ll(): missing data frames of localizations"))
    return(NULL)
  }
  if (is.null(names(df.lst))){
    print(sprintf("mapLayers.ll(): missing names for data frames in the list - use their order as default"))
    names(df.lst)<-sprintf("df%d",c(1:length(df.lst)))
  }
  df.desc<-names(df.lst)
  if (is.null(tags.list)){
    tags.list<-unique(sapply(df.lst, 
                             function(x){unique(x$TAG)}))
  }

  if(is.null(startTime)){
    startTime=min(sapply(df.lst, 
                         function(x){min(x$TIME)}))
  } else{
    if (is.POSIXct(startTime)){
      startTime=as.double(as.POSIXct(startTime, tz="UTC"))*1000
    }
  }
  endTime=startTime+timeDuration
  color_map<-genColorMap()
 
  ll.maps<-list()
  map.ix=0
  # for each data frame in the list p
  for (tg in tags.list){
    shortTag<- tg%%10000
   
     # Generate leaflet map for each tag
    llix<-leaflet() %>% addTiles(options=tileOptions(opacity=0.5))
    layersList<-c()
    for (i in 1:length(df.lst)){
      tmp.df<-df.lst[[i]]
     
      tmp.df<-as.data.frame(tmp.df%>%
                              filter(TAG==tg)%>%
                              filter(between(TIME, startTime, endTime))%>%
                              arrange(TIME))
    
      if (nrow(tmp.df)==0){
        next
      }
     
      tmp.spdf<-convertSpatial.ITM2WGS84(tmp.df)
      
      cntRows<-nrow(tmp.spdf@data)
      layerName.crcl<-sprintf("%s tag %d [%d rows]",
                         df.desc[i],
                         shortTag,
                         cntRows)
      layerName.ln<-sprintf("%s tag %d lines",
                              df.desc[i],
                              shortTag)
      llix<-llix%>%addCircleMarkers(
        data = tmp.spdf,  #  points data from spdf
        radius = 3,  # cycle radios
        opacity = 1,  # cycle transparency
        fillOpacity = 0.7,
        color = color_map[i],  # point colour
        stroke = FALSE,
        label = ~(sprintf("DateTime=%s",dateTime)),
        popup = ~(sprintf("%s<br>DateTime=%s",
                          df.desc[i],dateTime)),
        group=layerName.crcl)
      llix<-llix%>%addPolylines(
        data=tmp.spdf@coords,
        color=color_map[i],
        weight = 0.7,
        dashArray = "1,1",
        opacity = 0.7,
        group=sprintf("%s tag %d lines",
                              df.desc[i],
                              shortTag))
      layersList<-c(layersList, layerName.crcl,layerName.ln)
    }
    
    # add layers controller
    llix<-llix%>%addLayersControl(
      overlayGroups = layersList,
      options = layersControlOptions(collapsed = FALSE, autoZIndex=TRUE)
    ) %>%
    # add scale bar
    addScaleBar(position = c("bottomleft"),
                  options = scaleBarOptions(imperial=FALSE,maxWidth=200)) 
    map.ix=map.ix+1
    ll.maps[[map.ix]]<-llix
    names(ll.maps)[[map.ix]]=shortTag
   }
   return(ll.maps)
  
}





# generate leaflet map to view single localizations set with ADP position
map.ADP.ll<-function(df=NULL,
                 adp.df=NULL,
                 startTime=NULL,
                 timeDuration=3600*1000){
  
  if (is.null(df)){
    print("map.ADP.ll(): missing data frames of localizations")
    return(NULL)
  }


  # if (is.null(adp.df)){
  #   print("map.ADP.ll(): missing data frames for ADP locations")
  #   
  # }
  tg<-unique(df$TAG)[1]
  
  if(is.null(startTime)){
    startTime=min(df$TIME)
  } else{
    if (is.POSIXct(startTime)){
      startTime=as.double(as.POSIXct(startTime, tz="UTC"))*1000
    }
  }
  endTime=startTime+timeDuration
  start.txt<-as.character(as.POSIXct(startTime/1000, tz="UTC", origin="1970-01-01"),
                          format="%Y-%B-%d %H:%M")
  end.txt<-as.character(as.POSIXct(endTime/1000, tz="UTC", origin="1970-01-01"),
                          format="%H:%M")
  color_map<-genColorMap()
  
  
  shortTag<- tg%%10000
  # Generate leaflet map for tg
  llix<-leaflet() %>% addTiles(options=tileOptions(opacity=0.5))
  layersList<-c()
  
  # atlas localization layer   
  tmp.df<-as.data.frame(df%>%
                        filter(TAG==tg)%>%
                        filter(between(TIME, startTime, endTime))%>%
                        arrange(TIME))
  
  adp.df<-tmp.df%>%
    filter(!is.na(ADP.ID))%>%
    group_by(TAG, ADP.ID, ADP.X, ADP.Y, ADP.quality)%>%
    summarise(start.dateTime=as.POSIXct(min(TIME)/1000, tz="UTC", origin="1970-01-01"),
              end.dateTime=as.POSIXct(max(TIME)/1000, tz="UTC", origin="1970-01-01"),
              start=min(TIME),
              end=max(TIME),
              duration=max(TIME)-min(TIME),
              num_loc=n())%>%
    arrange(start)
  
  cntRows<-nrow(tmp.df)
  if (cntRows==0){
      print ("map.ADP.ll(): tag %d does not have localizatiosn between %s-%s",
             shortTag, start.txt, end.txt)
  } else{
        
    tmp.spdf<-convertSpatial.ITM2WGS84(tmp.df)
        
    
    layerName.crcl.move<-sprintf("tag %d move [%d rows]",
                            shortTag,
                            nrow(subset(tmp.df, is.na(ADP.ID))))
    layerName.crcl.adp<-sprintf("tag %d stopped [%d rows]",
                            shortTag,
                            nrow(subset(tmp.df, !is.na(ADP.ID))))
  
    layerName.ln<-sprintf("tag %d lines",
                          shortTag)
    llix<-llix%>%addCircleMarkers(
      data = subset(tmp.spdf,is.na(ADP.ID)),  #  points data from spdf
      radius = 3,  # cycle radios
      opacity = 1,  # cycle transparency
      fillOpacity = 0.5,
      color="#8856a7",
      fillColor = color_map[1],  # point colour
      stroke = TRUE,
      weight =0.8,
      popup = ~(sprintf("DateTime=%s",dateTime)),
      group=layerName.crcl.move)
    
    llix<-llix%>%addCircleMarkers(
      data = subset(tmp.spdf, !is.na(ADP.ID)),  #  points data from spdf
      radius = 3,  # cycle radios
      opacity = 1,  # cycle transparency
      fillOpacity = 0.5,
      fill=TRUE,
      fillColor =color_map[1], 
      color="#2c7fb8",
      stroke = TRUE,
      weight =0.8,
      popup = ~(sprintf("DateTime=%s<br>ADP#%d<br>ADP quality=%.2f",
                        dateTime, ADP.ID, ADP.quality)),
      group=layerName.crcl.adp)
    
    
    
    llix<-llix%>%addPolylines(
      data=tmp.spdf@coords,
      color=color_map[1],
      weight = 1,
     # dashArray = "1,1",
      opacity = 0.7,
      group=sprintf("tag %d lines",
                    shortTag))
   
    
    
    layersList<-c(layerName.crcl.move,layerName.crcl.adp,layerName.ln )
  
  }

  # add layer for ADP 
  
  # adp.df$start.dateTime<-as.POSIXct(adp.df$start/1000, tz="UTC", origin="1970-01-01")
  # adp.df$end.dateTime<-as.POSIXct(adp.df$end/1000, tz="UTC", origin="1970-01-01")
  # tmp.adp.df<-as.data.frame(adp.df%>%
  #                           filter(TAG==tg)%>%
  #                           filter(between(start, startTime, endTime))%>%
  #                           arrange(start))
  cntRows<-nrow(adp.df)
  if (cntRows==0){
    print ("map.ADP.ll(): tag %d does not have ADP locations between %s-%s",
           shortTag, start.txt, end.txt)
  } else{
    
    tmp.adp.spdf<-convertSpatial.ITM2WGS84(adp.df, xyColNames = c("ADP.X", "ADP.Y"))
    
    
    layerName.crcl<-sprintf("ADP [%d rows]",
                            cntRows)
    llix<-llix%>%addCircleMarkers(
      data = tmp.adp.spdf,  #  points data from spdf
      radius = 7,  # cycle radios
      opacity = 1,  # cycle transparency
      fillOpacity = 0.9,
      color = "red",  # point colour
      stroke = FALSE,
      label = ~(sprintf("ADP start=%s",start.dateTime)),
      popup = ~(sprintf("ADP#%d<br>from %s to %s<br>duration=%.2f min<br>quality=%.2f<br>#loc=%d",
                        ADP.ID,
                        as.character(start.dateTime, format="%H:%M:%S"),
                        as.character(end.dateTime, format="%H:%M:%S"),
                        (duration/60000),
                        ADP.quality,
                        num_loc)),
      group=layerName.crcl)
    layersList<-c(layersList, layerName.crcl)
  }
  

    # add layers controller
    llix<-llix%>%addLayersControl(
      overlayGroups = layersList,
      options = layersControlOptions(collapsed = FALSE, autoZIndex=TRUE)
    ) %>%
      # add scale bar
      addScaleBar(position = c("bottomleft"),
                  options = scaleBarOptions(imperial=FALSE,maxWidth=200))
  return(llix)
  
}



