#' @export
#' @author Nick McKay
#' @title Read varve thicknesses from shape file
#' @description Pulls varve thicknesses and error codes from a shapefile counted in GIS
#' @importFrom sf st_read
#' @import dplyr tibble magrittr
#' @param filename the path to a shape file
#' @param codeCol name of the column of varve confidence ratings. Can be a single and or multiple strings in a vector. Default = c("conf","VarveID")
#' @param markerCol name of the column of marker layer names. Can be a single and or multiple strings in a vector. Default = c("Marker","Markers")
#' @param varveTop where is the top of the varve sequence on the shape file. Options are "top","bottom","left","right". Default = "top".
#' @param scaleToThickness Value to scale the sum of the measurements. NA retains the scale in the .shp file.
#' @return A data.frame with the varve data
readVarveShapefile <- function(filename,codeCol = c("conf","VarveID"),markerCol = c("Marker","Markers", "Markers2", "marker", "markers"),tiePointCol = c("TiePoint", "TiePoints", "tiePoint", "tiepoint"), varveTop = "top",scaleToThickness = NA){

  shape <-  sf::st_read(filename,quiet = TRUE,stringsAsFactors = TRUE)

  #old way
  #varveCoords <- tibble::tibble(x1 = unlist(lapply(shape$geometry,"[[",1)),
                            # x2 = unlist(lapply(shape$geometry,"[[",2)),
                            # y1 = unlist(lapply(shape$geometry,"[[",3)),
                            # y2 = unlist(lapply(shape$geometry,"[[",4)))

  varveCoords <- as.data.frame(t(purrr::map_dfc(shape$geometry,magrittr::extract,1:4)))
  names(varveCoords) <- c("x1","x2","y1","y2")

  goodRows <- which(rowSums(is.finite(as.matrix(varveCoords)))==4)#remove nonfinite values




  #pull out varve code metadata
  hasCode <- FALSE
  if(any(names(shape) %in% codeCol)){
    gc <- which(names(shape) %in% codeCol)
    varveCoords$varveCode = shape[[gc]]
    hasCode <- TRUE
  }else{
    warning(paste(filename,"- no varve codes. Other names include:"))
    print(names(shape))
  }

  #pull out marker layers
  hasMarker <- FALSE
  if(any(names(shape) %in% markerCol)){
    gc <- which(names(shape) %in% markerCol)
    varveCoords$markers = shape[[gc[1]]]
    hasMarker <- TRUE
  }else{
    warning(paste(filename,"- no marker layers. Other names include:"))
    print(names(shape))
  }

  #pull out tie points (markers between cores)
  hasTiePoint <- FALSE
  if(any(names(shape) %in% tiePointCol)){
    gc <- which(names(shape) %in% tiePointCol)
    varveCoords$tiePoint = shape[[gc[1]]]
    hasTiePoint <- TRUE
  }else{
    warning(paste(filename,"- no tie points layers. Other names include:"))
    print(names(shape))
  }

  if(varveTop == "top"){
    varves <- varveCoords %>%
      dplyr::arrange(desc(y1)) %>%
      dplyr::mutate(thick = abs(y1-y2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(y1))
  }else if(varveTop == "right"){
    varves <- varveCoords %>%
      dplyr::arrange(desc(x1)) %>%
      dplyr::mutate(thick = abs(x1-x2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(x1))
  }else if(varveTop == "left"){
    varves <- varveCoords %>%
      dplyr::arrange(x1) %>%
      dplyr::mutate(thick = abs(x1-x2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(x1))
  }else if(varveTop == "bottom"){
    varves <- varveCoords %>%
      dplyr::arrange(y1) %>%
      dplyr::mutate(thick = abs(y1-y2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(y1))
}

  if(hasCode & hasMarker){
    out <- dplyr::select(varves,count,thick,varveCode, markers)
  }else if(!hasCode & hasMarker){
    out <- dplyr::select(varves,count,thick, markers)
  }else if(hasCode & !hasMarker){
    out <- dplyr::select(varves,count,thick,varveCode)
  }else{
    out <- dplyr::select(varves,count,thick)
  }

  if(hasTiePoint){
    out$tiePoint <- varves$tiePoint
  }


  #Scale measurements to match section thickness?
  if(!is.na(scaleToThickness)){
    out$thick <- out$thick*(scaleToThickness/sum(out$thick))
  }

  return(out)
}

#' @export
#' @author Nick McKay
#' @title Read multiple varve sequences thicknesses from a directory of shape files
#' @description Reads in a directory of varve thicknesses files, using readVarveShapefile.
#' @param directory the path to to the directory. Leaving it blank will open a file chooser.
#' @param codeCol name of the column of varve confidence ratings. Can be a single and or multiple strings in a vector. Default = c("conf","VarveID")
#' @param markerCol name of the column of marker layer names. Can be a single and or multiple strings in a vector. Default = c("Marker","Markers")
#' @param varveTop where is the top of the varve sequence on the shape file. Options are "top","bottom","left","right". Default = "top".
#' @param scaleToThickness Value to scale the sum of the measurements. NA retains the scale in the .shp file.
#' @return A list of data.frames with the varve data
readVarveDirectory <- function(directory = NULL,codeCol = c("conf","VarveID"),varveTop = "top",markerCol = c("Marker","Markers", "Markers2", "marker", "markers"),scaleToThickness = NA){
if(is.null(directory)){
  print("Select a file in the desired directory")
  directory <- dirname(file.choose())
}

  f2g <- list.files(directory,pattern = "*.shp")
  vOut <- vector(mode = "list",length = length(f2g))
  for(i in 1:length(f2g)){
    vOut[[i]] <- readVarveShapefile(file.path(directory,f2g[i]),codeCol = codeCol,varveTop = varveTop,markerCol = markerCol,scaleToThickness = scaleToThickness)

  }
  names(vOut) <- f2g
return(vOut)
}

#' @export
#' @author Nick McKay
#' @title Heuristically determine varve sequence order by marker layers
#' @description Tries to determine the order of the sequence based on the marker layers
#' @import dplyr
#' @param sequence a list of varve data.frames in any arbitrary order
#' @return A vector of indices with the sequence order
determineMarkerLayerOrder <- function(sequence){

  #take a guess at top and bottom of sequence:

  #pull out marker layers for each
  mls <- sapply(sapply(sequence,"[[","markers"),levels)

  #which ones might be top and bottom?
  tb <- which(sapply(mls,length)==1)
  top <- c()
  bot <- c()
  for(i in tb){
    #find a match
    mms <- unlist(sapply(mls,function(x){which(x==mls[[i]])}))
    if(max(mms)==1){
    top <- c(top,i)
    }
    if(max(mms)>1){
      bot <- c(bot,i)
    }
  }

  #do we have top and bottom?
  if(length(top)==1 & length(bot)==1){#YAY!
    secOrder <- rep(NA,length(sequence))
    secOrder[1] <- top
    secOrder[length(sequence)] <- bot
    secRemaining <- dplyr::setdiff(seq_along(sequence),c(top,bot))
  }else{
    stop("cant identify a top and bottom section")
  }

  #loop through and figure out the rest of the order...

  while(length(secRemaining)>0){
    remSpots <- which(is.na(secOrder))
    i <- min(remSpots)
    mms <- sapply(mls,function(x){ifelse(length(which(mls[[secOrder[i-1]]] %in% x))>0,which(mls[[secOrder[i-1]]] %in% x),0)})
    hasMatch <- setdiff(which(mms>0),secOrder)
    if(length(hasMatch)==1){#only one new match!
      secOrder[i] <- hasMatch
    }else if(length(hasMatch)>1){
      stop("multiple matches of marker layers. Should be a pair.")
    }else{
      stop("cannot find a matching marker layer")
    }
    secRemaining <- dplyr::setdiff(seq_along(sequence),na.omit(secOrder))
  }

return(secOrder)

}

#' @export
#' @author Nick McKay
#' @title Heuristically determine varve sequence order by marker layers
#' @description combines sections into sequence
#' @import dplyr
#' @param sortedSequence a list of varve data.frames in stratigraphic order
#' @return A data.frame of the full varve sequence.
combineSectionsByMarkerLayer <- function(sortedSequence){
#presently tries to handle multiple marker layer options stochastically.


#build a giant, concatenated dataframe
  allCounts <- dplyr::bind_rows(sortedSequence, .id = "section")
  sections <- unique(allCounts$section)

#get all the markerlayer names.
  mln <- na.omit(unique(allCounts$markers))

  tryCatch({
  #see if there are multiple marker layers that connect sections
  nMark <- allCounts %>%
    group_by(section) %>%
    summarize(markerLayers = sum(!is.na(markers))) %>%
    arrange(match(section,sections))

  nMark$expected <- 2
  nMark$expected[c(1,nrow(nMark))] <- 1
  #find sections with "extra" markers
  extra <- which(nMark$markerLayers-nMark$expected > 0)

  if(length(extra)>0){
    #look for similar pairs.
    extraCounts <- dplyr::bind_rows(sortedSequence[extra], .id = "section") %>%
      group_by(markers) %>%
      summarize(nMark = n())

    #summarize into a 2 column matrix, only one is needed.
    extraMarks <- matrix(na.omit(extraCounts$markers[extraCounts$nMark>1]),ncol = 2,byrow = T)


    for(i in 1:nrow(extraMarks)){
      mln <- setdiff(mln,sample(extraMarks[i,],1)) #remove one at random.
    }
  }

  },  error = function(e){cat("Error:", conditionMessage(e), " did not have multiple marker layers in one section. Carry on")})
#start at the top to the first markerlayer
  sp <- min(which(allCounts$markers==mln[1]))
  combSeq <- allCounts[1:sp,]

  #Deal with intermediate sections
  for(n in 2:length(mln)){
    #find second instance of previous markerlayer,
    st <- max(which(allCounts$markers==mln[n-1]))

    #find first instance of next marker layer.
    sp <- min(which(allCounts$markers==mln[n]))

    if(sp <= st){#that's bad, markers out of order
      stop(paste("marker layers",mln[n-1],"&",mln[n],"out of order"))
      }

    combSeq <- bind_rows(combSeq,allCounts[(st+1):sp,])
  }


  #Then append the bottom section
  st <- max(which(allCounts$markers==mln[n]))
  combSeq <- bind_rows(combSeq,allCounts[(st+1):nrow(allCounts),])

#rename and extend tibble
combSeq <- combSeq %>%
  mutate(sectionCounts = count) %>%
  mutate(count = seq_along(count))


return(combSeq)



}


#' @export
#' @author Stephanie Arcusa
#' @title Calibrate cumulative core length to real depth
#' @description Calculates real depth. Will set the new top of core depth to 0.
#' @param CombSeq A data.frame of the full varve sequence, output of combineSectionsByMarkerLayer
#' @param topOfVarvedSequence an integer of the depth of first lamination in mm e.g. 0 mm
#' @param botOfVarvedSequence an integer of the depth of the last lamination in mm
#' @return A data.frame of the full calibrated varve sequence.
calibrateCummulativeDepth <- function(CombSeq, topOfVarvedSequence, botOfVarvedSequence){
  
  test3 <- vector() 
  OldMin = topOfVarvedSequence
  NewMax = botOfVarvedSequence
  NewMin = 0
  
  for(i in 1:nrow(CombSeq)){
    
    OldMax = max(cumsum(CombSeq$thick), na.rm = T)
    OldValue = CombSeq$thick[i]
    
    OldRange = (OldMax - OldMin)  
    NewRange = (NewMax - NewMin)  
    NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
    
    test3[i] <- NewValue
    
  }
  CombSeq$thick <- test3
  CalSeq = CombSeq
  
  return(CalSeq)
  
}

#' @export
#' @author Stephanie Arcusa
#' @title Calibrate cumulative core length to real depth
#' @description Calculates real depth. Will set the new top of core depth to 0.
#' @param data A list. Output of varveModel.
#' @param botOfVarvedSequence an integer of the depth of the last lamination in mm. 
#' @param ML a vector of characters indicating all marker layers in the core
#' @param core an integer indicating the core to use as saved in input "data". Usually the dated core.
#' @param n the number of iterations used for varveModel
#' @return A list of varveModel output with additional items, including calibrated cumulative depth and age. Output can be used to plot age-depth.

calibrateVarveModelDepth <- function(data, botOfVarvedSequence, ML, n, core){
  
  #data[[1]] <- data[[1]]+2 #make values positive 
  find.first.na <- apply(data[[2]][[core]]$thicks,MARGIN = 2, function(x) which(is.na(x))[1])
  find.last.year <- find.first.na-1
  
  #lastTP <- data[[3]][which(ML == tail(ML, 1)),] # the last tie point must also be the end of the varved sequence
  
  data[[4]] <- matrix(NA, nrow = nrow(data[[2]][[core]]$thicks), ncol = n)
  
  OldMin = 0
  NewMax = botOfVarvedSequence
  NewMin = 0
  # Calibrate the last tie point depth to the bottom of the entire varve sequence
  for(j in 1:n){
    
   # if(is.na(lastTP[j]) | lastTP[j] < mean(lastTP, na.rm = T) - sd(lastTP, na.rm = T)){ 
    #  next #if the last TP is not reached, the iteration will be discarded from future analyses; similarly, if the last TP is much younger than the mean, the iteration is also discarded
    #}
    if(find.last.year == 0){
      next #if the last year is 0 this means the iteration is empty and should be skipped
    }
    
    for(i in 1:find.last.year[j]){
    #for(i in 1:lastTP[j]){ used to be this
     
      #tp <- lastTP[j] # SHA added this 23-09-2020 and 1:tp below
      tp <- find.last.year[j]
      OldMax = max(cumsum(data[[2]][[core]]$thicks[1:tp,j]), na.rm = T)
      OldValue = data[[2]][[core]]$thicks[i,j]
      
      OldRange = (OldMax - OldMin)  
      NewRange = (NewMax - NewMin)  
      NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
      
      data[[4]][i,j] <- NewValue
      
    }
    
  }
  
  # Add matrix with the cumulatitve depth
  data[[5]] <- matrix(NA, nrow = nrow(data[[2]][[core]]$thicks), ncol = n)
  for(i in 1:n){
    data[[5]][,i] <- cumsum(data[[4]][,i])
  }
  # Add matrix with age
  data[[6]] <- matrix(NA, nrow = nrow(data[[2]][[core]]$thicks), ncol = n)
  for(i in 1:n){
    
    if(find.last.year[i] == 0 | is.na(find.last.year[i])){
    #if(is.na(lastTP[i])){
      next
    }
    
    #data[[6]][,i] <- c(seq(1,lastTP[i],1), rep(NA,nrow(data[[2]][[core]]$thicks)-length(seq(1,lastTP[i],1))))
    data[[6]][,i] <- c(seq(1,find.last.year[i],1), rep(NA,nrow(data[[2]][[core]]$thicks)-length(seq(1,find.last.year[i],1))))
    
    
  }
  
  names(data) <- c("ModelOutput", "PostInfo","MLpostAge","VarveThick", "CumDepth", "Age")
  cal.data = data
  
  return(cal.data)
  
}
