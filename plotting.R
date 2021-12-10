#' Plot a bunch of varve sections
#'
#' @param varveSequenceList
#' @param output A pdf file to concatenate the output into
#' @import purrr pdftools
#'
#' @return a
#' @export
plotSections <- function(varveSequenceList,output = NA){
 
  for(i in 1:length(varveSequenceList)){
    
    if(anyNA(varveSequenceList[[i]]$varveCode)){
      print(paste0("check slide ", names(varveSequenceList[i])," because varveID looks wrong around varve ", which(is.na(varveSequenceList[[i]]$varveCode))))
    }
    
  }
  
  allPlots <- purrr::map2(varveSequenceList,names(varveSequenceList),plotVarveSection)
 if(!is.na(output)){
   purrr::map2(paste0(file.path(tempdir(),names(varveSequenceList)),".pdf"),allPlots,ggsave)
   apdf <- list.files(tempdir(),pattern = "*.shp.pdf",full.names = T)
   pdftools::pdf_combine(apdf,output = output)
 }
 return(allPlots)
}


#' Plot varve section
#'
#' @param section
#' @param sectionName
#' @param varveCodeScale
#' @import ggplot2 dplyr pracma
#' @return a ggplot list
#' @export
#'
#' @examples
plotVarveSection <- function(section,sectionName, varveCodeScale = 0.8){
  #scale varve codes to counts
  section$scaledVarveCode <-  varveCodeScale*section$varveCode*max(section$thick)/max(section$varveCode)
  svcLabels <- c(as.character(unique(sort(section$varveCode))),"Varve Codes")
  svcLabelY <- varveCodeScale*unique(sort(section$varveCode))*max(section$thick)/max(section$varveCode)
  svcLabelY <- c(svcLabelY,max(svcLabelY)*1.1)
  svcLabelX <- quantile(section$count,probs = 1)

# nearest neighbor varve codes for plotting.
  iCount <- seq(from = min(section$count),to = max(section$count),by = 0.1)
  nn <- data.frame(count = iCount,
                   varveCode = pracma::interp1(section$count,section$scaledVarveCode,xi = iCount,method = "nearest" ))


  #get marker layers for plotting
  ml <- filter(section, !is.na(markers))


  secPlot <- ggplot()+
    #plot varve codes
    geom_area(data = nn, aes(x = count, y = varveCode),fill = "Gray",alpha = 0.5)+
    #add varve code scale
    geom_label(aes(x= svcLabelX,y = svcLabelY, label = svcLabels),colour = "Gray")+
    geom_line(data = section, aes(x = count, y = thick))+
    ggtitle(sectionName)+
    ylab("Varve Thickness (arbitary units)")+
    xlab("Count")+
    theme_bw()

  if(nrow(ml)){
    secPlot <- secPlot+
    #add marker layers
    geom_vline(data = ml, aes(xintercept = count),color = "red")+
    geom_label(data = ml, aes(x= count,y = 0, label = markers),colour = "Red")
  }


  return(secPlot)


}

#' plot a composite sequence
#'
#' @param compSeq
#' @import RColorBrewer
#' @return a plot
#' @export
plotCompositeSequence <- function(compSeq){
  colourCount = length(unique(compSeq$section))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

  return(ggplot(compSeq)+geom_line(aes(x = count, y = thick,color=section))+theme_bw()+scale_colour_manual(values = getPalette(colourCount)))
}

#' @title Plots a varve age depth model with uncertainty
#' @author Stephanie Arcusa
#' @description The function will plot the varve age depth model with uncertainty  
#' @param data The output of calibrateVarveModel.
#' @return An age-depth ggplot object

plotVarveModel <- function(data){
  
  #Remove rows with all NAs
  data$VarveThick <- data$VarveThick[rowSums(is.na(data$VarveThick)) != ncol(data$VarveThick),colSums(is.na(data$VarveThick))<nrow(data$VarveThick)]
  count <- seq_len(nrow(data$VarveThick))
  
  ageEns <- createEnsembleAgeDepthModel(data$VarveThick)
  plotTable <- as.data.frame(ageEns$summaryTable)
  plotTable <- plotTable[-nrow(plotTable),]
  names(plotTable) <- c("depth", "c2.5", "c25", "c50", "c75", "c97.5")
  
  VarvePlot <- ggplot(plotTable)+
    geom_ribbon(aes(x = depth,ymin = c2.5  ,ymax = c97.5),
                fill = "gray80")+
    geom_ribbon(aes(x = depth,ymin = c25  ,ymax = c75), fill = "gray50")+ 
    geom_line(aes(x = depth, y = c50))+
    scale_x_reverse()+
    ylab(label = "Varve years")+
    xlab(label = "Depth (mm)")+
    theme_bw()+
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    coord_flip()
    
  
  return(VarvePlot)
  
}

