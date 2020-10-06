#

#' create cores data.frame for varve model
#'
#' @param ensOut output of generateThicknessEnsemble()
#' @param ocPrior
#' @param ucPrior
#' @param mlPrior
#'
#' @return
#' @export
createCoreInputForModel <- function(ensOut,ocPrior = 0.05, ucPrior = 0.05, mlPrior = 0.80){

  #pick a random column
  if(NCOL(ensOut$ensThick) > 1){ #09-19-20  should it be < SHA
    r <- sample.int(NCOL(ensOut$ensThick),1)
    #create a data.frame
    core <- data.frame(thickness = ensOut$ensThick[,r],markerLayers = as.character(ensOut$tiePoints[,r]))
  }else{
    core <- data.frame(thickness = ensOut$ensThick,markerLayers = as.character(ensOut$tiePoints))

  }
  #assign in prioris
  core$ocPrior <- ocPrior
  core$ucPrior <- ucPrior
  core$mlPrior <- mlPrior

  #create a  simple year vector
  core$year <- seq_along(core$thickness)


  #FOR NOW, CORES have to have a hard top! Top marker layer must have a probability of 1
  #this can be fixed, but needs some work.
  #make sure the top year is a marker layer
  core$markerLayers <- as.character(core$markerLayers)
  core$markerLayers[1] <- "topOfCore"
  core$mlPrior[1] <- 1

  core <- dplyr::filter(core,!is.na(thickness))

  return(core)
}




#' model varve thickness
#'
#' @param ensList
#' @param nSim
#' @param allMarkerLayers
#' @import matrixStats purrr
#' @return
#' @export
#'
#' @examples
varveModel <- function(ensList,nSim = 100,  allMarkerLayers = NA){


  #get one core to initialize
  cores <- purrr::map(ensList,createCoreInputForModel)


  nCores <- length(cores)

  #need a list of all marker layers in order
  if(all(is.na(allMarkerLayers))){
    allMarkerLayers <- as.character(na.omit(unique(c(purrr::flatten(ensList)$tiePoints))))
  }

  #setup Post of markerLayers
  MLpostAge = matrix(NA,nrow=length(allMarkerLayers),ncol=nSim)

  #setup posterior matrix storage.
  corePost = vector("list",length = nCores)
  for(c in 1:nCores){
    corePost[[c]]$mlPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    corePost[[c]]$ucPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    corePost[[c]]$ocPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    corePost[[c]]$meanCorr=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    corePost[[c]]$thicks=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    corePost[[c]]$thick2=c()
    corePost[[c]]$mlCum=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
  }


  totYears=sapply(cores,FUN = function(x)length(x$thickness))

  #Setup composites matrix
  compMat = matrix(NA,nrow = round(2*max(totYears)),ncol = nSim)
  #Simulate segments, update priors and create composites.

  for(i in 1:nSim){
    #sample from ensemble
  max.attempts = 5 # SHA
  
  for(iii in seq_len(max.attempts)) {
    
    OK <- tryCatch({
    
    cores <- purrr::map(ensList,createCoreInputForModel)

    #get some info about the thicknesses
    startYears = sapply(cores,FUN = function(x)min(x$year))
    endYears = sapply(cores,FUN = function(x)max(x$year))

    #make a matrix
    overlapMat = sapply(cores,FUN = function(x)x$thickness[which(x$year==max(startYears)):which(x$year==min(endYears))])

    #any NAs
    overlapMat = overlapMat[which(!apply(is.na(overlapMat),MARGIN = 1,any)),]


    meanThick=colMeans(overlapMat)
    sdThick=colSds(overlapMat)

    totYears=sapply(cores,FUN = function(x)length(x$thickness))

    #Update priors based on years between marker layers
    MLsim=list()

    #Setup storage for new prior information
    for(c in 1:nCores){
      cores[[c]]$bigOCA=matrix(data = 0,nrow = length(cores[[c]]$ocPrior),ncol = 1)
      cores[[c]]$bigUCA=matrix(data = 0,nrow = length(cores[[c]]$ucPrior),ncol = 1)
    }


    #run ML simulation
    for(c in 1:nCores){
      isML=which(!is.na(cores[[c]]$markerLayers))
      realTest=runif(length(isML))<cores[[c]]$mlPrior[isML]
      MLsim[[c]]=cores[[c]]$markerLayers
      MLsim[[c]][isML[which(!realTest)]]=NA
    }

    #identify which tie points are found in all cores in this simulation
    allML = lapply(cores,FUN = function(x)x$markerLayers[which(!is.na(x$markerLayers))])
    lastML = sapply(cores,FUN = function(x)x$markerLayers[max(which(!is.na(x$markerLayers)))])
    firstML = sapply(cores,FUN = function(x)x$markerLayers[min(which(!is.na(x$markerLayers)))])

    #create matrix of where we have marker layers
    simML = matrix(FALSE,ncol = nCores,nrow= length(allMarkerLayers))
    for(a in 1:length(allMarkerLayers)){
      for(c in 1:nCores){
        if(length(which(allML[[c]]==allMarkerLayers[a])) > 0 ){
          simML[a,c]=TRUE
        }
      }
    }

    notEnough <- allMarkerLayers[rowSums(simML)<2]

    #now fill in FALSEs that are surrounded by TRUES
    for(c in 1:nCores){
      #which are false
      isF=which(!simML[,c])
      for(f in isF){
        if(f<max(which(simML[,c])) &  f>min(which(simML[,c]))){
          simML[f,c]=TRUE
        }
      }
    }

    nMLsim = cbind(allMarkerLayers,apply(cbind(rowSums(simML),rep(2,times = length(rowSums(simML)))),1,max))



    TPlist = lapply(MLsim,unique)
    allTP =unique(unlist(TPlist))
    allTP=allTP[which(!is.na(allTP))]
    inAll=c()
    inTwo=c()
    for (t in 1:length(allTP)){
      inAll[t]=all(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))
      inTwo[t]=sum(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))==nMLsim[which(allTP[t]==nMLsim[,1]),2]

    }

    #which Marker Layers are we simulating here
    ML2sim=as.vector(allTP[which(inTwo)])
    ML2sim=na.exclude(ML2sim)
    MLyears = matrix(NA,nrow = length(ML2sim),ncol = nCores)
    for(m in 1:nrow(MLyears)){
      allMLList = (sapply(cores,FUN = function(x)which(x$markerLayers==ML2sim[m])))
      for(c in 1:nCores){
        thisML = allMLList[[c]]
        #remove MLs that didn't pass
        if(length(thisML==1)){
          MLyears[m,c] = thisML
        }
      }
    }

    #make sure that the MLyears is in the correct order!!!
    ody = apply(MLyears,2,diff)
    if(any(ody[is.finite(ody)]<0)){
      stop("There are reversals in your marker layers. Check your input files.")
    }
    #print(MLyears)


    #WARNING: AS CURRENTLY CONSTRUCTED, THIS ONLY MAKES SENSE FOR SYMMETRIC OC UC PRIORS!
    #calculate diffs
    dML=diff(MLyears)

    #ASSUME A NORMAL DISTRIBUTION FOR THE DURATION of the YEARS and DRAW FROM IT
    #THIS MAY NOT BE IDEAL!!!!!!!!!
    diffMeans = rowMeans(diff(MLyears),na.rm=TRUE)
    diffSds = rowSds(diff(MLyears),na.rm=TRUE)
    #simTrueDiff = rnorm(n = length(dML),mean = diffMeans,sd = diffSds)
    simTrueDiff = diffMeans
    OCUC = dML/simTrueDiff
    #rowSds(diff(MLyears))
    OCadd = OCUC-1
    OCadd[OCadd<=0]=0
    UCadd = -1*(OCUC-1)
    UCadd[UCadd<=0]=0


    for (c in 1:nCores){
      for (mi in 1:(nrow(MLyears)-1)){
        mlstart=which(ML2sim[mi]==cores[[c]]$markerLayers )
        mlstop=which(ML2sim[mi+1]==cores[[c]]$markerLayers )
        if(length(mlstart) ==1 & length(mlstop)==1){
          #find index over which to update priors
          pri=which(ML2sim[mi]==cores[[c]]$markerLayers ):which(ML2sim[mi+1]==cores[[c]]$markerLayers )
          cores[[c]]$bigOCA[pri,1]=OCadd[mi,c]
          cores[[c]]$bigUCA[pri,1]=UCadd[mi,c]
        }
      }
    }




    #use simulations on marker layers to get distributions on yearly priors.
    for (c in 1:nCores){
      #this way makes single estimates based on the mean. A distribuiton is better.
      # cores[[c]]$ocPrior=rowSums(cbind(cores[[c]]$ocPrior,rowMeans(cores[[c]]$bigOCA,na.rm=TRUE)),na.rm=TRUE)
      # cores[[c]]$ucPrior=rowSums(cbind(cores[[c]]$ucPrior,rowMeans(cores[[c]]$bigUCA,na.rm=TRUE)),na.rm=TRUE)

      #create a matrix of priors
      cores[[c]]$ocPriorMat = cores[[c]]$ocPrior+cores[[c]]$bigOCA
      cores[[c]]$ucPriorMat = cores[[c]]$ucPrior+cores[[c]]$bigUCA



    }

    #empty out this year and last year
    lastMLyear = 0
    thisMLyear = 0

    #make sure that there are no reversals!!!
    #MLyears
    ody = apply(MLyears,2,diff)

    #store Marker Layer posterior data
    for(c in 1:nCores){
      for(ml in 1:length(ML2sim)){
        wmlp=which(ML2sim[ml]==cores[[c]]$markerLayers)
        corePost[[c]]$mlPost[wmlp,i]=ML2sim[ml]
      }
    }

    #create a vector that gives instructions about matching seg lengths
    mustMatchLengths = matrix(TRUE,nrow=(nrow(MLyears)-1))

    #add to ML years if need be.
    #does the first segment start at 1? (is the surface a marker layer?)
    ML1test = MLyears[1,!is.na(MLyears[1,])]==1
    if (!all(ML1test) & any(ML1test)){
      stop("You cant have a marker layer that starts at 1 for some cores and not all of them")
    }
    if (!all(ML1test) ){
      #if no top seg, add one
      MLyears=rbind(matrix(1,ncol=nCores),MLyears)
      mustMatchLengths=rbind(TRUE,mustMatchLengths)#ALERT!!!!! THIS MEANS THAT WE'VE GIVEN A prob of 1 to the top layer!
    }

    lastMLtest = MLyears[nrow(MLyears),!is.na(MLyears[nrow(MLyears),])]
    totYearTest = totYears[!is.na(MLyears[nrow(MLyears),])]
    if(any((totYearTest - lastMLtest)==0) & !all((totYearTest - lastMLtest)==0)){
      stop("Your last marker layer must either be the oldest years, or at least 2 years a way from your oldest years")
    }
    if(!any((totYearTest - lastMLtest)==0)){
      #add in the bottom MLyears layer
      MLyears=rbind(MLyears,totYears)
      ML2sim=append(ML2sim,"bottom")

      mustMatchLengths=rbind(mustMatchLengths,FALSE)

    }

    #set up segments
    seg=vector("list",nrow(MLyears)-1)


    for(s in 1:(nrow(MLyears)-1)){
      for(c in 1:nCores){
        if(s==(nrow(MLyears)-1) | !is.na(MLyears[s,c]) & !is.na(MLyears[s+1,c])){
          if(s==1){
            seg[[s]][[c]]=cores[[c]][(MLyears[s,c]):MLyears[s+1,c],]#include the first marker layer the first year

          }else if(s< (nrow(MLyears)-1)){
            seg[[s]][[c]]=cores[[c]][(MLyears[s,c]+1):MLyears[s+1,c],]
          }else{#IF ITS THE LAST ROW, make the segment from the last, non-NA marker layer
            lastRealMarker = max(MLyears[-nrow(MLyears),c],na.rm=T)
            seg[[s]][[c]]=cores[[c]][(lastRealMarker+1):MLyears[s+1,c],]
          }
        }
      }
    }
    segTotalThicks = matrix(NA,nrow = 5,ncol=length(seg))
    segTotalYears = matrix(NA,nrow = 5,ncol=length(seg))

    ###CHECK THAT SEGMENT THICKNESS ARE CORRECT
    for(c in 1:nCores){
      for(s in 1:(nrow(MLyears)-1)){
        if(length(seg[[s]])>=c){
          segTotalThicks[c,s] = sum(seg[[s]][[c]]$thickness)
          segTotalYears[c,s] = length(seg[[s]][[c]]$year)
        }
      }
    }


    #   if(!all(sapply(seg[[1]],FUN = function(x)min(x$year)) == 1)){
    #     #make a new marker layer that goes to the top, that doesn't have the same length requirement.
    #     topSeg[[c]] = cores[[c]][MLyears[s,c]:MLyears[s+1,c],]
    #   }
    #segment

    newSeg=vector("list",nrow(MLyears)-1)
    #simulate between segments
    tSim=vector("list",nCores)

    composite = c()
    for(s in 1:length(seg)){#for each segment...
      #which cores are we simulating?
      core2sim=which(!sapply(seg[[s]],is.null))
      orig.core2sim = core2sim

      sameLength = FALSE
      t=0
      while(!sameLength){#
        t=t+1

        #build in failsafe
        if(t>1e4){
          if(length(core2sim)>2){
            print(paste("segment",as.character(s),"hit count limit.", as.character(length(core2sim)),"cores. Reducing number of cores in simulation"))
            core2sim = core2sim[-sample(core2sim)[1]]#remove one of the cores.
            t=0
          }else{
            stop("Too many iterations..no more cores to remove.. stopped.")
          }
        }
        #DO YAHTZEEE HERE
        for(c in core2sim){
          #simulate new segments given the UC and OC priors.
          soacDf <- data.frame(thick = seg[[s]][[c]]$thickness, ucProb = seg[[s]][[c]]$ucPriorMat, ocProb = seg[[s]][[c]]$ocPriorMat)
          tSim[[c]]= simulateOverAndUndercounting(soacDf)
        }

        newSegLengths = sapply(tSim,FUN = function(x)length(x$newThicks))[core2sim]
        sameLength= (all(newSegLengths==newSegLengths[1]) | !mustMatchLengths[s])


      }

      #see if any didn't run
      didnt.sim = orig.core2sim[!orig.core2sim %in% core2sim]
      if(length(didnt.sim)>0){#if any didn;t run...
        for(dd in didnt.sim){#one by one, simulate the length for any that got excluded

          ddt=0
          while(TRUE){#and simulate until you match the length
            ddt=ddt+1
            soacDf <- data.frame(thick = seg[[s]][[dd]]$thickness, ucProb = seg[[s]][[dd]]$ucPriorMat, ocProb = seg[[s]][[dd]]$ocPriorMat)
            tSim[[dd]]= simulateOverAndUndercounting(soacDf)
            #check to see if it's the correct length
            thisSegLength = length(tSim[[dd]]$newThicks)
            if(thisSegLength == newSegLengths[1]){#found it!!!
              core2sim = append(core2sim,dd)
              newSegLengths = append(newSegLengths,thisSegLength)
              break
            }
            if(ddt>1e5){#give up
              print(paste("gave up trying to get core",as.character(dd),"to be",as.character(newSegLengths[1]),"years long for segment",as.character(s)))
              break
            }
          }
        }
      }





      #Once same length is matched build new data.frame
      #And assign results into posterior.
      segthickmat=matrix(NA,nrow=max(newSegLengths),ncol=length(core2sim))
      for(c in 1:length(core2sim)){
        #build Thickmat
        segthickmat[1:newSegLengths[c],c] = tSim[[core2sim[c]]]$newThicks
      }
      cm=cor(segthickmat,use="complete.obs")
      cm[which(cm<0)]=0
      meanCorr=mean(cm[lower.tri(cm)])

      ###ASSIGN IN POSTERIOR INFORMATION###
      for(c in core2sim){
        #assign posterior info
        #Undercounted
        ucIndex = tSim[[c]]$UCi+min(seg[[s]][[c]]$year)-1
        corePost[[c]]$ucPost[ucIndex,i]=1
        #Overcounted
        ocIndex = tSim[[c]]$OCi+min(seg[[s]][[c]]$year)-1
        corePost[[c]]$ocPost[ocIndex,i]=1

        #mean correlation
        corePost[[c]]$meanCorr[min(seg[[s]][[c]]$year) : max(seg[[s]][[c]]$year) ,i] = meanCorr

        #store thicknesses for each core
        if(mustMatchLengths[s]){
          if(c==core2sim[1]){#only update if its the first one.
            thisMLyear = lastMLyear+newSegLengths[1]
          }
          corePost[[c]]$mlCum[thisMLyear,i]=ML2sim[s+1]
          #test if sum thicknesses are the same as the original (THEY SHOULD BE!)
          if(abs(sum(tSim[[c]]$newThicks)-sum(seg[[s]][[c]]$thickness))>.001){
            stop("The thicknesses don't sum to the same number!")
          }
          corePost[[c]]$thick2=append(corePost[[c]]$thick2,tSim[[c]]$newThicks)
          corePost[[c]]$thicks[(lastMLyear+1) : thisMLyear ,i] = tSim[[c]]$newThicks


        }else{
          lastSimYear = max(which(!is.na(corePost[[c]]$thicks[,i])))
          corePost[[c]]$thicks[(lastSimYear+1) : (lastSimYear+newSegLengths[which(c==core2sim)]) ,i] = tSim[[c]]$newThicks
        }
      }


      #is the total thickness the same as the total thickness should be?
      curThick=sum(corePost[[c]]$thicks[1 : sum(!is.na(corePost[[c]]$thicks[,i])) ,i] )
      lowerMLname= seg[[s]][[c]]$markerLayers[length(seg[[s]][[c]]$markerLayers)]

      #assign age into MLpostAge
      if(!is.na(lowerMLname)){
        wML = which(lowerMLname==allMarkerLayers)
        MLpostAge[wML,i] = thisMLyear
      }

      origMLi = which(cores[[c]]$markerLayers==lowerMLname)
      if(length(origMLi)<1){
        origMLi = length(cores[[c]]$thickness)
      }
      origThick = sum(cores[[c]]$thickness[1 : origMLi])
      # if(abs(origThick-curThick)> .001){
      #   stop("The current cumulative thickness isn't right")
      # }


      lastMLyear=thisMLyear#update with new last year

      #build Composite
      segCom= rowMeans(scale(segthickmat,center=meanThick[core2sim],scale =sdThick[core2sim]),na.rm = TRUE)
      composite = append(composite,segCom)
    }
    print(i)
    compMat[1:length(composite),i]=composite

    }, #Added by SHA
    error = function(e) {
      FALSE
    }
    )
    if(OK) {
      
      next
      
   }
 }

    
    
  }
  
  #return the output
  out <- list(compMat,corePost,MLpostAge)

  return(out)


}






#' model varve thickness, 1-iteration
#'
#' @param ensList
#' @param allMarkerLayers
#' @import matrixStats purrr
#' @description this is designed for bayesian applications.
#' @return
#' @export
#'
#' @examples
fastVarveModel <- function(ensList, allMarkerLayers = NA){
  nSim <- i <- 1

  #get one core to initialize
  cores <- purrr::map(ensList,createCoreInputForModel)


  nCores <- length(cores)

  #need a list of all marker layers in order
  if(all(is.na(allMarkerLayers))){
    allMarkerLayers <- as.character(na.omit(unique(c(purrr::flatten(ensList)$tiePoints))))
  }

  #setup Post of markerLayers
  MLpostAge = matrix(NA,nrow=length(allMarkerLayers),ncol=nSim)

  # #setup posterior matrix storage.
  corePost = vector("list",length = nCores)
  for(c in 1:nCores){
    #   corePost[[c]]$mlPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    #   corePost[[c]]$ucPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    #   corePost[[c]]$ocPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    #   corePost[[c]]$meanCorr=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    corePost[[c]]$thicks=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
    #corePost[[c]]$thick2=c()
    #   corePost[[c]]$mlCum=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim)
  }
  #

  totYears=sapply(cores,FUN = function(x)length(x$thickness))

  #Setup composites matrix
  compMat = matrix(NA,nrow = round(2*max(totYears)),ncol = nSim)
  #Simulate segments, update priors and create composites.

  #get some info about the thicknesses
  startYears = sapply(cores,FUN = function(x)min(x$year))
  endYears = sapply(cores,FUN = function(x)max(x$year))

  #make a matrix
  overlapMat = sapply(cores,FUN = function(x)x$thickness[which(x$year==max(startYears)):which(x$year==min(endYears))])

  #any NAs
  overlapMat = overlapMat[which(!apply(is.na(overlapMat),MARGIN = 1,any)),]


  meanThick=colMeans(overlapMat)
  sdThick=colSds(overlapMat)

  totYears=sapply(cores,FUN = function(x)length(x$thickness))

  #Update priors based on years between marker layers
  MLsim=list()

  #Setup storage for new prior information
  for(c in 1:nCores){
    cores[[c]]$bigOCA=matrix(data = 0,nrow = length(cores[[c]]$ocPrior),ncol = 1)
    cores[[c]]$bigUCA=matrix(data = 0,nrow = length(cores[[c]]$ucPrior),ncol = 1)
  }


  #run ML simulation
  for(c in 1:nCores){
    isML=which(!is.na(cores[[c]]$markerLayers))
    realTest=runif(length(isML))<cores[[c]]$mlPrior[isML]
    MLsim[[c]]=cores[[c]]$markerLayers
    MLsim[[c]][isML[which(!realTest)]]=NA
  }

  #identify which tie points are found in all cores in this simulation
  allML = lapply(cores,FUN = function(x)x$markerLayers[which(!is.na(x$markerLayers))])
  lastML = sapply(cores,FUN = function(x)x$markerLayers[max(which(!is.na(x$markerLayers)))])
  firstML = sapply(cores,FUN = function(x)x$markerLayers[min(which(!is.na(x$markerLayers)))])

  #create matrix of where we have marker layers
  simML = matrix(FALSE,ncol = nCores,nrow= length(allMarkerLayers))
  for(a in 1:length(allMarkerLayers)){
    for(c in 1:nCores){
      if(length(which(allML[[c]]==allMarkerLayers[a])) > 0 ){
        simML[a,c]=TRUE
      }
    }
  }

  notEnough <- allMarkerLayers[rowSums(simML)<2]

  #now fill in FALSEs that are surrounded by TRUES
  for(c in 1:nCores){
    #which are false
    isF=which(!simML[,c])
    for(f in isF){
      if(f<max(which(simML[,c])) &  f>min(which(simML[,c]))){
        simML[f,c]=TRUE
      }
    }
  }

  nMLsim = cbind(allMarkerLayers,apply(cbind(rowSums(simML),rep(2,times = length(rowSums(simML)))),1,max))



  TPlist = lapply(MLsim,unique)
  allTP =unique(unlist(TPlist))
  allTP=allTP[which(!is.na(allTP))]
  inAll=c()
  inTwo=c()
  for (t in 1:length(allTP)){
    inAll[t]=all(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))
    if(any(allTP[t]==nMLsim[,1])){
    inTwo[t]=sum(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))==nMLsim[which(allTP[t]==nMLsim[,1]),2]}else{
      inTwo[t] <- FALSE
    }

  }

  #which Marker Layers are we simulating here
  ML2sim=as.vector(allTP[which(inTwo)])
  ML2sim=na.exclude(ML2sim)
  MLyears = matrix(NA,nrow = length(ML2sim),ncol = nCores)
  for(m in 1:nrow(MLyears)){
    allMLList = (sapply(cores,FUN = function(x)which(x$markerLayers==ML2sim[m])))
    for(c in 1:nCores){
      thisML = allMLList[[c]]
      #remove MLs that didn't pass
      if(length(thisML==1)){
        MLyears[m,c] = thisML
      }
    }
  }

  #make sure that the MLyears is in the correct order!!!
  ody = apply(MLyears,2,diff)
  if(any(ody[is.finite(ody)]<0)){
    stop("There are reversals in your marker layers. Check your input files.")
  }
  #print(MLyears)


  #WARNING: AS CURRENTLY CONSTRUCTED, THIS ONLY MAKES SENSE FOR SYMMETRIC OC UC PRIORS!
  #calculate diffs
  dML=diff(MLyears)

  #ASSUME A NORMAL DISTRIBUTION FOR THE DURATION of the YEARS and DRAW FROM IT
  #THIS MAY NOT BE IDEAL!!!!!!!!!
  diffMeans = rowMeans(diff(MLyears),na.rm=TRUE)#diffMeans = diff(MLyears) #diffMeans = rowMeans(diff(MLyears),na.rm=TRUE)
  rowSds(diff(MLyears),na.rm=TRUE)#diffSds = diff(MLyears) #rowSds(diff(MLyears),na.rm=TRUE)
  #simTrueDiff = rnorm(n = length(dML),mean = diffMeans,sd = diffSds)
  simTrueDiff = diffMeans
  OCUC = dML/simTrueDiff
  #rowSds(diff(MLyears))
  OCadd = OCUC-1
  OCadd[OCadd<=0]=0
  UCadd = -1*(OCUC-1)
  UCadd[UCadd<=0]=0


  for (c in 1:nCores){
    for (mi in 1:(nrow(MLyears)-1)){
      mlstart=which(ML2sim[mi]==cores[[c]]$markerLayers )
      mlstop=which(ML2sim[mi+1]==cores[[c]]$markerLayers )
      if(length(mlstart) ==1 & length(mlstop)==1){
        #find index over which to update priors
        pri=which(ML2sim[mi]==cores[[c]]$markerLayers ):which(ML2sim[mi+1]==cores[[c]]$markerLayers )
        cores[[c]]$bigOCA[pri,1]=OCadd[mi,c]
        cores[[c]]$bigUCA[pri,1]=UCadd[mi,c]
      }
    }
  }




  #use simulations on marker layers to get distributions on yearly priors.
  for (c in 1:nCores){
    #this way makes single estimates based on the mean. A distribuiton is better.
    # cores[[c]]$ocPrior=rowSums(cbind(cores[[c]]$ocPrior,rowMeans(cores[[c]]$bigOCA,na.rm=TRUE)),na.rm=TRUE)
    # cores[[c]]$ucPrior=rowSums(cbind(cores[[c]]$ucPrior,rowMeans(cores[[c]]$bigUCA,na.rm=TRUE)),na.rm=TRUE)

    #create a matrix of priors
    cores[[c]]$ocPriorMat = cores[[c]]$ocPrior+cores[[c]]$bigOCA
    cores[[c]]$ucPriorMat = cores[[c]]$ucPrior+cores[[c]]$bigUCA



  }

  #empty out this year and last year
  lastMLyear = 0
  thisMLyear = 0

  #make sure that there are no reversals!!!
  #MLyears
  ody = apply(MLyears,2,diff)

  # #store Marker Layer posterior data
  # for(c in 1:nCores){
  #   for(ml in 1:length(ML2sim)){
  #     wmlp=which(ML2sim[ml]==cores[[c]]$markerLayers)
  #     corePost[[c]]$mlPost[wmlp,i]=ML2sim[ml]
  #   }
  # }

  #create a vector that gives instructions about matching seg lengths
  mustMatchLengths = matrix(TRUE,nrow=(nrow(MLyears)-1))

  #add to ML years if need be.
  #does the first segment start at 1? (is the surface a marker layer?)
  ML1test = MLyears[1,!is.na(MLyears[1,])]==1
  if (!all(ML1test) & any(ML1test)){
    stop("You cant have a marker layer that starts at 1 for some cores and not all of them")
  }
  if (!all(ML1test) ){
    #if no top seg, add one
    MLyears=rbind(matrix(1,ncol=nCores),MLyears)
    mustMatchLengths=rbind(TRUE,mustMatchLengths)#ALERT!!!!! THIS MEANS THAT WE'VE GIVEN A prob of 1 to the top layer!
  }

  lastMLtest = MLyears[nrow(MLyears),!is.na(MLyears[nrow(MLyears),])]
  totYearTest = totYears[!is.na(MLyears[nrow(MLyears),])]
  if(any((totYearTest - lastMLtest)==0) & !all((totYearTest - lastMLtest)==0)){
    stop("Your last marker layer must either be the oldest years, or at least 2 years a way from your oldest years")
  }
  if(!any((totYearTest - lastMLtest)==0)){
    #add in the bottom MLyears layer
    MLyears=rbind(MLyears,totYears)
    ML2sim=append(ML2sim,"bottom")

    mustMatchLengths=rbind(mustMatchLengths,FALSE)

  }

  #set up segments
  seg=vector("list",nrow(MLyears)-1)


  for(s in 1:(nrow(MLyears)-1)){
    for(c in 1:nCores){
      if(s==(nrow(MLyears)-1) | !is.na(MLyears[s,c]) & !is.na(MLyears[s+1,c])){
        if(s==1){
          seg[[s]][[c]]=cores[[c]][(MLyears[s,c]):MLyears[s+1,c],]#include the first marker layer the first year

        }else if(s< (nrow(MLyears)-1)){
          seg[[s]][[c]]=cores[[c]][(MLyears[s,c]+1):MLyears[s+1,c],]
        }else{#IF ITS THE LAST ROW, make the segment from the last, non-NA marker layer
          lastRealMarker = max(MLyears[-nrow(MLyears),c],na.rm=T)
          seg[[s]][[c]]=cores[[c]][(lastRealMarker+1):MLyears[s+1,c],]
        }
      }
    }
  }
  segTotalThicks = matrix(NA,nrow = 5,ncol=length(seg))
  segTotalYears = matrix(NA,nrow = 5,ncol=length(seg))

  ###CHECK THAT SEGMENT THICKNESS ARE CORRECT
  for(c in 1:nCores){
    for(s in 1:(nrow(MLyears)-1)){
      if(length(seg[[s]])>=c){
        segTotalThicks[c,s] = sum(seg[[s]][[c]]$thickness)
        segTotalYears[c,s] = length(seg[[s]][[c]]$year)
      }
    }
  }


  #   if(!all(sapply(seg[[1]],FUN = function(x)min(x$year)) == 1)){
  #     #make a new marker layer that goes to the top, that doesn't have the same length requirement.
  #     topSeg[[c]] = cores[[c]][MLyears[s,c]:MLyears[s+1,c],]
  #   }
  #segment

  newSeg=vector("list",nrow(MLyears)-1)
  #simulate between segments
  tSim=vector("list",nCores)

  composite = c()
  for(s in 1:length(seg)){#for each segment...
    #which cores are we simulating?
    core2sim=which(!sapply(seg[[s]],is.null))
    orig.core2sim = core2sim

    sameLength = FALSE
    t=0
    while(!sameLength){#
      t=t+1

      #build in failsafe
      if(t>1e4){
        if(length(core2sim)>2){
          print(paste("segment",as.character(s),"hit count limit.", as.character(length(core2sim)),"cores. Reducing number of cores in simulation"))
          core2sim = core2sim[-sample(core2sim)[1]]#remove one of the cores.
          t=0
        }else{
          stop("Too many iterations..no more cores to remove.. stopped.")
        }
      }
      #DO YAHTZEEE HERE
      for(c in core2sim){
        #simulate new segments given the UC and OC priors.
        soacDf <- data.frame(thick = seg[[s]][[c]]$thickness, ucProb = seg[[s]][[c]]$ucPriorMat, ocProb = seg[[s]][[c]]$ocPriorMat)
        tSim[[c]]= simulateOverAndUndercounting(soacDf)
      }

      newSegLengths = sapply(tSim,FUN = function(x)length(x$newThicks))[core2sim]
      sameLength= (all(newSegLengths==newSegLengths[1]) | !mustMatchLengths[s])


    }

    #see if any didn't run
    didnt.sim = orig.core2sim[!orig.core2sim %in% core2sim]
    if(length(didnt.sim)>0){#if any didn;t run...
      for(dd in didnt.sim){#one by one, simulate the length for any that got excluded

        ddt=0
        while(TRUE){#and simulate until you match the length
          ddt=ddt+1
          soacDf <- data.frame(thick = seg[[s]][[dd]]$thickness, ucProb = seg[[s]][[dd]]$ucPriorMat, ocProb = seg[[s]][[dd]]$ocPriorMat)
          tSim[[dd]]= simulateOverAndUndercounting(soacDf)
          #check to see if it's the correct length
          thisSegLength = length(tSim[[dd]]$newThicks)
          if(thisSegLength == newSegLengths[1]){#found it!!!
            core2sim = append(core2sim,dd)
            newSegLengths = append(newSegLengths,thisSegLength)
            break
          }
          if(ddt>1e5){#give up
            print(paste("gave up trying to get core",as.character(dd),"to be",as.character(newSegLengths[1]),"years long for segment",as.character(s)))
            break
          }
        }
      }
    }





    #Once same length is matched build new data.frame
    #And assign results into posterior.
    segthickmat=matrix(NA,nrow=max(newSegLengths),ncol=length(core2sim))
    for(c in 1:length(core2sim)){
      #build Thickmat
      segthickmat[1:newSegLengths[c],c] = tSim[[core2sim[c]]]$newThicks
    }
    # cm=cor(segthickmat,use="complete.obs")
    # cm[which(cm<0)]=0
    # meanCorr=mean(cm[lower.tri(cm)])

    ###ASSIGN IN POSTERIOR INFORMATION###
    for(c in core2sim){
      # #assign posterior info
      # #Undercounted
      # ucIndex = tSim[[c]]$UCi+min(seg[[s]][[c]]$year)-1
      # corePost[[c]]$ucPost[ucIndex,i]=1
      # #Overcounted
      # ocIndex = tSim[[c]]$OCi+min(seg[[s]][[c]]$year)-1
      # corePost[[c]]$ocPost[ocIndex,i]=1
      #
      # #mean correlation
      # corePost[[c]]$meanCorr[min(seg[[s]][[c]]$year) : max(seg[[s]][[c]]$year) ,i] = meanCorr

      #store thicknesses for each core
      if(mustMatchLengths[s]){
        if(c==core2sim[1]){#only update if its the first one.
          thisMLyear = lastMLyear+newSegLengths[1]
        }
        #          corePost[[c]]$mlCum[thisMLyear,i]=ML2sim[s+1]
        #test if sum thicknesses are the same as the original (THEY SHOULD BE!)
        if(abs(sum(tSim[[c]]$newThicks)-sum(seg[[s]][[c]]$thickness))>.001){
          stop("The thicknesses don't sum to the same number!")
        }
       # corePost[[c]]$thick2=append(corePost[[c]]$thick2,tSim[[c]]$newThicks)
        corePost[[c]]$thicks[(lastMLyear+1) : thisMLyear ,i] = tSim[[c]]$newThicks


      }else{
        lastSimYear = max(which(!is.na(corePost[[c]]$thicks[,i])))
        corePost[[c]]$thicks[(lastSimYear+1) : (lastSimYear+newSegLengths[which(c==core2sim)]) ,i] = tSim[[c]]$newThicks
      }
    }


    #is the total thickness the same as the total thickness should be?
    curThick=sum(corePost[[c]]$thicks[1 : sum(!is.na(corePost[[c]]$thicks[,i])) ,i] )
    lowerMLname= seg[[s]][[c]]$markerLayers[length(seg[[s]][[c]]$markerLayers)]

    #assign age into MLpostAge
    if(!is.na(lowerMLname)){
      wML = which(lowerMLname==allMarkerLayers)
      MLpostAge[wML,i] = thisMLyear
    }

    origMLi = which(cores[[c]]$markerLayers==lowerMLname)
    if(length(origMLi)<1){
      origMLi = length(cores[[c]]$thickness)
    }
    origThick = sum(cores[[c]]$thickness[1 : origMLi])
    # if(abs(origThick-curThick)> .001){
    #   stop("The current cumulative thickness isn't right")
    # }


    lastMLyear=thisMLyear#update with new last year

    #build Composite
    segCom= rowMeans(scale(segthickmat,center=meanThick[core2sim],scale =sdThick[core2sim]),na.rm = TRUE)
    composite = append(composite,segCom)
  }


#return the output
out <- list(composite,corePost,MLpostAge)

return(out)


}
