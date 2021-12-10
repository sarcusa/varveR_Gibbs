## Objective function

# Arguments
#
# ov1..n : overcounting probabilities for varve ID 1..n
# uc1..n : undercounting probabilities for varve ID 1..n
# C14 : a nested list with each sublist representing a radiocarbon depth and containing [1] the probability densities of each year, [2] the age value for each probability, and [3] the actual sample depth
# core 1 and core 2 should be col172_sequence and col173_sequence

objectiveFunction = function(ov1,uc1,ov2,uc2,ov3,uc3,C14,core1,core2){
  
  # creates a table of codes
  tt <- data.frame(codes = 1:6,
                   overcount = c(ov1,ov2,ov3,0.5,-2,-3), 
                   undercount =  c(uc1,uc2,uc3,0,-2,-3)) 
  # converts the codes to probabilities
  seq1 <- translateCodesToProbabilities(core1,translationTable = tt,baselineProb = .05)
  seq2 <- translateCodesToProbabilities(core2,translationTable = tt,baselineProb = .05)
  # finds tie points and calculates varve thicknesses in each sequence
  sseq1 <- simulateOverAndUndercounting(seq1)
  g <- which(!is.na(sseq1$newThicks))
  sseq1$ensThick <- sseq1$newThicks[g]
  sseq1$tiePoints <- sseq1$newTiePoint[g]
  
  sseq2 <- simulateOverAndUndercounting(seq2)
  g <- which(!is.na(sseq2$newThicks))
  sseq2$ensThick <- sseq2$newThicks[g]
  sseq2$tiePoints <- sseq2$newTiePoint[g]
  
  #Keep the tie points as characters
  sseq1$newTiePoint = as.character(sseq1$newTiePoint)
  sseq2$newTiePoint = as.character(sseq2$newTiePoint)
  
  # creates a list of cores
  ensList <- list(sseq1,sseq2)
  #finds the common tiepoints
  allTiePoints <- c("topOfCore", 
                    na.omit(unique(c(unique(sseq1$newTiePoint)),
                                   na.omit(unique(sseq2$newTiePoint)))))
  #runs the model
  MD <- fastVarveModel(ensList, allMarkerLayers = allTiePoints) 
  # Add cumulatitve depth in real unit
  MD[[4]] <- cumsum(abs(MD[[2]][[2]]$thicks[which(!is.na(MD[[2]][[2]]$thicks))])) 
  # Add year
  MD[[5]] <- seq(1,length(MD[[4]]),1)
  # Rename each item in the varve model output 
  names(MD) <- c("ModelOutput", "PostInfo","MLpostAge","CumDepth", "Age")
  
  y = vector() # vector of varve age at each depth
  t = vector() # vector of positions of closest year in the each 14C sample age distribution
  d = vector() # vector of age probability density for each 14C sample
  
  for(j in 1:length(C14)){
    
    # This finds the varve age at the same depth as the radiocarbon/lead sample. 
    # Note: the age scales must be the same so added calculation to do that
    y[j] = 1950-(2018- MD$Age[which(abs(MD$CumDepth-as.numeric(C14[[j]][3])) == 
                                      min(abs(MD$CumDepth-as.numeric(C14[[j]][3])), 
                                          na.rm = T))]) 
    # Now looks for the closest year in the pdf of 14C 
    t[j] <- which(abs(unlist(C14[[j]][2])-y[j]) == min(abs(unlist(C14[[j]][2])-y[j]), 
                                                       na.rm = T))
    # and find its density
    d[j] = unlist(C14[[j]][1])[which(abs(unlist(C14[[j]][2])-y[j]) ==
                                       min(abs(unlist(C14[[j]][2])-y[j]),
                                           na.rm = T))]
  }
  #print(d)
  prob = prod(d) # multiply the density of all radiocarbon samples at that age
  
  return(prob)
  
}


logObjFun <-  function(param,C14, core1, core2){
  
  ov1 <- param[1]
  ov2 <- param[3]
  ov3 <- param[5]
  uc1 <- param[2]
  uc2 <- param[4]
  uc3 <- param[6]
  
  return(log(objectiveFunction(ov1,uc1,ov2,uc2,ov3,uc3,C14,core1,core2))) 
}

# this function sets the priors for varve ID 1 and 2 as gamma curves, ID 3 as flat line priors
pV<- function(priors){
  
  pout <- vector()
  pout[1] <- dgamma(priors[1],shape = 2,rate = 10)/100
  pout[2] <- dgamma(priors[2],shape = 2,rate = 10)/100
  pout[3] <- dgamma(priors[3],shape = 10,rate = 20)/100
  pout[4] <- dgamma(priors[4],shape = 6,rate = 20)/100
  pout[5] <- dunif(priors[5])/100
  pout[6] <- dunif(priors[6]) /100
  
  return(prod(pout))
  
}

## Plotting varve model results to check

#' Plots a varve age depth model with uncertainty using the last 1000 sets of parameters from the Gibbs sampler
#'
#' @description The function will loop through the last 1000 sets of parameters from the Gibbs sampler using the 
#' fastVarveModel and plot the varve age depth model with uncertainty against the independent age depth model 
#' 
#' @param param the paramChain output from the Gibbs sampler
#' @param core1 The first core sequence output from the function combineSectionsByMarkerLayer
#' @param core2 The second core sequence output from the function combineSectionsByMarkerLayer
#' @param N The last possible iteration number of the Gibbs sampler. For example, if the sampler was run one 
#' million times, then N would be set to one million
#' @param ind.plot The independent age depth model to compare the varve model to. This should be the output of 
#' GeoChronR function plotChron
#' @param somedepths A vector of depths at which to calculate the age spread. This could be the depths of 
#' radiocarbon samples
#' @param it The number of iterations to run the fastVarveModel for
#'
#' @return A plot comparing the independent age depth model against the varve age depth model
#' @return A 
#' 

plot_varve_model_Gibbs <- function(param, core1, core2, N, ind.plot, somedepths, it){
  
  ens.thicks <- list()
  for(i in 1:it){
    print(i)
    tt <- data.frame(codes = 1:6,
                     overcount = c(param[N-i,1],param[N-i,3],param[N-i,5],0,-2,-3), 
                     undercount =  c(param[N-i,2],param[N-i,4],param[N-i,6],.50,-2,-3)) 
    
    seq1 <- translateCodesToProbabilities(core1,translationTable = tt,baselineProb = .05)
    seq2 <- translateCodesToProbabilities(core2,translationTable = tt,baselineProb = .05)
    
    sseq1 <- simulateOverAndUndercounting(seq1)
    g <- which(!is.na(sseq1$newThicks))
    sseq1$ensThick <- sseq1$newThicks[g]
    sseq1$tiePoints <- sseq1$newTiePoint[g]
    
    sseq2 <- simulateOverAndUndercounting(seq2)
    g <- which(!is.na(sseq2$newThicks))
    sseq2$ensThick <- sseq2$newThicks[g]
    sseq2$tiePoints <- sseq2$newTiePoint[g]
    
    ensList <- list(sseq1,sseq2)
    
    allTiePoints <- c("topOfCore", 
                      na.omit(unique(c(unique(sseq1$newTiePoint)),
                                     na.omit(unique(sseq2$newTiePoint)))))
    
    MD <- fastVarveModel(ensList, allMarkerLayers = allTiePoints) 
    
    ens.thicks[[i]] <- MD[[2]][[2]]$thicks
    
  } 
  
  ens.cumSum <- lapply(ens.thicks,  function(x) cumsum(abs(x[which(!is.na(x))])))
  ens.age <- lapply(ens.cumSum, function(x) seq(1,length(x),1))
  indx <- sapply(ens.thicks, length)
  res <- as.data.frame(do.call(cbind,lapply(ens.thicks, `length<-`,max(indx))))
  
  ageEns <- createEnsembleAgeDepthModel(res)
  ageModel <- createDepth2AgeFunction(ageDepthMat = ageEns$ageDepthEns)
  varveAges <- ageModel(somedepths*10) # somedepths must be in mm
  
  compare.models <- ind.plot +
    geom_ribbon(aes(y = ageEns$summaryTable[,1]/10,xmin = -67+ageEns$summaryTable[,2], 
                    xmax = -67+ageEns$summaryTable[,6]),fill = "gray80")+
    geom_line(aes(x = -67+ageEns$summaryTable[,4], y = ageEns$summaryTable[,1]/10), color = "red")
  
  print(compare.models)
  
  return(varveAges)
  return(ageEns)
  
  
}

plot_Gibbs_results  <- function(obj,param,n, plot.save){
  
  # Plots the objective chain
  if(plot.save == TRUE){
    
    pdf(file = "~/Varves/Gibbs_sampler/ObjChain_plot.pdf")
    plot((obj[,1]),type = "l", xlim = c(0,n), 
         xlab = "# iterations", ylab = "Objective Function (log)")
    dev.off()
    
    b  <- readline("What burn in value to use?")
    b  <- as.numeric(unlist(strsplit(b, ",")))
    
    names_param <- c("Overcounting prob.varve ID1", 
                     "Undercounting prob. varve ID1",
                     "Overcounting prob. varve ID2",
                     "Undercounting prob. varve ID2", 
                     "Overcounting prob. varve ID3",
                     "Undercounting prob. varve ID3")
    
    pdf(file = "~/Varves/Gibbs_sampler/ParamChain_plot.pdf")
    par(mfrow = c(3,2))
    for(i in 1:6){
      plot(param[,i],type = "l",xlim = c(0,n), 
           xlab = "# iterations", ylab = names_param[i])
    }
    dev.off()
    
    tki <- seq(b,n,by = 1)
    
    pdf(file = "~/Varves/Gibbs_sampler/HistParam_plot.pdf")
    par(mfrow = c(3,2))
    for(i in 1:6){
      hist(param[tki,i], xlab = names_param[i], main = "")
    }
    dev.off()
  }else{
    
    plot((obj[,1]),type = "l", xlim = c(0,n), 
         xlab = "# iterations", ylab = "Objective Function (log)")
    
    b  <- readline("What burn in value to use?")
    b  <- as.numeric(unlist(strsplit(b, ",")))
    
    names_param <- c("Overcounting prob.varve ID1", 
                     "Undercounting prob. varve ID1",
                     "Overcounting prob. varve ID2",
                     "Undercounting prob. varve ID2", 
                     "Overcounting prob. varve ID3",
                     "Undercounting prob. varve ID3")
    par(mfrow = c(3,2))
    for(i in 1:6){
      plot(param[,i],type = "l",xlim = c(0,n), 
           xlab = "# iterations", ylab = names_param[i])
    }
    
    tki <- seq(b,n,by = 1)
    
    dev.new()
    par(mfrow = c(3,2))
    for(i in 1:6){
      hist(param[tki,i], xlab = names_param[i], main = "")
    }
    
  }
}
