######
# Original gibbs model algorithm for one observer. Best to run this not too many iterations.
#
# S. Arcusa 
# 17-04-2020

lpdfile = "" # LipD file for COlumbine Lake. Available at https://doi.org/10.6084/m9.figshare.14417999
Obs3dir_172 = "" # Location folder containing data from observer 3 core COL17-2. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
Obs3dir_173 = "" # Location folder containing data from observer 3 core COL17-3. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
outputdir = "" # Folder location for output figures and output workstation

library("devtools")
library("varveR")
library("ggplot2")
library(geoChronR)
library("dplyr")
library(plyr)
library("lipdR")
library(matrixStats)
library(purrr)
library(tictoc)


# Load radiocarbon and lead data contained in the LiPD file
D <- readLipd(lpdfile)

# Run the BACON model
D <- runBacon(D,cutoff = 1e-10, lab.id.var = NULL, age.14c.var = "age14c", 
              age.14c.uncertainty.var = "age14cuncertainty", 
              age.uncertainty.var = "ageuncertainty", 
              reservoir.age.14c.var = NULL, rejected.ages.var = NULL,
              reservoir.age.14c.uncertainty.var = NULL)

# Map the age ensemble to the data
core = mapAgeEnsembleToPaleoData(D, paleo.meas.table.num = 1, age.var = "ageensemble") 

# Retrieving the age density distributions
# Prepare the 210Pb data
PB_distribution <- core$chronData[[1]]$model[[1]]$distributionTable
PB_probdens <- purrr::flatten(lapply(PB_distribution, function (x) x[c('probabilityDensity')]))
PB_dens <- purrr::flatten(lapply(PB_probdens, function (x) x[c('values')]))
PB_age <- purrr::flatten(lapply(PB_distribution, function (x) x[c('age')]))
PB_age_values <- purrr::flatten(lapply(PB_age, function (x) x[c('values')]))
# Combine 210Pb with 14C depths
depths_14C <- as.numeric(unlist(purrr::flatten(lapply(PB_distribution, function (x) x[c('depth')]))))
names(PB_dens) <- paste("depth_",depths_14C, sep = "")
names(PB_age_values) <- paste("depth_",depths_14C, sep = "")
# Prepare the 237Cs data
Cs_age <- seq(-15,-10,length.out = 250)
Cs_dens <-dnorm(Cs_age,-13,.25)
Cs_depth <- 3.25*10
Cs <- list(list(prop.table(Cs_dens),Cs_age,Cs_depth))
# Subset the 210Pb data 
PB_dens <- PB_dens[c(13:19)]
PB_age_values <- PB_age_values[c(13:19)]
depths_14C <- depths_14C[c(13:19)]
# Prepare object containing 210Pb,237Cs,14Cs depths, densities and ages
C14 <- list()
l <- list()
for(i in 1: length(depths_14C)){
  
  l[[1]] <- PB_dens[[i]]
  l[[2]] <- PB_age_values[[i]]
  l[[3]] <- depths_14C[i]*10
  
  C14[[i]] <- l
}

C14 <- c(Cs,C14)

#Plot age-depth model
start = -70
end = 4000
chron.plot = plotChron(core, model.num = 1, age.var = "ageensemble",
                       x.bin = seq(start,end, by= 100), 
                       dist.color = "#1E88E5",
                       probs = c(.1,.5,.9),y.bin = seq(0,300, by = 1))+
  labs(title = "", x = "Age (yr BP)", y = "")+
  xlim(end,start)+ 
  theme_linedraw()+
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
plot(chron.plot)


## Load varve data for each core

col172 <- readVarveDirectory(Obs3dir_172,varveTop = "left")
col173 <- readVarveDirectory(Obs3dir_173,varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

col172_sequence <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence <- combineSectionsByMarkerLayer(col173[o3])

# the following will calibrate the cumulative depth to the real depth at the last slide max depth
test3 <- vector() 
OldMin = 0
NewMax = 1230
NewMin = 0

for(i in 1:nrow(col173_sequence)){
  
  OldMax = max(cumsum(col173_sequence$thick), na.rm = T)
  OldValue = col173_sequence$thick[i]
  
  OldRange = (OldMax - OldMin)  
  NewRange = (NewMax - NewMin)  
  NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
  
  test3[i] <- NewValue
  
}
col173_sequence$thick <- test3

test2 <- vector() 
OldMin = 0
NewMax = 1270
NewMin = 0

for(i in 1:nrow(col172_sequence)){
  
  OldMax = max(cumsum(col172_sequence$thick), na.rm = T)
  OldValue = col172_sequence$thick[i]
  
  OldRange = (OldMax - OldMin)  
  NewRange = (NewMax - NewMin)  
  NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
  
  test2[i] <- NewValue
  
}
col172_sequence$thick <- test2


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
                           overcount = c(ov1,ov2,ov3,0,-2,-3), 
                           undercount =  c(uc1,uc2,uc3,.50,-2,-3)) 
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

tic("myobjF")
objectiveFunction(0,0,0,0,0,0,C14 = C14,col172_sequence, col173_sequence)
objectiveFunction(0.01,0.01,0,0.1,0,0.25,C14 = C14,col172_sequence, col173_sequence)
objectiveFunction(0.05,0.05,0,0.2,0,0.6,C14 = C14,col172_sequence, col173_sequence)
toc()

logObjFun <-  function(param,C14, core1, core2){
  
  ov1 <- param[1]
  ov2 <- param[3]
  ov3 <- param[5]
  uc1 <- param[2]
  uc2 <- param[4]
  uc3 <- param[6]
  
  return(log(objectiveFunction(ov1,uc1,ov2,uc2,ov3,uc3,C14,core1,core2))) 
}

tic("mylogF")
logObjFun(param = c(0,0,0,0,0,0),C14,col172_sequence, col173_sequence)
logObjFun(param = c(0.01,0.01,0,0.1,0,0.25),C14,col172_sequence, col173_sequence)
logObjFun(param = c(0.05,0.05,0,0.2,0,0.6),C14,col172_sequence, col173_sequence)
toc()

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

tic("ourGibbs")
nIt <- 1e6
paramChain = matrix(NA,nrow = nIt,6)
objChain = matrix(NA,nrow = nIt,6)
paramChain[1,] = c(0.05,0.05,0.1,0.2,0.3,0.6)

objChain[1,1] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3],
                    paramChain[1,4],paramChain[1,5],paramChain[1,6]),
                    C14,col172_sequence,col173_sequence)
#objChain[1,2] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3],
#                    paramChain[1,4],paramChain[1,5],paramChain[1,6]),
#                    C14,col172_sequence, col173_sequence)
#objChain[1,3] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3],
#                    paramChain[1,4],paramChain[1,5],paramChain[1,6]),
#                    C14,col172_sequence, col173_sequence)
#objChain[1,4] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3],
#                    paramChain[1,4],paramChain[1,5],paramChain[1,6]),
#                    C14,col172_sequence, col173_sequence)
#objChain[1,5] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3],
#                    paramChain[1,4],paramChain[1,5],paramChain[1,6]),
#                    C14,col172_sequence, col173_sequence)
#objChain[1,6] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3],
#                    paramChain[1,4],paramChain[1,5],paramChain[1,6]),
 #                  C14,col172_sequence, col173_sequence)

oldProb <- exp(objChain[1,1]) # changed from objChain[1,2]
objChain[1,] <- objChain[1,1]



tic("chain")
# loop through each step in the chain
for(i in 2:2){ #choose here if running loop in segments
#for(i in 2:nrow(paramChain)){ #choose here if running loop in one go
  print(i)
  #each time through initialize with the last time
  paramChain[i,]=paramChain[i-1,] 
  
  # loop through each parameter
  for(j in seq(1,ncol(paramChain),2)){ 
    # This limits the parameters to between 0 and 1
    repeat{
      a = abs(rnorm(1,sd = 0.1)) # use 0.1 to start, then 0.05
      b = abs(rnorm(1,sd = 0.1))
      if(runif(1) > 0.5){
        a = a*-1
      }else{
        b = b*-1
      }
      test1 = paramChain[i-1,j]+a 
      test2 = paramChain[i-1,j+1]+b 
      if(test1 >= 0 & test1 <=1 & test2 >= 0 & test2 <=1) break 
      
    }
   
    #innovate on the parameter by adjusting by a random normal number
    paramChain[i,j]=test1
    paramChain[i,j+1]=test2 # this was added
    
    #this tries the varve model and if an error arrises, it repeats the previous row
    # Calculates the log objective function
    #newProb <- try(logObjFun(param = c(paramChain[i,1],paramChain[i,2],paramChain[i,3],
    #                    paramChain[i,4],paramChain[i,5],paramChain[i,6]),
    #                    C14,col172_sequence, col173_sequence)) #removed on 22-06-20
    
    # Takes into account the priors through the function pv
    newProb <- try(exp(logObjFun(param = c(paramChain[i,1],paramChain[i,2],
                                           paramChain[i,3],paramChain[i,4],
                                           paramChain[i,5],paramChain[i,6]),
                                 C14,col172_sequence,
                             col173_sequence))*pV(paramChain[i,]))
    
    
    if(!class(newProb)=="try-error"){
      
      #trying to get the number bigger
      if( newProb/oldProb >= runif(1)){
        #if it passes, update old prob  
        oldProb <- newProb       
      }else{ 
        #if it fails, keep it the same
        paramChain[i,j] <- paramChain[i-1,j]
        paramChain[i,j+1] <- paramChain[i-1,j+1]
      } 
      #update the objective tracker with the updated, or previous prob
      objChain[i,j] <- log(oldProb)
      
      #if there is an error, repeats
    } else {
      print("got error")
     
        paramChain[i,j] <- paramChain[i-1,j]
        paramChain[i,j+1] <- paramChain[i-1,j+1]
    
      objChain[i,j] <- log(oldProb)
      
      }
      
   #print(oldProb)
    }
  
  p <- seq(0,nIt,10)
  if(i %in% p){
    
    print(objChain[(i-5):i,]) 
    print(paramChain[(i-5):i,])
    
  }
  
  q <- seq(0,nIt,1000)
  if(i %in% q){
    
    save.image(paste0(outputdir,"gibbs_sampling_workstation.RData"))
    
  }
  
}

save.image(paste0(outputdir,"gibbs_sampling_workstation.RData"))

toc()

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

comp.plot <- plot_varve_model(param = paramChain, core1 = col172_sequence, core2 = col173_sequence, N = 30000, ind.plot = chron.plot, somedepths = depths_14C, it = 10)


## Plotting results of Gibbs

plot_Gibbs_results  <- function(obj,param,n, plot.save){
  
  # Plots the objective chain
  if(plot.save == TRUE){
    
    pdf(file = paste0(outputdir,"ObjChain_plot.pdf"))
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
    
    pdf(file = paste0(outputdir,"ParamChain_plot.pdf"))
    par(mfrow = c(3,2))
    for(i in 1:6){
      plot(param[,i],type = "l",xlim = c(0,n), 
           xlab = "# iterations", ylab = names_param[i])
    }
    dev.off()
    
    tki <- seq(b,n,by = 1)
    
    pdf(file = paste0(outputdir,"HistParam_plot.pdf"))
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

plot_Gibbs_results(obj = objChain, param = paramChain, n = 40000, plot.save = F)

min(which(is.na(objChain[,1])))
