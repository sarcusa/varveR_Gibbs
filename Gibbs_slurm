load("~/Varves/Gibbs_sampler/Gibbs_observer1.RData")

library("dplyr")
library(matrixStats)
library(purrr)
library(tictoc)


for(i in 2:nIt){ #choose here if running loop in segments
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
    newProb <- try(exp(logObjFun(param = c(paramChain[i,1],paramChain[i,2],paramChain[i,3],
                                           paramChain[i,4],paramChain[i,5],paramChain[i,6]),C14,col172_sequence,
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
    
    save.image("~/Varves/Gibbs_sampler/Gibbs_observer1.RData")
    
  }
  
}

save.image("~/Varves/Gibbs_sampler/Gibbs_observer1.RData")
