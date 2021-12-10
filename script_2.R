# Arcusa et al Columbine Lake varve
# Script 2
# To be run after Gibbs slurm
# Once Gibbs has run 50,000 times, resume here for individual observer results. Will combine observers and cores.

Rsourcedir = "" # Location of the source code files. These include file GibbsrRelatedFunctions.R, ImportFiles.R, plotting.R, simulateOverAndUnderCounting.R, varveModel.R
Gibbsresultsdir = "" # Folder location for the Gibbs results
resultsdir = "" # Folder location for the results

lpdfile = "" # Change to point to LipD file. LipD file contains the data in a reuseable format. Available at https://doi.org/10.6084/m9.figshare.14417999.v1 
figdir = "" # Folder location to save figures of age depth models
CRSfile = "" # Point to CRS lead model file (txt), datasets available at https://doi.org/10.6084/m9.figshare.17156702
CFCSfile = "" # Point to CFCS lead model file (txt),datasets available at https://doi.org/10.6084/m9.figshare.17156702
CICfile = "" # Point to CIC lead model file (txt), datasets available at https://doi.org/10.6084/m9.figshare.17156702
Csfile = "" # Point to raw cesium data (csv) available at https://doi.org/10.6084/m9.figshare.17157245

# For observer 1
Obs1dir_172 = "" # Location folder containing data from observer 1 core COL17-2. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
Obs1dir_173 = "" # Location folder containing data from observer 1 core COL17-3. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
output_obs1_172 = "" # PDF output name for observer 1 core COL17-2. Must end with .pdf
output_obs1_173 = "" # PDF output name for observer 1 core COL17-3. Must end with .pdf

# For observer 2
Obs2dir_172 = "" # Location folder containing data from observer 2 core COL17-2. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
Obs2dir_173 = "" # Location folder containing data from observer 2 core COL17-3. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
output_obs2_172 = "" # PDF output name for observer 2 core COL17-2. Must end with .pdf
output_obs2_173 = "" # PDF output name for observer 2 core COL17-3. Must end with .pdf

# For observer 3
Obs3dir_172 = "" # Location folder containing data from observer 3 core COL17-2. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
Obs3dir_173 = "" # Location folder containing data from observer 3 core COL17-3. Data available at https://doi.org/10.6084/m9.figshare.14251400.v1
output_obs3_172 = "" # PDF output name for observer 3 core COL17-2. Must end with .pdf
output_obs3_173 = "" # PDF output name for observer 3 core COL17-3. Must end with .pdf

library("devtools")
#install_github("nickmckay/varveR")
#library("varveR")
library("ggplot2")
library(geoChronR)
library("dplyr")
library("matrixStats")
library("lipdR")
library(rplum)
library(cowplot)
library(LaplacesDemon)

# Observer 1

loadOneName(objName = "Bacon.plot", file = paste0(resultsdir,"multiple_cores_observers_workstation_1.RData"))

load(paste0(resultsdir,"Gibbs_observer1.RData"))

plot_Gibbs_results(obj = objChain, param = paramChain, n = 25000, plot.save = F, obs = 1, dir = Gibbsresultsdir)

tic()
ages1 <- plot_varve_model_Gibbs(param = paramChain, obj = objChain, core1 = col172_sequence, core2 = col173_sequence, ind.plot = Bacon.plot, somedepths = depths_14C, it = 300,dir = Gibbsresultsdir,obs = 1,PlotOpt = F,burnOpt = F,burn = 100)
toc()

loadOneName(objName = "col173_sequence1", file = paste0(resultsdir,"multiple_cores_observers_workstation_1.RData"))

count <- seq_len(nrow(ages1$ThicksEns[[2]]))
I1 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(ages1$ThicksEns[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(ages1$ThicksEns[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

#I2O1 <- sedRatesFun(depth = I2O1.depth, age = I1.age,tocumulate = T)
I1.age2 <- as.data.frame(ages1$ageEns$summaryTable)

#Observer 2

load(paste0(resultsdir,"Gibbs_observer2.RData"))

plot_Gibbs_results(obj = objChain, param = paramChain, n = 50000, plot.save = T, obs = 2, dir = Gibbsresultsdir)

tic()
ages2 <- plot_varve_model_Gibbs(param = paramChain, obj = objChain, core1 = col172_sequence, core2 = col173_sequence, ind.plot = Bacon.plot, somedepths = depths_14C, it = 300,dir = Gibbsresultsdir,obs = 2,PlotOpt = T,burnOpt = T)
toc()

loadOneName(objName = "col173_sequence2", file = paste0(resultsdir,"multiple_cores_observers_workstation_1.RData"))

count <- seq_len(nrow(ages2$ThicksEns[[2]]))
I2 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(ages2$ThicksEns[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(ages2$ThicksEns[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)


#I2O3 <- sedRatesFun(depth = I2O2.depth, age = I2.age,tocumulate = T)
I2.age2 <- as.data.frame(ages2$ageEns$summaryTable)

# Observer 3

load(paste0(resultsdir,"Gibbs_observer3.RData"))

plot_Gibbs_results(obj = objChain, param = paramChain, n = 50000, plot.save = T, obs = 3, dir = Gibbsresultsdir)

tic()
ages3 <- plot_varve_model_Gibbs(param = paramChain, obj = objChain, core1 = col172_sequence, core2 = col173_sequence, ind.plot = Bacon.plot, somedepths = depths_14C, it = 300,dir = Gibbsresultsdir,obs = 3,PlotOpt = T,burnOpt = F,burn = 50)
toc()

loadOneName(objName = "col173_sequence3", file = paste0(resultsdir,"multiple_cores_observers_workstation_1.RData"))

count <- seq_len(nrow(ages3$ThicksEns[[2]]))
I3 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(ages3$ThicksEns[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(ages3$ThicksEns[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

#I2O3 <- sedRatesFun(depth = I2O3.depth, age = I3.age,tocumulate = T)
I3.age2 <- as.data.frame(ages3$ageEns$summaryTable)

## SCRIPT NOW INTEGRATES ALL MODELS

library("devtools")
library("ggplot2")
library(geoChronR)
library("dplyr")
library("matrixStats")
library("lipdR")
library(rplum)
library(cowplot)
library(LaplacesDemon)
library(purrr)
library(HistogramTools)

# Observer 1

load(paste0(resultsdir,"Gibbs_observer1.RData"))

obs1 <- prepareObserverGibbs(param = paramChain, obj = objChain, core1 = col172_sequence, core2 = col173_sequence,it = 300,burnOpt = T)

# Observer 2

load(paste0(resultsdir,"Gibbs_observer2.RData"))

obs2 <- prepareObserverGibbs(param = paramChain, obj = objChain, core1 = col172_sequence, core2 = col173_sequence,it = 300,burnOpt = T)

# Observer 3

load(paste0(resultsdir,"Gibbs_observer3.RData"))

obs3 <- prepareObserverGibbs(param = paramChain, obj = objChain, core1 = col172_sequence, core2 = col173_sequence,it = 300,burnOpt = F, burn = 50)

int <- CombineObserverGibbs(Thick = c(obs1[[1]],obs2[[1]],obs3[[1]]), Param = list(obs1[[2]],obs2[[2]],obs3[[2]][1:1000,]), high = list(obs1[[3]],obs2[[3]],obs3[[3]]), somedepths = depths_14C, ind.plot = Bacon.plot, PlotOpt = T, dir = Gibbsresultsdir,original = c(0.05,0.05,0.1,0.2,0.3,0.6))

# Stats
# Integrated model 
round(mean(apply(int$ageEns[[1]][,2:901], FUN = function(x) max(x,na.rm = T), MARGIN = 2)),0)
round(quantile(x = apply(int$ageEns[[1]][,2:901], FUN = function(x) max(x,na.rm = T), MARGIN = 2),probs = c(0.025,0.975)),0)

PDC(3375,3137)

# Sed rates
Int.thick = c(obs1[[1]],obs2[[1]],obs3[[1]])
ens.cumSum <- lapply(Int.thick,  function(x) cumsum(abs(x[which(!is.na(x))])))
ens.age <- lapply(ens.cumSum, function(x) seq(1,length(x),1))
indx <- sapply(Int.thick, length)
Int.age <- as.data.frame(do.call(cbind,lapply(ens.age, `length<-`,max(indx))))
Int.depth <- as.data.frame(do.call(cbind,lapply(Int.thick, `length<-`,max(indx))))

count <- seq_len(nrow(Int.depth))
Int.rates.q <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(Int.depth), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(Int.depth, FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

loadOneName(objName = "AD.data", file = paste0(resultsdir,"multiple_cores_observers_workstation_1.RData"))

AD.data$Model <- c(rep("VarveR", nrow(AD.data)))

AD.data2 <- data.frame(Obs = c(rep('Obs1', nrow(I1.age2)), rep('Obs2', nrow(I2.age2)), rep('Obs3', nrow(I3.age2))), Depth = c(I1.age2$depthVec, I2.age2$depthVec, I3.age2$depthVec), c2.5 = c(I1.age2$`2.5%`,I2.age2$`2.5%`,I3.age2$`2.5%`), c25 = c(I1.age2$`25%`,I2.age2$`25%`,I3.age2$`25%`), c50 = c(I1.age2$`50%`,I2.age2$`50%`,I3.age2$`50%`), c75 = c(I1.age2$`75%`,I2.age2$`75%`,I3.age2$`75%`), c97.5 = c(I1.age2$`97.5%`,I2.age2$`97.5%`,I3.age2$`97.5%`), c2.5BP = c(I1.age2$`2.5%`,I2.age2$`2.5%`,I3.age2$`2.5%`)-67, c25BP = c(I1.age2$`25%`,I2.age2$`25%`,I3.age2$`25%`)-67, c50BP = c(I1.age2$`50%`,I2.age2$`50%`,I3.age2$`50%`)-67, c75BP = c(I1.age2$`75%`,I2.age2$`75%`,I3.age2$`75%`)-67, c97.5BP = c(I1.age2$`97.5%`,I2.age2$`97.5%`,I3.age2$`97.5%`)-67)

AD.data2$Model = c(rep("Integrated", nrow(AD.data2)))

all.data <- rbind(AD.data,AD.data2)

Obs1.plot <- AD.data2 %>%
  filter(Obs == "Obs1") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5), fill = "#999999", alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75), fill = "#999999", alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50), color = "#999999")+
  #geom_line(data =col173_sequence1, aes(x = col173_sequence1$cumsum, y = col173_sequence1$count ))+
  scale_y_continuous(breaks = seq(0,3000,500))+
  ylab(label = "Varve years")+
  xlab(label = "Depth (mm)")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        plot.margin = unit(c(0,0.1,0,0), "cm"))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse()

Obs1.plot

Obs2.plot <- AD.data2 %>%
  filter(Obs == "Obs2") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5), fill = "#56B4E9", alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75), fill = "#56B4E9", alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50), color = "#56B4E9")+
  #geom_line(data =col173_sequence2, aes(x = col173_sequence2$cumsum, y = col173_sequence2$count ))+
  scale_y_continuous(breaks = seq(0,3000,500))+
  ylab(label = "Varve years")+
  xlab("")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        axis.text.y = element_blank(), 
        plot.margin = unit(c(0,0.1,0,0), "cm")
  )+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse()

Obs2.plot

Obs3.plot <- AD.data2 %>%
  filter(Obs == "Obs3") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5), fill = "#E69F00", alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75), fill = "#E69F00", alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50), color = "#E69F00")+
  #geom_line(data =col173_sequence3, aes(x = col173_sequence3$cumsum, y = col173_sequence3$count ))+
  scale_y_continuous(breaks = seq(0,3000,500))+
  ylab(label = "Varve years")+
  xlab("")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(), 
        axis.text.y = element_blank(), 
        plot.margin = unit(c(0,0.1,0,0), "cm"))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse()

Obs3.plot

IntegratedObs123.plot <- cowplot::plot_grid(Obs1.plot,Obs2.plot, Obs3.plot,labels = c("A","B","C"), nrow = 1, align = "hv", vjust = 5.5)
IntegratedObs123.plot

# Integrated model stats

### Modeled varves stats
# Obs 1
round(mean(apply(ages1$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2)),0)
round(quantile(x = apply(ages1$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)*is.finite(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)), na.rm = T),0)
round(quantile(x = apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)*is.finite(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),probs = c(0.025,0.975), na.rm = T),0)

#Obs 2
round(mean(apply(ages2$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2)),0)
round(quantile(x = apply(ages2$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(ages2$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(ages2$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#Obs 3
round(mean(apply(ages3$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2)),0)
round(quantile(x = apply(ages3$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#COL17-3
round(mean(c(apply(ages1$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),apply(ages2$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),apply(ages3$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2))),0)

round(quantile(c(apply(ages1$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),apply(ages2$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2),apply(ages3$ageEns[[1]][,2:301], FUN = function(x) dplyr::last(x), MARGIN = 2)), probs = c(0.025,0.975)),0)

# COL17-2

round(mean(c(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)*is.finite(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)), na.rm = T),0)

round(quantile(c(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)*is.finite(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)), na.rm = T, probs = c(0.025,0.975)),0)


round(PDC(3138,3137),1)
round(PDC(3375,3137),0)

# Sed rates

mean(c(as.matrix(ages1$ThicksEns[[2]]),as.matrix(ages2$ThicksEns[[2]]),as.matrix(ages3$ThicksEns[[2]])), na.rm = T)
quantile(c(as.matrix(ages1$ThicksEns[[2]]),as.matrix(ages2$ThicksEns[[2]]),as.matrix(ages3$ThicksEns[[2]])), na.rm = T, probs = c(0.05,0.95))

I1.plot <- ggplot()+
  geom_ribbon(aes(x = -67+I1$x.bin,ymin = I1$quants[,1],ymax = I1$quants[,5]), fill = "#999999",alpha = 0.25)+
  geom_ribbon(aes(x = -67+I1$x.bin,ymin = I1$quants[,2],ymax = I1$quants[,4]), fill = "#999999", alpha = 0.5)+
  geom_line(aes(x = -67+I1$x.bin, y = I1$quants[,3] ), color = "#999999", size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(0,3)+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

I2.plot <- ggplot()+
  geom_ribbon(aes(x = -67+I2$x.bin,ymin = I2$quants[,1],ymax = I2$quants[,5]), fill = "#56B4E9",alpha = 0.25)+
  geom_ribbon(aes(x = -67+I2$x.bin,ymin = I2$quants[,2],ymax = I2$quants[,4]), fill = "#56B4E9", alpha = 0.5)+
  geom_line(aes(x = -67+I2$x.bin, y = I2$quants[,3] ), color = "#56B4E9",size = 1)+
  scale_x_reverse()+
  #ylab("")+
  ylab(label = expression("Sedimentation rate (mm yr"^{-1}*")"))+
  xlim(c(3300,-100))+
  ylim(0,3)+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

I3.plot <- ggplot()+
  geom_ribbon(aes(x = -67+I3$x.bin,ymin = I3$quants[,1],ymax = I3$quants[,5]), fill = "#E69F00",alpha = 0.25)+
  geom_ribbon(aes(x = -67+I3$x.bin,ymin = I3$quants[,2],ymax = I3$quants[,4]), fill = "#E69F00", alpha = 0.5)+
  geom_line(aes(x = -67+I3$x.bin, y = I3$quants[,3] ), color = "#E69F00",size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(0,3)+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
        
all.on.one.plot <- ggplot()+
  geom_ribbon(aes(x = -67+I1$x.bin,ymin = I1$quants[,1],ymax = I1$quants[,5]), fill = "#999999",alpha = 0.25)+
  #geom_ribbon(aes(x = -67+I1$x.bin,ymin = I1$quants[,2],ymax = I1$quants[,4]), fill = "#999999", alpha = 0.5)+
  #geom_line(aes(x = -67+I1$x.bin, y = I1$quants[,3] ), color = "#999999", size = 1)+
  geom_ribbon(aes(x = -67+I2$x.bin,ymin = I2$quants[,1],ymax = I2$quants[,5]), fill = "#56B4E9",alpha = 0.25)+
  #geom_ribbon(aes(x = -67+I2$x.bin,ymin = I2$quants[,2],ymax = I2$quants[,4]), fill = "#56B4E9", alpha = 0.5)+
  #geom_line(aes(x = -67+I2$x.bin, y = I2$quants[,3] ), color = "#56B4E9",size = 1)+
  geom_ribbon(aes(x = -67+I3$x.bin,ymin = I3$quants[,1],ymax = I3$quants[,5]), fill = "#E69F00",alpha = 0.25)+
  #geom_ribbon(aes(x = -67+I3$x.bin[1:133],ymin = I3$quants[1:133,2],ymax = I3$quants[1:133,4]), fill = "#E69F00", alpha = 0.5)+
  #geom_line(aes(x = -67+I3$x.bin[1:133], y = I3$quants[1:133,3] ), color = "#E69F00",size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(0,3)+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

all.on.one.plot

I1.r <- as.vector(I1$quants) #[,5]-I1$quants[,1]
I2.r <- as.vector(I2$quants) #[,5]-I2$quants[,1]
I3.r <- as.vector(I3$quants) #[,5]-I3$quants[,1]

mean(c(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)*is.finite(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)), na.rm = T)

quantile(c(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)*is.finite(apply(ages1$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), apply(ages3$ThicksEns[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)), na.rm = T, probs = c(0.025,0.975))

# Stats

stats <- rbind(I1.age2,I2.age2,I3.age2)
stats$obs <- c(rep("obs1",nrow(I1.age2)), rep("obs2",nrow(I2.age2)), rep("obs3", nrow(I3.age2)))
colnames(stats) <- c("depth","p2.5","p25","p50","p75","p97.5","core")

# Analysis of varve type shift

# Mid-core shift from type 1 to type 2
varve.shift <- stats %>%
  filter(between(depth, 500, 700 )) %>%
  dplyr::summarise(#age.lowVY = min(p2.5),
                   #age.highVY = max(p97.5),
                   #age.meanVY = mean(p50),
                   age.lowCE = 2017-min(p50),
                   age.highCE = 2017-max(p50),
                   age.meanCE = 2017-mean(p50),
                   #age.lowBP = 1950-(2017-min(p2.5)),
                   #age.highBP = 1950-(2017-max(p97.5)),
                   #age.meanBP = 1950-(2017-mean(p50))
                   ) 
varve.shift

# Varve shift to type 3 @ 7 cm in COL17-3
stats %>%
  filter(between(depth, 69.5, 70 )) %>%
  dplyr::summarise(
    age.lowCE = 2017-min(p50),
    age.highCE = 2017-max(p50),
    age.meanCE = 2017-mean(p50),
    ) 
# Deep massive layer start @ 8 cm in COL17-3
stats %>%
  filter(between(depth, 80, 80.5 )) %>%
  dplyr::summarise(
    age.lowCE = 2017-min(p50),
    age.highCE = 2017-max(p50),
    age.meanCE = 2017-mean(p50),
  ) 
# Shallowest massive layer start @ 2.5 cm in COL17-3
stats %>%
  filter(between(depth, 24, 25.5 )) %>%
  dplyr::summarise(
    age.lowCE = 2017-min(p50),
    age.highCE = 2017-max(p50),
    age.meanCE = 2017-mean(p50),
  ) 

Int.plot <- ggplot()+
  geom_ribbon(aes(x = -67+Int.rates.q$x.bin,ymin = Int.rates.q$quants[,1],ymax = Int.rates.q$quants[,5]), fill = "darkolivegreen4",alpha = 0.25)+
  geom_ribbon(aes(x = -67+Int.rates.q$x.bin,ymin = Int.rates.q$quants[,2],ymax = Int.rates.q$quants[,4]), fill = "darkolivegreen4", alpha = 0.5)+
  geom_line(aes(x = -67+Int.rates.q$x.bin, y = Int.rates.q$quants[,3] ), color = "darkolivegreen4", size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(0,3)+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
Int.plot

Int.with.obs.plot <- ggplot()+
  geom_ribbon(aes(x = -67+Int.rates.q$x.bin,ymin = Int.rates.q$quants[,1],ymax = Int.rates.q$quants[,5]), fill = "darkolivegreen4",alpha = 0.25)+
  geom_ribbon(aes(x = -67+Int.rates.q$x.bin,ymin = Int.rates.q$quants[,2],ymax = Int.rates.q$quants[,4]), fill = "darkolivegreen4", alpha = 0.5)+
  geom_line(aes(x = -67+Int.rates.q$x.bin, y = Int.rates.q$quants[,3] ), color = "darkolivegreen4", size = 0.1)+
  geom_line(aes(x = -67+I1$x.bin, y = I1$quants[,3] ), color = "#999999", size = 0.1)+
  geom_line(aes(x = -67+I2$x.bin, y = I2$quants[,3] ), color = "#56B4E9", size = 0.1)+
  geom_line(aes(x = -67+I3$x.bin, y = I3$quants[,3] ), color = "#E69F00", size = 0.1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(0,3)+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
Int.with.obs.plot

Int.r <- as.vector(Int.rates.q$quants) #[,5]-Int.rates.q$quants[,1]
mean(Int.r, na.rm = T)
quantile(Int.r, na.rm = T, probs = c(0.025,0.975))

Int.age <- as.data.frame(int$ageEns$summaryTable)
colnames(Int.age) <- c("depth","p2.5","p25","p50","p75","p97.5")
# Analysis of varve type shift

# Mid-core shift from type 1 to type 2
Int.age %>%
  filter(between(depth, 500, 600 )) %>% # use 500 and 600
  dplyr::summarise(age.lowVY = min(p2.5),
    age.highVY = max(p97.5),
    age.meanVY = mean(p50),
    age.lowCE = 2017-min(p2.5),
    age.highCE = 2017-max(p97.5),
    age.meanCE = 2017-mean(p50),
    age.lowBP = 1950-(2017-min(p2.5)),
    age.highBP = 1950-(2017-max(p97.5)),
    age.meanBP = 1950-(2017-mean(p50))) 


# Varve shift to type 3 @ 7 cm in COL17-3
Int.age %>%
  filter(between(depth, 69.9, 70 )) %>%
  dplyr::summarise(
    age.lowCE = 2017-min(p2.5),
    age.highCE = 2017-max(p97.5),
    age.meanCE = 2017-mean(p50),
  ) 
# Deep massive layer start @ 8 cm in COL17-3
Int.age %>%
  filter(between(depth, 80, 80.1 )) %>%
  dplyr::summarise(
    age.lowCE = 2017-min(p2.5),
    age.highCE = 2017-max(p97.5),
    age.meanCE = 2017-mean(p50),
  ) 
# Shallowest massive layer start @ 2.5 cm in COL17-3
Int.age %>%
  filter(between(depth, 25.4, 25.5 )) %>%
  dplyr::summarise(
    age.lowCE = 2017-min(p2.5),
    age.highCE = 2017-max(p97.5),
    age.meanCE = 2017-mean(p50),
  ) 

# Sed rates during the DACP
Int.rates <- data.frame(cbind(convertBP2AD(Int.rates.q$x.bin),Int.rates.q$quants))
colnames(Int.rates) <- c("age","p2.5","p25","p50","p75","p97.5")
Int.rates %>%
  filter(between(age,419,882)) %>%
  dplyr::summarise(
    sr.min = min(p50),
    sr.max = max(p50),
    sr.mean = mean(p50),
  ) 
Int.rates %>%
  filter(between(age,882,1300)) %>%
  dplyr::summarise(
    sr.min = min(p50),
    sr.max = max(p50),
    sr.mean = mean(p50),
  ) 
Int.rates %>%
  filter(between(age,1,419)) %>%
  dplyr::summarise(
    sr.min = min(p50),
    sr.max = max(p50),
    sr.mean = mean(p50),
  )

sed.rate.error <- data.frame(Observer = c(rep("Radiometric",length(b.r)), rep("1", length(O1.r)+length(I1.r)), rep("2",length(O2.r)+length(I2.r)), rep("3",length(O3.r)+length(I3.r)), rep("All",length(Int.r))), Rates = c(b.r, O1.r, O2.r,O3.r,I1.r,I2.r,I3.r, Int.r), Model = c(rep("Bacon",length(b.r)), rep("VarveR", length(O1.r)), rep("Integrated",length(I1.r)), rep("VarveR",length(O2.r)),rep("Integrated",length(I2.r)), rep("VarveR",length(O3.r)), rep("Integrated",length(I3.r)),rep("Integrated",length(Int.r))))
                             
ggplot(sed.rate.error, aes(Rates, color = Observer))+
  facet_wrap(~Model)+
  geom_density()
  theme_bw()+
  xlab("95% HDR sedimentation rates")

violin.sed.rates <- ggplot(sed.rate.error, aes(x = Model,y =Rates, fill = Observer))+
  geom_violin()+
 # facet_wrap(~Model)
  scale_fill_manual(values = c("Radiometric" = "purple", "1" = "#999999",
                                "2" = "#56B4E9", "3" = "#E69F00", 
                               "All" = "darkolivegreen4"))+
  theme_bw()+
  ylab(label = expression("Sedimentation rate (mm yr"^{-1}*")"))+
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
violin.sed.rates
