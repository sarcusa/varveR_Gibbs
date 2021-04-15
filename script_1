# Arcusa et al. Columbine Lake varve record
# To be used in conjunction with data files
# Script 1

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

rm(list = ls(all.names = TRUE))
# Info useful for all analyses

TP_depth <- c(0,8.5,20,35,52,58,63.5,66.5,114.5,124.5) # from the core surface of COL17-3

# aGE-DEPTH MODEL: using the .out from PLUM & BACON
# Steps:
# start with clean workspace
# All of this has to be run mannually

# Read lipd file. LiPD file available on Figshare: 10.6084/m9.figshare.14417999
# Lipd file must contain up to date chrondata table

C <- readLipd("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Columbine.Arcusa.2020.lpd")

# Run bacon from the lipd file and up to date chrondata table

C <- runBacon(C,lab.id.var = NULL, age.14c.var = "age14c", 
              age.14c.uncertainty.var = "age14cuncertainty", 
              age.var = "age",age.uncertainty.var = "ageuncertainty",
              depth.var = "depth", reservoir.age.14c.var = NULL, 
              reservoir.age.14c.uncertainty.var = NULL, 
              rejected.ages.var = NULL, d.max = 130, cutoff = 0.0001)

# Map the age ensemble to the paleo data, use chron table 1, measurement table 1
C = mapAgeEnsembleToPaleoData(C)

# Now go to the folder that contains the bacon output and swap the .out file with the PLUM & BACOn .out. Make sure the name replacement is the same. Then run loadbaconoutput and select which table to use. Make sure the .bacon file is there too. For example, in this instance, use the _66 files.
# .out file provided in Figshare at 10.6084/m9.figshare.14417999

C <- loadBaconOutput(C,site.name = C$dataSetName) #using chrontable 2 but measurement table 1

# Map the age ensemble to the paleo data, use chron table 2 but measurement table 1
C = mapAgeEnsembleToPaleoData(C)

# Now two chron tables exist: (1) is the IntCal20 Bacon, (2) is the Bacon and Plum combo using IntCal13

# Plot the data with geochronR. Make sure to select the correct model table. 
start = -70
end = 3500

# this plot bacon intcal20
Bacon.plot = plotChron(C, model.num = 1, age.var = "ageensemble",
                       x.bin = seq(start,end, by= 100), 
                       dist.color = "#1E88E5", truncate.dist = NA,
                       probs = c(.1,.5,.9))+
  labs(title = "", x = "Age (yr BP)", y = "Depth (cm)")+
  xlim(end,start)+ 
  theme_linedraw()+
  theme(aspect.ratio = 1, legend.position="none",
        panel.grid.major = element_line(color = "grey70"), 
        panel.grid.minor = element_line(color = "grey70"))
plot(Bacon.plot)

# this plots bacon and plum combo with intcal13
BaconPlum.plot = plotChron(C, model.num = 2, 
                           age.var = "ageensemble", 
                           x.bin = seq(-70,100, by= 10), 
                           dist.color = "#1E88E5", truncate.dist = NA,
                                      probs = c(.1,.5,.9))+
  #labs(title = "", x = "Age (yr BP)", y = "Depth (cm)")+
  labs(title = "", x = "", y = "")+
  coord_cartesian(xlim = c(100, -70), ylim = c(12, 0), expand = c(0,0))+
  #scale_x_reverse(limits = c(100, -70), expand = c(0,0))+ 
  #scale_y_reverse(limits = c(12, 0), expand = c(0,0))+
  theme_linedraw()+
  theme(aspect.ratio = 1, legend.position="none", 
        panel.grid.major = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))
plot(BaconPlum.plot)

# Combining both plots

plot.with.inset <-
  ggdraw() +
  draw_plot(Bacon.plot) +
  draw_plot(BaconPlum.plot, x = 0.42, y = .097, width = .45, height = .45)
plot(plot.with.inset)

# To extract the plum model only
# Extract Lead data from Plum model

D <- readLipd("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Columbine.Arcusa.2020.lpd")
D <- loadBaconOutput(D,site.name = D$dataSetName) # create a new model and select _101 (remember that the .out file provided by Marco needs to have the first column removed)
D = mapAgeEnsembleToPaleoData(D)

# Plot the lead with geochronR. Need to have run "multiple_cores_observers" for this to work. This does not work well...
start = -70
end = 50

Plum.plot = plotChron(D, model.num = 1, age.var = "ageensemble",
                       x.bin = seq(start,end, by= 1), 
                       dist.color = "#1E88E5", truncate.dist = NA,
                       probs = c(.1,.5,.9))+
  labs(title = "", x = "Age (yr BP)", y = "")+
  coord_cartesian(xlim = c(end, start), ylim = c(12, 0), expand = c(0,0))+
  #scale_x_reverse(limits = c(end, start), expand = c(0,0))+ 
  #scale_y_reverse(limits = c(10, 0), expand = c(0,0))+
  theme_light()+
  theme(legend.position="none", aspect.ratio = 1)
plot(Plum.plot)

plum.varve.plot <- Plum.plot+ #ggplot()+
  geom_ribbon(aes(y = AD.data$Depth/10,xmin = AD.data$c2.5-67,xmax = -67+AD.data$c97.5, fill = AD.data$Obs), alpha = 0.25)+
  geom_ribbon(aes(y = AD.data$Depth/10,xmin = -67+AD.data$c25,xmax = -67+AD.data$c75, fill = AD.data$Obs), alpha = 0.5)+ 
  #geom_line(aes(y = AD.data$Depth/10, x = AD.data$c50-67, color = AD.data$Obs))+
  #scale_x_reverse(limits = c(10,0))+
  #scale_y_reverse(limits = c(end, start), expand = c(0,0))+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))#+
  #scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  #ylab(label = "Varve years")+
  #xlab(label = "Depth (mm)")+
  #theme_bw(aspect.ratio = 1)+
  #theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
  #      panel.grid.minor = element_blank(), 
  #      legend.title = element_blank(), legend.position = c(.8,.8))+
  #coord_flip()
plum.varve.plot

## NOTE: SYMMETRICAL PRIORS IN THE FOLLOWING SECTION

library("devtools")
#library("varveR")
library("ggplot2")
library(geoChronR)
library("dplyr")
library(plyr)
library("lipdR")
library(matrixStats)
library(purrr)
library(tictoc)
library('patchwork')
library("LaplacesDemon")

# Need the following functions uploaded in the environment, available on the repo:
# GibbsRelatedFunctions.R
# importFiles.R
# plotting.R
# simulateOverAndUndercounting.R
# varveModel.R

n = 500
PlotOp = F
x <- seq(0,1,length.out = 100)

plotSymmPriors(x, dir = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/")

# Observer 1 (repeat for each observer, changing the input data)

# Data available here: 10.6084/m9.figshare.14251400
# Download and change directory in code below

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Charlotte_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Charlotte_counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2_charlotte.pdf")
plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3_charlotte.pdf")
}

col172_sequence1 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence1 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence1 <- calibrateCummulativeDepth(CombSeq = col172_sequence1, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence1 <- calibrateCummulativeDepth(CombSeq = col173_sequence1, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
plotCompositeSequence(col172_sequence1)
plotCompositeSequence(col173_sequence1)
}

modeledVarves1 <- list()
for(i in 1:n){
  print(i)
  
  modeledVarves1[[i]] <- SYMgammaFastVarveModel(core1 = col172_sequence1, 
                         core2 = col173_sequence1)
  
}

modeledVarves1 <- formatOutput(modeledVarves1)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves1[[1]][,2])), y = modeledVarves1[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves1[[1]][,1])), y = modeledVarves1[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model
#choose the dated core (core = 2, COL17-3)
calModeledVarves1 <- calibrateVarveModelDepth(data = modeledVarves1, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2) 

Vplot1 <- plotVarveModel(calModeledVarves1)
print(Vplot1)

# Keep just in case
#probabilities based on observations
#transtable <- data.frame(codes = 1:6,overcount = c(0.05,0.10,0.30,0.5,-2,-3), undercount =  c(0.05,0.10,0.30,0,-2,-3)) 
#transtable <- data.frame(codes = 1:6,overcount = c(0.01,0.05,0.10,0.5,-2,-3), undercount =  c(0.02,0.10,0.20,0,-2,-3)) 

transtable <- data.frame(codes = 1:6,overcount = c(0.01,0.05,0.10,0.5,-2,-3), undercount =  c(0.01,0.05,0.1,0,-2,-3)) 
col172_sequence1 <- translateCodesToProbabilities(col172_sequence1,translationTable = transtable,baselineProb = .05)
col173_sequence1 <- translateCodesToProbabilities(col173_sequence1,translationTable = transtable,baselineProb = .05)

# Convert factors to characters
col172_sequence1$tiePoint = as.character(col172_sequence1$tiePoint)
col173_sequence1$tiePoint = as.character(col173_sequence1$tiePoint)

col172Ens1 <- generateThicknessEnsemble(col172_sequence1,nEns = 100)
col173Ens1 <- generateThicknessEnsemble(col173_sequence1,nEns = 100)

ensList1 <- list(col172Ens1,col173Ens1)

allTiePoints <- apply(cbind(na.omit(unique(col172_sequence1$tiePoint)),na.omit(unique(col173_sequence1$tiePoint))),1,unique)
allMarkerLayers <- allTiePoints
allMarkerLayers <- c("topOfCore",allMarkerLayers)

modeledVarves1 <- varveModel(ensList1, nSim = n, allMarkerLayers)

# Observer 2 (repeat for each observer, changing the input data)

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Sela counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Sela counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2_Sela.pdf")
plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3_Sela.pdf")
}

col172_sequence2 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence2 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence2 <- calibrateCummulativeDepth(CombSeq = col172_sequence2, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence2 <- calibrateCummulativeDepth(CombSeq = col173_sequence2, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
plotCompositeSequence(col172_sequence2)
plotCompositeSequence(col173_sequence2)
}

modeledVarves2 <- list()
for(i in 1:n){
  print(i)
  
  modeledVarves2[[i]] <- SYMgammaFastVarveModel(core1 = col172_sequence2, 
                                                core2 = col173_sequence2)
  
}

modeledVarves2 <- formatOutput(modeledVarves2)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves2[[1]][,2])), y = modeledVarves2[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves2[[1]][,1])), y = modeledVarves2[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model

calModeledVarves2 <- calibrateVarveModelDepth(data = modeledVarves2, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2)

Vplot2 <- plotVarveModel(calModeledVarves2)
print(Vplot2)

# Keep

#probabilities based on observations
#transtable <- data.frame(codes = 1:6,overcount = c(0.01,0.05,0.10,0.5,-2,-3), undercount =  c(0.02,0.10,0.20,0,-2,-3)) 
transtable <- data.frame(codes = 1:6,overcount = c(0.01,0.05,0.10,0.5,-2,-3), undercount =  c(0.01,0.05,0.1,0,-2,-3)) 
col172_sequence2 <- translateCodesToProbabilities(col172_sequence2,translationTable = transtable,baselineProb = .05)
col173_sequence2 <- translateCodesToProbabilities(col173_sequence2,translationTable = transtable,baselineProb = .05)

# Convert factors to characters
col172_sequence2$tiePoint = as.character(col172_sequence2$tiePoint)
col173_sequence2$tiePoint = as.character(col173_sequence2$tiePoint)

col172Ens2 <- generateThicknessEnsemble(col172_sequence2,nEns = 100)
col173Ens2 <- generateThicknessEnsemble(col173_sequence2,nEns = 100)

ensList2 <- list(col172Ens2,col173Ens2)

allTiePoints <- apply(cbind(na.omit(unique(col172_sequence2$tiePoint)),na.omit(unique(col173_sequence2$tiePoint))),1,unique)
allMarkerLayers <- allTiePoints
allMarkerLayers <- c("topOfCore",allMarkerLayers)

modeledVarves2 <- varveModel(ensList2, nSim = n, allMarkerLayers)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves2[[1]][,2])), y = modeledVarves2[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves2[[1]][,1])), y = modeledVarves2[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model

calModeledVarves2 <- calibrateVarveModelDepth(data = modeledVarves2, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2)

Vplot2 <- plotVarveModel(calModeledVarves2)
print(Vplot2)

# Observer 3 (repeat for each observer, changing the input data)

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2_Steph.pdf")
plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3_Steph.pdf")
}

col172_sequence3 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence3 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence3 <- calibrateCummulativeDepth(CombSeq = col172_sequence3, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence3 <- calibrateCummulativeDepth(CombSeq = col173_sequence3, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
plotCompositeSequence(col172_sequence3)
plotCompositeSequence(col173_sequence3)
}

modeledVarves3 <- list()
for(i in 1:n){
  print(i)
  
  modeledVarves3[[i]] <- SYMgammaFastVarveModel(core1 = col172_sequence3, 
                                                core2 = col173_sequence3)
  
}

modeledVarves3 <- formatOutput(modeledVarves3)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves3[[1]][,2])), y = modeledVarves3[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves3[[1]][,1])), y = modeledVarves3[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model

calModeledVarves3 <- calibrateVarveModelDepth(data = modeledVarves3, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2)

Vplot3 <- plotVarveModel(calModeledVarves3)
print(Vplot3)

##### Comparisons
## Plot raw observer counts for each core

COL172_all <- data.frame(Obs = c(rep('Obs1', length(col172_sequence1$section)), rep('Obs2', length(col172_sequence2$section)), rep('Obs3', length(col172_sequence3$section))), Count = c(col172_sequence1$count, col172_sequence2$count, col172_sequence3$count), Thick = c(col172_sequence1$thick, col172_sequence2$thick, col172_sequence3$thick), QC = c(col172_sequence1$varveCode, col172_sequence2$varveCode, col172_sequence3$varveCode), core = c(rep('COL17-2', length(sum(length(col172_sequence1$section), length(col172_sequence2$section), length(col172_sequence3$section))))), y = rep(1,length(col172_sequence1$section)+length(col172_sequence2$section)+length(col172_sequence3$section)))

COL173_all <- data.frame(Obs = c(rep('Obs1', length(col173_sequence1$section)), rep('Obs2', length(col173_sequence2$section)), rep('Obs3', length(col173_sequence3$section))), Count = c(col173_sequence1$count, col173_sequence2$count, col173_sequence3$count), Thick = c(col173_sequence1$thick, col173_sequence2$thick, col173_sequence3$thick), QC = c(col173_sequence1$varveCode, col173_sequence2$varveCode, col173_sequence3$varveCode), core = c(rep('COL17-3', length(sum(length(col173_sequence1$section), length(col173_sequence2$section), length(col173_sequence3$section))))),y = rep(1,length(col173_sequence1$section)+length(col173_sequence2$section)+length(col173_sequence3$section)))

core1 <- ggplot(COL172_all)+ 
  geom_line(aes(x = Count, y = Thick, color = Obs))+
  facet_grid(rows = vars(Obs), scales = "free_y")+
  ggtitle('COL17-2')+
  scale_y_continuous(limits = c(0,20))+
  #scale_x_continuous(labels = NULL)+
  #xlab('')+
  ylab('Thickness (mm)')+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  scale_x_reverse()+ #added this, may not work
  theme_light()+
  theme(legend.position="none", panel.spacing = unit(2, "lines"))
  
core2 <- ggplot(COL173_all)+
  geom_line(aes(x = Count, y = Thick, color = Obs))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  facet_grid(rows = vars(Obs), scales = "free_y")+
  theme_light()+
  scale_y_continuous(limits = c(0,20))+
  theme(legend.position="none", panel.spacing = unit(2, "lines"))+
  ggtitle('COL17-3')+
  scale_x_reverse(limits = c(2500,0))+ #added this, may not work
  ylab('')

both_counts <- core1 + core2 + plot_layout(nrow = 1)
both_counts

# Varve quality plots

greens <- c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","#edf8fb")

core1 <- ggplot(COL172_all,aes(x = Count, y =1, fill = factor(QC)))+ 
  geom_tile()+
  facet_grid(rows = vars(Obs))+
  ggtitle('COL17-2')+
  scale_fill_manual(values=greens)+
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank())
core1

core2 <- ggplot(COL173_all,aes(x = Count, y =1, fill = factor(QC)))+ 
  geom_tile()+
  facet_grid(rows = vars(Obs))+
  ggtitle('COL17-3')+
  scale_fill_manual(values=greens)+
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank())
core2

both_QC <-core1 + core2 + plot_layout(nrow = 2)
both_QC

# Compare age-depth models

AD.data <- data.frame(Obs = c(rep('Obs1', nrow(Vplot1$data)), rep('Obs2', nrow(Vplot2$data)), rep('Obs3', nrow(Vplot3$data))), Depth = c( Vplot1$data$depth, Vplot2$data$depth, Vplot3$data$depth), c2.5 = c(Vplot1$data$c2.5,Vplot2$data$c2.5,Vplot3$data$c2.5), c25 = c(Vplot1$data$c25,Vplot2$data$c25,Vplot3$data$c25), c50 = c(Vplot1$data$c50,Vplot2$data$c50,Vplot3$data$c50), c75 = c(Vplot1$data$c75,Vplot2$data$c75,Vplot3$data$c75), c97.5 = c(Vplot1$data$c97.5,Vplot2$data$c97.5,Vplot3$data$c97.5), c2.5BP = c(Vplot1$data$c2.5,Vplot2$data$c2.5,Vplot3$data$c2.5)-67, c25BP = c(Vplot1$data$c25,Vplot2$data$c25,Vplot3$data$c25)-67, c50BP = c(Vplot1$data$c50,Vplot2$data$c50,Vplot3$data$c50)-67, c75BP = c(Vplot1$data$c75,Vplot2$data$c75,Vplot3$data$c75)-67, c97.5BP = c(Vplot1$data$c97.5,Vplot2$data$c97.5,Vplot3$data$c97.5)-67)

output1 <- data.frame(Depth = round(Vplot1$data$depth[1:260],1), Age = round(2018-Vplot1$data$c50[1:260],0))
output2 <- data.frame(Depth = round(Vplot2$data$depth[1:260],1), Age = round(2018-Vplot2$data$c50[1:260],0))
output3 <- data.frame(Depth = round(Vplot3$data$depth[1:260],1), Age = round(2018-Vplot1$data$c50[1:260],0))

CompareAgeDepthPlot <- ggplot(AD.data)+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5, fill = Obs), alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75, fill = Obs), alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50, color = Obs))+
  scale_y_continuous(breaks = seq(0,3000,500))+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  ylab(label = "Varve years")+
  xlab(label = "Depth (mm)")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(), legend.position = c(.8,.2))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse()

CompareAgeDepthPlot

compareVarveBacon <- Bacon.plot+
  geom_ribbon(aes(y = AD.data$Depth/10, xmin = AD.data$c2.5BP,xmax = AD.data$c97.5BP, fill = AD.data$Obs), alpha = .25)+
  geom_ribbon(aes(y = AD.data$Depth/10, xmin = AD.data$c25BP,xmax = AD.data$c75BP, fill = AD.data$Obs), alpha = .5)+ 
  geom_line(aes(x = AD.data$c50BP, y = AD.data$Depth/10, color = AD.data$Obs))+
  scale_y_reverse()+
  scale_x_reverse()+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00","red","red", "red", "red", "red"))

compareVarveBacon

stats <- rbind(col172_sequence1, col172_sequence2, col172_sequence3,
               col173_sequence1[,1:9],col173_sequence2[,1:9],
               col173_sequence3[,1:9])

stats$core <- c(rep('COL17-2', nrow(col172_sequence1)+nrow(col172_sequence2)+nrow(col172_sequence3)),rep('COL17-3', nrow(col173_sequence1)+nrow(col173_sequence2)+nrow(col173_sequence3)))

stats$obs <- c(rep('obs1', nrow(col172_sequence1)),
               rep('obs2', nrow(col172_sequence2)), 
               rep('obs3',nrow(col172_sequence3)),
               rep('obs1', nrow(col173_sequence1)),
               rep('obs2', nrow(col173_sequence2)), 
               rep('obs3',nrow(col173_sequence3)))

stats %>%
  filter(!varveCode %in% c(4,5,6)) %>%
  #group_by(.dots = c("obs", "core")) %>%
  group_by(obs, core, .drop = F) %>%
  dplyr::summarise(total = round(length(thick),2),
                   minThick = round(min(thick),2),
                   maxThick = round(max(thick),2),
                   meanThick = round(mean(thick),2),
                   medianThick = round(median(thick),2),
                   sdThick = round(sd(thick),2)) %>%
  group_by(core, .drop = F) %>%
  dplyr::summarise(core_total = mean(total))


stats %>%
  filter(varveCode %in% 4) %>%
  group_by(obs, core, .drop = F) %>%
  dplyr::summarise(cumulative = length(thick))

ML.172.1 <- diff(c(1,which(!is.na(col172_sequence1$tiePoint))))
ML.172.2 <- diff(c(1,which(!is.na(col172_sequence2$tiePoint))))
ML.172.3 <- diff(c(1,which(!is.na(col172_sequence3$tiePoint))))

ML.173.1 <- diff(c(1,which(!is.na(col173_sequence1$tiePoint))))
ML.173.2 <- diff(c(1,which(!is.na(col173_sequence2$tiePoint))))
ML.173.3 <- diff(c(1,which(!is.na(col173_sequence3$tiePoint))))

ML <- data.frame(COL172 = c(ML.172.1,ML.172.2,ML.172.3), 
           COL173 = c(ML.173.1,ML.173.2,ML.173.3), 
           Obs = c(rep(1,length(ML.172.1)),rep(2,length(ML.172.2)),rep(3,length(ML.172.3)))) 
 
ML <- ML %>%
  mutate(Difference = COL172 - COL173)
 
ML$PercentDiff <- PDC(ML$COL172, ML$COL173)

write.csv(ML,file = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/marker_layer_201020.csv")

stable <- stats %>%
  filter(!varveCode %in% c(4,5,6))%>%
  group_by(core,.drop = F) %>%
  select(thick)

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
a <- min(c(stable$thick[1:6056], stable$thick[6057:12395]))
b <- max(c(stable$thick[1:6056], stable$thick[6057:12395]))
ax <- pretty(a:b, n = 12)
col2 <- hist(stable$thick[1:6056], plot = F, breaks = ax)
col3 <- hist(stable$thick[6057:12395], plot = F, breaks = ax)
plot(col2, col = c1, ylim = c(0,3800), main = "", xlab = "Thickness (mm)")
plot(col3, col = c2, add = T, ylim = c(0,3800), xlim = c(0,3))
legend("topright", legend = c("COL17-2", "COL17-3"), col = c(c1,c2),pch = 15)

### Varve code 6 analysis

stats %>%
  filter(varveCode %in% 6) %>%
  group_by(core,obs, .drop = F) %>%
  dplyr::summarise(cumulative = sum(thick)) %>%
  group_by(core,.drop = F) %>%
  dplyr::summarise(av.extra = mean(cumulative))

### Varve code 5 analysis

stats %>%
  select(section, thick, varveCode, obs, core) %>%
  filter(varveCode %in% 5) %>%
  group_by(core,obs, .drop = F) %>%
  dplyr::summarise(cumulative = sum(thick),
                   maxi = max(thick),
                   av = mean(thick)) %>%
  group_by(core,.drop = F) %>%
  dplyr::summarise(av.extra = mean(cumulative),
                   av.mean = mean(av))

249/1270*100 #core 2
107/1230*100 #core 3

### Varve code 2, 3 and 4 analysis

stats %>%
  select(section, thick, varveCode, obs, core) %>%
  filter(varveCode %in% c(2,3,4)) %>%
  group_by(core,obs, .drop = F) %>%
  dplyr::summarise(cumulative = sum(thick)) %>%
  group_by(core,.drop = F) %>%
  dplyr::summarise(av.extra = mean(cumulative))
  
994/1270*100 #core 2
964/1230*100 #core 3

# Mid-core shift from type 1 to type 2
varve.shift <- stats %>%
  select(section, thick, varveCode, obs, core) %>%
  group_by(core,obs,.drop = F) %>%
  dplyr::mutate(cumulative = cumsum(thick)) %>%
  dplyr::mutate(number = 1) %>%
  dplyr::mutate(age = cumsum(number)) %>%
  filter(between(cumulative, 500, 700 )) %>%
  ungroup() %>%
  dplyr::summarise(age.lowCE = 2017-min(age),
                   age.highCE = 2017-max(age)) 

varve.shift

varve.shift <- stats %>%
  select(section, thick, varveCode, obs, core) %>%
  group_by(core,obs,.drop = F) %>%
  dplyr::mutate(cumulative = cumsum(thick)) %>%
  dplyr::mutate(number = 1) %>%
  dplyr::mutate(age = cumsum(number)) %>%
  filter(between(age, 30, 50 )) %>%
  ungroup() %>%
  dplyr::summarise(min.depth = min(cumulative),
                   max.depth = max(age)) 
  

varve.shift

# Start of deepest massive layer
stats %>%
  select(section, thick, varveCode, obs, core) %>%
  group_by(core,obs,.drop = F) %>%
  filter(core == "COL17-2") %>%
  dplyr::mutate(cumulative = cumsum(thick)) %>%
  dplyr::mutate(number = 1) %>%
  dplyr::mutate(age = cumsum(number)) %>%
  filter(between(cumulative, 80, 95 )) %>%
  ungroup() %>%
  dplyr::summarise(age.lowCE = 2017-min(age),
                   age.highCE = 2017-max(age),
                   age.meanCE = 2017-mean(age)) 

# End the deepest massive layer
stats %>%
  select(section, thick, varveCode, obs, core) %>%
  group_by(core,obs,.drop = F) %>%
  filter(core == "COL17-2") %>%
  dplyr::mutate(cumulative = cumsum(thick)) %>%
  dplyr::mutate(number = 1) %>%
  dplyr::mutate(age = cumsum(number)) %>%
  filter(between(cumulative, 79, 80)) %>%
  ungroup() %>%
  dplyr::summarise(age.lowCE = 2017-min(age),
                   age.highCE = 2017-max(age),
                   age.meanCE = 2017-mean(age)) 

### Modeled varves stats
# Obs 1
round(mean(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#Obs 2
round(mean(apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#Obs 3
round(mean(apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#COL17-2
round(mean(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) )),0)

round(quantile(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) ), probs = c(0.025,0.975)),0)

#COL17-3
round(mean(c(apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) )),0)

round(quantile(c(apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) ), probs = c(0.025,0.975)),0)

# All
mean(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)  ))

quantile(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)  ), probs = c(0.025,0.975))

# Percentage diff calculator

round(PDC(v1 = 2419, v2 = 2380),1)
round(PDC(v1 = 2847, v2 = 2423),0)

bacon.ages <- selectData(L = C) # pick 1, then 4
bacon.depths <- selectData(L= C) # pick 1, then 1
ageDiff = diff(bacon.ages$values)
depthDiff = diff(bacon.depths$values*10)
replicated.depth.matrix = matrix(data=depthDiff,nrow = length(depthDiff),ncol = ncol(ageDiff))
srEns = replicated.depth.matrix/ageDiff
empty.values = as.matrix(t(rep(NaN, 1000)))
srEns = rbind(srEns, empty.values)
bacon.rates <- quantile2d(bacon.ages$values,srEns,probs = c(0.025,0.25,0.5,0.75,0.975))

mean(bacon.rates$quants, na.rm = T)
quantile(bacon.rates$quants, na.rm = T, probs = c(0.05,0.95))

count <- seq_len(nrow(modeledVarves1$PostInfo[[2]]))
M1 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(modeledVarves1$PostInfo[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(modeledVarves1$PostInfo[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

count <- seq_len(nrow(modeledVarves2$PostInfo[[2]]))
M2 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(modeledVarves2$PostInfo[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(modeledVarves2$PostInfo[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

count <- seq_len(nrow(modeledVarves3$PostInfo[[2]]))
M3 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(modeledVarves3$PostInfo[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(modeledVarves3$PostInfo[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

mean(c(M2O1,M2O2,M2O3), na.rm = T)
quantile(c(M2O1,M2O2,M2O3), na.rm = T, probs = c(0.025,0.975))

Bacon.sed.plot <- plotTimeseriesEnsRibbons(X = bacon.ages,Y = srEns,probs = c(.05,0.25,.5,.75,.95),color.low = "plum1", color.high = "blue",alp = 0.25, color.line = "purple")+
  scale_x_reverse()+xlab("Age (years BP)")+
  xlim(c(3300,-100))+
  labs(title = "",y = "")+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
Bacon.sed.plot
  
M1.plot <- ggplot()+
  geom_ribbon(aes(x = -67+M1$x.bin,ymin = M1$quants[,1],ymax = M1$quants[,5]), fill = "#999999",alpha = 0.25)+
  geom_ribbon(aes(x = -67+M1$x.bin,ymin = M1$quants[,2],ymax = M1$quants[,4]), fill = "#999999", alpha = 0.5)+
  geom_line(aes(x = -67+M1$x.bin, y = M1$quants[,3] ), color = "#999999", size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(c(0,3))+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
  
M2.plot <- ggplot()+
  geom_ribbon(aes(x = -67+M2$x.bin,ymin = M2$quants[,1],ymax = M2$quants[,5]), fill = "#56B4E9",alpha = 0.25)+
  geom_ribbon(aes(x = -67+M2$x.bin,ymin = M2$quants[,2],ymax = M2$quants[,4]), fill = "#56B4E9", alpha = 0.5)+
  geom_line(aes(x = -67+M2$x.bin, y = M2$quants[,3] ), color = "#56B4E9",size = 1)+
  scale_x_reverse()+
  #ylab("")+
  ylab(label = expression("Sedimentation rate (mm yr"^{-1}*")"))+
  xlim(c(3300,-100))+
  ylim(c(0,3))+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

M3.plot <- ggplot()+
  geom_ribbon(aes(x = -67+M3$x.bin,ymin = M3$quants[,1],ymax = M3$quants[,5]), fill = "#E69F00",alpha = 0.25)+
  geom_ribbon(aes(x = -67+M3$x.bin,ymin = M3$quants[,2],ymax = M3$quants[,4]), fill = "#E69F00", alpha = 0.5)+
  geom_line(aes(x = -67+M3$x.bin, y = M3$quants[,3] ), color = "#E69F00",size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(c(0,3))+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
        
## NOTE: ASYMMETRICAL PRIORS IN THIS SECTION

library("devtools")
#library("varveR")
library("ggplot2")
library(geoChronR)
library("dplyr")
library(plyr)
library("lipdR")
library(matrixStats)
library(purrr)
library(tictoc)
library('patchwork')
library("LaplacesDemon")

# Observer 1 (use same data as for the symmetrical priors)

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Charlotte_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Charlotte_counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
  plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2_charlotte.pdf")
  plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3_charlotte.pdf")
}

col172_sequence1 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence1 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence1 <- calibrateCummulativeDepth(CombSeq = col172_sequence1, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence1 <- calibrateCummulativeDepth(CombSeq = col173_sequence1, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
  plotCompositeSequence(col172_sequence1)
  plotCompositeSequence(col173_sequence1)
}

modeledVarves1 <- list()
for(i in 1:n){
  print(i)
  
  modeledVarves1[[i]] <- ASYMgammaFastVarveModel(core1 = col172_sequence1, core2 = col173_sequence1)
  
}

modeledVarves1 <- formatOutput(modeledVarves1)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves1[[1]][,2])), y = modeledVarves1[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves1[[1]][,1])), y = modeledVarves1[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model
#choose the dated core (core = 2, COL17-3)
#calModeledVarves1 <- calibrateVarveModelDepth(data = modeledVarves1, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2) 

VAplot1 <- plotVarveModel(modeledVarves1)
print(VAplot1)

# Observer 2

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Sela counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Sela counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
  plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2_Sela.pdf")
  plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3_Sela.pdf")
}

col172_sequence2 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence2 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence2 <- calibrateCummulativeDepth(CombSeq = col172_sequence2, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence2 <- calibrateCummulativeDepth(CombSeq = col173_sequence2, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
  plotCompositeSequence(col172_sequence2)
  plotCompositeSequence(col173_sequence2)
}

modeledVarves2 <- list()
for(i in 1:n){
  print(i)
  
  modeledVarves2[[i]] <- ASYMgammaFastVarveModel(core1 = col172_sequence2,  core2 = col173_sequence2)
  
}

modeledVarves2 <- formatOutput(modeledVarves2)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves2[[1]][,2])), y = modeledVarves2[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves2[[1]][,1])), y = modeledVarves2[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model

#calModeledVarves2 <- calibrateVarveModelDepth(data = modeledVarves2, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2)

VAplot2 <- plotVarveModel(modeledVarves2)
print(VAplot2)

# Observer 3

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
  plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2_Steph.pdf")
  plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3_Steph.pdf")
}

col172_sequence3 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence3 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence3 <- calibrateCummulativeDepth(CombSeq = col172_sequence3, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence3 <- calibrateCummulativeDepth(CombSeq = col173_sequence3, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
  plotCompositeSequence(col172_sequence3)
  plotCompositeSequence(col173_sequence3)
}

modeledVarves3 <- list()
for(i in 1:n){
  print(i)
  
  modeledVarves3[[i]] <- ASYMgammaFastVarveModel(core1 = col172_sequence3, core2 = col173_sequence3)
  
}

modeledVarves3 <- formatOutput(modeledVarves3)

if(PlotOp == TRUE){
  ggplot() +
    geom_line(aes(x = 2018-seq_len(length(modeledVarves3[[1]][,2])), y = modeledVarves3[[1]][,2]))+
    geom_line(aes(x = 2018-seq_len(length(modeledVarves3[[1]][,1])), y = modeledVarves3[[1]][,1]+25),color = "red")+
    theme_bw()
}

# Creating varve age depth model

#calModeledVarves3 <- calibrateVarveModelDepth(data = modeledVarves3, botOfVarvedSequence = 1230, ML = allMarkerLayers, n = n, core = 2)

VAplot3 <- plotVarveModel(modeledVarves3)
print(VAplot3)

# Compare to symmetrical VarveR

VA.AD.data <- data.frame(Obs = c(rep('Obs1', nrow(VAplot1$data)), rep('Obs2', nrow(VAplot2$data)), rep('Obs3', nrow(VAplot3$data))), Depth = c( VAplot1$data$depth, VAplot2$data$depth, VAplot3$data$depth), c2.5 = c(VAplot1$data$c2.5,VAplot2$data$c2.5,VAplot3$data$c2.5), c25 = c(VAplot1$data$c25,VAplot2$data$c25,VAplot3$data$c25), c50 = c(VAplot1$data$c50,VAplot2$data$c50,VAplot3$data$c50), c75 = c(VAplot1$data$c75,VAplot2$data$c75,VAplot3$data$c75), c97.5 = c(VAplot1$data$c97.5,VAplot2$data$c97.5,VAplot3$data$c97.5), c2.5BP = c(VAplot1$data$c2.5,VAplot2$data$c2.5,VAplot3$data$c2.5)-67, c25BP = c(VAplot1$data$c25,VAplot2$data$c25,VAplot3$data$c25)-67, c50BP = c(VAplot1$data$c50,VAplot2$data$c50,VAplot3$data$c50)-67, c75BP = c(VAplot1$data$c75,VAplot2$data$c75,VAplot3$data$c75)-67, c97.5BP = c(VAplot1$data$c97.5,VAplot2$data$c97.5,VAplot3$data$c97.5)-67)

VA.CompareAgeDepthPlot <- ggplot(VA.AD.data)+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5, fill = Obs), alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75, fill = Obs), alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50, color = Obs))+
  scale_y_continuous(breaks = seq(0,3000,500))+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  ylab(label = "Varve years")+
  xlab(label = "Depth (mm)")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(), legend.position = c(.8,.2))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse()

VA.CompareAgeDepthPlot

varve.data <- rbind(AD.data, VA.AD.data)
varve.data$Level <- c(rep("Symmetrical",nrow(AD.data)),rep("Asymmetrical",nrow(VA.AD.data)))

LINES <- c("c2.5" = "NA", "c97.5" = "NA","c25" = "NA","c75" = "NA","c50" = "solid")

Obs1.plot <- varve.data %>%
  filter(Obs == "Obs1") %>%
  #filter(Level == 'Asymmetrical') %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5, fill = Level, color = Level), alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75, fill = Level, color = Level), alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50, color = Level, linetype = Level))+
  scale_fill_manual(values = c("#999999","NA"))+
  scale_linetype_manual(values = c('solid','longdash'))+
  scale_color_manual(values = c("#999999","black"))+
  #scale_y_continuous()+
  ylab(label = "Varve years")+
  xlab(label = "Depth (mm)")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        plot.margin = unit(c(0,0.1,0,0), "cm"))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse(breaks = seq(0,4000,1000))
Obs1.plot

Obs2.plot <- varve.data %>%
  filter(Obs == "Obs2") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5, fill = Level, color = Level), alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75, fill = Level, color = Level), alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50, color = Level, linetype = Level))+
  scale_fill_manual(values = c("#56B4E9","NA"))+
  scale_linetype_manual(values = c('solid','longdash'))+
  scale_color_manual(values = c("#56B4E9","black"))+
  #scale_y_continuous(breaks = seq(0,4000,500), )+
  ylab(label = "Varve years")+
  xlab(label = "Depth (mm)")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        plot.margin = unit(c(0,0.1,0,0), "cm"))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse(breaks = seq(0,4000,1000), limits = c(4000,0))
Obs2.plot

Obs3.plot <- varve.data %>%
  filter(Obs == "Obs3") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5, fill = Level, color = Level), alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75, fill = Level, color = Level), alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50, color = Level, linetype = Level))+
  scale_fill_manual(values = c("#E69F00","NA"))+
  scale_linetype_manual(values = c('solid','longdash'))+
  scale_color_manual(values = c("#E69F00","black"))+
  #scale_y_continuous(breaks = seq(0,4000,500))+
  ylab(label = "Varve years")+
  xlab(label = "Depth (mm)")+
  theme_bw()+
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        plot.margin = unit(c(0,0.1,0,0), "cm"))+
  coord_flip()+
  scale_x_reverse()+
  scale_y_reverse(breaks = seq(0,4000,1000))
Obs3.plot

Obs123.plot <- cowplot::plot_grid(Obs1.plot,Obs2.plot, Obs3.plot,labels = c("A","B","C"), nrow = 2, align = "hv", vjust = 5.5)
Obs123.plot

compareVarveBacon <- Bacon.plot+
  geom_ribbon(aes(y = VA.AD.data$Depth/10, xmin = VA.AD.data$c2.5BP,xmax = VA.AD.data$c97.5BP, fill = VA.AD.data$Obs), alpha = .25)+
  geom_ribbon(aes(y = VA.AD.data$Depth/10, xmin = VA.AD.data$c25BP,xmax = VA.AD.data$c75BP, fill = VA.AD.data$Obs), alpha = .5)+ 
  geom_line(aes(x = VA.AD.data$c50BP, y = VA.AD.data$Depth/10, color = VA.AD.data$Obs))+
  scale_y_reverse()+
  scale_x_reverse()+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00","red","red", "red", "red", "red"))
  
### Modeled varves stats
# Obs 1
round(mean(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#Obs 2
round(mean(apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#Obs 3
round(mean(apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

round(mean(apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0)
round(quantile(x = apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),probs = c(0.025,0.975)),0)

#COL17-2
round(mean(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) )),0)

round(quantile(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) ), probs = c(0.025,0.975)),0)

#COL17-3
round(mean(c(apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) )),0)

round(quantile(c(apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2) ), probs = c(0.025,0.975)),0)

# All
mean(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)  ))

quantile(c(apply(modeledVarves1[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves1[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves2[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2),apply(modeledVarves3[[2]][[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)  ), probs = c(0.025,0.975))

round(PDC(v1 = 2813, v2 = 2740),0)
round(PDC(v1 = 2813, v2 = 2740),1)

# Sedimentation rates

bacon.ages <- selectData(L = C) # pick 1, then 4
bacon.depths <- selectData(L= C) # pick 1, then 1
ageDiff = diff(bacon.ages$values)
depthDiff = diff(bacon.depths$values*10)
replicated.depth.matrix = matrix(data=depthDiff,nrow = length(depthDiff),ncol = ncol(ageDiff))
srEns = replicated.depth.matrix/ageDiff
empty.values = as.matrix(t(rep(NaN, 1000)))
srEns = rbind(srEns, empty.values)
bacon.rates <- quantile2d(bacon.ages$values,srEns,probs = c(0.025,0.25,0.5,0.75,0.975))

mean(bacon.rates$quants, na.rm = T)
quantile(bacon.rates$quants, na.rm = T, probs = c(0.05,0.95))

count <- seq_len(nrow(modeledVarves1$PostInfo[[2]]))
M1 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(modeledVarves1$PostInfo[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(modeledVarves1$PostInfo[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

count <- seq_len(nrow(modeledVarves2$PostInfo[[2]]))
M2 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(modeledVarves2$PostInfo[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(modeledVarves2$PostInfo[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

count <- seq_len(nrow(modeledVarves3$PostInfo[[2]]))
M3 <- plotTimeseriesEnsRibbons(X = -67+count, Y = as.matrix(modeledVarves3$PostInfo[[2]]), probs = c(0.05,0.25,0.5,0.75,0.95), color.low = "#999999", color.high = "#999999", alp = 0.25,x.bin = -67:round(mean(apply(modeledVarves3$PostInfo[[2]], FUN = function(x) min(which(is.na(x))), MARGIN = 2)),0), limit.outliers.x = T, export.quantiles = T)

mean(c(as.matrix(modeledVarves1$PostInfo[[2]]),as.matrix(modeledVarves2$PostInfo[[2]]),as.matrix(modeledVarves3$PostInfo[[2]])), na.rm = T)
quantile(c(as.matrix(calModeledVarves1$VarveThick),as.matrix(calModeledVarves2$VarveThick),as.matrix(calModeledVarves3$VarveThick)), na.rm = T, probs = c(0.025,0.975))

M1.plot <- ggplot()+
  geom_ribbon(aes(x = -67+M1$x.bin,ymin = M1$quants[,1],ymax = M1$quants[,5]), fill = "#999999",alpha = 0.25)+
  geom_ribbon(aes(x = -67+M1$x.bin,ymin = M1$quants[,2],ymax = M1$quants[,4]), fill = "#999999", alpha = 0.5)+
  geom_line(aes(x = -67+M1$x.bin, y = M1$quants[,3] ), color = "#999999", size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(c(0,3))+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

M2.plot <- ggplot()+
  geom_ribbon(aes(x = -67+M2$x.bin,ymin = M2$quants[,1],ymax = M2$quants[,5]), fill = "#56B4E9",alpha = 0.25)+
  geom_ribbon(aes(x = -67+M2$x.bin,ymin = M2$quants[,2],ymax = M2$quants[,4]), fill = "#56B4E9", alpha = 0.5)+
  geom_line(aes(x = -67+M2$x.bin, y = M2$quants[,3] ), color = "#56B4E9",size = 1)+
  scale_x_reverse()+
  #ylab("")+
  ylab(label = expression("Sedimentation rate (mm yr"^{-1}*")"))+
  xlim(c(3300,-100))+
  ylim(c(0,3))+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

M3.plot <- ggplot()+
  geom_ribbon(aes(x = -67+M3$x.bin,ymin = M3$quants[,1],ymax = M3$quants[,5]), fill = "#E69F00",alpha = 0.25)+
  geom_ribbon(aes(x = -67+M3$x.bin,ymin = M3$quants[,2],ymax = M3$quants[,4]), fill = "#E69F00", alpha = 0.5)+
  geom_line(aes(x = -67+M3$x.bin, y = M3$quants[,3] ), color = "#E69F00",size = 1)+
  scale_x_reverse()+
  ylab("")+
  xlim(c(3300,-100))+
  ylim(c(0,3))+
  coord_flip()+
  theme_light()+
  theme(aspect.ratio = 1.5,axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
        
b.r <- as.vector(bacon.rates$quants) #[,5]-bacon.rates$quants[,1]
AS.O1.r <- as.vector(M1$quants) #[,5]-M1$quants[,1]
AS.O2.r <- as.vector(M2$quants) #[,5]-M2$quants[,1]
AS.O3.r <- as.vector(M3$quants) #[,5]-M3$quants[,1]

save.image(file = "C:/Users/steph/Documents/varveR_steph/multiple_cores_observers_workstation_1.RData")

# Now the script switches to producing the files necessary to run the gibbs model for each observer.

library("devtools")
library("ggplot2")
library(geoChronR)
library("dplyr")
library(plyr)
library("lipdR")
library(matrixStats)
library(purrr)
library(tictoc)

# Prepare data used by all observers
# Note this calculates the bacon model using the IntCal20 curve.

# Load radiocarbon and lead data contained in the LiPD file
D <- readLipd("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Columbine.Arcusa.2020.lpd")

# Run the BACON model
D <- runBacon(D,cutoff = 1e-10, lab.id.var = NULL, age.14c.var = "age14c",               age.14c.uncertainty.var = "age14cuncertainty", 
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

# Observer data (same as before)

# Observer 1

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Charlotte_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Charlotte_counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

col172_sequence <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence <- calibrateCummulativeDepth(CombSeq = col172_sequence, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence <- calibrateCummulativeDepth(CombSeq = col173_sequence, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

# Use this to check if the algorithm works
objectiveFunction(0,0,0,0,0,0,C14 = C14,col172_sequence, col173_sequence)

nIt <- 1e6
paramChain = matrix(NA,nrow = nIt,6)
objChain = matrix(NA,nrow = nIt,6)
paramChain[1,] = c(0.05,0.05,0.1,0.1,0.3,0.3) # was c(0.05,0.05,0.1,0.2,0.3,0.6)

objChain[1,1] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3], paramChain[1,4],paramChain[1,5],paramChain[1,6]),C14,col172_sequence,col173_sequence)
oldProb <- exp(objChain[1,1])
objChain[1,] <- objChain[1,1]

save.image(file = "C:/Users/steph/Documents/varveR_steph/Gibbs_observer1.RData")

# Observer 2

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Sela counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Sela counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

col172_sequence <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence <- calibrateCummulativeDepth(CombSeq = col172_sequence, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence <- calibrateCummulativeDepth(CombSeq = col173_sequence, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

# Use this to check if the algorithm works
objectiveFunction(0,0,0,0,0,0,C14 = C14,col172_sequence, col173_sequence)
logObjFun(param = c(0.01,0.01,0,0.1,0,0.25),C14,col172_sequence, col173_sequence)


nIt <- 1e6
paramChain = matrix(NA,nrow = nIt,6)
objChain = matrix(NA,nrow = nIt,6)
paramChain[1,] = c(0.05,0.05,0.1,0.2,0.3,0.6)

objChain[1,1] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3], paramChain[1,4],paramChain[1,5],paramChain[1,6]),C14,col172_sequence,col173_sequence)
oldProb <- exp(objChain[1,1])
objChain[1,] <- objChain[1,1]

save.image(file = "C:/Users/steph/Documents/varveR_steph/Gibbs_observer2.RData")

# Observer 3

col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-3/",varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

col172_sequence <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence <- calibrateCummulativeDepth(CombSeq = col172_sequence, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence <- calibrateCummulativeDepth(CombSeq = col173_sequence, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

# Use this to check if the algorithm works
objectiveFunction(0,0,0,0,0,0,C14 = C14,col172_sequence, col173_sequence)
logObjFun(param = c(0.01,0.01,0,0.1,0,0.25),C14,col172_sequence, col173_sequence)


nIt <- 1e6
paramChain = matrix(NA,nrow = nIt,6)
objChain = matrix(NA,nrow = nIt,6)
paramChain[1,] = c(0.05,0.05,0.1,0.2,0.3,0.6)

objChain[1,1] <-logObjFun(param = c(paramChain[1,1],paramChain[1,2],paramChain[1,3], paramChain[1,4],paramChain[1,5],paramChain[1,6]),C14,col172_sequence,col173_sequence)
oldProb <- exp(objChain[1,1])
objChain[1,] <- objChain[1,1]

save.image(file = "C:/Users/steph/Documents/varveR_steph/Gibbs_observer3.RData")



### END of script ##
# Now run Gibbs_slurm script 
# Once completed use script_2
