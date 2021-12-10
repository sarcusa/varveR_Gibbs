####
# Script that prepares the data frames of laminations from multiple cores and observers. To be used in conjunction with "Gibbs_sampling_model.R" and "COL_script.R" for now. 

Rsourcedir = "" # Location of the source code files. These include file GibbsrRelatedFunctions.R, ImportFiles.R, plotting.R, simulateOverAndUnderCounting.R, varveModel.R
resultsdir = "" # Folder location for the results

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

file.sources = list.files(Rsourcedir, 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)

sapply(file.sources,source,.GlobalEnv)

n = 100
PlotOp = F

# Obs 1

col172 <- readVarveDirectory(Obs1dir_172,varveTop = "left")
col173 <- readVarveDirectory(Obs1dir_173,varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
plotSections(col172,output = output_obs1_172)
plotSections(col173,output = output_obs1_173)
}

col172_sequence1 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence1 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence1 <- calibrateCummulativeDepth(CombSeq = col172_sequence1, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence1 <- calibrateCummulativeDepth(CombSeq = col173_sequence1, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
plotCompositeSequence(col172_sequence1)
plotCompositeSequence(col173_sequence1)
}

#probabilities based on observations
transtable <- data.frame(codes = 1:6,overcount = c(0.05,0.10,0.30,0.5,-2,-3), undercount =  c(0.05,0.10,0.30,0,-2,-3)) 
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

### Observer 2

col172 <- readVarveDirectory(Obs2dir_172,varveTop = "left")
col173 <- readVarveDirectory(Obs2dir_173,varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
plotSections(col172,output = output_obs2_172)
plotSections(col173,output = output_obs2_173)
}

col172_sequence2 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence2 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence2 <- calibrateCummulativeDepth(CombSeq = col172_sequence2, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence2 <- calibrateCummulativeDepth(CombSeq = col173_sequence2, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
plotCompositeSequence(col172_sequence2)
plotCompositeSequence(col173_sequence2)
}

#probabilities based on observations
transtable <- data.frame(codes = 1:6,overcount = c(0.05,0.10,0.30,0.5,-2,-3), undercount =  c(0.05,0.10,0.30,0,-2,-3)) 
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

### Observer 3

col172 <- readVarveDirectory(Obs3dir_172,varveTop = "left")
col173 <- readVarveDirectory(Obs3dir_173,varveTop = "left")

o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
if(PlotOp == TRUE){
plotSections(col172,output = output_obs3_172)
plotSections(col173,output = output_obs3_173)
}

col172_sequence3 <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence3 <- combineSectionsByMarkerLayer(col173[o3])

col172_sequence3 <- calibrateCummulativeDepth(CombSeq = col172_sequence3, topOfVarvedSequence = 0, botOfVarvedSequence = 1270)
col173_sequence3 <- calibrateCummulativeDepth(CombSeq = col173_sequence3, topOfVarvedSequence = 0, botOfVarvedSequence = 1230)

if(PlotOp == TRUE){
plotCompositeSequence(col172_sequence3)
plotCompositeSequence(col173_sequence3)
}

#probabilities based on observations
transtable <- data.frame(codes = 1:6,overcount = c(0.05,0.10,0.30,0.5,-2,-3), undercount =  c(0.05,0.10,0.30,0,-2,-3)) 
col172_sequence3 <- translateCodesToProbabilities(col172_sequence3,translationTable = transtable,baselineProb = .05)
col173_sequence3 <- translateCodesToProbabilities(col173_sequence3,translationTable = transtable,baselineProb = .05)

# Convert factors to characters
col172_sequence3$tiePoint = as.character(col172_sequence3$tiePoint)
col173_sequence3$tiePoint = as.character(col173_sequence3$tiePoint)

col172Ens3 <- generateThicknessEnsemble(col172_sequence3,nEns = 100)
col173Ens3 <- generateThicknessEnsemble(col173_sequence3,nEns = 100)

ensList3 <- list(col172Ens3,col173Ens3)

allTiePoints <- apply(cbind(na.omit(unique(col172_sequence3$tiePoint)),na.omit(unique(col173_sequence3$tiePoint))),1,unique)
allMarkerLayers <- allTiePoints
allMarkerLayers <- c("topOfCore",allMarkerLayers)

modeledVarves3 <- varveModel(ensList3, nSim = n, allMarkerLayers)

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
  theme(legend.position="none", panel.spacing = unit(2, "lines"))+
  ggtitle('COL17-3')+
  scale_x_reverse(limits = c(2500,0))+ #added this, may not work
  ylab('')

both_counts <- core1 + core2 + plot_layout(nrow = 1)
both_counts

ggsave(filename = paste0(resultsdir,"all_cores_observer_counts.pdf"), plot = both_counts)
ggsave(filename = paste0(resultsdir,"all_cores_observer_counts.png"), plot = both_counts)

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

ggsave(filename = paste0(resultsdir,"all_cores_observer_QC.pdf"), plot = both_QC)
ggsave(filename = paste0(resultsdir,"all_cores_observer_QC.png"), plot = both_QC)

# Compare age-depth models

AD.data <- data.frame(Obs = c(rep('Obs1', nrow(Vplot1$data)), rep('Obs2', nrow(Vplot2$data)), rep('Obs3', nrow(Vplot3$data))), Depth = c( Vplot1$data$depth, Vplot2$data$depth, Vplot3$data$depth), c2.5 = c(Vplot1$data$c2.5,Vplot2$data$c2.5,Vplot3$data$c2.5), c25 = c(Vplot1$data$c25,Vplot2$data$c25,Vplot3$data$c25), c50 = c(Vplot1$data$c50,Vplot2$data$c50,Vplot3$data$c50), c75 = c(Vplot1$data$c75,Vplot2$data$c75,Vplot3$data$c75), c97.5 = c(Vplot1$data$c97.5,Vplot2$data$c97.5,Vplot3$data$c97.5), c2.5BP = c(Vplot1$data$c2.5,Vplot2$data$c2.5,Vplot3$data$c2.5)-67, c25BP = c(Vplot1$data$c25,Vplot2$data$c25,Vplot3$data$c25)-67, c50BP = c(Vplot1$data$c50,Vplot2$data$c50,Vplot3$data$c50)-67, c75BP = c(Vplot1$data$c75,Vplot2$data$c75,Vplot3$data$c75)-67, c97.5BP = c(Vplot1$data$c97.5,Vplot2$data$c97.5,Vplot3$data$c97.5)-67)

output1 <- data.frame(Depth = round(Vplot1$data$depth[1:260],1), Age = round(2018-Vplot1$data$c50[1:260],0))
output2 <- data.frame(Depth = round(Vplot2$data$depth[1:260],1), Age = round(2018-Vplot2$data$c50[1:260],0))
output3 <- data.frame(Depth = round(Vplot3$data$depth[1:260],1), Age = round(2018-Vplot1$data$c50[1:260],0))

write.table(x = output, file = paste0(resultsdir,"Columbine_210Pb_varves.txt"), col.names = T, row.names = F)

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

ggsave(filename = paste0(resultsdir,"all_age_depth_models.pdf"), plot = CompareAgeDepthPlot)
ggsave(filename = paste0(resultsdir,"all_age_depth_models.png"), plot = CompareAgeDepthPlot)

## Compare varveR model to raw counts

col173_sequence1$cumsum <- cumsum(col173_sequence1$thick)
col173_sequence2$cumsum <- cumsum(col173_sequence2$thick)
col173_sequence3$cumsum <- cumsum(col173_sequence3$thick)

Obs1.plot <- AD.data %>%
  filter(Obs == "Obs1") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5), fill = "#999999", alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75), fill = "#999999", alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50), color = "#999999")+
  geom_line(data =col173_sequence1, aes(x = col173_sequence1$cumsum, y = col173_sequence1$count ))+
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

Obs2.plot <- AD.data %>%
  filter(Obs == "Obs2") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5), fill = "#56B4E9", alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75), fill = "#56B4E9", alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50), color = "#56B4E9")+
  geom_line(data =col173_sequence2, aes(x = col173_sequence2$cumsum, y = col173_sequence2$count ))+
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

Obs3.plot <- AD.data %>%
  filter(Obs == "Obs3") %>%
  ggplot()+
  geom_ribbon(aes(x = Depth,ymin = c2.5,ymax = c97.5), fill = "#E69F00", alpha = 0.25)+
  geom_ribbon(aes(x = Depth,ymin = c25,ymax = c75), fill = "#E69F00", alpha = 0.5)+ 
  geom_line(aes(x = Depth, y = c50), color = "#E69F00")+
  geom_line(data =col173_sequence3, aes(x = col173_sequence3$cumsum, y = col173_sequence3$count ))+
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

Obs123.plot <- cowplot::plot_grid(Obs1.plot,Obs2.plot, Obs3.plot,labels = c("A","B","C"), nrow = 1, align = "hv", vjust = 5.5)
Obs123.plot

ggsave(plot = Obs123.plot, paste0(resultsdir,"Obs123.pdf"))

#save.image(file = "C:/Users/steph/Documents/varveR_steph/multiple_cores_observers_workstation.RData")

# Must have Run COL_script for this to work

compareVarveBacon <- Bacon.plot+
  geom_ribbon(aes(y = AD.data$Depth/10, xmin = AD.data$c2.5BP,xmax = AD.data$c97.5BP, fill = AD.data$Obs), alpha = .25)+
  geom_ribbon(aes(y = AD.data$Depth/10, xmin = AD.data$c25BP,xmax = AD.data$c75BP, fill = AD.data$Obs), alpha = .5)+ 
  geom_line(aes(x = AD.data$c50BP, y = AD.data$Depth/10, color = AD.data$Obs))+
  scale_y_reverse()+
  scale_x_reverse()+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00"))+
  scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00","red","red", "red", "red", "red"))

compareVarveBacon

ggsave(filename = paste0(resultsdir,"all_age_depth_models_plus_Bacon.pdf"), plot = compareVarveBacon)
ggsave(filename = paste0(resultsdir,"all_age_depth_models_plus_Bacon.png"), plot = compareVarveBacon)
  
## Statistics

# column bind the core sequences instead, then add a column for the core name. then use dplyr to calculate the stats. That way can exclude ID 4,5,6

stats <- rbind(col172_sequence1, col172_sequence2, col172_sequence3,
               col173_sequence1,col173_sequence2,col173_sequence3)

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
  dplyr::summarise(basalAge = round(length(thick),2),
                   minThick = round(min(thick),2),
                   maxThick = round(max(thick),2),
                   meanThick = round(mean(thick),2),
                   medianThick = round(median(thick),2),
                   sdThick = round(sd(thick),2))

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
 
write.csv(ML,file = paste0(resultsdir,"marker_layer.csv"))

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
  dplyr::summarise(cumulative = sum(thick)) %>%
  group_by(core,.drop = F) %>%
  dplyr::summarise(av.extra = mean(cumulative))

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

#### Missing laminations analysis

pdf(file = paste0(resultsdir,"basal_age_colours.pdf"))
par(mfrow = c(3,1))

hist(apply(modeledVarves1[[2]][[2]]$thicks, FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(1,0,0,0.5), main = "Observer 1", ylim = c(0,40), xlab = "",xlim = c(2000,4500))
hist(apply(modeledVarves1[[2]][[1]]$thicks, FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(0,0,1,0.5), add = T)
hist(apply(modeledVarves1[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(0.1,0.1,0.1,0.5), add = T)
box()

hist(apply(modeledVarves2[[2]][[2]]$thicks, FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(1,0,0,0.5), main = "Observer 2", ylim = c(0,40), xlab = "",xlim = c(2000,4500))
hist(apply(modeledVarves2[[2]][[1]]$thicks, FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(0,0,1,0.5), add = T)
hist(apply(modeledVarves2[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(0.1,0.1,0.1,0.5), add = T)
box()

hist(apply(modeledVarves3[[2]][[2]]$thicks, FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(1,0,0,0.5), main = "Observer 2", ylim = c(0,40), xlab = "",xlim = c(2000,4500))
hist(apply(modeledVarves3[[2]][[1]]$thicks, FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(0,0,1,0.5), add = T)
hist(apply(modeledVarves3[[1]], FUN = function(x) min(which(is.na(x))), MARGIN = 2), col=rgb(0.1,0.1,0.1,0.5), add = T)
box()

dev.off()

# Analysis of varve type shift

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

varve.shift
