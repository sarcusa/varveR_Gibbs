###############
# Script for Columbine Paper. Plots age-depth models.

lpdfile = "" # Change to point to LipD file. LipD file contains the data in a reuseable format. Available at https://doi.org/10.6084/m9.figshare.14417999.v1 
figdir = "" # Folder location to save figures of age depth models
CRSfile = "" # Point to CRS lead model file (txt), datasets available at https://doi.org/10.6084/m9.figshare.17156702
CFCSfile = "" # Point to CFCS lead model file (txt),datasets available at https://doi.org/10.6084/m9.figshare.17156702
CICfile = "" # Point to CIC lead model file (txt), datasets available at https://doi.org/10.6084/m9.figshare.17156702
Csfile = "" # Point to raw cesium data (csv) available at https://doi.org/10.6084/m9.figshare.17157245

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

# Read lipd file. Lipd file must contain up to date chrondata table

C <- readLipd(lpdfile)

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

ggsave(filename = paste0("COL_bacon_", Sys.Date(), ".pdf"), plot = Bacon.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_bacon_", Sys.Date(), ".png"), plot = Bacon.plot, device = "png", path = figdir)

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

ggsave(filename = paste0("COL_BaconPlum_", Sys.Date(), ".pdf"), plot = BaconPlum.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_BaconPlum_", Sys.Date(), ".png"), plot = BaconPlum.plot, device = "png", path = figdir)

# Combining both plots

plot.with.inset <-
  ggdraw() +
  draw_plot(Bacon.plot) +
  draw_plot(BaconPlum.plot, x = 0.42, y = .097, width = .45, height = .45)
plot(plot.with.inset)

ggsave(filename = paste0("COL_BaconPlum_inset_", Sys.Date(), ".pdf"), plot = plot.with.inset, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_BaconPlum_inset_", Sys.Date(), ".png"), plot = plot.with.inset, device = "png", path = figdir)

# To extract the plum model only
# Extract Lead data from Plum model

D <- readLipd(lpdfile)

#Plum model output available at https://doi.org/10.6084/m9.figshare.14417999.v1 with file name COL17-3_66.out

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

ggsave(filename = paste0("COL_plum_", Sys.Date(), ".pdf"), plot = Plum.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_plum_", Sys.Date(), ".png"), plot = Plum.plot, device = "png", path = figdir)

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

ggsave(filename = paste0("COL_plumVarve_", Sys.Date(), ".pdf"), plot = Plum.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_plumVarve_", Sys.Date(), ".png"), plot = Plum.plot, device = "png", path = figdir)

## From the serac package, calculate lead models, datasets available at https://doi.org/10.6084/m9.figshare.17156702

CRS <- read.table(file = CRSfile, header = T, col.names = c("Depth_avg","BestAD","MinAD","MaxAD","SAR","SAR_err"))

CIC <- read.table(file = CICfile, header = T, col.names = c("Depth_avg","BestAD","MinAD","MaxAD","SAR","SAR_err"))

CFCS <- read.table(file = CFCSfile, header = T, col.names = c("Depth_avg","BestAD","MinAD","MaxAD","SAR","SAR_err"))

all.lead <- rbind(CRS,CFCS) #rbind(CRS,CIC,CFCS)
all.lead$Model <- c(rep("CRS",nrow(CRS)), rep("CFCS",nrow(CFCS))) #rep("CIC",nrow(CIC))
all.lead$MinBP <- 1950-all.lead$MinAD
all.lead$MaxBP <- 1950-all.lead$MaxAD
all.lead$BestBP <- 1950-all.lead$BestAD
rownames(all.lead) <- NULL

Pb.plots <- ggplot(all.lead)+
  geom_ribbon(aes(x = Depth_avg,ymin = MinBP,ymax = MaxBP, fill = Model), alpha = 0.25)+
  geom_line(aes(x = Depth_avg, y = BestBP, color = Model))+
  scale_x_reverse()+
  scale_y_reverse()+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086", "black", "black", "black"))+
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086", "black", "black", "black"))+
  coord_flip(xlim = c(125,0), ylim = c(50,-70) , expand = c(0,0))+
  theme_light()+
  xlab(label = "year BP")+
  ylab(label = "Depth (mm)")
Pb.plots

ggsave(filename = paste0("COL_Pb_", Sys.Date(), ".pdf"), plot = Pb.plots, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_Pb_", Sys.Date(), ".png"), plot = Pb.plots, device = "png", path = figdir)

combine_models <- data.frame(Model_Obs = c(AD.data$Obs, all.lead$Model), Depth = c(AD.data$Depth,all.lead$Depth_avg), MinBP = c(AD.data$c2.5-67, all.lead$MinBP), MaxBP = c(AD.data$c97.5-67, all.lead$MaxBP), BestBP = c(AD.data$c50-67, all.lead$BestBP))

VarvesLeadPlot <- ggplot(combine_models)+
  geom_ribbon(aes(y = Depth,xmin = MinBP,xmax = MaxBP, fill = Model_Obs), alpha = 0.25)+
  geom_line(aes(y = Depth, x = BestBP, color = Model_Obs))+
  geom_point(aes(x = -13, y = 32.5))+
  geom_errorbar(aes(ymin = 37.5, ymax = 22.5, x = -13))+
  scale_x_reverse()+
  scale_y_reverse()+
  coord_cartesian(xlim = c(50,-70), ylim = c(125,0))+
  theme_light()+
  ylab("Depth (mm)")+
  scale_color_manual(values = c("#7fc97f","#beaed4","#999999", "#56B4E9", "#E69F00"))+
  scale_fill_manual(values = c("#7fc97f","#beaed4","#999999", "#56B4E9", "#E69F00"))
VarvesLeadPlot

ggsave(filename = paste0("COL_Varves_Pb_", Sys.Date(), ".pdf"), plot = VarvesLeadPlot, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_Varves_Pb_", Sys.Date(), ".png"), plot = VarvesLeadPlot, device = "png", path = figdir)

plum.lead.plot <- Plum.plot +
  geom_ribbon(aes(y = all.lead$Depth_avg/10,xmin = all.lead$MinBP,xmax = all.lead$MaxBP, fill = all.lead$Model), alpha = 0.25)+
  geom_line(aes(y = all.lead$Depth_avg/10, x = all.lead$BestBP, color = all.lead$Model))+
  geom_point(aes(x = -13, y = 3.25))+
  geom_errorbar(aes(ymin = 3.75, ymax = 2.25, x = -13))+
  #scale_x_reverse()+
  #scale_y_reverse()+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086", "grey10", "grey10", "grey10", "grey10", "grey10", "grey10"))+
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086", "black", "black", "black"), name = "Model or Observer")
  #coord_flip(xlim = c(125,0), ylim = c(50,-70) , expand = c(0,0))+
  #theme_light()+
  #xlab(label = "year BP")+
  #ylab(label = "Depth (mm)")
plum.lead.plot

ggsave(filename = paste0("COL_Plum_lead_", Sys.Date(), ".pdf"), plot = plum.lead.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_Plum_lead_", Sys.Date(), ".png"), plot = plum.lead.plot, device = "png", path = figdir)


all.age_models <- Plum.plot +
  geom_ribbon(aes(y = combine_models$Depth/10,xmin = combine_models$MinBP,xmax = combine_models$MaxBP, fill = combine_models$Model_Obs), alpha = 0.25)+
  geom_line(aes(y = combine_models$Depth/10, x = combine_models$BestBP, color = combine_models$Model_Obs))+
  geom_point(aes(x = -13, y = 3.25))+
  geom_errorbar(aes(ymin = 3.75, ymax = 2.25, x = -13))+
  scale_x_reverse()+
  scale_y_reverse()+
  #coord_cartesian(xlim = c(50,-70), ylim = c(125,0))+
  theme_light()+
  #coord_flip(xlim = c(50,-70), ylim = c(125,0))+
  ylab("Depth (mm)")+
  scale_color_manual(values = c("#7fc97f","#beaed4", "#999999", "#56B4E9", "#E69F00", "red", "red", "red", "red", "red", "red"), guide = F)+
  scale_fill_manual(values = c("#7fc97f","#beaed4","#999999", "#56B4E9", "#E69F00"), name = "Model or Observer")+
  theme(legend.position = "bottom")+
  theme(aspect.ratio = 1)
all.age_models

ggsave(filename = paste0("COL_Plum_varve_lead_", Sys.Date(), ".pdf"), plot = all.age_models, device = "pdf", path = figdir)
ggsave(filename = paste0("COL_Plum_varve_lead_", Sys.Date(), ".png"), plot = all.age_models, device = "png", path = figdir)

# Plot raw Cs and Pb 

raw.Pb <- read.csv(file = Csfile, header = T)

raw.Pb.m <- reshape2::melt(raw.Pb)

lead.raw.plot <- ggplot(raw.Pb)+
  geom_point(aes(x = Pb210, y = Depth))+
  scale_y_reverse()+
  theme_light()+
  geom_errorbarh(aes(y = Depth,xmin=Pb210-Pb210E, xmax=Pb210+Pb210E))+
  geom_point(aes(x = Pb214, y = Depth), shape = 0)+
  geom_errorbarh(aes(y = Depth,xmin=Pb214-Pb214E, xmax=Pb214+Pb214E))+
  xlab("214Pb and 210Pb activity 
       [Bq kg-1]")+
  ylab("Depth (mm)")+
  theme(aspect.ratio = 3,
        plot.margin = unit(c(1,0,0.5,0), "cm"))

cesium.raw.plot <- ggplot(raw.Pb)+
  geom_point(aes(x = Cs137, y = Depth))+
  scale_y_reverse()+
  theme_light()+
  geom_errorbarh(aes(y = Depth,xmin=Cs137-Cs137E, xmax=Cs137+Cs137E))+
  ylab("")+
  xlab("137Cs activity 
       [Bq kg-1]")+
  theme(aspect.ratio = 3,
        axis.text.y = element_blank(),
        plot.margin = unit(c(1,0,0.5,0), "cm"))
cesium.raw.plot

raw.plot <- lead.raw.plot + cesium.raw.plot + plot_layout(nrow = 1)
raw.plot

ggsave(filename = paste0("raw_lead_", Sys.Date(), ".pdf"), plot = raw.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("raw_lead_", Sys.Date(), ".png"), plot = raw.plot, device = "png", path = figdir)

# Final figure

final.plota <- cowplot::plot_grid(lead.raw.plot,cesium.raw.plot, plot.with.inset, rel_widths = c(1,1,2), nrow = 1, labels = c("A","B","C"))
final.plota
final.plotb <- cowplot::plot_grid(all.age_models,compareVarveBacon, rel_widths = c(1,1), nrow = 1, labels = c("E","F"))
final.plotb

final.chronology.plot <- cowplot::plot_grid(final.plota ,final.plotb, nrow = 2)

ggsave(filename = paste0("final_chrono_plot_", Sys.Date(), ".pdf"), plot = final.chronology.plot, device = "pdf", path = figdir)
ggsave(filename = paste0("final_chrono_plot_", Sys.Date(), ".png"), plot = final.chronology.plot, device = "png", path = figdir)

#### Now switch to multiple-cores-observer script
#####


## Extra analyzes. Comment out.


ind.AD <- convertBP2AD(core$chronData[[1]]$model[[2]]$ensembleTable[[1]]$ageEnsemble$values)
ind.depth <- core$chronData[[1]]$model[[2]]$ensembleTable[[1]]$depth$values
ind.max.age <- colMins(ind.AD)

# Add a table to the LiPD 

core$paleoData[[1]]$measurementTable[[1]]$depth$variableName <- "Depth"
core$paleoData[[1]]$measurementTable[[1]]$depth$units <- "cm"
core$paleoData[[1]]$measurementTable[[1]]$depth$values <- TP_depth

core = mapAgeEnsembleToPaleoData(core) # choose 2 for chronData model, and 1 for paleodata measurement table
Depth_forInd <- selectData(core, varName = "depth", where = "paleoData", model.num = 1, tableType = "ensemble") # choose 1
Age_ens_for_Ind <- selectData(L = core, varName = "ageensemble", where = "paleoData") # choose
ageModel_Ind <- Age_ens_for_Ind$values+68 # To convert to the same age scale as the varves, will add 68 so the top is 1

PB_distribution <- core$chronData[[1]]$model[[2]]$distributionTable
PB_probdens <- purrr::flatten(lapply(PB_distribution, function (x) x[c('probabilityDensity')]))
PB_dens <- purrr::flatten(lapply(PB_probdens, function (x) x[c('values')]))
PB_age <- purrr::flatten(lapply(PB_distribution, function (x) x[c('age')]))
PB_age_values <- purrr::flatten(lapply(PB_age, function (x) x[c('values')]))
depths_14C <- as.numeric(unlist(purrr::flatten(lapply(PB_distribution, function (x) x[c('depth')]))))
names(PB_dens) <- paste("depth_",depths_14C, sep = "")
names(PB_age_values) <- paste("depth_",depths_14C, sep = "")

### VarveR model

#step 1, load in shape files
col172 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-2/",varveTop = "left")
col173 <- readVarveDirectory("D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Steph_counts/COL17-3/",varveTop = "left")

#order by marker layers
o2 <- determineMarkerLayerOrder(col172)
names(col172)[o2]
o3 <- determineMarkerLayerOrder(col173)
names(col173)[o3]

#plot each section
plotSections(col172,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-2.pdf")
plotSections(col173,output = "D:/OneDrive for Business/OneDrive - Northern Arizona University/PhD thesis/Varves/Results/COL17-3.pdf")

# combine sections into sequence
col172_sequence <- combineSectionsByMarkerLayer(col172[o2])
col173_sequence <- combineSectionsByMarkerLayer(col173[o3])

#make a colorful plot
plotCompositeSequence(col172_sequence)
plotCompositeSequence(col173_sequence)

#overcounting and uncercounting probs
unique(col172_sequence$varveCode)
unique(col173_sequence$varveCode)
#transtable <- data.frame(codes = 1:5,overcount = c(0.05,0.10,0.30,.01,0), undercount =  c(0.05,0.10,0.30,.50,.90))
transtable <- data.frame(codes = 1:6,overcount = c(0.05,0.10,0.30,0,-2,-3), undercount =  c(0.05,0.10,0.30,.50,-2,-3))
col172_sequence <- translateCodesToProbabilities(col172_sequence,translationTable = transtable,baselineProb = .05)
col173_sequence <- translateCodesToProbabilities(col173_sequence,translationTable = transtable,baselineProb = .05)

#generate thickness ensemble
col172Ens <- generateThicknessEnsemble(col172_sequence,nEns = 100)
col173Ens <- generateThicknessEnsemble(col173_sequence,nEns = 100)

plotTimeseriesEnsLines(X = 2018-seq_len(nrow(col172Ens$ensThick)),Y = col172Ens$ensThick,maxPlotN = 2,color = "black") %>%
  plotTimeseriesEnsLines(X = 2018-seq_len(nrow(col173Ens$ensThick)),Y = col173Ens$ensThick,maxPlotN = 2,color = "red")+
  ylim(c(0,4000)) + xlim(c(1000,2017))

ggplot()+geom_line(aes(x = 2018-seq_len(nrow(col172Ens$ensThick)),y = col172Ens$ensThick[,2]))+
  geom_line(aes(x = 2018-seq_len(nrow(col173Ens$ensThick)),y = col173Ens$ensThick[,2]+2000),color = "red")+theme_bw()

ensList <- list(col172Ens,col173Ens)

allTiePoints <- apply(cbind(na.omit(unique(col172_sequence$tiePoint)),na.omit(unique(col173_sequence$tiePoint))),1,unique)

allMarkerLayers <- allTiePoints
#allMarkerLayers = as.character(na.omit(unique(c(purrr::flatten(ensList)$tiePoints))))
allMarkerLayers <- c("topOfCore",allMarkerLayers)
#allTiePoints <- c("topOfCore",apply(cbind(na.omit(unique(col172_sequence$tiePoint)),na.omit(unique(col173_sequence$tiePoint))),1,unique))

for(i in 1:100){
  print(tail(unique(ensList[[1]]$tiePoints[,i]),1))
}

n = 10

modeledVarves <- varveR:::varveModel(ensList, n, allMarkerLayers)


ggplot() +
  geom_line(aes(x = 2018-seq_len(length(modeledVarves[[1]][,2])), y = modeledVarves[[1]][,2]))+
  geom_line(aes(x = 2018-seq_len(length(modeledVarves[[1]][,1])), y = modeledVarves[[1]][,1]+25),color = "red")+
  theme_bw()

# Creating varve age depth model

# Converting arbitrary depth to actual depth

modeledVarves[[1]] <- modeledVarves[[1]]+1 #making all values positive by adding arbitrary number

lastTP <- modeledVarves[[3]][which(allMarkerLayers == tail(allMarkerLayers, 1)),] # Find the year/index for last tiepoint, if the core doesn't reach the last TP then will be NA and be disqualified from the analysis (for now)

modeledVarves[[4]] <- matrix(NA, nrow = nrow(modeledVarves[[1]]), ncol = n)
OldMin = 0
NewMax = 124.5
NewMin = 0

for(j in 1:n){
  
  if(is.na(lastTP[j])){ 
    next
  }
  for(i in 1:lastTP[j]){
    
    OldMax = max(cumsum(modeledVarves[[1]][,j]), na.rm = T)
    OldValue = modeledVarves[[1]][i,j]
    
    OldRange = (OldMax - OldMin)  
    NewRange = (NewMax - NewMin)  
    NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
    
    modeledVarves[[4]][i,j] <- NewValue
    
  }
  
}
  
# Add matrix with the cumulatitve depth

modeledVarves[[5]] <- matrix(NA, nrow = nrow(modeledVarves[[1]]), ncol = n)
for(i in 1:n){
  
  modeledVarves[[5]][,i] <- cumsum(modeledVarves[[4]][,i])
  
}

# Add matrix with age

modeledVarves[[6]] <- matrix(NA, nrow = nrow(modeledVarves[[1]]), ncol = n)
for(i in 1:n){
  
  if(is.na(lastTP[i])){
    next
  }
  
  modeledVarves[[6]][,i] <- c(seq(1,lastTP[i],1), rep(NA,nrow(modeledVarves[[1]])-length(seq(1,lastTP[i],1))))
  
}
names(modeledVarves) <- c("ModelOutput", "PostInfo","MLpostAge","VarveThick", "CumDepth", "Age")

# Age depth model from varves

modeledVarves$VarveThick <- modeledVarves$VarveThick[rowSums(is.na(modeledVarves$VarveThick)) != ncol(modeledVarves$VarveThick),colSums(is.na(modeledVarves$VarveThick))<nrow(modeledVarves$VarveThick)] # this removes rows and columns where ALL are NAs

count <- seq_len(nrow(modeledVarves$VarveThick))

geoChronR::plotTimeseriesEnsRibbons(X = count, Y = modeledVarves$VarveThick) %>%
  geoChronR::plotTimeseriesEnsLines(X = count, Y = modeledVarves$VarveThick, maxPlotN = 2,color = "red")

ageEns <- createEnsembleAgeDepthModel(modeledVarves$VarveThick) #calculate an age/depth matrix

ageModel_OriV <- createDepth2AgeFunction(ageDepthMat = ageEns$ageDepthEns)
someDepthSequence <- ind.depth # Maybe choose the same values as the 210Pb and 14C? Or Same depth as independent age model
varveAges <- ageModel_OriV(someDepthSequence)

plot(varveAges$depth,varveAges$medianAge)
hist(varveAges$ageEnsemble[5,]) # to check the spread in the selected age

# Plot the varve model alone

plotTable <- as.data.frame(ageEns$summaryTable)

names(plotTable) <- c("depth", "c2.5", "c25", "c50", "c75", "c97.5")

VarvePlot <- ggplot(plotTable)+
  geom_ribbon(aes(x = depth,ymin = c2.5  ,ymax = c97.5),fill = "gray80")+
  geom_ribbon(aes(x = depth,ymin = c25  ,ymax = c75), fill = "gray50")+ 
  geom_line(aes(x = depth, y = c50))+
  scale_x_reverse()+
  theme_bw()+
  coord_flip()


geoChronR::plotTimeseriesEnsRibbons(X = 2018-count, Y = modeledVarves$VarveThick, x.bin = 1:2017) +
  xlim(c(1981,2017))+ 
  ylim(c(0,0.25))+
  xlab("Year AD") + 
  ylab("Varve thickness") +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, 
                                   hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =20))


# Plot Varve and independent model

chron.plot + 
  geom_ribbon(aes(y = plotTable$depth, xmin = -67+plotTable$c2.5, xmax = -67+plotTable$c97.5), fill = "gray80")+
  geom_ribbon(aes(y = plotTable$depth, xmin = -67+plotTable$c25, xmax = -67+plotTable$c75), fill = "gray50")+
  geom_line(aes(x = -67+plotTable$c50, y = plotTable$depth))

test <- quantile2d(x = varveAges$ageEnsemble, y = varveAges$depth)
  
chron.plot +
  geom_line(aes(y = varveAges$depth, x = -67+varveAges$medianAge ))
  
#### Other stuff 

#What does the following do? Is it necessary?
cv <- ageEns$ageDepthEns[,-1]
endRow <- cv[nrow(cv),]
centEns <- which(near(endRow,median(endRow),tol = .6))
centEns <- sample(centEns,size = 1)

#plot with central
geoChronR::plotTimeseriesEnsRibbons(X = 2018-count, Y = modeledVarves$VarveThick) +
  geom_line(aes(x = 2018-count, y = varveEnsemble[,centEns], color = "Red"), size = .2)


##########
## Testing for the varve model that most resembles the independent model

# Setting the info

n = 10 # number of simulations

id1 <- seq(0,0.1,0.01) # the opposite would be 0.1-id1 i.e. the reverse
id2 <- seq(0.0,0.2,0.02)
id3 <- seq(0.0,0.6,0.05)

maxage <- list() # list of vectors of max ages for each simulation
settings <- list() # list of priors used in the simulation


# Model test

#pb <- txtProgressBar(min = 1, max = length(id1), style = 3)
counter = 0
for(i in 1:length(id1)){
  for(j in 1:length(id2)){
    for(k in 1:length(id3)){
      
      counter = counter + 1
      print(c("Starting simulation",i,j,k))
      # iterating through the possible probabilities
      ttable <- data.frame(codes = 1:6,
                           overcount = c(id1[i],id2[j],id3[k],0,-2,-3), 
                           undercount = c(rev(id1)[i],rev(id2)[j],rev(id3)[k],.50,-2,-3))
      
      col172_sequence <- translateCodesToProbabilities(col172_sequence,translationTable = ttable,
                                                       baselineProb = .05)
      col173_sequence <- translateCodesToProbabilities(col173_sequence,translationTable = ttable,
                                                       baselineProb = .05)
      
      col172Ens <- generateThicknessEnsemble(col172_sequence,nEns = 100)
      col173Ens <- generateThicknessEnsemble(col173_sequence,nEns = 100)
      
      
      ensList <- list(col172Ens,col173Ens)
      allTiePoints <- apply(cbind(na.omit(unique(col172_sequence$tiePoint)),
                                  na.omit(unique(col173_sequence$tiePoint))),1,unique)
      allMarkerLayers <- allTiePoints
      allMarkerLayers <- c("topOfCore",allMarkerLayers)
      
      
      BestModeledVarves <- varveR:::varveModel(ensList, n, allMarkerLayers)
      
  
      # saving all max ages for each simulation but removing NAs
      maxage[[counter]] <- tail(BestModeledVarves[[3]],1)[!is.na(tail(BestModeledVarves[[3]],1))] 
        
      # saving all the settings for each simulation
      settings[[counter]] <- c(ttable$overcount,ttable$undercount)
      
      print(c("Completed",i,j,k))
      #setTxtProgressBar(pb, i)
    }
  }
}

# Analysis of test results

hist(unlist(maxage)) # how to convert to the same timescale? Here removing 67 years so the start year is 2

n.obs <- sapply(maxage, length)
seq.max <- seq_len(max(n.obs))
maxage_mat <- t(sapply(maxage, "[", i = seq.max)) # Each simulation is a row

mean_age <- lapply(maxage, function(x) mean(x))
hist(unlist(mean_age))
mean_age_vector <- unlist(mean_age)
mean_age_vector[which.max(mean_age_vector)] # finding which simulation has the oldest mean age

n.obs <- sapply(settings, length)
seq.max <- seq_len(max(n.obs))
settings_mat <- t(sapply(settings, "[", i = seq.max)) # Each simulation is a row

settings_mat[which.max(mean_age_vector),] # Finding the settings of that simulation

####### Re-running the model with the best settings

Besttable <- data.frame(codes = 1:6,overcount = c(0.05,0,0,0,-2,-3), undercount =  c(0.05,0.20,0.60,.50,-2,-3))
col172_sequence <- translateCodesToProbabilities(col172_sequence,translationTable = Besttable,baselineProb = .05)
col173_sequence <- translateCodesToProbabilities(col173_sequence,translationTable = Besttable,baselineProb = .05)
col172Ens <- generateThicknessEnsemble(col172_sequence,nEns = 100)
col173Ens <- generateThicknessEnsemble(col173_sequence,nEns = 100)
ensList <- list(col172Ens,col173Ens)
allTiePoints <- apply(cbind(na.omit(unique(col172_sequence$tiePoint)),na.omit(unique(col173_sequence$tiePoint))),1,unique)
allMarkerLayers <- allTiePoints
allMarkerLayers <- c("topOfCore",allMarkerLayers)

n = 10

AdjustedVarves <- varveR:::varveModel(ensList, n, allMarkerLayers)

AdjustedVarves[[1]] <- AdjustedVarves[[1]]+1 #making all values positive by adding arbitrary number

lastTP <- AdjustedVarves[[3]][which(allMarkerLayers == tail(allMarkerLayers, 1)),] # Find the year/index for last tiepoint, if the core doesn't reach the last TP then will be NA and be disqualified from the analysis (for now)

AdjustedVarves[[4]] <- matrix(NA, nrow = nrow(AdjustedVarves[[1]]), ncol = n)
OldMin = 0
NewMax = 124.5
NewMin = 0

for(j in 1:n){
  
  if(is.na(lastTP[j])){ 
    next
  }
  for(i in 1:lastTP[j]){
    
    OldMax = max(cumsum(AdjustedVarves[[1]][,j]), na.rm = T)
    OldValue = AdjustedVarves[[1]][i,j]
    
    OldRange = (OldMax - OldMin)  
    NewRange = (NewMax - NewMin)  
    NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
    
    AdjustedVarves[[4]][i,j] <- NewValue
    
  }
  
}

# Add matrix with the cumulatitve depth

AdjustedVarves[[5]] <- matrix(NA, nrow = nrow(AdjustedVarves[[1]]), ncol = n)
for(i in 1:n){
  
  AdjustedVarves[[5]][,i] <- cumsum(AdjustedVarves[[4]][,i])
  
}

# Add matrix with age

AdjustedVarves[[6]] <- matrix(NA, nrow = nrow(AdjustedVarves[[1]]), ncol = n)
for(i in 1:n){
  
  if(is.na(lastTP[i])){
    next
  }
  
  AdjustedVarves[[6]][,i] <- c(seq(1,lastTP[i],1), rep(NA,nrow(AdjustedVarves[[1]])-length(seq(1,lastTP[i],1))))
  
}
names(AdjustedVarves) <- c("ModelOutput", "PostInfo","MLpostAge","VarveThick", "CumDepth", "Age")

# Age depth model from varves

AdjustedVarves$VarveThick <- AdjustedVarves$VarveThick[rowSums(is.na(AdjustedVarves$VarveThick)) != ncol(AdjustedVarves$VarveThick),colSums(is.na(AdjustedVarves$VarveThick))<nrow(AdjustedVarves$VarveThick)] # this removes rows and columns where ALL are NAs

count <- seq_len(nrow(AdjustedVarves$VarveThick))

geoChronR::plotTimeseriesEnsRibbons(X = count, Y = AdjustedVarves$VarveThick) %>%
  geoChronR::plotTimeseriesEnsLines(X = count, Y = AdjustedVarves$VarveThick, maxPlotN = 2,color = "red")

ageEns <- createEnsembleAgeDepthModel(AdjustedVarves$VarveThick) #calculate an age/depth matrix

ageModel_AdjV <- createDepth2AgeFunction(ageDepthMat = ageEns$ageDepthEns)
someDepthSequence <- ind.depth # Maybe choose the same values as the 210Pb and 14C? Or Same depth as independent age model
varveAges <- ageModel_AdjV(someDepthSequence)

plot(varveAges$depth,varveAges$medianAge)
hist(varveAges$ageEnsemble[5,]) # to check the spread in the selected age

# Plot the varve model alone
plotTable <- as.data.frame(ageEns$summaryTable)

names(plotTable) <- c("depth", "c2.5", "c25", "c50", "c75", "c97.5")

# Plot Varve and independent model

chron.plot + 
  geom_ribbon(aes(y = plotTable$depth, xmin = -67+plotTable$c2.5, xmax = -67+plotTable$c97.5), fill = "gray80")+
  geom_ribbon(aes(y = plotTable$depth, xmin = -67+plotTable$c25, xmax = -67+plotTable$c75), fill = "gray50")+
  geom_line(aes(x = -67+plotTable$c50, y = plotTable$depth))

##### Age of tie points

TPAges_OriV <- ageModel_OriV(TP_depth) # From the original varve model
TPAges_AdjV <- ageModel_AdjV(TP_depth) # From the adjusted varve model
#TPAges_Ind <- ageModel_Ind(TP_depth) # From the independent age depth model (PLUM&BACON), this doesn't work
TPAges_Ind <- ageModel_Ind[, sample(ncol(ageModel_Ind), n)]

compare <- data.frame(Age = c(as.vector(TPAges_Ind), as.vector(TPAges_OriV$ageEnsemble), as.vector(TPAges_AdjV$ageEnsemble)), Model = c(rep("Independent",length(as.vector(TPAges_Ind))), rep("Original",length(as.vector(TPAges_OriV$ageEnsemble))), rep("Adjusted",length(as.vector(TPAges_AdjV$ageEnsemble)))), TiePoint = c(rep(TP_depth,ncol(TPAges_Ind)+ncol(TPAges_OriV$ageEnsemble)+ncol(TPAges_AdjV$ageEnsemble))))


ggplot(compare, aes(x = Age, fill = Model))+
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 20)+
  facet_wrap(. ~ TiePoint)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.1))

