
#####################
###clear all lists###
#####################

rm(list=ls())

###################
####Load packages
###################
library(dplyr)
library(qgraph)
library(bootnet)
library(mgm)
library(skimr)
library(tidyr)
library(reshape)
library(ggplot2)
library(mgm)
library(optmatch)
library(EstimateGroupNetwork)
library(MatchIt)
library(corrplot)
library(networktools)
library(readxl)
library(qgraph)
library(EGAnet)


###################
####Read in data 
###################
setwd("")

data <- read_excel("SchizophreniaImputed.xlsx")

dim(data) #3006 obserations of 47 variables 

###################
####Sample characteristics & data exploration:
###################
#This section was omitted from the open access script and datasets to safeguard anonimity of participants. 

###################
####AIM 1: STRUCTURE AND CENTRALITY OF THE SCHIZOPHRENIA SYMPTOM NETWORK MODEL
###################

MyData <- select(data, p1, p2, p3, p4, p5, p6, p7, n1, n2, n3, n4, n5, n6, n7, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16)



#Create labels for the items:
group.item <- list("PANSS-positive"=c(1:7), 
                   "PANSS-negative"=c(8:14), 
                   "PANSS-GeneralPsychopathology"=c(15:30))

longnames <- c("Delusions",
               "Conceptual Disorganization",
               "Hallucinations",
               "Excitement",
               "Grandiosity",
               "Suspiciousness/persecution",
               "Hostility",
               "Blunted Affect",
               "Emotional Withdrawal",
               "Poor Rapport",
               "Passive/Apathetic",
               "Difficulty in Abstract Thinking",
               "Lack of Spontaneity and Flow of Conversation",
               "Stereotyped Thinking",
               "Somatic Concern",
               "Anxiety",
               "Guilt Feelings",
               "Tension",
               "Mannerisms and Posturing",
               "Depression",
               "Motor Retardation",
               "Uncooperativeness",
               "Unusual Thought Content",
               "Disorientation",
               "Poor Attention",
               "Lack of Judgment and Insight",
               "Disturbance of Volition",
               "Poor Impulse Control",
               "Preoccupation",
               "Active Social Avoidance")

nodenames <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G14", "G15", "G16")
myvars <- c("p1","p2", "p3", "p4", "p5", "p6", "p7", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10", "g11", "g12", "g13", "g14", "g15", "g16")
#Note: The myvars-object is used below to specify the part of the datasets used for the analyses (PCL-C items). As such, this object is redundant for 
#      the open access datasets which only consist of the PCL-C items and is only of impact for those running the script on more extensive datasets.


#test the distribution:
hist(MyData[,1])
hist(MyData[,2])
hist(MyData[,3])
hist(MyData[,4])
hist(MyData[,5])
hist(MyData[,6])
hist(MyData[,7])
hist(MyData[,8])
hist(MyData[,9])
hist(MyData[,10])
hist(MyData[,11])
hist(MyData[,12])
hist(MyData[,13])
hist(MyData[,14])
hist(MyData[,15])
hist(MyData[,16])
hist(MyData[,17])
hist(MyData[,18])
hist(MyData[,19])
hist(MyData[,20])
hist(MyData[,21])
hist(MyData[,22])
hist(MyData[,23])
hist(MyData[,24])
hist(MyData[,25])
hist(MyData[,26])
hist(MyData[,27])
hist(MyData[,28])
hist(MyData[,29])
hist(MyData[,30])

library("huge")
MyData_norm <- huge.npn(MyData) #transformation in function of normality
hist(MyData_norm[,1])
hist(MyData_norm[,2])
hist(MyData_norm[,3])
hist(MyData_norm[,4])
hist(MyData_norm[,5])
hist(MyData_norm[,6])
hist(MyData_norm[,7])
hist(MyData_norm[,8])
hist(MyData_norm[,9])
hist(MyData_norm[,10])
hist(MyData_norm[,11])
hist(MyData_norm[,12])
hist(MyData_norm[,13])
hist(MyData_norm[,14])
hist(MyData_norm[,15])
hist(MyData_norm[,16])
hist(MyData_norm[,17])
hist(MyData_norm[,18])
hist(MyData_norm[,19])
hist(MyData_norm[,20])
hist(MyData_norm[,21])
hist(MyData_norm[,22])
hist(MyData_norm[,23])
hist(MyData_norm[,24])
hist(MyData_norm[,25])
hist(MyData_norm[,26])
hist(MyData_norm[,27])
hist(MyData_norm[,28])
hist(MyData_norm[,29])
hist(MyData_norm[,30])

#check for reduntant nodes in the model: 
library(networktools)
gb <- goldbricker(MyData_norm, p = 0.05, method = "hittner2003", threshold = 0.25,
                  corMin = 0.5, progressbar = TRUE)
gb #no redundant nodes

# Obtain correlation matrix:
cor_matrix <- cor_auto(MyData_norm)

# Obtain GGM 
Network <- qgraph(cor_matrix, layout = "spring", 
                  groups=group.item, labels=nodenames, nodeNames = longnames, 
                  graph = "glasso", tuning = 0.5, sampleSize = 3006, 
                  legend.cex = 0.4, palette="pastel", borders=FALSE, theme="colorblind", usePCH=TRUE, minimum=0, details=TRUE, edge.labels=TRUE)

#Corresponding weightmatrix:
x <- getWmat(Network)
table(x)
Network$Edgelist


#Expected influence calculates the expected value of a node's influence on other nodes in the network, based on its own centrality and the centrality of its neighbors. 
#It assumes that a node's influence depends not only on its own characteristics, but also on the characteristics of its neighbors. 
#Nodes with high expected influence are those that are highly central in the network and are surrounded by other highly central nodes, 
#making them more likely to have a strong impact on other nodes in the network.

#Expected influence:
centralityPlot(Network, include = "ExpectedInfluence", orderBy = "ExpectedInfluence")

#Add predictability to the plot
#Pred :
set.seed(148218)
fit_obj1 <- mgm(data=MyData[,myvars],
                type = rep('g', 30),
                level = rep('1', 30),
                lambdSel = 'EBIC',
                ruleReg = 'OR')

pred_obj1 <- predict(object = fit_obj1, 
                     data = MyData[,myvars], 
                     errorCon = 'R2')

print(pred_obj1)
table(pred_obj1$errors)
mean(pred_obj1$errors[,2]) #54.38%

#Plot schizophrenia network model including predictability (with or without labels for items):
R2_network_Spring_withoutlabels <- qgraph(getWmat(Network),
                                          layout = "spring",
                                          pie = as.numeric(as.character(pred_obj1$error[,2])), 
                                          pieColor = rep('#377EB8', 30),
                                          groups=group.item, 
                                          legend.cex = 0.4, palette="pastel", theme="colorblind", minimum=0, details=TRUE)


R2_network_Spring <- qgraph(getWmat(Network),
                            layout = "spring",
                            pie = as.numeric(as.character(pred_obj1$error[,2])), 
                            pieColor = rep('#377EB8', 30),
                            groups=group.item, labels=nodenames, nodeNames = longnames, 
                            legend.cex = 0.4, palette="pastel", theme="colorblind", minimum=0, details=TRUE)


###################
####AIM 2: COMMUNITY DETECTION:
###################

# Walktrap community detection:
comm1<- EGA(MyData[,myvars], plot.EGA = TRUE) #4 dimensions detected:
#'Group 1'=c(p1,p2,p3,p4,p5,p6,p7), 'Group 2'=c(n1,n2,n3,n4,n5,n6,n7), 
#'Group 3'=c(g4,g5,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16) , 'Group 4' = c(g1,g2,g3,g6)

#Create new group item based on the observed communities:
group.item2 <- list("Positive Scale"= c(1,2,3,4,5,6,7), 
                    "Negative Scale"= c(8,9,10,11,12,13,14),
                    "General Psychopathology Scale part 1"= c(18,19,21,22,23,24,25,26,27,28,29,30),
                    "General Psychopathology Scale part 2"= c(15,16,17,20))

#Create network with adjusted communities:                   
Network <- qgraph(cor_matrix, layout = "spring", 
                      groups=group.item2, labels=nodenames, nodeNames = longnames, 
                      graph = "glasso", tuning = 0.5, sampleSize = 3006, 
                      legend.cex = 0.4, palette="pastel", borders=FALSE, theme="colorblind", usePCH=TRUE, minimum=0, details=TRUE)


#The bridge function in the qgraph package calculates bridge centrality, which is a measure of how important a node is for connecting different clusters or communities within the network. 
#Bridges are nodes that connect two or more clusters, and their removal can result in the fragmentation of the network into disconnected components.
#The plot function with the "Bridge Expected Influence (1-step)" argument shows the expected influence of each node on other nodes in the network, taking into account their role as bridges. 
#The z-score and order options allow you to standardize the values and sort them by magnitude, respectively.

# Bridge centrality:
b <- bridge(Network, communities = group.item)
plot(b, include = "Bridge Expected Influence (1-step)", zscore=TRUE, order="value")

#Plot Figures:
dpi=600    #pixels per square inch
jpeg("Figure 1.jpeg", width=10*dpi, height=7*dpi, res=dpi)
plot(Network)
dev.off()

write.csv(x, "weights table.csv")

###############################################
### AIM 3: Evaluation of network properties ###
###############################################

library("bootnet")  
par(mfrow = c(2, 1))
Bootnet_Network <- estimateNetwork(MyData,default="EBICglasso", labels=longnames)

set.seed(148218)
Network_boot1 <- bootnet(Bootnet_Network, nBoots = 1000, nCores = 6)
set.seed(148218)
Network_boot2 <- bootnet(Bootnet_Network, nBoots = 1000, type = "case")

Accuracy <- plot(Network_boot1, labels = F, order = "sample") 
plot(Accuracy)
STR_Stability <- plot(Network_boot2)
plot(STR_Stability)
corStability(Network_boot2)

library(bootnet)


sessionInfo()
