#install.packages("geosphere")
library(geosphere)
library(dplyr)
library(ggplot2)
library(gridExtra)


#########################
#BATHYACMAEA
#########################

setwd(choose.dir())
#Calculate geographic distances
b_latlon <- read.csv("Bathyacmaea_latlon.csv", header = FALSE)
str(b_latlon)
b_dist_mat <- distm(b_latlon, fun = distGeo)
summary(b_dist_mat)
#Remove upper triangle of pairwise comparisons
b_dist_mat[upper.tri(b_dist_mat,diag=TRUE)] <- "NA"
#Stack matrix values into one column
b_geo_dist <- as.data.frame(as.table(b_dist_mat))
#Convert values to numeric
b_geo_dist$Freq <- as.numeric(b_geo_dist$Freq)
#c_geo_dist <- na.omit(c_geo_dist)
summary(b_geo_dist)

#Load in depth data
b_depth_mat <- read.csv("Bathyacmaea_DepthDifferences.csv", header = TRUE)
b_depth_mat <- as.matrix(b_depth_mat)
b_depth_mat[upper.tri(b_depth_mat,diag=TRUE)] <- "NA"
#Remove first column of names
b_depth_mat <- b_depth_mat[,-1]
#Stack matrix values into one column
b_depth_dist <- as.data.frame(as.table(b_depth_mat))
#Convert values to numeric
b_depth_dist$Freq <- as.numeric(b_depth_dist$Freq)
summary(b_depth_dist)

#Load in genetic data
b_genetic_mat <- read.csv("Bathyacmaea_GeneticDistances.csv", header = TRUE)
b_genetic_mat <- as.matrix(b_genetic_mat)
#Remove first column of names
b_genetic_mat <- b_genetic_mat[,-1]
#Stack matrix values into one column
b_genetic_dist <- as.data.frame(as.table(b_genetic_mat))
#Convert values to numeric
b_genetic_dist$Freq <- as.numeric(b_genetic_dist$Freq)
summary(b_genetic_dist)

#Create dataframe
bathy <- as.data.frame(b_genetic_dist$Freq)
bathy$depth <- b_depth_dist$Freq
bathy$geo <- b_geo_dist$Freq
#Make sure all the column names are good
colnames(bathy)
colnames(bathy) [1] <- c("gene")
colnames(bathy)
#Remove missing data
bathy.na <- na.omit(bathy)
summary(bathy.na)
#Test
plot(depth ~ geo, data = bathy.na)

#MLR with interaction
bathy.lm <- lm(gene ~ depth * geo, data = bathy.na)
summary(bathy.lm)
#depth: ***
#geo: ***
#interaction: *
#Check model
step(bathy.lm, direction = "backward")
#All stays in

#Z-score depth & geography distances
bathy.na$z_depth <- scale(bathy.na$depth, center = TRUE, scale = TRUE)
bathy.na$z_geo <- scale(bathy.na$geo, center = TRUE, scale = TRUE)
summary(bathy.na)

#Plots
bathy_d <- ggplot(data = bathy.na, mapping = aes(x = depth, y = gene)) +
  geom_point(color = "dodgerblue2", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  labs(title = "Bathyacmaea", x = "Difference in Depth (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

bathy_g <- ggplot(data = bathy.na, mapping = aes(x = geo, y = gene)) +
  geom_point(color = "black", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "black", size = 2, se = FALSE) +
  labs(title = "Bathyacmaea", x = "Geographic Distance (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

bathy <- ggplot(data = bathy.na) +
  geom_point(mapping = aes(x = z_geo, y = gene), color = "black",shape = 21, size = 2, alpha = .66)+
  geom_smooth(mapping = aes(x = z_geo, y = gene), method = "lm", color = "black", size = 2, se = FALSE) +
  geom_point(mapping = aes(x = z_depth, y = gene), color = "dodgerblue2", shape = 21, size = 2, alpha = .66) +
  geom_smooth(mapping = aes(x = z_depth, y = gene), method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  labs(title = "(a) Bathyacmaea", x = "z-distance", y = "Pariwise Genetic Distance") +
  theme(
    plot.title = element_text(face = "italic", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))


######################
#COCCULINA
######################

setwd(choose.dir())
#Calculate geographic distances
c_latlon <- read.csv("Cocculina_latlon.csv", header = FALSE)
str(c_latlon)
c_dist_mat <- distm(c_latlon, fun = distGeo)
summary(c_dist_mat)
#Remove upper triangle of pairwise comparisons
c_dist_mat[upper.tri(c_dist_mat,diag=TRUE)] <- "NA"
#Stack matrix values into one column
c_geo_dist <- as.data.frame(as.table(c_dist_mat))
#Convert values to numeric
c_geo_dist$Freq <- as.numeric(c_geo_dist$Freq)
#c_geo_dist <- na.omit(c_geo_dist)
summary(c_geo_dist)

#Load in depth data
c_depth_mat <- read.csv("Cocculina_DepthDifferences.csv", header = TRUE)
c_depth_mat <- as.matrix(c_depth_mat)
#Remove first column of names
c_depth_mat <- c_depth_mat[,-1]
#Stack matrix values into one column
c_depth_dist <- as.data.frame(as.table(c_depth_mat))
#Convert values to numeric
c_depth_dist$Freq <- as.numeric(c_depth_dist$Freq)
summary(c_depth_dist)

#Load in genetic data
c_genetic_mat <- read.csv("Cocculina_GeneticDistances.csv", header = TRUE)
c_genetic_mat <- as.matrix(c_genetic_mat)
#Remove first column of names
c_genetic_mat <- c_genetic_mat[,-1]
#Stack matrix values into one column
c_genetic_dist <- as.data.frame(as.table(c_genetic_mat))
#Convert values to numeric
c_genetic_dist$Freq <- as.numeric(c_genetic_dist$Freq)
summary(c_genetic_dist)

#Create dataframe
cocculina <- as.data.frame(c_genetic_dist$Freq)
cocculina$depth <- c_depth_dist$Freq
cocculina$geo <- c_geo_dist$Freq
#Make sure all the column names are good
colnames(cocculina)
colnames(cocculina) [1] <- c("gene")
colnames(cocculina)
#Remove missing data
cocculina.na <- na.omit(cocculina)
summary(cocculina.na)
#Test
plot(depth ~ geo, data = cocculina.na)

#MLR with interaction
cocculina.lm <- lm(gene ~ depth * geo, data = cocculina.na)
summary(cocculina.lm)
#depth: ***
#geo: ***
#interaction: ***
#Check model
step(cocculina.lm, direction = "backward")
#All stay in

#Z-score depth & geography distances
cocculina.na$z_depth <- scale(cocculina.na$depth, center = TRUE, scale = TRUE)
cocculina.na$z_geo <- scale(cocculina.na$geo, center = TRUE, scale = TRUE)
summary(cocculina.na)

#Plots
cc_d <- ggplot(data = cocculina.na, mapping = aes(x = depth, y = gene)) +
  geom_point(color = "dodgerblue2", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  labs(title = "Cocculina", x = "Difference in Depth (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

cc_g <- ggplot(data = cocculina.na, mapping = aes(x = geo, y = gene)) +
  geom_point(color = "black", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "black", size = 2, se = FALSE) +
  labs(title = "Cocculina", x = "Geographic Distance (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

cc <-
ggplot(data = cocculina.na) +
  geom_point(mapping = aes(x = z_geo, y = gene), color = "black",shape = 21, size = 2, alpha = .66)+
  geom_smooth(mapping = aes(x = z_geo, y = gene), method = "lm", color = "black", size = 2, se = FALSE) +
  geom_point(mapping = aes(x = z_depth, y = gene), color = "dodgerblue2", shape = 21, size = 2, alpha = .5) +
  geom_smooth(mapping = aes(x = z_depth, y = gene), method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  labs(title = "(b) Cocculina", x = "z-distance", y = "Pairwise Genetic Distance") +
  theme(
    plot.title = element_text(face = "italic", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#########################
#PARALEPETOPSIS
########################

setwd(choose.dir())
#Calculate geographic distances
para_latlon <- read.csv("Paralepetopsis_latlon.csv", header = FALSE)
summary(para_latlon)
para_dist_mat <- distm(para_latlon, fun = distGeo)
summary(para_dist_mat)
#Remove upper triangle of pairwise comparisons
para_dist_mat[upper.tri(para_dist_mat,diag=TRUE)] <- "NA"
#Stack matrix values into one column
para_geo_dist <- as.data.frame(as.table(para_dist_mat))
#Convert values to numeric
para_geo_dist$Freq <- as.numeric(para_geo_dist$Freq)
#c_geo_dist <- na.omit(c_geo_dist)
summary(para_geo_dist)

#Load in depth data
para_depth_mat <- read.csv("Paralepetopsis_DepthDifferences.csv", header = TRUE)
para_depth_mat <- as.matrix(para_depth_mat)
#Remove first column of names
para_depth_mat <- para_depth_mat[,-1]
#Stack matrix values into one column
para_depth_dist <- as.data.frame(as.table(para_depth_mat))
#Convert values to numeric
para_depth_dist$Freq <- as.numeric(para_depth_dist$Freq)
summary(para_depth_dist)

#Load in genetic data
para_genetic_mat <- read.csv("Paralepetopsis_GeneticDistances.csv", header = TRUE)
para_genetic_mat <- as.matrix(para_genetic_mat)
#Remove first column of names
para_genetic_mat <- para_genetic_mat[,-1]
#Stack matrix values into one column
para_genetic_dist <- as.data.frame(as.table(para_genetic_mat))
#Convert values to numeric
para_genetic_dist$Freq <- as.numeric(para_genetic_dist$Freq)
summary(para_genetic_dist)

#Create dataframe
para <- as.data.frame(para_genetic_dist$Freq)
para$depth <- para_depth_dist$Freq
para$geo <- para_geo_dist$Freq
#Make sure all the column names are good
colnames(para)
colnames(para) [1] <- c("gene")
colnames(para)
#Remove missing data
para.na <- na.omit(para)
summary(para.na)
#Test
plot(depth ~ geo, data = para.na)

#MLR with interaction
para.lm <- lm(gene ~ depth * geo, data = para.na)
summary(para.lm)
#Depth: **
#geo: ***
#interaction: NS
#Check model
step(para.lm, direction = "backward")
#take out interaction (depth + geo)

#Z-score depth & geography distances
para.na$z_depth <- scale(para.na$depth, center = TRUE, scale = TRUE)
para.na$z_geo <- scale(para.na$geo, center = TRUE, scale = TRUE)
summary(para.na)

#Plots
para_d <- ggplot(data = para.na, mapping = aes(x = depth, y = gene)) +
  geom_point(color = "dodgerblue2", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  labs(title = "Paralepetopsis", x = "Difference in Depth (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

para_g <- ggplot(data = para.na, mapping = aes(x = geo, y = gene)) +
  geom_point(color = "black", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "black", size = 2, se = FALSE) +
  labs(title = "Paralepetopsis", x = "Geographic Distance (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

para <- 
ggplot(data = para.na) +
  geom_point(mapping = aes(x = z_geo, y = gene), color = "black",shape = 21, size = 2, alpha = .66)+
  geom_smooth(mapping = aes(x = z_geo, y = gene), method = "lm", color = "black", size = 2, se = FALSE) +
  geom_point(mapping = aes(x = z_depth, y = gene), color = "dodgerblue2", shape = 21, size = 2, alpha = .66) +
  geom_smooth(mapping = aes(x = z_depth, y = gene), method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  labs(title = "(c) Paralepetopsis", x = "z-distance", y = "Pairwise Genetic Distance") +
  theme(
    plot.title = element_text(face = "italic", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

#########################
#PYROPELTA
########################

setwd(choose.dir())
#Calculate geographic distances
pyro_latlon <- read.csv("Pyropelta_latlon.csv", header = FALSE)
summary(pyro_latlon)
pyro_dist_mat <- distm(pyro_latlon, fun = distGeo)
summary(pyro_dist_mat)
#Remove upper triangle of pairwise comparisons
pyro_dist_mat[upper.tri(pyro_dist_mat,diag=TRUE)] <- "NA"
#Stack matrix values into one column
pyro_geo_dist <- as.data.frame(as.table(pyro_dist_mat))
#Convert values to numeric
pyro_geo_dist$Freq <- as.numeric(pyro_geo_dist$Freq)
#c_geo_dist <- na.omit(c_geo_dist)
summary(pyro_geo_dist)

#Load in depth data
pyro_depth_mat <- read.csv("Pyropelta_DepthDifferences.csv", header = TRUE)
pyro_depth_mat <- as.matrix(pyro_depth_mat)
#Remove first column of names
pyro_depth_mat <- pyro_depth_mat[,-1]
#Stack matrix values into one column
pyro_depth_dist <- as.data.frame(as.table(pyro_depth_mat))
#Convert values to numeric
pyro_depth_dist$Freq <- as.numeric(pyro_depth_dist$Freq)
summary(pyro_depth_dist)

#Load in genetic data
pyro_genetic_mat <- read.csv("Pyropelta_GeneticDistances.csv", header = TRUE)
pyro_genetic_mat <- as.matrix(pyro_genetic_mat)
#Remove first column of names
pyro_genetic_mat <- pyro_genetic_mat[,-1]
#Stack matrix values into one column
pyro_genetic_dist <- as.data.frame(as.table(pyro_genetic_mat))
#Convert values to numeric
pyro_genetic_dist$Freq <- as.numeric(pyro_genetic_dist$Freq)
summary(pyro_genetic_dist)

#Create dataframe
pyro <- as.data.frame(pyro_genetic_dist$Freq)
pyro$depth <- pyro_depth_dist$Freq
pyro$geo <- pyro_geo_dist$Freq
#Make sure all the column names are good
colnames(pyro)
colnames(pyro) [1] <- c("gene")
colnames(pyro)
#Remove missing data
pyro.na <- na.omit(pyro)
summary(pyro.na)
#Test
plot(depth ~ geo, data = pyro.na)

#MLR with interaction
pyro.lm <- lm(gene ~ depth * geo, data = pyro.na)
summary(pyro.lm)
#depth: ***
#geo: ***
#interaction: **
#Check model
step(pyro.lm, direction = "backward")
#All stays in

#Z-score depth & geography distances
pyro.na$z_depth <- scale(pyro.na$depth, center = TRUE, scale = TRUE)
pyro.na$z_geo <- scale(pyro.na$geo, center = TRUE, scale = TRUE)
summary(pyro.na)

#Plots
pyro_d <- ggplot(data = pyro.na, mapping = aes(x = depth, y = gene)) +
  geom_point(color = "dodgerblue2", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  labs(title = "Pyropelta", x = "Difference in Depth (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

pyro_g <- ggplot(data = pyro.na, mapping = aes(x = geo, y = gene)) +
  geom_point(color = "black", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "black", size = 2, se = FALSE) +
  labs(title = "Pyropelta", x = "Geographic Distance (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

pyro <- 
ggplot(data = pyro.na) +
  geom_point(mapping = aes(x = z_geo, y = gene), color = "black",shape = 21, size = 2, alpha = .66)+
  geom_smooth(mapping = aes(x = z_geo, y = gene), method = "lm", color = "black", size = 2, se = FALSE) +
  geom_point(mapping = aes(x = z_depth, y = gene), color = "dodgerblue2", shape = 21, size = 2, alpha = .66) +
  geom_smooth(mapping = aes(x = z_depth, y = gene), method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  labs(title = "(d) Pyropelta", x = "z-distance", y = "Pairwise Genetic Distance") +
  theme(
    plot.title = element_text(face = "italic", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))


#########################
#LEPETODRILUS
########################

setwd(choose.dir())
#Calculate geographic distances
lep_latlon <- read.csv("Lepetodrilus_latlon.csv", header = FALSE)
summary(lep_latlon)
lep_dist_mat <- distm(lep_latlon, fun = distGeo)
summary(lep_dist_mat)
#Remove upper triangle of pairwise comparisons
lep_dist_mat[upper.tri(lep_dist_mat,diag=TRUE)] <- "NA"
#Stack matrix values into one column
lep_geo_dist <- as.data.frame(as.table(lep_dist_mat))
#Convert values to numeric
lep_geo_dist$Freq <- as.numeric(lep_geo_dist$Freq)
#c_geo_dist <- na.omit(c_geo_dist)
summary(lep_geo_dist)

#Load in depth data
lep_depth <- read.csv("Lepetodrilus_depth.csv", header = FALSE)
summary(lep_depth)
#lowercase v:
v1 <- lep_depth$V1
v2 <- lep_depth$V1
lep_depth_mat <- outer(v1, v2, FUN = "-")
lep_depth_mat <- abs(lep_depth_mat)
summary(lep_depth_mat)
#Stack matrix values into one column
lep_depth_dist <- as.data.frame(as.table(lep_depth_mat))
#Convert values to numeric
lep_depth_dist$Freq <- as.numeric(lep_depth_dist$Freq)
summary(lep_depth_dist)

#Load in genetic data
lep_genetic_mat <- read.csv("Lepetodrilus_GeneticDistances.csv", header = TRUE)
lep_genetic_mat <- as.matrix(lep_genetic_mat)
#Remove first column of names
lep_genetic_mat <- lep_genetic_mat[,-1]
#Stack matrix values into one column
lep_genetic_dist <- as.data.frame(as.table(lep_genetic_mat))
#Convert values to numeric
lep_genetic_dist$Freq <- as.numeric(lep_genetic_dist$Freq)
summary(lep_genetic_dist)

#Create dataframe
lep <- as.data.frame(lep_genetic_dist$Freq)
lep$depth <- lep_depth_dist$Freq
lep$geo <- lep_geo_dist$Freq
#Make sure all the column names are good
colnames(lep)
colnames(lep) [1] <- c("gene")
colnames(lep)
#Remove missing data
lep.na <- na.omit(lep)
summary(lep.na)
#Test
plot(depth ~ geo, data = lep.na)

#MLR with interaction
lep.lm <- lm(gene ~ depth * geo, data = lep.na)
summary(lep.lm)
#depth: ***
#geo: ***
#interaction: **
#Check model
step(lep.lm, direction = "backward")
# ***

#Z-score depth & geography distances
lep.na$z_depth <- scale(lep.na$depth, center = TRUE, scale = TRUE)
lep.na$z_geo <- scale(lep.na$geo, center = TRUE, scale = TRUE)
summary(lep.na)

#Plots
lep_d <- ggplot(data = lep.na, mapping = aes(x = depth, y = gene)) +
  geom_point(color = "dodgerblue2", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  labs(title = "Lepetodrilus", x = "Difference in Depth (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

lep_g <- ggplot(data = lep.na, mapping = aes(x = geo, y = gene)) +
  geom_point(color = "black", shape = 21, size = 2, alpha = .75) +
  geom_smooth(method = "lm", color = "black", size = 2, se = FALSE) +
  labs(title = "Lepetodrilus", x = "Geographic Distance (m)", y = "Genetic Distance") +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

lep <- 
ggplot(data = lep.na) +
  geom_point(mapping = aes(x = z_geo, y = gene), color = "black",shape = 21, size = 2, alpha = .66)+
  geom_smooth(mapping = aes(x = z_geo, y = gene), method = "lm", color = "black", size = 2, se = FALSE) +
  geom_point(mapping = aes(x = z_depth, y = gene), color = "dodgerblue2", shape = 21, size = 2, alpha = .66) +
  geom_smooth(mapping = aes(x = z_depth, y = gene), method = "lm", color = "dodgerblue2", size = 2, se = FALSE) +
  geom_hline(yintercept=.03, linetype="dashed", color = "red") +
  labs(title = "(e) Lepetodrilus", x = "z-distance", y = "Pairwise Genetic Distance") +
  theme(
    plot.title = element_text(face = "italic", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size =12, color = "black", face = "italic"),
    legend.key = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA))

########################
#Plotting all together
#######################
grid.arrange(bathy, cc, para, pyro, lep, nrow = 2, ncol = 3)
grid.arrange(bathy_d, cc_d, para_d, pyro_d, lep_d, nrow = 2, ncol = 3)
grid.arrange(bathy_g, cc_g, para_g, pyro_g, lep_g, nrow = 2, ncol = 3)

