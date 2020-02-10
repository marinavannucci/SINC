library(huge)
library(ggplot2)
library(reshape2)
library(reticulate)

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(igraph)
library(ggraph)
library(colormap)

# Load data
OTU_names = as.character(c(as.matrix(read.table("/home/nathan/EMGS/mLDM/RealData/otu_names.csv",header = F))))
otu_order = c(as.matrix(read.table("/home/nathan/EMGS/mLDM/RealData/otu_order.txt",header = F)))
Cyt_names = read.table("/home/nathan/EMGS/mLDM/RealData/Cytokine_names.txt")
load("/home/nathan/EMGS/mLDM/NewData/momspi16S_tax.rda")
tax = momspi16S_tax[OTU_names,]
omega = read.csv("/home/nathan/EMGS/mLDM/NewData/Results/EZ_t.csv",header = F)
network = (abs(omega) > 0.50)*1

dataUU = cbind(OTU_names,network*lower.tri(network))
dataUU[dataUU == 0] = NA
colnames(dataUU) = c("from",OTU_names)
rownames(dataUU) = 1:90
dataUU = as.data.frame(dataUU)

# Transform the adjacency matrix in a long format
connect <- dataUU %>% 
  gather(key="to", value="value", -1) %>%
  mutate(to = gsub("\\.", " ",to)) %>%
  na.omit() 

# Number of connection per person
c( as.character(connect$from), as.character(connect$to)) %>%
  as.tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> coauth
colnames(coauth) <- c("name", "n")
#dim(coauth)

# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = TRUE )

# Find community
# com <- walktrap.community(mygraph)
un = unique(tax[,2])
genus = tax[coauth$name,2]
com = c()
for(i in 1:90) {
  com = c(com,which(un == genus[i]))
}

taxx = tax
taxx[is.na(tax)] = "none"
un = unique(taxx[,6])
genus = taxx[coauth$name,6]
com2 = c()
for(i in 1:90) {
  com2 = c(com2,which(un == genus[i]))
}
#max(com$membership)

coauth <- coauth %>% 
  # mutate( grp = com) 
  mutate( grp = com) %>%
  mutate(subgrp = com2) %>%
  arrange(grp,subgrp) 
# mutate(name=factor(name, name))

# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = TRUE )

OTU_order = coauth$name

labels_6 = tax[coauth$name,6]
labels_6[is.na(labels_6)] = "No Genus"

####

phylum = tax[,2]
unique_phylums = unique(phylum)

order = c()
col = c()
group_size = c()
for(i in 1:length(unique_phylums)) {
  order = c(order,which(phylum == unique_phylums[i]))
  col = c(col,rep(i,length(which(phylum == unique_phylums[i]))))
  group_size = c(group_size,length(which(phylum == unique_phylums[i])))
}
cols = col

#################
### cytokine plot
B = read.csv("/home/nathan/EMGS/mLDM/NewData/Results/phi_t.csv",header = F)
## B = B[,-1]

mean(abs(B) > 0.025)
B_edges = (abs(B) > 0.50)*1
mean(B_edges)

col = c()
for(i in 1:90) {
  col = c(col,cols[i,i] - 1)
}

col_B = matrix(0,nrow = 90,ncol = 29)
for(i in 1:90) {
  col_B[i,] = cols[i]
}

names = c()
for(i in 1:29) {
  names = c(names,paste("C",i,sep = ""))
}

colnames(col_B) = c(as.matrix(Cyt_names))

ppi.hpHC.hm <- melt((col_B * as.matrix(B_edges[otu_order,])))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

colors = c('white','blue','red','orange','yellow','green','purple')

ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = c(28,48,88), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_text(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "none") + 
  ggtitle("Microbiome-Cytokine")

#######################
########### microbe - microbe plots
cols_network = matrix(1,nrow = 90,ncol = 90)
for(i in 1:90){
  cols_network[which(cols == i),which(cols == i)] = i + 1 
}

omega = read.csv("/home/nathan/EMGS/mLDM/NewData/Results/EZ_t.csv",header = F)

network = (abs(omega) > 0.50)*1

ppi.hpHC.hm <- melt((cols_network * as.matrix(network[otu_order,otu_order])))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

colors = c('white', 'grey','blue','red','orange','yellow','green','purple')

ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") + 
  ggtitle("Microbe - Microbe")

edges = apply(network[otu_order,otu_order],2,sum)
edges = as.data.frame(cbind(edges,1:90))

colors = c('blue','red','orange','yellow','green','purple')

ggplot(edges,aes(x=V2,xend=V2,y=0,yend=edges,colour = factor(cols))) + 
  geom_segment() +
  scale_color_manual(values = colors) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") + ylab("Node Degree")

angles = c((90 - 4*seq(1:90))[1:43],(270 - 4.08*seq(1:90))[44:87])

ggraph(mygraph,layout = 'linear', circular = TRUE) + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=FALSE)+
  geom_node_point(aes(size=1.1^(n), color=factor(grp), fill=grp), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=colors) +
  geom_node_text(aes(label = labels_6), angle = angles, hjust=c(rep(0,43),rep(1,44)), size=2.3)+
  theme_void() + 
  xlim(-1.2,1.2) + 
  ylim(-1.2,1.2) +
  theme(
    legend.position="none",
    plot.margin=unit(c(0.0,0.0,0.0,0.0), "null"),
    panel.spacing=unit(c(0.0,0.0,0.0,0.0), "null")
  ) 

####### tables #####
OTUs_connected = apply(network,1,sum)
most_otu_cyt = colnames(otu_data)[order(OTUs_connected,decreasing = TRUE)[1:10]]
most_otu_cyt_table = cbind(sort(OTUs_connected,decreasing = TRUE)[1:10],tax[as.character(most_otu_cyt),2:6])

least_otu_cyt = colnames(otu_data)[order(OTUs_connected,decreasing = FALSE)[1:10]]
least_otu_cyt_table = cbind(sort(OTUs_connected,decreasing = FALSE)[1:10],tax[as.character(least_otu_cyt),2:6])

Cyt_connected = apply(B_edges,2,sum)
Cyt_keep = names(sort(Cyt_connected,decreasing = TRUE))[1:6]
Cyt_names = colnames(cyt_data)[sapply(1:6,function(x) which(order(Cyt_connected,decreasing = TRUE) == x))]

OTU_cyt_connected = apply(B_edges,1,sum)
most_cyt_otu = colnames(otu_data)[order(OTU_cyt_connected,decreasing = TRUE)[1:10]]
most_cyt_otu_table = cbind(sort(OTU_cyt_connected,decreasing = TRUE)[1:10],tax[as.character(most_cyt_otu),2:6])

#################
######### cytokine and microbes
OTUs = c()
for(i in 1:26) {
  OTUs = c(OTUs,paste("OTU",i))
}

Cyt_connected = apply(B_edges[1:26,],2,sum)
Cyt_keep = names(sort(Cyt_connected,decreasing = TRUE))[1:5]
Cyt_names = colnames(cyt_data)[sapply(1:5,function(x) which(order(Cyt_connected,decreasing = TRUE) == x))]

n1 = network[otu_order[1:26],otu_order[1:26]]
colnames(n1) = rownames(n1) = OTUs

g1 = graph_from_adjacency_matrix(n1)

n2 = matrix(0,ncol = 26+5,nrow = 26+5)
n2[1:26,27:31] = B_edges[otu_order[1:26],Cyt_keep]
colnames(n2) = rownames(n2) = c(OTUs,Cyt_names)
g2 = graph_from_adjacency_matrix(t(n2))

g3 = g1 + g2
g3$layout = layout_with_fr
g3$layout = layout_with_graphopt
plot(g3,edge.arrow.mode = c(rep(2,length(E(g2))),rep(0,length(E(g1)))),vertex.color=c(rep("darkseagreen",26),rep("lightcoral",6)),
     edge.arrow.size=0.4,vertex.shape =c(rep('circle',26),rep('square',6)),
     vertex.label.cex = .75,edge.width = c(rep(1,length(E(g2))),rep(4,length(E(g1)))),
     edge.color = 'black',vertex.label.color = 'black')

####### supplemental plots
ggraph(mygraph,layout = 'linear', circular = TRUE) + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.1, edge_width=0.3, fold=TRUE)+
  geom_node_point(aes(size=1.1^(n), color=factor(grp), fill=grp), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  # # scale_color_manual(values=mycolor) +
  geom_node_text(aes(label = labels_6), angle = (90 - 4*seq(1:90)), hjust=0, size=2.3)+
  theme_void() + 
  theme(
    legend.position="none",
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "null")
    # panel.spacing=unit(c(0,0,0,0), "null")
  ) 
