B = read.csv("/home/nathan/SINC/DataApplication/Estimated_CytoMicro_Edge_Prob.csv",header = F)

B_edges = (abs(B) > 0.50)*1
B_edges = B_edges[otu_order,]
B_edges = rbind(matrix(0,nrow = 5,ncol = 29),B_edges)
B_edges = rbind(B_edges[1:59,],matrix(0,nrow = 5,ncol = 29),B_edges[60:95,] + B_edges[60:95,])
B_edges = rbind(B_edges[1:78,],matrix(0,nrow = 5,ncol = 29),B_edges[79:100,] + B_edges[79:100,])
B_edges = rbind(B_edges[1:92,],matrix(0,nrow = 5,ncol = 29),B_edges[93:105,] + B_edges[93:105,])
B_edges = rbind(B_edges[1:104,],matrix(0,nrow = 5,ncol = 29),B_edges[105:110,] + B_edges[105:110,])
B_edges = rbind(B_edges[1:114,],matrix(0,nrow = 5,ncol = 29),B_edges[115,] + B_edges[115,])
# B_edges = rbind(B_edges[1:114,],matrix(0,nrow = 5,ncol = 29),B_edges[114:115,] + B_edges[114:115,])

col_B = matrix(0,nrow = 115,ncol = 29)
for(i in 1:115) {
  col_B[i,] = cols[i]
}

names = c()
for(i in 1:29) {
  names = c(names,paste("C",i,sep = ""))
}

colnames(B_edges) = c(as.matrix(Cyt_names))

ppi.hpHC.hm <- melt((1 * as.matrix(B_edges)))
colnames(ppi.hpHC.hm) <- c("Microbes","Row","PPI")

colors = c('white','blue','red','orange','yellow','green','purple')

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/CytokineMicrobe.pdf")
ggplot(data = ppi.hpHC.hm, aes(x=Microbes, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey') + 
  scale_fill_manual(values = colors) + 
  geom_rect(xmin = 0, xmax = 5.5,ymin = 0.5,ymax = 29.5,color = 'grey',size = .2,fill = 'white') +
  geom_rect(xmin = 59.5, xmax = 64.5,ymin = 0.5,ymax = 29.5,color = 'grey',size = .2,fill = 'white') +
  geom_rect(xmin = 78.5, xmax = 83.5,ymin = 0.5,ymax = 29.5,color = 'grey',size = .2,fill = 'white') +
  geom_rect(xmin = 92.5, xmax = 97.5,ymin = 0.5,ymax = 29.5,color = 'grey',size = .2,fill = 'white') +
  geom_rect(xmin = 104.5, xmax = 109.5,ymin = 0.5,ymax = 29.5,color = 'grey',size = .2,fill = 'white') +
  geom_rect(xmin = 114.5, xmax = 119.5,ymin = 0.5,ymax = 29.5,color = 'grey',size = .2,fill = 'white') +

  geom_text(x = 2.5, y = 14.5, label = "F i r m i c u t e s",angle = 90,hjust = 0.5,size = 6) +
  geom_text(x = 62, y = 14.5, label = "A c t i n o b a c t e r i a",angle = 90,hjust = 0.5,size = 6) +
  geom_text(x = 81, y = 14.5, label = "B a c t e r o i d e t e s",angle = 90,hjust = 0.5,size = 6) +
  geom_text(x = 95, y = 14.5, label = "P r o t e o b a c t e r i a",angle = 90,hjust = 0.5,size = 6) +
  geom_text(x = 107, y = 14.5, label = "F u s o b a c t e r i a",angle = 90,hjust = 0.5,size = 6) +
  geom_text(x = 117, y = 14.5, label = "T M 7",angle = 90,hjust = 0.5,size = 6) +

  scale_x_continuous(breaks = c(27.0,61.5,73,81.0,87.0,90), labels = labels) +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(face = "bold",size = 14),
        axis.title.x=element_text(size = 15,face = 'bold'),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "none",axis.ticks.x = element_blank())

dev.off()
##########
omega = read.csv("/home/nathan/SINC/DataApplication/Estimated_MicroMicro_Edge_Prob.csv",header = F)
network = (abs(omega) > 0.50)*1
network = network[otu_order,otu_order]
network[upper.tri(network)] = network[upper.tri(network)]*100

# colors = c('white','grey','blue','red','orange','yellow','green','purple')
colors = c('white','grey','blue','red','orange','yellow','green','grey90', '#ccccff','#ffb2b2','#ffcc99','#ffff99','#ccff99')
ppi.hpHC.hm <- melt((cols_network * as.matrix(network)))
ppi.hpHC.hm[,2]  = sort(rep(1:90,90))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/MicrobeMicrobe_2.pdf")
ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_text(size = 16,face = 'bold'),axis.title.y=element_text(size = 16,face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +

  ## base of boxes
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[2], y = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[3], y = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[4], y = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[5], y = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
  geom_segment(aes(x = 89.5, y = 89.5, xend = 90.5, yend = 89.5)) +
  ## right side of boxes
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[1], x = c(54.5,68.5,77.5,84.5,89.5)[1], yend = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[2], x = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[3], x = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[4], x = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[5], x = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[5], xend = c(54.5,68.5,77.5,84.5,89.5)[5]),color = 'black') +
  geom_segment(aes(y = 0.5, x = 0.5, yend = 90.5, xend = 90.5)) +
  ##
  geom_text(y = 27.5, x = 26, label = "F i r m i c u t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 61.5, x = 60, label = "A c t i n o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 73, x = 71.5, label = "B a c t e r o i d e t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 81.5, x = 80, label = "P r o t e o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 87.5, x = 86, label = "F u s o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 90.0, x = 89, label = "T M 7",angle = 0,hjust = 1,size = 6,vjust = 0.0) +
  xlab("Microbe") + ylab("Microbe")
dev.off()

##############
network = (abs(omega) > 0.50)*1
network = network[otu_order,otu_order]
col = rep(1,90)
edges = apply(network,2,sum)
col = c(rep(0,5),col)
edges = c(rep(0,5),edges)
col = c(col[1:59],rep(0,5),col[60:95]+col[60:95])
edges = c(edges[1:59],rep(0,5),edges[60:95])
col = c(col[1:78],rep(0,5),col[79:100]+col[79:100])
edges = c(edges[1:78],rep(0,5),edges[79:100])
col = c(col[1:92],rep(0,5),col[93:105]+col[93:105])
edges = c(edges[1:92],rep(0,5),edges[93:105])
col = c(col[1:104],rep(0,5),col[105:110]+col[105:110])
edges = c(edges[1:104],rep(0,5),edges[105:110])
col = c(col[1:114],rep(0,5),col[115]+col[115])
edges = c(edges[1:114],rep(0,5),edges[115])

edges = as.data.frame(cbind(edges,1:120))


colors = c('white','blue','red','orange','yellow','green','purple')
# colors = rep('black',6)

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/NodeDegree.pdf")
ggplot(edges,aes(x=V2,xend=V2,y=0,yend=edges,colour = factor(col))) + 
  geom_segment(size = 0.30) +
  # geom_vline(xintercept = end_phylum) + 
  scale_color_manual(values = colors) + 
  geom_text(x = 2.5, y = 0.05, label = "F i r m i c u t e s",angle = 90,hjust = 0.0,col = 'black',size = 6) +
  geom_text(x = 62, y = 0.05, label = "A c t i n o b a c t e r i a",angle = 90,hjust = 0.0,col = 'black',size = 6) +
  geom_text(x = 81, y = 0.05, label = "B a c t e r o i d e t e s",angle = 90,hjust = 0.0,col = 'black',size = 6) +
  geom_text(x = 95, y = 0.05, label = "P r o t e o b a c t e r i a",angle = 90,hjust = 0.0,col = 'black',size = 6) +
  geom_text(x = 107, y = 0.05, label = "F u s o b a c t e r i a",angle = 90,hjust = 0.0,col = 'black',size = 6) +
  geom_text(x = 117, y = 0.05, label = "T M 7",angle = 90,hjust = 0.0,col = 'black',size = 6) +
  
  # annotate("rect", xmin = -0.5, xmax = 54.5, ymin = -Inf, ymax = Inf, fill = "gray 80", alpha = 0.1) +
  # annotate("rect", xmin = 68.5, xmax = 77.5, ymin = -Inf, ymax = Inf, fill = "gray 80", alpha = 0.1) +
  # annotate("rect", xmin = 84.5, xmax = 89.5, ymin = -Inf, ymax = Inf, fill = "gray 80", alpha = 0.1) +
  theme(axis.text.x = element_blank(),axis.title.y=element_text(size = 18,face = 'bold'),axis.text.y=element_text(size = 14),
        plot.title = element_text(hjust = 0.5),panel.grid.major.y = element_line(colour = 'grey'),
        panel.grid.minor.y = element_line(colour = 'light grey'),
        panel.grid.major.x = element_blank(), axis.title.x=element_text(size = 18,face = 'bold'),
        panel.background = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none") + ylab("Node Degree") + xlab("Microbe")
dev.off()
############
omega = read.csv("/home/nathan/SINC/DataApplication/Estimated_MicroMicro_Edge_Prob.csv",header = F)
network = (abs(omega) > 0.50)*1

angles = c((90 - 4*seq(1:90))[1:45],(270 - 4.0*seq(1:90))[46:90])

colors = c('blue','red','orange','yellow','green','purple')
grp = diag(cols_network)
shapes = grp
shapes[grp == 2] = 0
shapes[grp == 3] = 1
shapes[grp == 4] = 2
shapes[grp == 5] = 3
shapes[grp == 6] = 8
shapes[grp == 7] = 9
with_edges = apply(network[otu_order,otu_order],2,sum) != 0
# shapes = c(0,1,2,3,8,9,10)

colors = c('blue','red','orange','yellow','green','purple')
node_labels = c()
for(i in 1:length(labels_6)){
  print(i)
  node_labels = c(node_labels,paste("   ",labels_6[i],"   "))
}
labels_6 = tax[coauth$name,6]

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/Circle_network.pdf")
ggraph(mygraph,layout = 'linear', circular = TRUE) + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=FALSE)+
  geom_node_point(aes(size=1.2^n, color=factor(grp), fill=grp), alpha=1.0,shape = shapes[with_edges]) +
  # geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 2*3.14,start = 0, end = 2*3.14*54/87, fill = 'gray')) +
  # scale_size_continuous(range=c(0.5,5)) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 1.37, r = 1.5, start = 0, end = 3.14*54*2/90),color = colors[1],alpha = 0.2,fill = NA) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 1.37, r = 1.5, start = 3.14*68*2/90, end = 3.14*77*2/90),color = colors[3],alpha = 0.2,fill = NA) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 1.37, r = 1.5, start = 3.14*84*2/90, end = 3.14*89*2/90),color = colors[5],alpha = 0.2,fill = NA) +

  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 1.37, r = 1.5, start = 3.14*54*2/90, end = 3.14*68*2/90),color = colors[2],alpha = 0.2,fill = NA) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 1.37, r = 1.5,start = 3.14*77*2/90, end = 3.14*84*2/90),color = colors[4],alpha = 0.2,fill = NA) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 1.37, r = 1.5, start = 3.14*89*2/90, end = 3.14*90*2/90),color = colors[6],alpha = 0.2,fill = NA) +
  scale_color_manual(values=colors) +
  geom_node_text(aes(label = node_labels), angle = angles, hjust=c(rep(0,45),rep(1,45)), size=2.0)+
  theme_void() + 
  xlim(-1.5,1.5) + 
  ylim(-1.5,1.5) +
  theme(
    legend.position="none",
    plot.margin=unit(c(0.0,0.0,0.0,0.0), "null"),
    panel.spacing=unit(c(0.0,0.0,0.0,0.0), "null")
  )
dev.off()


##########
## methods comparison plot
B = read.csv("/home/nathan/SINC/DataApplication/Estimated_CytoMicro_Edge_Prob.csv",header = F)

B_edges = (abs(B) > 0.50)*1
B_edges = B_edges[otu_order,]

OTUs = c()
for(i in c(14,19,20)) {
  OTUs = c(OTUs,paste("OTU",i))
}

B_edges_named = B_edges
colnames(B_edges_named) = Cyt_names$x
Cyt_names_small = names(sort(Cyt_connected_named,decreasing = TRUE))[1:5]
Cyt_names_small = as.character(Cyt_names$x[c(10,18,20)])

Cyt_connected = apply(B_edges[otu_order[c(14,19,20)],],2,sum)
Cyt_keep = names(sort(Cyt_connected,decreasing = TRUE))[1:5]

n1 = SE[c(14,19,20),c(14,19,20)]
colnames(n1) = rownames(n1) = OTUs
g1 = graph_from_adjacency_matrix(n1)
# E(g1)$lty = c(3,3,3,1,3,1)

n2 = matrix(0,ncol = 6,nrow = 6)
n2[1:3,4:6] = B_edges[c(14,19,20),c(10,18,20)]
colnames(n2) = rownames(n2) = c(OTUs,Cyt_names_small)
g2 = graph_from_adjacency_matrix(t(n2))

g3 = g1 + g2
g3$layout = layout_with_fr
g3$layout = layout_with_graphopt
# g3$lty = 3
pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/Methods_Comparison.pdf")
plot(g3,edge.arrow.mode = c(rep(2,length(E(g2))),rep(0,length(E(g1))),rep(0,length(E(g4)))),vertex.color=c(rep("darkseagreen",3),rep("lightcoral",3)),
     edge.arrow.size=0.6,vertex.shape =c(rep('circle',3),rep('square',3)),
     vertex.label.cex = 2.0,edge.width = c(rep(1,length(E(g2))),3,3,4,3),
     edge.color = 'black',vertex.label.color = 'black',edge.lty = c(rep(1,length(E(g2))),c(0,0,1,0,0,1)),vertex.size = 48)
dev.off()


plot(g1,edge.arrow.mode = c(rep(0,length(E(g2))),rep(0,length(E(g1))),rep(0,length(E(g4)))),vertex.color=c(rep("darkseagreen",3),rep("lightcoral",3)),
     edge.arrow.size=0.6,vertex.shape =c(rep('circle',3),rep('square',3)),
     vertex.label.cex = 2.0,edge.width = c(rep(4,length(E(g2))),3,3,4,3),
     edge.color = 'black',vertex.label.color = 'black',edge.lty = c(rep(1,length(E(g2))),c(3,3,1,3,3,1)),vertex.size = 48)

############
##
#################
######### cytokine and microbes
OTUs = c()
for(i in 1:26) {
  OTUs = c(OTUs,paste("OTU",i))
}

B_edges_named = B_edges
colnames(B_edges_named) = Cyt_names$x
Cyt_connected_named = apply(B_edges_named[otu_order[1:26],],2,sum)
Cyt_names_small = names(sort(Cyt_connected_named,decreasing = TRUE))[1:5]

Cyt_connected = apply(B_edges[otu_order[1:26],],2,sum)
Cyt_keep = names(sort(Cyt_connected,decreasing = TRUE))[1:5]

n1 = network[otu_order[1:26],otu_order[1:26]]
colnames(n1) = rownames(n1) = OTUs

g1 = graph_from_adjacency_matrix(n1)

n2 = matrix(0,ncol = 26+5,nrow = 26+5)
n2[1:26,27:31] = B_edges[otu_order[1:26],Cyt_keep]
colnames(n2) = rownames(n2) = c(OTUs,Cyt_names_small)
g2 = graph_from_adjacency_matrix(t(n2))

g3 = g1 + g2
g3$layout = layout_with_fr
g3$layout = layout_with_graphopt
plot(g3,edge.arrow.mode = c(rep(2,length(E(g2))),rep(0,length(E(g1)))),vertex.color=c(rep("darkseagreen",26),rep("lightcoral",6)),
     edge.arrow.size=0.4,vertex.shape =c(rep('circle',26),rep('square',6)),
     vertex.label.cex = .75,edge.width = c(rep(1,length(E(g2))),rep(4,length(E(g1)))),
     edge.color = 'black',vertex.label.color = 'black')


##############
cn = c()
for(i in 1:100){
  cn = c(cn,paste("V",i,sep = ""))
}
rn = c()
for(i in 1:100){
  rn = c(rn,paste("v",i,sep = ""))
}
omega = read.table("/home/nathan/Documents/SINC_from_Marina2/SINC/Simulations/SimulatedData/Band/adj1.csv",header = F)
colnames(omega)  = cn
rownames(omega) = rn

colors = c('white','black','blue','red','orange','yellow','green','purple')
ppi.hpHC.hm <- melt(as.matrix(omega))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/Sample_Band.pdf")
ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey93') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 40,face = 'bold'),
        legend.position = "none") + ggtitle("Band")
dev.off()

omega = read.table("/home/nathan/Documents/SINC_from_Marina2/SINC/Simulations/SimulatedData/Cluster/adj1.csv",header = F)
colnames(omega)  = cn
rownames(omega) = rn

colors = c('white','black','blue','red','orange','yellow','green','purple')
ppi.hpHC.hm <- melt(as.matrix(omega))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/Sample_Cluster.pdf")
ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey93') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 40,face = 'bold'),
        legend.position = "none") + ggtitle("Cluster")
dev.off()


omega = read.table("/home/nathan/Documents/SINC_from_Marina2/SINC/Simulations/SimulatedData/Hub/adj1.csv",header = F)
colnames(omega)  = cn
rownames(omega) = rn

colors = c('white','black','blue','red','orange','yellow','green','purple')
ppi.hpHC.hm <- melt(as.matrix(omega))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/Sample_Hub.pdf")
ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey93') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 40,face = 'bold'),
        legend.position = "none") + ggtitle("Hub")
dev.off()

omega = read.table("/home/nathan/Documents/SINC_from_Marina2/SINC/Simulations/SimulatedData/Random/adj1.csv",header = F)
colnames(omega)  = cn
rownames(omega) = rn

colors = c('white','black','blue','red','orange','yellow','green','purple')
ppi.hpHC.hm <- melt(as.matrix(omega))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")

pdf("/home/nathan/EMGS/mLDM/Revision/NewPlots/Sample_Random.pdf")
ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey93') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 40,face = 'bold'),
        legend.position = "none") + ggtitle("Random")
dev.off()
