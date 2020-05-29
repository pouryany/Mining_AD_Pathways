rm(list = ls())
print("Mapping to PCxN")
options(StringAsFactors = F)

library(RColorBrewer)
library(igraph)
library("GGally")
library(intergraph)
library(intergraph)

output.dir  <- paste0("Output")
fdr.paths   <- 0.1
edge.fdr    <- 0.05
cor.thresh  <- 0

# Loading PCxN
pcxn.dir   <- paste0("data/improved_PCxN_KEGG_PCxN.RDS")
PCxN       <- readRDS(pcxn.dir)
temp       <- PCxN

temp$p.Adjust <- p.adjust(temp$p.value, method="fdr")
res <- temp[temp$p.Adjust < edge.fdr,]
res <- res[abs(res$PathCor) > cor.thresh,]
writeLines(paste0("Parameter Setting: \n"))
writeLines("Number of pathways examined: \n")
print(length(union(unique(res$Pathway.A),unique(res$Pathway.B))))
writeLines("Parameter Setting: \n")


if(nrow(res) < 1) next

path.dir <- read.csv("new_data/Overall_rank_minusDUPandADDs.txt", sep = "\t") 
path.dir <- path.dir$Pathway

San.paths  <- path.dir
San.paths  <- paste0("Pathway.KEGG_",San.paths)
all.paths  <- union(unique(PCxN$Pathway.B),unique(PCxN$Pathway.A))
San.paths  <- San.paths[which(((San.paths) %in% all.paths))]
San.paths  <- San.paths[1:50]
all.paths  <- union(unique(res$Pathway.B),unique(res$Pathway.A))
San.paths  <- San.paths[which((San.paths) %in% all.paths)]
San.paths <- gsub("Pathway.KEGG_","", San.paths)


if(length(San.paths) <=1) next

net <- res[,c(1,2,4)]
net <- as.matrix(net)
net <- igraph::graph_from_data_frame(net, directed = F)

V(net)$name <- gsub("Pathway.KEGG_","", V(net)$name)

cols <- brewer.pal(8,"Set2")
sub.mods <- induced_subgraph(net,(San.paths))
comps    <- components(sub.mods)
color.it <- which(table(comps$membership) > 1)

node.cols   <- ifelse(comps$membership %in% color.it,cols[comps$membership],"#D0D0D0")
 
V(sub.mods)$color <- node.cols

set.seed(2)

E(sub.mods)$color    <- ifelse(as.numeric(E(sub.mods)$PathCor)>0,
                               "#E41A1C","#377EB8")

set.seed(2)
plot(sub.mods, vertex.size = 5,
     vertex.color = node.cols, edge.curved = 0.1)

legend(x = "topleft", legend=c("Positive Cor", "Negative Cor"),
       col=c("#E41A1C","#377EB8"), lty=1, lwd =2, cex = 1.6 , bty = "n")


apop.network <- intergraph::asNetwork(sub.mods)

network::network.vertex.names(apop.network)
#7
set.seed(7)
ggnet2(apop.network,layout.exp = 0.2,
       arrow.size = 4,
       label = T, edge.size = 0.5,
       layout.par = list("cell.jitter" = 0.1),
       edge.color = "color",
       node.color = "color",
       node.size = 6,
       label.size = 3,label.trim = T,arrow.gap = 0.015,
       mode = "fruchtermanreingold") +
    theme(legend.title=element_blank())


network::get.edge.attribute(apop.network,"color")





## Neighbors of AD

neighbors <- neighborhood(net,1 , nodes =San.paths[1:5])



cols <- brewer.pal(8,"Set2")
sub.mods <- induced_subgraph(net,unlist(neighbors))
comps    <- components(sub.mods)
color.it <- which(table(comps$membership) > 1)

node.cols   <- ifelse(V(sub.mods)$name %in% San.paths[1:5],cols[2],"#D0D0D0")

V(sub.mods)$color <- node.cols

set.seed(2)

E(sub.mods)$color    <- ifelse(as.numeric(E(sub.mods)$PathCor)>0,
                               "#E41A1C","#377EB8")

apop.network <- intergraph::asNetwork(sub.mods)

network::network.vertex.names(apop.network)
#7
set.seed(3)
ggnet2(apop.network,layout.exp = 0.2,
       arrow.size = 4,
       label = T,
       edge.size = 0.5,
       layout.par = list("cell.jitter" = 0.05),
       edge.color = "color",
       node.color = "color",
       node.size = 8,
       label.size = 4,label.trim = T,arrow.gap = 0.015,
       mode = "fruchtermanreingold") +
    theme(legend.title=element_blank())












