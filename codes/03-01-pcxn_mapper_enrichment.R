rm(list = ls())
print("Mapping to enrich.map")
options(StringAsFactors = F)

library(RColorBrewer)
library(igraph)
library("GGally")
library(intergraph)
library(intergraph)


output.dir  <- paste0("Output")
fdr.paths   <- 0.1
edge.fdr    <- 5e-2
cor.thresh  <- 0

# Loading enrich.map
enrich.map.dir   <- paste0("data/KEGG_KEGG_ENRICHMENT.RDS")
enrich.map       <- readRDS(enrich.map.dir)
temp       <- enrich.map

temp$p.Adjust <- p.adjust(temp$pval, method="fdr")
res <- temp[temp$p.Adjust < edge.fdr,]

writeLines("Number of pathways examined: \n")
print(length(union(unique(res$x),unique(res$y))))


if(nrow(res) < 1) next

path.dir <- read.csv("new_data/Overall_rank_minusDUPandADDs.txt", sep = "\t") 
path.dir <- path.dir$Pathway

San.paths  <- path.dir
San.paths  <- paste0("Pathway.KEGG_",San.paths)
all.paths  <- union(unique(enrich.map$x),unique(enrich.map$y))
San.paths  <- San.paths[which(((San.paths) %in% all.paths))]
San.paths  <- San.paths[1:50]
all.paths  <- union(unique(res$x),unique(res$y))
San.paths  <- San.paths[which((San.paths) %in% all.paths)]
San.paths <- gsub("Pathway.KEGG_","", San.paths)


if(length(San.paths) <=1) next

net <- res[,c(1,2)]
net$weight <- -log(res$p.Adjust)
net <- as.matrix(net)
net <- igraph::graph_from_data_frame(net, directed = F)

V(net)$name <- gsub("Pathway.KEGG_","", V(net)$name)

cols <- brewer.pal(8,"Set2")
sub.mods      <- induced_subgraph(net,(San.paths))
comps    <- components(sub.mods)
color.it <- which(table(comps$membership) > 1)

node.cols   <- ifelse(comps$membership %in% color.it,cols[comps$membership],"#D0D0D0")
 
V(sub.mods)$color <- node.cols
V(sub.mods)$name  <- ifelse(V(sub.mods)$name %in% 
                                union(San.paths[1:5],grep("disease",
                                                           V(sub.mods)$name,
                                                           value = T)),
                                      V(sub.mods)$name,"")

# 
# node.cols   <- cols[clusts$membership]
# small.clust <- which(table(clusts$membership) <= table(clusts$membership)[5],useNames = T)
# node.cols[clusts$membership %in% small.clust] <- NA


# pdf(paste0(save.dir,"/AllPathways_enrich.mapCor_",gsub('\\.', '', cor.thresh),
#            "_FDRDEP_",gsub('\\.', '', fdr.paths),ifelse(FLAG, "_FLAG200",""),".pdf"),
#     width = 18,height = 11)
# 
# l <- layout_nicely(sub.mods)
# l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

# set.seed(2)
# plot(sub.mods, vertex.size = 5,
#      vertex.color = node.cols)
# 
# 
# legend(x = "topleft", legend=c("Positive Cor", "Negative Cor"),
#        col=c("#E41A1C","#377EB8"), lty=1, lwd =2, cex = 1.6 , bty = "n")


apop.network <- intergraph::asNetwork(sub.mods)

network::network.vertex.names(apop.network)
#7

weight <- network::get.edge.attribute(apop.network,"weight")
weight <- 4 *as.numeric(weight)/max(as.numeric(weight))
network::set.edge.attribute(apop.network,"weight",(weight))


set.seed(7)
ggnet2(apop.network,
       arrow.size = 4,layout.exp = 0.2,
       label = T,
       edge.size = "weight",
       layout.par = list("cell.jitter" = 0.03),
       node.color = "color",
       node.size = 6,
       label.size =4,label.trim = T,arrow.gap = 0.015,
       mode = "fruchtermanreingold") +
    theme(legend.title=element_blank())







