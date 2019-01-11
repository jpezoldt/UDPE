# Author: Maria Litovchenko
# Date: 05.12.2018
# Function:
#       1) Process dynGENIE3 output to pipe to Cytoscape





library(data.table)
library(GENIE3)
library(igraph)


#####
#PATH & global
#####
subset_of_interest <- "Adv"
PATH_input <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/dynGENIE3/Seurat_selected_Split_Adv_NonAdv/",subset_of_interest,sep="")
PATH_output <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/scRNAseq/dynGENIE3/Seurat_selected_Split_Adv_NonAdv/",subset_of_interest,"/ToCytoscape",sep="")


# Functions -------------------------------------------------------------------
#' mytriangle
#' Defines triangle vertex shape for igraph
#' @param coords
mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}

# load Rds object which was saved from dynGENIE3
testNet <- readRDS(paste(PATH_input,"/dynGENIE3_Expr_thresh_0.05_inOne_withRegs.Rds",sep=""))
# select appropriate number of links. Recommended to test 100 - 1000
v_linkN <- c(100,200,500,1000,2000,4000,10000)
l_netLinks <- list()
for(i in 1:length(v_linkN)){
  v_linkN_i <- v_linkN[i]
  print(v_linkN_i)
  netLinks_i <- getLinkList(testNet$weight.matrix, v_linkN_i)
  l_netLinks[[i]] <- netLinks_i
}
netLinks <- getLinkList(testNet$weight.matrix, v_linkN)
saveRDS(l_netLinks, paste(PATH_output,"/l_netLinks_threshed_05_inOne.Rds",sep=""))
l_netLinks <- readRDS(paste(PATH_output,"/l_netLinks_threshed_05_inOne.Rds",sep=""))
netLinks <- l_netLinks[[6]]

#Subset table based on number of links per TF
tt_netLinks <- table(netLinks$regulatoryGene)
netLinks <- subset(netLinks, regulatoryGene %in% names(tt_netLinks[tt_netLinks > 8]))
levels(as.factor(as.character(netLinks$regulatoryGene)))
# Here I start to build a visualization with igraph, I don't put it in 
# function since I would like to leave it customizable for you
netLinks$regulatoryGene <- as.character(netLinks$regulatoryGene)
netLinks$targetGene <- as.character(netLinks$targetGene)

# create color and shape info of the nodes of the network
netInfo <- data.table(name = unique(c(netLinks$regulatoryGene,
                                      netLinks$targetGene)))
# put here color of the nodes you wish
netInfo[, color := rainbow(nrow(netInfo))]
# assign here shape of the node: currently it's random
netInfo[, shape := sample(c('triangle', 'square', 'circle'), nrow(netInfo),
                          replace = T)]
# in case you want only some names of the nodes to be dispayed, put the
# selected names here and set rest to NA
#coreCircGenes <- c("Hoxb3","Prrx1","Gata6","Mef2c","Pknox1","Hoxb5","Mybl2",
 #                  "Prrx2","Hoxb8","Sox18","Lhx8","Hoxc8","Irf3","Msx1","Hoxc6",
  #                 "Batf","Foxp1","Meis1","E2f7","Glis2","Mybl2","E2f1","Isl1",
   #                "Ebf3","Irf3","Zkscan3","Irf5","Irf8","Nfkb2","Meox1","Hoxc9",
    #               "Gbx2","Stat1","Pknox2","Hoxc2","Clock","Hlf")
#netInfo[, nameModified := sapply(name, function(x) if(x %in% coreCircGenes) {x}
       #                          else {NA})]
# create net itself
net <- graph_from_data_frame(netLinks, directed = T, vertices = netInfo)
layout <- layout_with_fr(net)
# plot
add_shape("triangle", clip = shapes("circle")$clip, plot = mytriangle)
plot(net, layout = layout, 
     vertex.shape = V(net)$shape, # shape of the node, comes from netInfo
     vertex.size = 3, # size of the nodes
     vertex.color = V(net)$color, # color of node, comes from netInfo
     vertex.label = V(net)$nameModified, # labels to dispay,comes from netInfo
     vertex.label.color = "black", # color of label
     vertex.label.cex = 1.2, # size of label
     edge.arrow.size = 0.1, # size of the arrow at the end
     #edge.width = 20 * E(net)$weight, # width of the link(edge)
     edge.color = 'black', # color of the link(edge)
     edge.curved = .4, # curvature of the edge
     main = 'iGraph test')

# Output for Cytoscape
write.table(netLinks, '~/Desktop/netLinksForCytoscape.txt', header = T, 
            sep = '\t', quote = F)
