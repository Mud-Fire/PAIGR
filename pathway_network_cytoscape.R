library(SPIA)
library(KEGGgraph)
library(igraph)
library(limma)
library(PADOG)
library(ROntoTools)
library(GSA)
library(CePa)
library(annotate)
library(KEGGREST)
library(EnrichmentBrowser)
library(XML)

path  <- "F:/branch2/kegg/"
files <- list.files(path = path, pattern = "*.xml")
netnode <- NULL
pathnames <- NULL
for(file in files){
    
    mapkpathway <- try(parseKGML(paste(path,file,sep = "")),TRUE)
    pathwayname <- mapkpathway@pathwayInfo@title
    pathnames <- c(pathnames,pathwayname)
    kgml <- xmlParse(paste(path,file,sep = ""))
    out_nodes = getNodeSet(kgml, "//entry[@type='map']", fun = xmlToList)
    links_summary <- data.frame(do.call(rbind, out_nodes))
    for(names in links_summary$graphics){
        nextNode <- names[1]
        if (grepl(pattern = "TITLE:",nextNode)){
            nextNode <- substring(nextNode,7)
        }
        link <- c(pathwayname,nextNode)
        netnode <- rbind(netnode,link)
    }
}
col_name <- c("node1","node2")
colnames(netnode) <- col_name

netnode_1 <- netnode[which(netnode[,2] %in% pathnames),]
netnode_2 <- netnode_1[which(as.character(netnode_1[,1]) != as.character(netnode_1[,2])),]
netnode_3 <- unique(netnode_2)
netnode_4 <- netnode_3[,c(2,1)]

netnode_5 <- as.vector(unlist(merge(netnode_3,netnode_4,c(1,2))[,1]))
netnode_6 <- as.vector(unlist(merge(netnode_3,netnode_4,c(1,2))[,2]))
netnode_7 <- cbind(netnode_5,netnode_6)
netnode_8 <- paste(netnode_7[,1],netnode_7[,2],"")
netnode_9 <- paste(netnode_3[,1],netnode_3[,2],"")
netnode_10 <- netnode_3[-which(netnode_9 %in% netnode_8),]
write.csv(netnode_10, file = "netnode_10.csv", quote = TRUE, na = "NA",  row.names = TRUE)


# 
# netnode_1 <- netnode[,]
# netnode_2 <- netnode_1[which(as.character(netnode_1[,1]) != as.character(netnode_1[,2])),]
# netnode_3 <- unique(netnode_2) 
# netnode_4 <- netnode_3[,c(2,1)]
# 
# netnode_5 <- as.vector(unlist(merge(netnode_3,netnode_4,c(1,2))[,1]))
# netnode_6 <- as.vector(unlist(merge(netnode_3,netnode_4,c(1,2))[,2]))
# netnode_7 <- cbind(netnode_5,netnode_6)
# netnode_8 <- paste(netnode_7[,1],netnode_7[,2],"")
# netnode_9 <- paste(netnode_3[,1],netnode_3[,2],"")
# netnode_10 <- netnode_3[-which(netnode_9 %in% netnode_8),]
# write.csv(netnode_10, file = "netnode_10_1.csv", quote = TRUE, na = "NA",  row.names = TRUE)
