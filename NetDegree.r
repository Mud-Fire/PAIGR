kgml.path = system.file("extdata/keggxml/hsa",package="SPIA")
mdir = kgml.path
paths <- dir(mdir,pattern=".xml")

options(digits=4)

library(SPIA)
library(KEGGgraph)
library(igraph)
library(limma)
library(PADOG)
library(EnrichmentBrowser)
library(ggplot2)


pathnames <- NULL
pathIDs <- NULL
for (ij in 1:length(paths)){
  mapkpathway <- try(parseKGML(paste(mdir,paths[[ij]],sep="/")),TRUE)
  pathwayname <- mapkpathway@pathwayInfo@title
  pathID <- mapkpathway@pathwayInfo@name
  pathnames <- c(pathnames,pathwayname)
  pathIDs <- c(pathIDs,pathID)
}
pathIDsp <- strsplit(as.character(pathIDs),":")
pathIDsps <- do.call(rbind,pathIDsp)
pathIDs <- pathIDsps[,2]
names(pathnames) <- pathIDs

pathnetlist <- NULL
pathwaynames <- NULL
for (i in 1:length(paths)){
  mapkpathway <- try(parseKGML(paste(mdir,paths[[i]],sep="/")),TRUE)
  pathnodes <- nodes(mapkpathway)
  pathwayname <- mapkpathway@pathwayInfo@title
  #pathwaynames <- c(pathwaynames,pathwayname)
  for(j in 1:length(pathnodes)){
    if(pathnodes[[j]]@type == 'map'){
      pathnetnode <- pathnodes[[j]]@graphics@name
      nodesp <- strsplit(as.character(pathnetnode),":")
      pathnetnode <- do.call(rbind,nodesp)
      if(length(pathnetnode) == 2){
        pathnetnode <- pathnetnode[2]
      }
      else{
        pathnetnode <- pathnetnode[1]
      }
      pathnetnode <- intersect(pathnetnode,pathnames)
      if(length(setdiff(pathnetnode,pathwayname)) !=0){
        pathnetlist <- rbind(pathnetlist,c(pathwayname,pathnetnode))
      }
    }
  }
  
}

edges1 <- paste(pathnetlist[,1],"|",pathnetlist[,2],sep="")
edges2 <- paste(pathnetlist[,2],"|",pathnetlist[,1],sep="")
coedge <- intersect(edges1,edges2)
alledges <- setdiff(edges1,coedge)
alledges <- alledges[!duplicated(alledges)]
circsp <- strsplit(as.character(alledges),"\\|")
alledgeslist <- do.call(rbind,circsp)
circsp <- strsplit(as.character(coedge),"\\|")
coedges <- do.call(rbind,circsp)
nums <- NULL
for(ii in 1:nrow(coedges)){
  for(jj in ii:nrow(coedges)){
    if(coedges[ii,1]==coedges[jj,2] & coedges[ii,2]==coedges[jj,1]){
      nums <- c(nums,ii)
    }
  }
}

coedges1 <- coedges[nums,]
alledgeslist <- rbind(alledgeslist,coedges1)
e <- unique.matrix(alledgeslist)
g<-graph_from_edgelist(e,directed = F)
netdegree <- igraph::degree(g)
netdegree <- -sort(-netdegree)
