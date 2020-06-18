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


getdataaslist = function(x) {
    x = get(x)
    exp = experimentData(x)
    dataset = exp@name
    disease = notes(exp)$disease
    dat.m = exprs(x)
    ano = pData(x)
    design = notes(exp)$design
    annotation = paste(x@annotation, ".db", sep = "")
    targetGeneSets = notes(exp)$targetGeneSets
    list = list(dataset, disease, dat.m, ano, design, annotation, targetGeneSets)
    names(list) = c("dataset", "disease", "dat.m", "ano", "design", "annotation", 
                    "targetGeneSets")
    return(list)
}

# parameter (set)
# return sample list(group, esetm)
getGSEdatamatrix = function(x){
    data(list = x, package = "KEGGdzPathwaysGEO")
    #==========================prepare-sample-data===================
    dlist = getdataaslist(x)
    block = dlist$ano$Block
    Block = block
    block = factor(Block)
    group = dlist$ano$Group
    G = factor(group)
    force(G)
    force(block)
    paired = dlist$design == "Paired"
    esetm = dlist$dat.m
    annotation = dlist$annotation
    #判断数据类型是否是矩阵
    stopifnot(class(esetm) == "matrix")
    #判断矩阵大小
    stopifnot(all(dim(esetm) > 4))
    #判断group类型是否为字符或者factor
    stopifnot(class(group)%in%c("factor","character"))
    #判断group的长度是否和esetm的长度吻合，即样本数目是否一致
    stopifnot(length(group) == dim(esetm)[2])

    #是否成对设计的，如果成对设计则满足
    if (paired) {
        #block数量是否与group相同
        stopifnot(length(block) == length(group))
        #是否所有block参数的数量都为2
        stopifnot(all(table(block) == 2))
    }

    #将基因编号换为extrez系统代号
    aT1 = filteranot(esetm, group, paired, block, annotation)
    esetm = esetm[rownames(esetm) %in% aT1$ID, ]
    rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]
    #计算样本数量
    topSigNum = dim(esetm)[1]
    return(list(group,esetm))
}

# parameter (fileName)
# return hsa list(hsa.Name,gMatrix, geneList)
getHsaGMatrix = function(x){    
    hsa.Data <- try(parseKGML(x),TRUE)
    hsa.Name <- hsa.Data@pathwayInfo@title
    hsa.pathway <- KEGGpathway2Graph(hsa.Data, expandGenes = T)
    hsa.Net <- igraph.from.graphNEL(hsa.pathway)
    hsa.Matrix <- as_adjacency_matrix(hsa.Net)
    gMatrix <- as.matrix(hsa.Matrix)
    hsa.GeneName <- rownames(gMatrix)
    genesList <- unlist(strsplit(hsa.GeneName,":"))[c(F,T)]
    rownames(gMatrix) <- genesList
    colnames(gMatrix) <- genesList
    return(list(hsa.Name, gMatrix, genesList))
}

# parameter ((group, esetm), geneList, gMatrix)
# return sample & hsa list(esetmC, esetmD, samgMatrix)
getHsaSamMatrix = function(x,y,z){
    group <- x[[1]]
    esetm <- x[[2]]
    genesList <- y
    g.Matrix <- z
    esetmSam <- subset(esetm, rownames(esetm) %in% c(genesList))
    samgMatrix <- g.Matrix[rownames(esetmSam), rownames(esetmSam)]
    names(group) = colnames(esetm)
    # numC = summary(as.factor(group))["c"]
    # numD = summary(as.factor(group))["d"]
    groupC = names(group[which(group == "c")])
    esetmC <-esetmSam[,groupC]
    groupD = names(group[which(group == "d")])
    esetmD <-esetmSam[,groupD]
    return(list(esetmC, esetmD, samgMatrix))
}


getHSADegree = function(x,y){
    
    path <- x
    files <- y
    pathnames <- NULL
    pathIDs <- NULL
    print(paste(path,files[1],sep = ""))
    for (file in files){
        
        mapkpathway <- try(parseKGML(paste(path,file,sep = "")),TRUE)
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
    for (file in files){
        mapkpathway <- try(parseKGML(paste(path,file,sep = "")),TRUE)
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
    return(netdegree)
}

options(scipen = 15)
# hsa prepare
path  <- "F:/branch2/kegg/"
files <- list.files(path = path, pattern = "*.xml")
# sample prepare
set = "GSE4107"
# set = "GSE23878"
esetmMatrix <- getGSEdatamatrix(set)
hsaNetDegree <- getHSADegree(path,files)
k <- 1
edgeName   <- NULL
corP.value <- NULL
for (file in files){
    fileName    <- paste(path,file, sep = "")
    gMatrixData <- getHsaGMatrix(fileName)
    gMatrix     <- gMatrixData[[2]]
    genesList   <- gMatrixData[[3]]
    esetmSam    <- getHsaSamMatrix(esetmMatrix, genesList,gMatrix)
    esetmC      <- esetmSam[[1]]
    esetmD      <- esetmSam[[2]]
    samgMatrix  <- esetmSam[[3]]
    if(sum(gMatrix) == 0){
        files[k] <- NA
        next
    }
    pC <- cor(t(esetmC))
    pT <- cor(t(esetmD))
    p.diff <- samgMatrix*(pT - pC)
    p.diff <- abs(p.diff)
    for(i in 1:nrow(gMatrix)){
        for (j in 1:nrow(gMatrix)){
            if (gMatrix[i,j] > 0){
                startName <- rownames(gMatrix)[i]
                endName <- colnames(gMatrix)[j]
                edgeName <- c(edgeName,paste(startName,endName,sep = "|"))
                if(startName %in% rownames(p.diff) && endName %in% rownames(p.diff)){
                    corP.value <- c(corP.value,p.diff[rownames(p.diff)==startName,colnames(p.diff) == endName])
                }else{
                    corP.value <- c(corP.value,0)
                }
                
            }
        }
    }
    remove(i, j, pC, pT, p.diff)
    k = k + 1
}
remove(k, fileName)
edgeP <- data.frame(edgeName, corP.value)
edgeP <- edgeP[!duplicated(edgeP[,1]),]
edgeP.sort <- edgeP[order(-(as.numeric(edgeP$corP.value))),]
edgeP.top <- edgeP.sort[1:(0.1*nrow(edgeP.sort)),1]
files.naomit <- na.omit(files)
print("step1: get edgeP and cor")


hsa.ID          <- NULL
hsa.Name        <- NULL
edge.intersect  <- NULL
line.numeber    <- NULL
P.value         <- NULL
for (file in files.naomit){
    
    fileName <- paste(path,file,sep = "")
    gMatrixData <- getHsaGMatrix(fileName)
    gMatrix <- gMatrixData[[2]]
    lineP <- NULL
    for(i in 1:nrow(gMatrix)){
        for (j in 1:nrow(gMatrix)){
            if (gMatrix[i,j] > 0){
                lineP <- c(lineP, paste(rownames(gMatrix)[i], colnames(gMatrix)[j], sep = "|"))
            }
        }
    }
    intcNumber <- length(intersect(lineP,edgeP.top))
    lineNumber <- length(intersect(lineP,edgeP.sort[,1]))
    if(lineNumber > 0){
        hsa.ID <- c(hsa.ID,file)
        hsa.Name <- c(hsa.Name,gMatrixData[[1]]) 
        edge.intersect <- c(edge.intersect,intcNumber)
        line.numeber <- c(line.numeber,lineNumber)
        PValue.value <- phyper(intcNumber-1,lineNumber,nrow(edgeP.sort)-length(intersect(lineP,edgeP.sort[,1])),length(edgeP.top),lower.tail = F)
        #PValue1.value <- (length(intersect(lineP,edgeP.top))/length(edgeP.top))/(length(lineP)/length(edgeP.sort[,1]))    
        P.value <- c(P.value,PValue.value)
    }else{
        next
    }
    remove(i, j)
}
P.just <- p.adjust(P.value,"fdr")
PResult <- data.frame(hsa.ID, hsa.Name, edge.intersect, line.numeber, P.value, P.just)
PResult.sort <- PResult[order(as.numeric(PResult$P.just)),]
rownames(PResult.sort) <- NULL


netdegree <- data.frame(hsa.Name = names(hsaNetDegree),netdegree = hsaNetDegree)
bindDegree = merge(PResult.sort,netdegree,by = 'hsa.Name',all.x=T)
bindDegree <- bindDegree[order(as.numeric(bindDegree$P.just)),]
Result.sort.degree <- subset(bindDegree ,select = c('hsa.ID', 'hsa.Name', 'edge.intersect', 'line.numeber',
                                           'P.value', 'P.just','netdegree'))

print("step2: get allp and cor")

write.csv(Result.sort.degree, file = "PResult.2.csv", quote = TRUE, na = "NA",  row.names = TRUE)

