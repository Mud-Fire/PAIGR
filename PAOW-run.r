library(SPIA)
library(KEGGgraph)
library(igraph)
library(limma)
library(PADOG)
library(ROntoTools)
library(CePa)
library(GSA)
library(annotate)
library(EnrichmentBrowser)
path  <- "F:/branch2/kegg/"
paths <- list.files(path = path, pattern = "*.xml")
# kgml.path=system.file("extdata/keggxml/hsa",package="SPIA")
# mdir=kgml.path
# paths<-dir(mdir,pattern=".xml")

#mapkG <- mergepathway(paths,switch = TRUE)
#wfg <- frequencyw(kgml.path)

datalist <- read.csv("F:/branch2/datalist.csv",header = T)
datasetnames <- datalist
set <- NULL

targetgs <- NULL
for(i in 1:33){
  set = as.character(datasetnames[[1]][i])
  #set = "GSE4107"
  print(set)
  data(list = set, package = "KEGGdzPathwaysGEO")
  #write a function to extract required info into a list
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
  dlist = getdataaslist(set)
  
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

  stopifnot(class(esetm) == "matrix")
  stopifnot(all(dim(esetm) > 4))

  stopifnot(class(group) %in% c("factor", "character"))
  stopifnot(length(group) == dim(esetm)[2])
  stopifnot(all(group %in% c("c", "d")))
  stopifnot(all(table(group) > 2))
  if (paired) {
    stopifnot(length(block) == length(group))
    stopifnot(all(table(block) == 2))
  }

  aT1 = filteranot(esetm, group, paired, block, annotation)
  esetm = esetm[rownames(esetm) %in% aT1$ID, ]
  rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]

  topSigNum = dim(esetm)[1]

  if (paired) {
    design <- model.matrix(~0 + G + block)
    colnames(design) <- substr(colnames(design), 2, 100)
  }
  if (!paired) {
    design <- model.matrix(~0 + G)
    colnames(design) <- levels(G)
  }
  fit <- lmFit(esetm, design)
  cont.matrix <- makeContrasts(contrasts = "d-c", levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  aT1 <- topTable(fit2, coef = 1, number = topSigNum)
  aT1$ID <- rownames(aT1)

  ########Ronto_tool#########
  # kpg <- NULL
  # kpg_name <- NULL
  # for (i in 1:length(paths)){
  #   fileName    <- paste(path,paths[i], sep = "")
  #   gMatrixData <- getHsaGMatrix(fileName)
  #   hsa.Data <- try(parseKGML(fileName),TRUE)
  #   hsa.Name <- hsa.Data@pathwayInfo@name
  #   hsa.pathway <- KEGGpathway2Graph(hsa.Data, expandGenes = T)
  #   kpg_name <- c(kpg_name,hsa.Name)
  #   kpg <- c(kpg,hsa.pathway)
  # }
  # names(kpg) <- kpg_name
  # #kpg1 <- keggPathwayGraphs("hsa", verbose = FALSE)
  # top <- aT1
  # top <- top[top$adj.P.Val <= 0.05,]
  # fc <- top$logFC
  # if(length(fc)<200){
  #     top <- aT1[which(aT1$P.Value<0.05),]
  #     top <- top[which(top$logFC>1.5),]
  #     fc = top$logFC
  # }
  # if(length(fc)<200){
  #     top <- aT1[1:(nrow(aT1)*0.01),]
  #     fc = top$logFC
  # }
  # # top <- aT1
  # # fc <- top$logFC[top$adj.P.Val <= 0.05]
  # pv <- top$P.Value[top$adj.P.Val <= 0.05]
  # names(fc) <- paste("hsa:",top$ID[top$adj.P.Val <= 0.05],sep="")
  # if(length(pv)<200){
  #     top <- aT1[which(aT1$P.Value<0.05),]
  #     top <- top[which(top$logFC>1.5),]
  #     pv = top$P.Value
  # }
  # if(length(pv)<200){
  #     top <- aT1[1:(nrow(aT1)*0.01),]
  #     pv = top$P.Value
  # }
  # 
  # names(pv) <- paste("hsa:",top$ID[top$adj.P.Val <= 0.05],sep="")
  # fcAll <- top$logFC
  # names(fcAll) <- paste("hsa:",top$ID,sep="")
  # pvAll <- top$P.Value
  # names(pvAll) <- paste("hsa:",top$ID,sep="")
  # ref <- paste("hsa:",top$ID,sep="")
  # kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
  # peRes <- pe(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = FALSE)
  # res_Ronto <- summary(peRes)
  # write.csv(res_Ronto,paste("F:/branch2/Compare/Ronto/",set,".csv",sep=""))

  #######SPIA########
  #makeSPIAdata(kgml.path=system.file("extdata/keggxml/hsa",package="SPIA"),organism="hsa",out.path="./")
  tg1 <- aT1[aT1$adj.P.Val<0.05,]
  DE = tg1$logFC
  names(DE) <- as.vector(tg1$ID)

  if(length(DE)<200){
      top <- aT1[which(aT1$P.Value<0.05),]
      top <- top[which(top$logFC>1.5),]
      DE = top$logFC
      names(DE) <- as.vector(top$ID)
  }
  if(length(DE)<200){
      top <- aT1[1:(nrow(aT1)*0.01),]
      DE = top$logFC
      names(DE) <- as.vector(top$ID)
  }

  ALL = aT1$ID
  res <- spia(de=DE,all=ALL,organism="hsa",data.dir=NULL,pathids=NULL,nB=2000,plots=FALSE,verbose=TRUE,beta=NULL,combine="fisher")
  write.csv(res,paste("F:/branch2/Compare/SPIA/",set,".csv",sep=""))

  # #######Fisher#####
  pNDE.fdr = p.adjust(res$pNDE,"fdr")
  res <- data.frame(Name=res$Name,ID=res$ID,pSize=res$pSize,NDE=res$NDE,pNDE=res$pNDE,pNDE.fdr=pNDE.fdr,KEGGLINK=res$KEGGLINK,stringsAsFactors=FALSE)
  res <- res[order(res$pNDE),]
  rownames(res) <- NULL
  write.csv(res,paste("F:/branch2/Compare/Fisher/",set,".csv",sep=""))
  #
  # #####GSA#######
  # library(KEGG.db)
  # # pw2id = as.list(KEGGPATHID2EXTID)
  # # gslist = pw2id[grep("hsa", names(pw2id))]
  # 
  # gslist <- NULL
  # for (i in 1:length(paths)){
  #     fileName <- paste(path,paths[i], sep = "")
  #     gMatrixData <- getHsaGMatrix(fileName)
  #     hsa.Name <- hsa.Data@pathwayInfo@title
  #     hsa.Data <- try(parseKGML(fileName),TRUE)
  #     hsa.pathway <- KEGGpathway2Graph(hsa.Data, expandGenes = T)
  #     hsa.Net <- igraph.from.graphNEL(hsa.pathway)
  #     hsa.Matrix <- as_adjacency_matrix(hsa.Net)
  #     gMatrix <- as.matrix(hsa.Matrix)
  #     hsa.GeneName <- rownames(gMatrix)
  #     genesList <- unlist(strsplit(hsa.GeneName,":"))[c(F,T)]
  #     genesList <- list(genesList)
  #     names(genesList) <- unlist(strsplit(paths[i],"[.]"))[c(T,F)]
  #     gslist <- c(gslist,genesList)
  # }
  # 
  # genesets <- gslist
  # y<-dlist$ano$Group
  # names(y)<-dlist$ano$Sample
  # y[which(y=="c")]<-1
  # y[which(y=="d")]<-2
  # x <- esetm
  # y <- y[colnames(x)]
  # y <- as.numeric(y)
  # genenames <- rownames(esetm)
  # if(paired){
  #   types <- "Two class paired"
  # }
  # if(!paired){
  #   types <- "Two class unpaired"
  # }
  # GSA.obj <- GSA(x,y, genenames=genenames, genesets=genesets, nperms=100)
  # # GSA.obj <- GSA(x,y, genenames=genenames, genesets=genesets, resp.type=types, nperms=100)
  # res <- rbind(GSA.listsets(GSA.obj, geneset.names=names(genesets),FDRcut=1)$negative,GSA.listsets(GSA.obj, geneset.names=names(genesets),FDRcut=1)$positive)
  # res <- res[order(res[,4]),]
  # write.csv(res,paste("F:/branch2/Compare/GSA/",set,".csv",sep=""))
  # 
  # ######MRGSE########
  # tv<-aT1$t
  # names(tv)<-aT1$ID
  # MRGSE <- function(gs){
  #   gsnodes <- intersect(gs,names(tv))
  #   geneSetTest(gsnodes,tv)
  # }
  # res <- lapply(gslist,MRGSE)
  # res <- do.call('c',res)
  # p.adj <- p.adjust(res,"fdr")
  # names(p.adj) <-NULL
  # p.value <- res
  # names(p.value) <- NULL
  # res <- data.frame(pathwayID=names(res),p.value=p.value,p.adj=p.adj)
  # res <- res[order(res$p.value),]
  # rownames(res) <- NULL
  # write.csv(res,paste("F:/branch2/Compare/MRGSE/",set,".csv",sep=""))
  # 
  # ######GSEA#########
  # # kegg.gs <- get.kegg.genesets("hsa")
  # # kegg.gs <- getGenesets("hsa",db = "kegg")
  # # kegg.gs<-kegg.gs[97:326]
  # 
  # # kegg.gs <- NULL
  # # for (i in 1:length(paths)){
  # #     print(i)
  # #     fileName <- paste(path,paths[i], sep = "")
  # #     gMatrixData <- getHsaGMatrix(fileName)
  # #     hsa.Data <- try(parseKGML(fileName),TRUE)
  # #     hsa.Name <- hsa.Data@pathwayInfo@title
  # #     hsa.pathway <- KEGGpathway2Graph(hsa.Data, expandGenes = T)
  # #     hsa.Net <- igraph.from.graphNEL(hsa.pathway)
  # #     hsa.Matrix <- as_adjacency_matrix(hsa.Net)
  # #     gMatrix <- as.matrix(hsa.Matrix)
  # #     hsa.GeneName <- rownames(gMatrix)
  # #     genesList <- unlist(strsplit(hsa.GeneName,":"))[c(F,T)]
  # #     genesList <- list(genesList)
  # #     names(genesList) <- unlist(strsplit(paths[i],"[.]"))[c(T,F)]
  # #     kegg.gs <- c(kegg.gs,genesList)
  # # }
  # # dlist = getdataaslist(set)
  # # group = dlist$ano$Group
  # # sets <- get(set)
  # # pData(sets)$Group <- ifelse(group == "d", 1, 0)
  # # eset <- probe2gene(sets)
  # # eset<- normalize(eset, norm.method="quantile")
  # # #table(eset$Group)
  # # eset<- deAna(eset)
  # # sbea.res <- sbea(method="gsea", eset=eset, gs=kegg.gs,alpha=1.001)
  # # res<-gs.ranking(sbea.res)
  # # write.csv(res,paste("F:/branch2/Compare/GSEA/",set,".csv",sep=""))
  # 
  # 
  # ######PADOG########
  # myr = padog(
  #   esetm = dlist$dat.m,
  #   group = dlist$ano$Group,
  #   paired = dlist$design == "Paired",
  #   block = dlist$ano$Block,
  #   targetgs = dlist$targetGeneSets,
  #   # targetgs = "05012",
  #   annotation = dlist$annotation,
  #   gslist = "KEGG.db",
  #   organism = "hsa",
  #   verbose = FALSE,
  #   Nmin = 3,
  #   NI = 50,
  #   plots = TRUE,
  #   dseed = 1
  # )
  # write.csv(myr,paste("F:/branch2/Compare/PADOG/",set,".csv",sep=""))
  #targetgs <- c(targetgs,dlist$targetGeneSets)
}

datasets<-data.frame(datasetnames,targetpathway=targetgs)
datasets$targetpathway[33] <- "05012"

taget <- NULL
rank <- NULL
path  <- "F:/branch2/Compare/Fisher/"
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- substr(as.character(datasets$targetpathway[i]),2,5)
    print(hsaId)
    if (hsaId %in% readlist$ID){
        ll <- readlist[which(readlist$ID == hsaId),]
        aa <- floor(order(readlist$pNDE.fdr))[which(readlist$ID == hsaId)]
    }
    else{
        ll <- ll
        aa <- aa
    }
    taget <- rbind(taget,ll)
    rank <- c(rank,aa)
}
fisherset<-data.frame(datasetnames,taget,rank)
write.csv(fisherset,paste("F:/branch2/Fishersum",".csv",sep=""))


taget <- NULL
path  <- "F:/branch2/Compare/GSA/"
rank <- NULL
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- paste("hsa",as.character(datasets$targetpathway[i]),sep = "")
    if (hsaId %in% readlist$Gene_set_name){
        print("i")
        ll <- readlist[which(readlist$Gene_set_name == hsaId),]
        aa <- floor(order(readlist$FDR))[which(readlist$Gene_set_name == hsaId)]
    }else{
        ll <- ll
        aa <- aa
    }
    taget <- rbind(taget,ll)
    rank <- c(rank,aa)
}
GSAset<-data.frame(datasetnames,taget,rank)
write.csv(GSAset,paste("F:/branch2/GSAsum",".csv",sep=""))


taget <- NULL
rank <- NULL
path  <- "F:/branch2/Compare/MRGSE/"
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- paste("hsa",as.character(datasets$targetpathway[i]),sep = "")
    if (hsaId %in% readlist$pathwayID){
        print("i")
        ll <- readlist[which(readlist$pathwayID == hsaId),]
        aa <- floor(order(readlist$p.adj))[which(readlist$pathwayID == hsaId)]
    }else{
        aa <- aa
        ll <- ll
    }
    taget <- rbind(taget,ll)
    rank <- c(rank,aa)
}
MRGSEset<-data.frame(datasetnames,taget,rank)
write.csv(MRGSEset,paste("F:/branch2/MRGSEsum",".csv",sep=""))


rank <- NULL
taget <- NULL
path  <- "F:/branch2/Compare/Ronto/"
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- substr(as.character(datasets$targetpathway[i]),2,5)
    if (hsaId %in% readlist$ID){
        ll <- readlist[which(readlist$ID == hsaId),]
        aa <- floor(order(readlist$p.adj))[which(readlist$ID == hsaId)]
    }else{
        print("i")
        ll <- ll
    }
    taget <- rbind(taget,ll)
}
Rontoset<-data.frame(datasetnames,taget)
write.csv(Rontoset,paste("F:/branch2/Rontosum",".csv",sep=""))

rank <- NULL
taget <- NULL
path  <- "F:/branch2/Compare/SPIA/"
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- substr(as.character(datasets$targetpathway[i]),2,5)
    if (hsaId %in% readlist$ID){
        ll <- readlist[which(readlist$ID == hsaId),]
        aa <- floor(order(readlist$pGFdr))[which(readlist$ID == hsaId)]
    }else{
        print(i)
        ll <- ll
    }
    taget <- rbind(taget,ll)
    rank <- c(rank,aa)
}
SPIAset<-data.frame(datasetnames,taget,rank)
write.csv(SPIAset,paste("F:/branch2/SPIAsum",".csv",sep=""))

rank <- NULL
taget <- NULL
path  <- "F:/branch2/Compare/PADOG/"
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- substr(as.character(datasets$targetpathway[i]),2,5)
    if (hsaId %in% readlist$ID){
        ll <- readlist[which(readlist$ID == hsaId),]
        aa <- floor(order(readlist$Ppadog))[which(readlist$ID == hsaId)]
    }else{
        print(i)
        ll <- ll
    }
    taget <- rbind(taget,ll)
    rank <- c(rank,aa)
}
PADOGset<-data.frame(datasetnames,taget,rank)
write.csv(PADOGset,paste("F:/branch2/PADOGsum",".csv",sep=""))

rank <- NULL
taget <- NULL
ll <- "a"
path  <- "F:/branch2/Compare/Ronto/"
for(i in 1:length(datasets$Dataset)){
    fileName <- paste(path,as.character(datasets$Dataset[i]),".csv", sep = "")
    readlist <- read.csv(fileName,header = T)
    hsaId <- paste("path:hsa",as.character(datasets$targetpathway[i]),sep = "")
    if (hsaId %in% readlist$X){
        ll <- readlist[which(readlist$X == hsaId),]
        aa <- floor(order(readlist$pComb.fdr))[which(readlist$X == hsaId)]
    }else{
        print(i)
        ll <- ll
    }
    taget <- rbind(taget,ll)
    rank <- c(rank,aa)
}
Rontoset<-data.frame(datasetnames,taget,rank)
write.csv(Rontoset,paste("F:/branch2/Rontosum",".csv",sep=""))

for (file in files){
    fileName <- paste(path,file, sep = "")
    datalist <- read.csv(fileName,header = T)
    
        
}




