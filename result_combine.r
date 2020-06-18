path1 <- "F:/branch2/Compare/Fisher/"
path2 <- "F:/branch2/Compare/GSA/"
path3 <- "F:/branch2/Compare/MRGSE/"
path4 <- "F:/branch2/Compare/PADOG/"
path5 <- "F:/branch2/Compare/SPIA/"
path6 <- "F:/branch2/RESULT/"
pathSum <- "F:/branch2/sum/gse/"

colnamelist <- c("hsaID","hsaName","p.just","Fisher.pNDE.far","GSA.FDR","MRGSE.p.adj","PADOG.Ppadog","SPIA.pGFWER")
datalist <- read.csv("F:/branch2/datalist.csv",header = T)
datasetnames <- datalist

for(i in 1:33){
    nameGSE <- as.character(datasetnames[[1]][i])

    nameGSE <- "GSE4107"
    hsaId <- NULL
    fileData <- read.csv(paste(path6,"result.",nameGSE,".csv",sep = ""))
    hsaResult <- fileData[which(fileData$P.just <= 0.05),c(2,3,7)]
    hsaId <- substr(hsaResult$hsa.ID,5,8)
    
    FisherData <- read.csv(paste(path1,nameGSE,".csv",sep = ""))
    fisherResult <- NULL
    for (i in 1:length(hsaId)){

        if (hsaId[i] %in% as.character(FisherData$ID)){
            FisherP <- FisherData[which(hsaId[i] == as.character(FisherData$ID)),6]
        }else{
            FisherP <- -1
        }
        fisherResult <- c(fisherResult,FisherP)
    }
    
    GSAData <- read.csv(paste(path2,nameGSE,".csv",sep = ""))
    GSAResult <- NULL
    for (i in 1:length(hsaId)){

        if (hsaId[i] %in% substr(as.character(GSAData$Gene_set_name),5,8)){
            GSAP <- GSAData[which(substr(as.character(GSAData$Gene_set_name),5,8) == hsaId[i]),6]
        }else{
            GSAP <- -1
        }
        GSAResult <- c(GSAResult,GSAP)
    }
    
    
    MRGSEData <- read.csv(paste(path3,nameGSE,".csv",sep = ""))
    MRGSEResult <- NULL
    for (i in 1:length(hsaId)){

        if (hsaId[i] %in% substr(as.character(MRGSEData$pathwayID),5,8)){
            MRGSEP <- MRGSEData[which(substr(as.character(MRGSEData$pathwayID),5,8) == hsaId[i]),4]
        }else{
            MRGSEP <- -1
        }
        MRGSEResult <- c(MRGSEResult,MRGSEP)
    }
    
    PADOGData <- read.csv(paste(path4,nameGSE,".csv",sep = ""))
    PADOGResult <- NULL
    for (i in 1:length(hsaId)){

        if (hsaId[i] %in% as.character(PADOGData$ID)){
            PADOGP <- PADOGData[which(as.character(PADOGData$ID) == hsaId[i]),8]
        }else{
            PADOGP <- -1
        }
        PADOGResult <- c(PADOGResult,PADOGP)
    }
    
    SPIAData <- read.csv(paste(path5,nameGSE,".csv",sep = ""))
    SPIAResult <- NULL
    for (i in 1:length(hsaId)){

        if (hsaId[i] %in% as.character(SPIAData$ID)){
            SPIAP <- SPIAData[which(as.character(SPIAData$ID) == hsaId[i]),11]
        }else{
            SPIAP <- -1
        }
        SPIAResult <- c(SPIAResult,SPIAP)
    }
    
    combineResult <- data.frame(hsaResult,fisherResult,GSAResult,MRGSEResult,PADOGResult,SPIAResult)
    colnames(combianResult) <- colnamelist
    write.csv(combineResult,paste(pathSum,nameGSE,".pathway.csv",sep = ""))
}
