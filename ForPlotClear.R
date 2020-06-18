
network <- read.table("F:/branch2/netnode_10.csv",header = T,sep = ",")
Cyt_result1 <- read.table("F:/branch2/netnode_10.csv default node.csv",header = T,sep = ",")

indexAndDegree <- NULL
index_degree <- paste(Cyt_result1$Index , Cyt_result1$Degree , sep = ":")
indexAndDegree <- c(indexAndDegree,index_degree )

Cyt_result2 <- cbind(Cyt_result1,indexAndDegree)

indexNetwork <- data.frame("node1" = character(length(network[,1])), "node2" = character(length(network[,1])),stringsAsFactors=FALSE)
for (i in 1:length(Cyt_result2[,1])){
    lable <- as.character(Cyt_result2[i,1])
    new_value <- as.character(Cyt_result2[i,4])
    indexNetwork[which(as.character(network[,1]) == lable),1] = new_value
    indexNetwork[which(as.character(network[,2]) == lable),2] = new_value
}

write.csv(indexNetwork,"F:/branch2/indexNetwork.csv")