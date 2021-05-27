load("C:/Users/golde/Desktop/CSCI_5461/HW3/Data1.RData")
#print top 200 genes from Data1 for GO analysis
paste(as.character(head(Data1$Gene_Name, 200)), sep="' '", collapse=", ")


#sometimes kmeans has to be run more than once due to bug in package, sometimes returns "Quick-TRANSfer stage steps exceeded maximum (= 1006550)" error
kmeans10 <- kmeans(Data1$SeqData, centers = 10, iter.max = 100, algorithm = "Hartigan-Wong")
kmeans20 <- kmeans(Data1$SeqData, centers = 20, iter.max = 100, algorithm = "Hartigan-Wong")
kmeans50 <- kmeans(Data1$SeqData, centers = 50, iter.max = 100, algorithm = "Hartigan-Wong")



hist(kmeans10$size, breaks = 100)
hist(kmeans20$size, breaks = 100)
hist(kmeans50$size, breaks = 100)


cluster_13_genes = character(0)
for (i in 1:length(kmeans20$cluster)){
  if (kmeans20$cluster[i] == 13){
    cluster_13_genes <- append(cluster_13_genes, Data1$Gene_Name[i])
  }
}

#print genes from cluster_13_genes for GO analysis
paste(as.character(head(cluster_13_genes, 200)), sep="' '", collapse=", ")
