setwd("C:/Users/golde/Desktop/CSCI_5461/HW2")

###pull in data###
ArrayData<-read.delim("ArrayData.txt",header=T,row.names = 1,sep="\t")
SeqData<-read.delim("SeqData.txt",header=T,row.names = 1,sep="\t")

ArrayData <- ArrayData[, !apply(ArrayData == 0, 2, all)]
SeqData <- SeqData[, !apply(SeqData == 0, 2, all)]
ArrayData <- na.omit(ArrayData)
ArrayData <- na.omit(ArrayData)

ArrayDatagrp1 <- ArrayData[ArrayData$group_ID == 1,]
ArrayDatagrp2 <- ArrayData[ArrayData$group_ID == 2,]
SeqDatagrp1<-SeqData[SeqData$group_ID == 1,]
SeqDatagrp2<-SeqData[SeqData$group_ID == 2,]
ArrayDatarows<-nrow(ArrayData)
SeqDatarows<-nrow(SeqData)
ArrayDatacols<-ncol(ArrayData)
SeqDatacols<-ncol(SeqData)

####
#ttest array
for (i in c(1:ncol(ArrayData)) ) {
  if (i == 1) {ArrayData_ttest <- "NA"}
  if (i > 1) {ArrayData_ttest <- append(ArrayData_ttest,t.test(ArrayDatagrp1[i],ArrayDatagrp2[i])$p.value)}
}
#ttest seq
for (i in c(1:ncol(SeqData)) ) {
  if (i == 1) {SeqData_ttest <- "NA"}
  if (i > 1) {SeqData_ttest <- append(SeqData_ttest,t.test(SeqDatagrp1[i],SeqDatagrp2[i])$p.value)}
}

#ttest array bon
for (i in c(1:ncol(ArrayData)) ) {
  if (i == 1) {ArrayData_ttest_bon <- "NA"}
  if (i > 1) {ArrayData_ttest_bon <- append(ArrayData_ttest_bon,t.test(ArrayDatagrp1[i],ArrayDatagrp2[i], p.adjust.method = "bonferroni")$p.value)}
}
#ttest seq bon
for (i in c(1:ncol(SeqData)) ) {
  if (i == 1) {SeqData_ttest_bon <- "NA"}
  if (i > 1) {SeqData_ttest_bon <- append(SeqData_ttest_bon,t.test(SeqDatagrp1[i],SeqDatagrp2[i], p.adjust.method = "bonferroni")$p.value)}
}

#ttest array fdr
for (i in c(1:ncol(ArrayData)) ) {
  if (i == 1) {ArrayData_ttest_fdr <- "NA"}
  if (i > 1) {ArrayData_ttest_fdr <- append(ArrayData_ttest_fdr,t.test(ArrayDatagrp1[i],ArrayDatagrp2[i], p.adjust.method = "fdr")$p.value)}
}
#ttest seq fdr
for (i in c(1:ncol(SeqData)) ) {
  if (i == 1) {SeqData_ttest_fdr <- "NA"}
  if (i > 1) {SeqData_ttest_fdr <- append(SeqData_ttest_fdr,t.test(SeqDatagrp1[i],SeqDatagrp2[i], p.adjust.method = "fdr")$p.value)}
}

#Wilcoxon array
for (i in c(1:ncol(ArrayData)) ) {
  if (i == 1) {ArrayDataWilcoxon <- "NA"}
  if (i > 1) {ArrayDataWilcoxon <- append(ArrayDataWilcoxon,wilcox.test(as.matrix(ArrayDatagrp1[i]),as.matrix(ArrayDatagrp2[i]))$p.value)}
}
#Wilcoxon seq
for (i in c(1:ncol(SeqData)) ) {
  if (i == 1) {SeqDataWilcoxon <- "NA"}
  if (i > 1) {SeqDataWilcoxon <- append(SeqDataWilcoxon,wilcox.test(as.matrix(SeqDatagrp1[i]),as.matrix(SeqDatagrp2[i]))$p.value)}
}


##array ttest
#append array ttest
ArrayData_ttest <- as.numeric(ArrayData_ttest)
ArrayData <- rbind(ArrayData,ArrayData_ttest)
rownames(ArrayData)[rownames(ArrayData) == "446"] = "ttest.p"

#order by array ttest
ttestorderedArrayData <- ArrayData[,order(ArrayData["ttest.p",])]
top10ttestArrayData <- colnames(ttestorderedArrayData[1:10])

#enumerate number of array data that make ttest
ttestorderedArrayData05.c.o. <- ttestorderedArrayData["ttest.p",]<0.05
ttestorderedArrayDataCOenumerate <-rowSums(ttestorderedArrayData05.c.o. == "TRUE", na.rm = TRUE)

##array wilcoxon
#append array wilcoxon
ArrayDataWilcoxon <- as.numeric(ArrayDataWilcoxon)
ArrayData <- rbind(ArrayData,ArrayDataWilcoxon)
rownames(ArrayData)[rownames(ArrayData) == "447"] = "Wilcoxon.p"

#order by array ttest
wilcoxonorderedArrayData <- ArrayData[,order(ArrayData["Wilcoxon.p",])]
top10wilcoxonArrayData <- colnames(wilcoxonorderedArrayData[1:10])

#enumerate number of seq data that make wilcoxon cut off
wilcoxonorderedArrayData05.c.o. <- wilcoxonorderedArrayData["Wilcoxon.p",]<0.05
wilcoxonorderedArrayDataCOenumerate <-rowSums(wilcoxonorderedArrayData05.c.o. == "TRUE", na.rm = TRUE)

##seq ttest
#append seq ttest
SeqData_ttest <- as.numeric(SeqData_ttest)
SeqData <- rbind(SeqData,SeqData_ttest)
rownames(SeqData)[rownames(SeqData) == "309"] = "ttest.p"

#order by seq ttest
ttestorderedSeqData <- SeqData[,order(SeqData["ttest.p",])]
top10ttestSeqData <- colnames(ttestorderedSeqData[1:10])

#enumerate number of seq data that make wilcoxon cut off
ttestorderedSeqData05.c.o. <- ttestorderedSeqData["ttest.p",]<0.05
ttestSeqDataCOenumerate <-rowSums(ttestorderedSeqData05.c.o. == "TRUE", na.rm = TRUE)

##seq wilcoxon
#append seq wilcoxon
SeqDataWilcoxon <- as.numeric(SeqDataWilcoxon)
SeqData <- rbind(SeqData,SeqDataWilcoxon)
rownames(SeqData)[rownames(SeqData) == "310"] = "Wilcoxon.p"

#order by seq wilcoxon
wilcoxonorderedSeqData <- SeqData[,order(SeqData["Wilcoxon.p",])]
top10wilcoxonSeqData <- colnames(wilcoxonorderedSeqData[1:10])

#enumerate number of seq data that make wilcoxon cut off
wilcoxonorderedSeqData05.c.o. <- wilcoxonorderedSeqData["Wilcoxon.p",]<0.05
wilcoxonSeqDataCOenumerate <-rowSums(wilcoxonorderedSeqData05.c.o. == "TRUE", na.rm = TRUE)

##make histograms
#Histogram of array data t-test, saved manually using Rstudio export function
hist(ArrayData_ttest,breaks=100)

#Histogram of seq data t-test, saved manually using Rstudio export function
hist(SeqData_ttest,breaks=100)

#Histogram of array data Wilcoxon, saved manually using Rstudio export function
hist(ArrayDataWilcoxon,breaks=100)

#Histogram of seq data Wilcoxon, saved manually using Rstudio export function
hist(SeqDataWilcoxon,breaks=100)


####3.1
##array ttest bonferroni
#append array ttest bonferroni
ArrayData_ttest_bon <- as.numeric(ArrayData_ttest_bon)
ArrayData <- rbind(ArrayData,ArrayData_ttest_bon)
rownames(ArrayData)[rownames(ArrayData) == "448"] = "ttest_bon"

#order array ttest bonferroni
ttest_bonorderedArrayData <- ArrayData[,order(ArrayData["ttest_bon",])]

#enumerate number of seq data that make make bonferroni cut off
ttest_bonorderedArrayData05.c.o. <- ttest_bonorderedArrayData["ttest_bon",]<0.05
ttest_bonArrayDataCOenumerate <-rowSums(ttest_bonorderedArrayData05.c.o. == "TRUE", na.rm = TRUE)

##seq ttest bonferroni
#append seq ttest bonferroni
SeqData_ttest_bon <- as.numeric(SeqData_ttest_bon)
SeqData <- rbind(SeqData,SeqData_ttest_bon)
rownames(SeqData)[rownames(SeqData) == "311"] = "ttest_bon"

#order seq ttest bonferroni
ttest_bonorderedSeqData <- SeqData[,order(SeqData["ttest_bon",])]

#enumerate number of seq data that make bonferroni cut off
ttest_bonorderedSeqData05.c.o. <- ttest_bonorderedSeqData["ttest_bon",]<0.05
ttest_bonSeqDataCOenumerate <-rowSums(ttest_bonorderedSeqData05.c.o. == "TRUE", na.rm = TRUE)
