## Mar, 2021
## Author: Sambhawa Priya, PhD candidate (BICB), Blekhman Lab. 

## This homework has been adapted based on the following papers:
## Code and data adapted from Zhou et al., 2019: https://www.frontiersin.org/articles/10.3389/fgene.2019.00579/full
## Dataset orginally published at:  
## Goodrich et al. 2014: http://dx.doi.org/10.1016/j.cell.2014.09.053
## Compiled by Duvallet et al.: https://www.ncbi.nlm.nih.gov/pubmed?Db=pubmed&Cmd=ShowDetailView&TermToSearch=29209090

## Preparation: Create a directory/folder named "ML_HW2" at a relevant location on your computer,
## e.g. if you have a course directory for BICB_8510, create this directory there. 
## Next, place this Rscript in the directory ML_HW1. 
## Next, create a directory within ML_HW1 called "input" and place the input files
## i.e. metadata.txt and otu_table.txt in the "input" directory. 

## Initialization
rm(list=ls()) ## Don't do this if other objects in workspace that you need later.
install.packages("rlang")
install.packages("generics")
install.packages("lifecycle")
install.packages("gower")
install.packages("caret",
                 repos = "http://cran.r-project.org", 
                 dependencies = c("Depends", "Imports", "Suggests"))

install.packages("randomForest", dependencies = c("Depends", "Suggests"))

install.packages("pROC", dependencies = c("Depends", "Suggests"))


library(caret) ## Most popular R package for machine learning
library(randomForest)
library(pROC)
#inherant
library(stringr)


## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

############### Input metadata and otu table #################
## Goodrich et al. data
metadata <- read.table(paste0(current_dir,"/input/ob_goodrich.metadata.txt"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE, check.names = F)
dim(metadata)


#Keep samples with n_sample == 0 and disease-states "OB" and "H"
metadata <- metadata[(metadata$n_sample == 0 & metadata$DiseaseState %in% c("H","OB")),]
dim(metadata)


#select one sample per individual and one individual per twin pair
metadata <- metadata[!duplicated(metadata$familyid),]
dim(metadata)


## drop unused level from metadata
metadata <- droplevels(metadata)

## How many samples per disease state?
table(metadata$DiseaseState)
# #H  OB 
# 279 135
## This will be our classification problem!

## Please be patient. This will take sometime (1-2 mins).
otu <- read.table(paste0(current_dir,"/input/ob_goodrich.otu_table.100.denovo.rdp_assigned"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE, check.names = F)
dim(otu)

## identify common samples between metadata and otu table
common_samples <- intersect(rownames(metadata), colnames(otu))
length(common_samples)

## Ensure order of samples in metadata and otu table are identical
otu <- otu[,common_samples]
dim(otu) 

metadata <- metadata[common_samples,]
dim(metadata) 

all(rownames(metadata) == colnames(otu)) 
# should be TRUE

##################### Do Preprocessing (see specific instructions where provided) ########################

## Q1. Remove OTUs with fewer than 10 reads (10 points)
otu <- otu[which(rowSums(otu)>=10),]
dim(otu)
## Q2. Remove OTUs which were present in fewer than 5% of samples (10 points)
otu <- otu[which(rowSums(otu>0) >= ncol(otu)*.05),]
dim(otu)

## Binning
# collapse OTUs to genus level by summing their abundance counts
## Split taxa names by ";" into 8 parts
tax_table <- str_split_fixed(rownames(otu),";",n=8)
## Append genus name (6th column in genus) to the otu table
otu <- data.frame(genus=tax_table[,6],otu)
## Filter out otus where genus is not characterized. 
otu <- otu[!(otu$genus=="g__"),]
dim(otu)

## summarize otu table by genus label
otu <- aggregate(otu[,-1], by=list(otu$genus), sum)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
dim(otu)


#calculate relative abundance of each genus by dividing its value by the total reads per sample
otu <- sweep(otu,2,colSums(otu),"/")

## Prepare for training
x <- data.matrix(otu)
## transpose to make rows are samples and feature (i.e. genera) as columns. 
x <- t(x)

## relevel disease state to allow "OB" to be case for ML model downstream
levels(metadata$DiseaseState) #[1] [1] "H"  "OB"
metadata$DiseaseState <- relevel( metadata$DiseaseState, "OB")
levels(metadata$DiseaseState) #[1] "OB" "H" 


################ Train and Test ###############
## Q3. Split data by 90% training and 10% test, and report the output of training (best mtry),
## and the output of testing (confusion matrix, sensitivity, specificity, precision, AUC).  
## Also show the ROC curve. (30 points)
set.seed(1000)
train_index <- createDataPartition(metadata$DiseaseState, ## outcome
                                   p = 0.9, ## percentage of training samples
                                   list = FALSE ## show subsamples as matrix, not list
                                   # times = 10 ## This will create 10 different 80% subsamples
) 
View(train_index)

x.train <- x[train_index,] 
y.train <- metadata$DiseaseState[train_index]

x.test <- x[-train_index,]
y.test <- metadata$DiseaseState[-train_index]

train_control <- trainControl(
  method = "cv",
  number = 5, ## also try 10
  summaryFunction=twoClassSummary, # computes area under the ROC curve
  classProbs = TRUE ## required for scoring models using ROC
)


set.seed(1000)
rf_train <- train( x = x.train, y = as.factor(y.train),
                   method='rf',
                   metric="ROC", ## default accuracy
                   trControl = train_control)
rf_train


## mtry: Number of variables randomly sampled as candidates at each split.
rf_train$resample

## We can modify tuning grid using tuneGrid param in train(). 
rf_test <- predict(rf_train, x.test) 
rf_test

# compare predicted outcome and true outcome
conf_matrix <- confusionMatrix(rf_test, y.test)
conf_matrix$table


## Compute precision
conf_matrix$byClass

## Can also spit out probability instead of predicted class
rf_test <- predict(rf_train, x.test, type = "prob")
rf_test
rf_test <- rf_test[,1]

############ Plot ROC curve ############
## ROC curve
rf <- roc(y.test,rf_test) ## pROC package
auc <- rf$auc
auc


## Plot ROC curve
pdf(paste0(current_dir,"/output/ob_goodrich90.pdf"))
plot(rf, col="blue",legacy.axes = TRUE)
dev.off()

## Q4. Split data by 70% training and 30% test, and report the output of training (best mtry),
## and the output of testing (confusion matrix, sensitivity, specificity, precision, AUC).  
## Also show the ROC curve. Did the AUC increase or decrease compared to result from Q3? (30 points) 
#####AUC decreased with %70 training data vs %90 training data#####
set.seed(1000)
train_index <- createDataPartition(metadata$DiseaseState, ## outcome
                                   p = 0.7, ## percentage of training samples
                                   list = FALSE ## show subsamples as matrix, not list
                                   # times = 10 ## This will create 10 different 80% subsamples
) 
View(train_index)

x.train <- x[train_index,] 
y.train <- metadata$DiseaseState[train_index]

x.test <- x[-train_index,]
y.test <- metadata$DiseaseState[-train_index]

train_control <- trainControl(
  method = "cv",
  number = 5, ## also try 10
  summaryFunction=twoClassSummary, # computes area under the ROC curve
  classProbs = TRUE ## required for scoring models using ROC
)


set.seed(1000)
rf_train <- train( x = x.train, y = as.factor(y.train),
                   method='rf',
                   metric="ROC", ## default accuracy
                   trControl = train_control)
rf_train


## mtry: Number of variables randomly sampled as candidates at each split.
rf_train$resample

## We can modify tuning grid using tuneGrid param in train(). 
rf_test <- predict(rf_train, x.test) 
rf_test

# compare predicted outcome and true outcome
conf_matrix <- confusionMatrix(rf_test, y.test)
conf_matrix$table


## Compute precision
conf_matrix$byClass

## Can also spit out probability instead of predicted class
rf_test <- predict(rf_train, x.test, type = "prob")
rf_test
rf_test <- rf_test[,1]

############ Plot ROC curve ############
## ROC curve
rf <- roc(y.test,rf_test) ## pROC package
auc <- rf$auc
auc


## Plot ROC curve
pdf(paste0(current_dir,"/output/ob_goodrich70.pdf"))
plot(rf, col="blue",legacy.axes = TRUE)
dev.off()


########## List feature importance in random forest ###########
## Q5. Show top 10 features ranked by their importance corresponding to Q4 output. (20 points) 
## Hint: Use varImp() as shown in class. Sort by importance and pick top-10 features
impVars <- varImp(rf_train)
ImpMeasure <- data.frame(order(impVars$importance))
ImpMeasureordered <- ImpMeasure[order(-ImpMeasure$Overall),,drop = FALSE]
ImpMeasuretop10 <- ImpMeasureordered[1:10,,drop = FALSE]
ImpMeasuretop10
