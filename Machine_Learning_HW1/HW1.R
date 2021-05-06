## Mar, 2021
## Author: Sambhawa Priya, PhD candidate (BICB), Blekhman Lab. 

## BICB 8510
## HW1 for microbiome analysis 
## Based on metadata and otu tables taken from https://github.com/LangilleLab

## Preparation: Create a directory/folder named "ML_HW1" at a relevant location on your computer,
## e.g. if you have a course directory for BICB_8510, create this directory there. 
## Next, place this Rscript in the directory ML_HW1. 
## Next, create a directory within ML_HW1 called "input" and place the input files
## i.e. metadata.txt and otu_table.txt in the "input" directory. 


##install libraries
install.packages("vegan") ## for distance computation for microbiome data
## vignette: https://cran.r-project.org/web/packages/vegan/vegan.pdf
##import libraries
library(vegan) 

## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## Input metadata
## This will be located in directory "ML_HW1/input/", if you followed the instruction above. 
metadata <- read.table(paste0(current_dir,"/input/metadata.txt"),sep="\t",header=T,row.names=1,stringsAsFactors=TRUE, comment.char="")
## Input otu table. Also should be located at "ML_HW1/input/" directory
otu <- read.table(paste0(current_dir,"/input/otu_table.txt"),sep="\t",header=T,row.names=1,stringsAsFactors=FALSE, comment.char="")


####### Q1: How many rows and columns do the metadata and otu have? (10 points) 
print("metadata: 2cols,40rows") 
print("otu: 40cols, 1000rows")

###### Q2. In metadata table, how many samples are categorized as inflammed? How many samples are control? Show script. (10 points)
## Hint: Use column "state" in the metadata table.
print("number inflamed:")
print(nrow(subset(metadata,state == "inflamed")))

print("number control:")
print(nrow(subset(metadata,state == "control")))

## Transpose otu table to make samples as rows and otus as columns.
otu <- as.data.frame(t(otu))
dim(otu)
# [1]   40 1000

####### Q3: Identify common samples between metadata table and otu table. Show script. (10 points)
common_samples <- intersect(row.names(metadata),rownames(otu))

## Ensure order of samples in metadata and otu table are identical
otu <- otu[common_samples,]
metadata <- metadata[common_samples,]

## Explore data distribution
## Pick an OTU (OTU_77)
otu_random <- otu[,grep("OTU_77$", colnames(otu))]

####### Q4: Plot the distribution of this OTU using a histogram. Show script and upload plot. (20 points)
hist(otu_random, br =110)

