load("C:/Users/golde/Desktop/CSCI_5461/HW3/Data2.RData")

install.packages("e1071")
library("e1071")

cl <- Data2$training_label
train1000 <- Data2$training_data[,1:1000]
test1000 <- Data2$testing_data[,1:1000]

svm.model <- svm(cl~., data = train1000,kernel = "linear", type = "C-classification", scale = FALSE)
summary(svm.model)

svm_predict <- predict(svm.model, test1000)

#make variable for classification accuracy
svmaccuracy_top1000 <- length(grep("TRUE", svm_predict == Data2$testing_label))/37