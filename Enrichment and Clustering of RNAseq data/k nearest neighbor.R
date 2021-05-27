load("C:/Users/golde/Desktop/CSCI_5461/HW3/Data2.RData")

install.packages("class")
library("class")

train <- Data2$training_data
test <- Data2$testing_data
cl <- Data2$training_label
knn1 <- knn(train, test, cl, k = 1)
knn3 <- knn(train, test, cl, k = 3)
knn5 <- knn(train, test, cl, k = 5)

knn1accuracy <- length(grep("TRUE", knn1 == Data2$testing_label))/37
knn3accuracy <- length(grep("TRUE", knn3 == Data2$testing_label))/37
knn5accuracy <- length(grep("TRUE", knn5 == Data2$testing_label))/37

train1000 <- Data2$training_data[,1:1000]
test1000 <- Data2$testing_data[,1:1000]
knn1_top1000 <- knn(train1000, test1000, cl, k = 1)
knn3_top1000 <- knn(train1000, test1000, cl, k = 3)
knn5_top1000 <- knn(train1000, test1000, cl, k = 5)

knn1accuracy_top1000 <- length(grep("TRUE", knn1_top1000 == Data2$testing_label))/37
knn3accuracy_top1000 <- length(grep("TRUE", knn3_top1000 == Data2$testing_label))/37
knn5accuracy_top1000 <- length(grep("TRUE", knn5_top1000 == Data2$testing_label))/37

32.2 
#accuracy is better because the genes with lower P-values are removed 
#providing little to no value for sorting