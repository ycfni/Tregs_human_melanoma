######################## fit a decision tree ############################

### The following code constructs decision tree models to identify combinations of features to prospectively classify cells and corresponds with data presented in Extended Data Figure 9 [Lu et al. Nat Immunol 2025]. The code reads in single-cell RNA expression matrix along with a binarized categorization for features of interest. The sample is then randomly divided into a training and testing cohort. Decision trees are then fit to the training dataset. Prediction statistics are then calculated on the testing dataset. 

#########################################################################

# import the library for building a decision tree
library(rpart)
library(rpart.plot)

# load the data from the path
load("/path/to/your/file.RData")

# filter data to contain the outcome of interest and all predictors. Examples provided included whether the sample is tumor antigen reactive ("ltsr.react"), KIR ("is.KIR"), KLRD1 ("is.KLRD1"), or TIGIT ("is.TIGIT"). 
data <- metadata[,c("ltsr.react","is.KIR","is.CX3CR1","is.KLRD1","is.TIGIT","is.KLRG1")] 

# indicator used to divide the training and testing set with sample size 1:1
ind <- sample(1:nrow(data), size = ceiling(nrow(data)/2), replace = FALSE)

################
# fit a decision tree using the training set
fit <- rpart(ltsr.react ~ .,
             data = data[ind,], method = "class", control = list(cp = 0.001))

# store 2-by-2 table showing the true and predicted results for the 
# training and testing set
class.pred.train <- table(predict(fit, data[ind,], type="class"), data[ind,]$ltsr.react)
class.pred.test <- table(predict(fit, data[-ind,], type="class"), data[-ind,]$ltsr.react)

# print the test error (i.e., 1-accuracy) for the training and testing set
paste("The test error on the training set is", round(1-sum(diag(class.pred.train))/sum(class.pred.train),4))
paste("The test error on the testing set is", round(1-sum(diag(class.pred.test))/sum(class.pred.test),4))

# plot the fitted tree
rpart.plot(fit, uniform=TRUE, main="Classification Tree for ltsr")

################
# further prune the fitted tree to get a simple tree structure
# "cp" controls the tree complexity
fit2 <- prune(fit, cp = 0.01)

# store 2-by-2 table showing the true and predicted results for the 
# training and testing set
class.pred.train2 <- table(predict(fit2, data[ind,], type="class"), data[ind,]$ltsr.react)
class.pred.test2 <- table(predict(fit2, data[-ind,], type="class"), data[-ind,]$ltsr.react)

# print the generated 2-by-2 tables, test error, and accuracy for the 
# training and testing set
rownames(class.pred.train2) <- c("predicted NO", "predicted Yes")
colnames(class.pred.train2) <- c("true NO", "true Yes")
class.pred.train2
paste("The test error on the training set is", round(1-sum(diag(class.pred.train2))/sum(class.pred.train2),4))

rownames(class.pred.test2) <- c("predicted NO", "predicted Yes")
colnames(class.pred.test2) <- c("true NO", "true Yes")
class.pred.test2
paste("The test error on the testing set is", round(1-sum(diag(class.pred.test2))/sum(class.pred.test2),4))
paste("The accuracy on the testing set is", round(sum(diag(class.pred.test2))/sum(class.pred.test2),4))

# plot the fitted tree
rpart.plot(fit2, uniform=TRUE, main="Classification Tree for ltsr")



