############### Demo of applying the Lasso Logistic Regression Model ###############

The following code demonstrates how to use the fitted Lasso Logistic Regression Classifer model for identifying tumor antigen-reactive CD8+ T cells in the blood of patients with melanoma, as described in Lu et al. [Nat Immunol 2025]. The code reads in the Lasso logistic model then randomly generates testing data for 100 subjects. The prediction performance is calculated and summarized as a receiver operating curve.

####################################################################################

######## set the working directory
#setwd("insert path name")

######## load the working environment, which includes the fitted model "Lasso_model" and the target surface markers "target_genes".
load("Lasso_logistic_model.RData")

######## load the following three packages
library(glmnet)
library(stringr)
library(pROC)

######## construct the input matrix, which is n\times 41 matrix and n is the number of patients. 
######## In this "input_mat", the first column corresponds to variable "Sex", which is a binary 
######## variable with 0 denoting female and 1 denoting male. The remaining 40 columns correspond 
######## to the surface markers in object "target_genes". "Y" denotes the response of n patients.

n = 100 #100 subject was chosen arbitrarily. Please change to the actual number of subjects
input_mat <- matrix(0, nrow = n, ncol = 41)
Y <- round(runif(n=n, min=0, max=1), 0) #Generates an random input for the clinical variable sex (F=0, M=1) for the purposes of the demo.

######## predict the probability of being active for the new dataset (input_mat).
predict_Y_prob <- predict(Lasso_model, input_mat, type = "response")

######## draw the roc curve
rocobj_test <- roc(Y, predict_Y_prob[,1])
plot(rocobj_test, print.auc = TRUE, auc.polygon = TRUE, grid = c(0.1,0.2), grid.col = c("green", "red"),
     max.auc.polygon = TRUE, auc.polygon.col = "skyblue", print.thres = TRUE)
