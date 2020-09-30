library(caret)
library(kernlab)
library(Seurat)

reference <- readRDS('../../citeseq/5k_v3.RDS')
reference <- doPCA(reference)

inTraining <- createDataPartition(reference$cell_type, p = .75, list = FALSE)

training = data.frame(as.matrix(Embeddings(reference, 'pca')[inTraining,1:5]))
training_pos = training[reference$cell_type == 'B cell',]
testing = data.frame(as.matrix(Embeddings(reference, 'pca')[-inTraining,1:5]))

cell_type = 'B cell'

labels = reference$cell_type == cell_type
levels(labels) <- c(TRUE, FALSE)

y_training = labels[inTraining]
y_testing = labels[-inTraining]

training$Class <- as.factor(y_training)

fitControl <- trainControl(method = "repeatedcv",
                           ## 10-fold CV...
                           number = 10,
                           ## repeated ten times
                           repeats = 10)


oc_model <- train(Class ~ ., data = training, 
                  method = ocSVM, 
                  tuneLength = 8,
                  trControl = fitControl)

oc_model
y_pred <- predict(oc_model, testing)



#### Test without caret
model = kernlab::ksvm(as.matrix(training_pos), kernel = "rbfdot", type = 'one-svc', nu = 0.1, scaled = FALSE)
model
y_pred2 = kernlab::predict(model, testing)
table(y_testing, y_pred2)

