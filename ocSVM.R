### Try to implement one-class SVM in caret
### Following this example: 
### https://topepo.github.io/caret/using-your-own-model-in-train.html
### Example 13.1, but use one-svc function from ksvm in kernlab

ocSVM <- list(type = "Classification",
              library = "kernlab",
              loop = NULL) 


prm <- data.frame(parameter = c("C", "sigma", "nu"),
                  class = rep("numeric", 3),
                  label = c("Cost", "Sigma", "nu"))

ocSVM$parameters <- prm

svmGrid <- function(x, y, len = NULL, search = "grid") {
  library(kernlab)
  ## This produces low, middle and high values for sigma 
  ## (i.e. a vector with 3 elements). 
  sigmas <- kernlab::sigest(as.matrix(x), na.action = na.omit, scaled = TRUE)  
  ## To use grid search:
  if(search == "grid") {
    out <- expand.grid(sigma = mean(as.vector(sigmas[-2])),
                       C = 2 ^((1:len) - 3), 
                       nu = seq(from = 0.1, to = 0.3, by = (0.3-0.1)/(len-1)))
  } else {
    ## For random search, define ranges for the parameters then
    ## generate random values for them
    rng <- extendrange(log(sigmas), f = .75)
    out <- data.frame(sigma = exp(runif(len, min = rng[1], max = rng[2])),
                      C = 2^runif(len, min = -5, max = 8),
                      nu = runif(len, min = 0.1, max = 0.3))
  }
  out
}

ocSVM$grid <- svmGrid

svmFit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) { 
  
  idx = which(y == TRUE)

  kernlab::ksvm(
    x = as.matrix(x[idx,]),
    kernel = "rbfdot",
    type = 'one-svc',
    kpar = list(sigma = param$sigma),
    C = param$C,
    nu = param$nu,
    prob.model = classProbs,
    scaled = FALSE,
    ...
  )
}

ocSVM$fit <- svmFit

svmPred <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
  
  kernlab::predict(modelFit, newdata)

}

ocSVM$predict <- svmPred


svmProb <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
  
  kernlab::predict(modelFit, newdata, type = "probabilities")

}

ocSVM$prob <- svmProb

ocSVM$levels <- function(x) kernlab::lev(x)

