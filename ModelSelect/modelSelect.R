modelSelect <- function(dat, folds,nnetSize, numPreds) {
  require(nnet)
  #Preprocess---------------------------------------------------------------------------------------------------
  allModelErrors <<- rep(0, numPreds)
  allModelFormulas = list()
  length(allModelFormulas) = numPreds
  names(dat)[2] = "DIAGNOSIS"
#   dat$DIAGNOSIS = gsub('PD',1,dat$DIAGNOSIS)
#   dat$DIAGNOSIS = gsub('HC',0,dat$DIAGNOSIS)
  remainingPreds= names(dat)[3:(length(names(dat))-6)]
  basePreds <<- c("CAUDATE_R", "CAUDATE_L", "PUTAMEN_R", "PUTAMEN_L", "CAUDATE_ASYMMETRY", "PUTAMEN_ASYMMETRY")
  
  #split data into test and kfold
  testInd = sample(1:dim(dat)[1],floor(.25*dim(dat)[1]))
  testdata = dat[testInd,]
  Kfolddata = dat[-testInd,]
  
  #KFold Cross-Validation---------------------------------------------------------------------------------------------------
  for (p in 1:numPreds) {
    errorList = rep(0,length(remainingPreds))
    list = 1:folds
    id = sample(1:folds,nrow(Kfolddata),replace = TRUE)
    for (i in 1:folds) {
      trainSet = subset(dat, id %in% list[-i])
      testSet = subset(dat, id %in% c(i))
      print(paste('Fold: ',i))
      for (m in 1:length(remainingPreds)) {
        print(m)
        modelNNet = nnet(as.formula(paste(paste("DIAGNOSIS~",paste(basePreds,collapse="+")),'+',remainingPreds[m])) ,
                         data = trainSet, size = nnetSize, rang = 0.4, decay = 5e-4, maxit = 300, trace = FALSE)
        modelPred = predict(object = modelNNet, newdata = testSet, type = 'raw')
        modelPred = round(modelPred, digits = 0)
        classTable = table(testSet$DIAGNOSIS, modelPred)
        errorList[m] = errorList[m] + (1-sum(diag(classTable))/sum(classTable))
      }
    }
    errorList = errorList / folds
    
    
    #find model with min error and add predictor to current model
    minIndex = match(min(errorList),errorList)
    predToAdd = remainingPreds[minIndex]
    remainingPreds = remainingPreds[-minIndex]
    basePreds <<- c(basePreds, predToAdd)
    print(basePreds)
    print(paste("Current model error rate: ",errorList[minIndex]))
    allModelErrors[p] <<- errorList[minIndex]
    allModelFormulas[[p]] <<- as.formula(paste("DIAGNOSIS~",paste(basePreds,collapse="+")))
    
    modelNNet = nnet(as.formula(paste(paste("DIAGNOSIS~",paste(basePreds,collapse="+")),'+',predToAdd)) ,
                     data = trainSet, size = nnetSize, rang = 0.4, decay = 5e-4, maxit = 300, trace = FALSE)
    modelPred = predict(object = modelNNet, newdata = testSet, type = 'raw')
    modelPred = round(modelPred, digits = 0)
    pred <- prediction(modelPred, testSet$DIAGNOSIS)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
    plot(perf, col=rainbow(10))
  }
  
  minIndex = match(min(allModelErrors),allModelErrors)
  modelToUse = allModelFormulas[minIndex]

  modelNNet = nnet(modelToUse, data = testdata, size = nnetSize, rang = 0.4, decay = 5e-4, maxit = 300, trace = FALSE)
  modelPred = predict(object = modelNNet, newdata = testdata, type = 'raw')
  modelPred = round(modelPred, digits = 0)
  classTable = table(testdata$DIAGNOSIS, modelPred)
  print(paste("FINAL ERROR RATE: ", (1-sum(diag(classTable))/sum(classTable))))
}