require(xgboost)
require(FNN)

XGBoost_train = function(train_object0, var.genes, do.scale=FALSE){
  
  predictor_Data = as.matrix(train_object0@data[var.genes,])
  if (do.scale) predictor_Data = t(scale(t(predictor_Data)))
  max.cells.per.ident = 700; train.frac = 0.8
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using mininum of ", 0.8*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
  for (i in as.character(levels(train_object0@ident))){
    cells.in.clust = which.cells(train_object0,i);
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(as.numeric(i)-1,length(train.temp))); validation.label = c(validation.label, rep(as.numeric(i)-1, length(validation.temp)));
  }
  train_matrix <- xgb.DMatrix(data = t(predictor_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(predictor_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- 200 # number of XGBoost rounds
  print(1)
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  plotConfusionMatrix(table(validation.label, valid_predlabels))
  
  return(bst_model)
}

RF_train = function(train_object0, var.genes, do.scale=FALSE){
  
  predictor_Data = as.matrix(train_object0@data[var.genes,])
  if (do.scale) predictor_Data = t(scale(t(predictor_Data)))
  max.cells.per.ident = 700; train.frac = 0.8
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using mininum of ", 0.8*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
  for (i in as.character(levels(train_object0@ident))){
    cells.in.clust = which.cells(train_object0,i);
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(as.numeric(i)-1,length(train.temp))); validation.label = c(validation.label, rep(as.numeric(i)-1, length(validation.temp)));
  }
  
  
  tmp = as.vector(table(training.label))
  sampsizes = rep(min(tmp),length(tmp))
  rf = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 301, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE) 
  validation_pred_rf = predict(rf,t(predictor_Data[,validation.set]))
  plotConfusionMatrix(table(validation.label, validation_pred_rf))
  
  return(rf)
}

GraphAssignment = function(object, k=30, knn=30, test.label=NULL){
  cells.ident0 = as.character(object@ident); names(cells.ident0) = names(object@ident)
  cells.ident = cells.ident0
  
  l=0
  iter=1
  while (l==0){
    print(iter)
    cells.to.assign = names(cells.ident)[cells.ident == "Unassigned"]
    cells.unassigned = cells.to.assign
    ref.cells = setdiff(rownames(object@pca.rot), cells.unassigned)
    frac.unassigned = length(cells.unassigned)/(length(rownames(object@pca.rot)))
    print(paste0(round(frac.unassigned*100), " percent cells are Unassigned"))
    Xquery = object@pca.rot[cells.to.assign,1:k]
    X = object@pca.rot[ref.cells,1:k]
    nearest = get.knnx(X,Xquery,k=knn+1, algorithm = "kd_tree" )
    new.assignments = apply(nearest$nn.index[,2:(knn+1)],1,function(x){
      Nvotes = table(object@ident[ref.cells][x]);
      if (max(Nvotes) > 0.33*knn){
        return(which.max(Nvotes))
      } else {
        return(46)
      }
    })
    
    cells.ident = as.character(cells.ident); names(cells.ident) = names(object@ident)
    cells.ident[cells.to.assign] = levels(object@ident)[new.assignments]
    
    A_bst = table(test.label,cells.ident)
    A_bst = A_bst[,levels(object@ident)]
    row.order = c(1,3,5,9,2,21,7,10,19,6,11,17,8,4,13,18,16,26,31,20,22,27,25,28,24,23,30,15,12,14,29,32)
    row.order = c(row.order,setdiff(1:32, row.order))
    col.order = c(1:4,29,5,8,18,6,7,9,26,10:11,16,35,12,14,13,15,24,17,20,22,19,25,21,23,30,27:28,31,40,32,45,33,41,34,36:39,42:44)
    col.order = c(col.order,setdiff(1:46, col.order))
    p2 = plotConfusionMatrix(A_bst[row.order,col.order],order=NULL,ylab.use = "ONC RGC", xlab.use="Atlas RGC", x.lab.rot = 45,col.low="white",col.high="darkblue") + ggtitle("XG Boost")
    
    frac.unassigned.new = table(cells.ident)[46] / length(cells.ident)
    percentage.reduction = (frac.unassigned - frac.unassigned.new)  / frac.unassigned
    if (percentage.reduction < 0.15) l = 1
    iter=iter+1
  }
  
  return(cells.ident)
  
  
  
}
