.cvFolds <- function(Y, V){
  #utility function copied from github/ledell/cvAUC
  # Create CV folds (stratify by outcome)	
  Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
  Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}  	
  return(folds)
}

.cvAUCDoFit <- function(v, folds, train, y){
  # utility function copied from github/ledell/cvAUC
  # Train & test a model; return predicted values on test samples
  set.seed(v)
  ycol <- which(names(train)==y)
  params <- list(x = train[-folds[[v]],-ycol,drop=FALSE],
                 y = as.factor(train[-folds[[v]],ycol]),
                 xtest = train[folds[[v]],-ycol,drop=FALSE])
  modelData <- as.data.frame(cbind(params$y,params$x))
  colnames(modelData) <- c("y",paste0("x",seq(ncol(modelData)-1)))
  glmfit <- glm(y~.,family=binomial(link="logit"),data=modelData)
  toPred <- params$xtest
  colnames(toPred) <- paste0("x",seq(ncol(modelData)-1))
  glmpredicted <- predict(glmfit,newdata=toPred,type="response")
  return(glmpredicted)
}

glmCVAUC <- function(train, y = "V1", V = 10) {
  # utility function copied from github/ledell/cvAUC
  # and then modified. Used to create folds
  folds <- .cvFolds(Y = train[,c(y)], V = V)
  
  # Generate CV predicted values
  glmpredictions <- foreach(v = 1:V, .combine="c", 
    .packages=c("stats")) %do% .cvAUCDoFit(v, folds, train, y)
  glmpredictions[unlist(folds)] <- glmpredictions

  # Get CV AUC and 95% confidence interval
  cvaucres <- ci.cvAUC(predictions = glmpredictions, 
                  labels = train[,c(y)],
                  folds = folds, 
                  confidence = 0.95)
  return(cvaucres)
}

computeAUCforBestCluster <- function(countMatrix,responderStatusDF,derivedARIs,derivedPhenotypes) {
    #
    #attach meta data to the count matrix derived from a clustering method.
    #by construction, the rownames of the countmatrix are the samples
    #and each row sums to the total number of observations in the sample.
    #
    countDF <- as.data.frame(countMatrix)
    cellPops <- colnames(countDF)
    countDF$parentCount <- apply(countDF,1,sum)
    countDF$sampleName <- rownames(countDF)
    countDF$subjectID <- rownames(countDF)
    modelDF <- inner_join(countDF,responderStatusDF,by="subjectID")
    #
    #screen clusters for assocaition with simulated outcome
    #using a biniomial GLMM 
    #
    pvalVec <- c()
    for (cellPop in cellPops) {
        mdf <- data.frame(
            resType=as.factor(modelDF[,"modelRT"]),
            subjectID=as.factor(modelDF[,"subjectID"]),
            parentCount=modelDF[,"parentCount"],
            childCount=modelDF[,cellPop]
        )
        suppressMessages(pv <- safeModForPerfectDisc(mdf))
        if (!is.na(pv)) {
            pvalVec <- append(pvalVec,pv[1])
            names(pvalVec)[length(pvalVec)] <- paste0(cellPop)
        }
    }
    #
    #identify the "best" cluster as that with the smallest pvalue
    #and then lookup up its derived ARI and associated phenotype
    #
    bestCluster <- names(pvalVec[which(pvalVec==min(pvalVec))[1]])
    bestClusterARI <- derivedARIs[bestCluster]
    bestClusterPhenotype <- derivedPhenotypes[bestCluster]
    #
    #finally, compute the per sample proportions of the "best" cluster
    #and then use those to compute cvAUC
    #
    selectionDF <- modelDF[,c("modelRT",bestCluster),drop=FALSE]
    for (sName in setdiff(colnames(selectionDF),"modelRT")) {
        selectionDF[,sName] <- (selectionDF[sName]/modelDF[,"parentCount"])
    }
    cvaucResults <- glmCVAUC(train = selectionDF, y = "modelRT", V = 5)
    cvaucResults$bestCluster <- bestCluster
    cvaucResults$bestClusterARI <- bestClusterARI
    cvaucResults$bestClusterPhenotype <- bestClusterPhenotype
    return(cvaucResults)
}

getLabelsForClustering <- function(clusterSource,expressionSource) {
    #
    #for a given clustering, derive labels by comparing cluster medians to 
    #sample medians
    #
    if (length(clusterSource) != nrow(expressionSource)) {
        stop("Labeling clustering of different length")
    }
    simExpMedians <- apply(expressionSource,2,median)
    methodClusterIDs <- names(table(clusterSource))
    methodLabels <- c()
    for (mClusterID in methodClusterIDs) {
        mClusterExprs <- expressionSource[which(clusterSource==mClusterID),,drop=FALSE]
        mClusterMeds <- apply(mClusterExprs,2,median)
        labelBase <- as.numeric(mClusterMeds >= simExpMedians)
        labelInterpret <- rep("-",ncol(expressionSource))
        labelInterpret[which(labelBase==1)] <- "+"
        mClusterLabel <- paste(paste0(paste0("V",seq(10)),labelInterpret),collapse="")
        methodLabels <- append(methodLabels,mClusterLabel)
        names(methodLabels)[length(methodLabels)] <- mClusterID
    }
    return(methodLabels)
}

getCluteringARIsRelativeSpike <- function(clusterSource,truthSource) {
    #
    #given clusterSource, an arbitrary clustering, and truthSource, a binary 
    #vector where 1 corresponds to the true spiked cluster,
    #compute the ARI of each cluster in clusterSource relative to the truthSource
    #
    if (length(clusterSource) != length(truthSource)) {
        stop("Scoring clusterings of different length")
    }
    methodClusterIDs <- names(table(clusterSource))
    methodARIs <- c()
    for (mClusterID in methodClusterIDs) {
        mClusterSpike <- rep(0,length(truthSource)) 
        mClusterSpike[which(clusterSource==mClusterID)] <- 1
        mClusterScore <- ClusterR::external_validation(
                                       true_labels=truthSource,
                                       clusters=mClusterSpike,
                                       method="adjusted_rand_index",
                                       summary_stats=FALSE
                                   )
        methodARIs <- append(methodARIs,mClusterScore)
        names(methodARIs)[length(methodARIs)] <- mClusterID
    }
    return(methodARIs)
}


getCountMatrixFromClustering <- function(clusterSource,numberOfCols,sampleNames,sampleLookups) {
    #
    #given clusterSource, an arbitrary clustering, compute a per-sample
    #count matrix with rows corresponding to samples and columns corresponding to clusterings.
    #assumes the clustering method generates a clusterSource consisting of numeric
    #cluster labels
    #
    countMatrix <- matrix(0,nrow=length(sampleNames),ncol=numberOfCols)
    rownames(countMatrix) <- sampleNames
    colnames(countMatrix) <- seq(numberOfCols)
    for (sampleName in sampleNames) {
        rowLookup <- which(rownames(countMatrix)==sampleName)
        dataLookup <- which(sampleLookups==sampleName)
        sampleCluster <- clusterSource[dataLookup]
        for (clusterNum in seq(numberOfCols)) {
            clusterCount <- length(which(sampleCluster==clusterNum))
            if (clusterCount > 0) {
                colLookup <- which(colnames(countMatrix)==as.character(clusterNum))
                countMatrix[rowLookup,colLookup] <- clusterCount
            }
        }
    }
    return(countMatrix)
}
