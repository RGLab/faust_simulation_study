Posdef <- function (n, ev = runif(n, 1, 2))
{
    #function written by Ravi Varadhan, from r-help mailing list
    #Thu Feb 7, 20:02:30 CET 2008
    #Generates a positive a positive definine matrix of dimension n
    #ev bounds the covariance by 2
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
}

labelMeanVector <- function(meanVector) {
    #
    #assign an interpretable label to the meanVector
    #assumes 0 maps to unexpressed, >0 expressed
    #
    baseStr <- paste0("V",seq(length(meanVector)))
    tokenStr <- rep("+",length(meanVector))
    tokenStr[which(meanVector==0)] <- "-"
    return(paste0(paste0(baseStr,tokenStr),collapse=""))
}

genClusterCentersWithLabels <- function(possibleCenterMat,
                                        nPop=5,
                                        seedVal=0)
{
    #given the possibleCenterMat, sample a location value from each column to determine 
    #the mean vector of a cluster. annotations are then determined by the sampled meanvector
    #Returns the fixedMeanMatrix, which locates each cluster's center
    #Returns the fixedLabelVector, which is the true population of the cluster.
    #if (seedVal > 0) set.seed(seedVal)    
    pcMat <- t(possibleCenterMat)
    allPops <- expand.grid(split(pcMat,rep(seq(nrow(pcMat)),ncol(pcMat))))
    if (nPop > nrow(allPops)) {
        print("Too many clusters relative to specMat/possibleCenterMat.")
        stop("Reduce number of clusters or increase possible mean vectors.")
    }
    clusterIndex <- sample(seq(nrow(allPops)),nPop)
    outMat <- allPops[clusterIndex,,drop=FALSE]
    outLab <- as.character(apply(outMat,1,labelMeanVector))
    outList <- list(fixedMeanMatrix=outMat,fixedLabelVector=outLab)
    return(outList)
}

simSample <- function(sampleDim,
                      sampleSize,
                      transformationType,
                      mixtureType="gaussianOnly",
                      fixedMeanMatrix=NA,
                      fixedLabelVector=NA,
                      noiseDim=0,
                      probVecSample,
                      isKnockout=FALSE,
                      isSpikeIn=FALSE,
                      sRegime=0,
                      targetRanks=c(length(probVecSample)-1,length(probVecSample)),
                      hasBatchEffect=FALSE,
                      batchEffect=0
                      )
{
    #
    #helper function, used to generate simulation samples according to a variety of
    #simulation assumptions.
    #
    numClusters <- nrow(fixedMeanMatrix)
    probVec <- probVecSample
    knockoutStatus <- "No knockout"
    if (isKnockout) {
        #select a cluster with median probability and zero it out
        targetProbs <- sort(probVec)
        targetProb <- targetProbs[ceiling(length(targetProbs)/2)]
        targetLookup <- which(probVec == targetProb)
        modProbLookups <- setdiff(seq(length(probVec)),targetLookup)
        knockoutStatus <- fixedLabelVector[targetLookup]
        probVec[modProbLookups] <- (probVec[modProbLookups] + (targetProb/length(modProbLookups)))
        probVec[targetLookup] <- 0
    }
    if ((sRegime) && (isSpikeIn)) {
        #
        #this logic is run if the sample being simulated is a spiked sample.
        #this means the sample has its 
        #
        #assumes the probVec is sorted in decreasing order.
        #this is guaranteed by the calling function.
        #
        #in the spiked in regime, we increase the prevalance of the target population
        #to the mass in the targetRank vector
        currentSpikeMass <- probVec[targetRanks[2]]
        targetSpikeMass <- probVec[targetRanks[1]]
        
        #spread existing mass equally across all populations.
        massIncrement <- currentSpikeMass/(length(probVec)-1)
        probVec <- probVec + massIncrement
        
        #zero out the spike
        probVec[targetRanks[2]] <- 0
        
        #proportionally decrement mass from other components 
        #so that response is correlated only with the spiked phenotype
        probVec <- probVec - (targetSpikeMass * probVec)
        probVec[targetRanks[2]] <- targetSpikeMass
        if (abs(sum(probVec) - 1) > 1e-9) {
            print(probVec)
            print(sum(probVec))
            print(abs(sum(probVec) - 1))
            stop("Error in probability reapportionment")
        }
        if (min(probVec) == probVec[targetRanks[2]]) {
            stop("Error in spiking in population")
        }
    }
    sampleSizeVec <- as.vector(t(rmultinom(1,sampleSize,probVec)))
    outData <- matrix(nrow=0,ncol=(sampleDim+noiseDim))
    outLabel <- c()
    for (clusterNumber in seq(numClusters)) {
        #if we are simulating a sample with a population knocked out, no draws
        if (sampleSizeVec[clusterNumber] == 0) {
            next
        }
        currentMu <- as.numeric(fixedMeanMatrix[clusterNumber,,drop=TRUE])
        if (hasBatchEffect) {
            currentMu <- currentMu + batchEffect #a constant batch effect translates the sample mean.
        }
        #simulate subject-to-subject variability by perturbing the mean vector
        subjectShift <- round(rnorm(length(currentMu),mean=0,sd=(1/sqrt(2))))
        currentMu <- currentMu + subjectShift #model sample specifc translations of the mean vector
        currentLabel <- fixedLabelVector[clusterNumber]
        outLabel <- append(outLabel,rep(currentLabel,sampleSizeVec[clusterNumber]))
        currentSigma<- Posdef(sampleDim)
        #the next logic determines if the sample is only contains samples from multivariate gaussian
        #multivariate T, or alternatig.
        if (mixtureType == "tPlusGauss") {
            if ((clusterNumber %% 2) == 0) {
                currentSample <- mvrnorm(sampleSizeVec[clusterNumber],mu=currentMu,Sigma=currentSigma)
            }
            else {
                currentSample <- rmvt(sampleSizeVec[clusterNumber],delta=currentMu,sigma=currentSigma,df=2)
            }
        }
        else if (mixtureType == "tOnly") {
            currentSample <- rmvt(sampleSizeVec[clusterNumber],delta=currentMu,sigma=currentSigma,df=2)
        }
        else  {
            currentSample <- mvrnorm(sampleSizeVec[clusterNumber],mu=currentMu,Sigma=currentSigma)
        }
        #if a cluster only has 1 observation in a sample, coerce it to a matrix.
        if (is.vector(currentSample)) {
            currentSample <- t(as.matrix(currentSample))
        }
        #if the simulation has nuisance variables, sample them and bind them to the sample.
        if (noiseDim > 0) {
            currentSigma<- Posdef(noiseDim)
            noiseSample <- mvrnorm(nrow(currentSample),mu=rep(5,noiseDim),Sigma=currentSigma)
            if (is.vector(noiseSample)) {
                noiseSample <- t(as.matrix(noiseSample))
            }
            currentSample <- cbind(currentSample,noiseSample)
        }
        #apply the transformation type {either identity or gamma(1+abs(x/4))} to the columns
        outData <- rbind(outData,t(apply(currentSample,1,transformationType)))
    }
    colnames(outData) <- paste0("V",seq(ncol(outData)))
    outList <- list(sampleMatrix=outData,sampleLabels=outLabel,knockoutInfo=knockoutStatus)
    return(outList)
}


getStrOrder <- function(faustLabel) {
    #helper function, used to figure out the depth-score order of faust labels    
    orderVec <- c()
    while (nchar(faustLabel) > 0) {
        currentToken <- substr(faustLabel,1,1)
        if (currentToken=="V") {
            orderStr <- c()
        }
        else if (currentToken %in% c("-","M","+")) {
            orderVec <- append(orderVec,as.numeric(orderStr))
        }
        else {
            orderStr <- paste0(append(orderStr,currentToken),collapse="")
        }
        if (nchar(faustLabel) == 1){
            faustLabel <- ""
        }
        else {
            faustLabel <- substr(faustLabel,2,nchar(faustLabel))
        }
    }
    return(orderVec)
}

truth2faust <- function(trueLabel,faustOrder) {
    #helper function, used to reorder a true label by the depth-score faustOrder.
    tokenStr <- labelVec <- c()
    vctr <- 0
    while (nchar(trueLabel) > 0) {
        currentToken <- substr(trueLabel,1,1)
        if (currentToken=="V") {
            vctr <- vctr + 1
            if (vctr > 1) {
                labelVec <- append(labelVec,tokenStr)
                tokenStr <- c()
            }
        }
        tokenStr <- paste0(append(tokenStr,currentToken),collapse="")
        if (nchar(trueLabel) == 1){
            trueLabel <- ""
        }
        else {
            trueLabel <- substr(trueLabel,2,nchar(trueLabel))
        }
    }
    labelVec <- append(labelVec,tokenStr)
    return(paste0(labelVec[faustOrder],collapse=""))
}

safeModForPerfectDisc <- function(dataSet) {
    out <- tryCatch(
    {
        m <- glmer(cbind(childCount,(parentCount-childCount)) ~ resType + (1|subjectID),
                   data=dataSet,family="binomial",
                   control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e5)))
        #turn off singularity as warning since methods can find perfect predictors.
        #check.conv.singular = .makeCC(action = "warning",  tol = 1e-4)))
        rv <- c(coefficients(summary(m))[2,4])
        return(rv)
    },
    error=function(cond){
        #message("Error!")
        #message(cond)
        return(NA)
    },
    warning=function(cond){
        #message("Warning!")
        #message(cond)
        return(NA)
    },
    finally={
        #message("done!")
    }
    )
    return(out)
}

simulateExperiment <- function(
                               meanVectorBoundsMatrix, #2xp matrix, first row determines "-" for the variable, second determines "+"
                               numSamples=100, #how many samples are in the simulated dataset
                               randomSeed=0,#reproducibility
                               transformationList=list(function(x){return(x)},
                                                       function(x){return(x)},
                                                       function(x){return(x)}), #column-wise transformations to apply to 
                               noiseDimension=5,#how many nuisance variables
                               probVecIn, #the mixture weights
                               minSampleSize=5000,
                               maxSampleSize=40000,
                               tncp=10000,
                               knockoutNumber=0, #should a component be knocked out, set to 1 if so
                               spikeInFlag=FALSE, #is there a spiked component
                               targetRanks=c(length(probVecIn)-1,length(probVecIn)),
                               useBatchFlag=FALSE, #is there a batch effect
                               batchEffectShift=0.75,
                               fixedResFlag=FALSE, #TRUE if a responder always has the spiked phenotype
                               responderStatusVec=c(NA)
                               )
{
    #
    #simulate an experimental dataset according to a variety of assumptions.
    #always sort the probability vector to guarantee smallest mass is the final element.    
    #
    probVecForSim <- sort(probVecIn,decreasing=TRUE)
    sampleSpecs <- genClusterCentersWithLabels(
        possibleCenterMat=meanVectorBoundsMatrix,
        nPop=length(probVecForSim),
        seedVal=randomSeed
    )
    currentTransformation <- transformationList[[1]]
    currentSample <- 1
    startKnockoutNum <- numSamples - knockoutNumber + 1
    isKnockoutSample <- FALSE
    #the regime determines if we spike in a population
    #regime==1 -> spike it in.
    regime <- rep(0,numSamples)
    regime[sample(seq(numSamples),(numSamples/2))] <- 1
    if (fixedResFlag) {
        #use this flag to always have responders exhibit the differentially abundant phenotype
        regime <- rep(0,numSamples)
        regime[which(responderStatusVec==1)] <- 1
    }
    nextBatchEffect <- rep((-1*batchEffectShift),ncol(sampleSpecs$fixedMeanMatrix))
    regimeList <- knockoutList <- labelsList <- flowList <- truthList <- list()
    while (currentSample <= numSamples) {
        #bound an experimental sample at 5000 observations
        nextSampleSize <- min(max(minSampleSize,round(rt(1,df=3,ncp=tncp))),maxSampleSize)
        nextRegime <- regime[currentSample]
        if (currentSample >= startKnockoutNum) {
            isKnockoutSample <- TRUE
        }
        else {
            isKnockoutSample <- FALSE
        }
        #if there is a batch effect, is affects samples mod 10.
        if ((useBatchFlag) && (((currentSample - 1) %% 10) == 0)) {
            nextBatchEffect <- nextBatchEffect + batchEffectShift
            #print(paste0("Batch: ",nextBatchEffect[1]))
            if ((currentSample - 1) > floor(numSamples/3)) {
                currentTransformation <- transformationList[[2]]
            }
            if ((currentSample - 1) > floor((2*(numSamples/3)))) {
                currentTransformation <- transformationList[[3]]
            }
        }

        sampleData <- simSample(
            sampleDim=ncol(sampleSpecs$fixedMeanMatrix),
            sampleSize=nextSampleSize,
            transformationType=currentTransformation,
            fixedMeanMatrix=sampleSpecs$fixedMeanMatrix,
            fixedLabelVector=sampleSpecs$fixedLabelVector,
            noiseDim=noiseDimension,
            probVecSample=probVecForSim,
            isKnockout=isKnockoutSample,
            isSpikeIn=spikeInFlag,
            sRegime=nextRegime,
            targetRanks=targetRanks,
            hasBatchEffect=useBatchFlag,
            batchEffect=nextBatchEffect
        )
        #record results
        if (currentSample < 10) {
            outName <- paste0("sample00",currentSample)
        }
        else if (currentSample < 100) {
            outName <- paste0("sample0",currentSample)
        }
        else {
            outName <- paste0("sample",currentSample)
        }
        ff <- flowFrame(sampleData$sampleMatrix)
        flowList <- append(flowList,ff)
        names(flowList)[length(flowList)] <- outName
        labelsList <- append(labelsList,list(sampleData$sampleLabels))
        names(labelsList)[length(labelsList)] <- outName
        truthList <- append(truthList,list(table(sampleData$sampleLabels)))
        names(truthList)[length(truthList)] <- outName
        knockoutList <- append(knockoutList,list(table(sampleData$knockoutInfo)))
        names(knockoutList)[length(knockoutList)] <- outName
        regimeList <- append(regimeList,list(nextRegime))
        names(regimeList)[length(regimeList)] <- outName
        currentSample <- currentSample + 1
    }
    truthMat <- getSimTrueCountMatrix(truthList)
    
    outputList <- list(
        "flowFrameList"=flowList,
        "truthValueList"=truthList,
        "truthMat"=truthMat,
        "labelsList"=labelsList,
        "knockoutList"=knockoutList,
        "spikedPop"=sampleSpecs$fixedLabelVector[targetRanks[2]],
        "spikedInPopList"=regimeList
    )
    return(outputList)
}

getSimTrueCountMatrix <- function(truthList)
{
    #helper function, get a sample by truth count matrix
    uniqueNames <- c()
    for (i in seq(length(truthList))) {
        uniqueNames <- append(uniqueNames,names(truthList[[i]]))
    }
    uniqueNames <- sort(unique(uniqueNames))
    truthMat <- matrix(0,nrow=length(truthList),ncol=length(uniqueNames))
    colnames(truthMat) <- uniqueNames
    rownames(truthMat) <- names(truthList)
    for (i in seq(length(truthList))) {
        cv <- truthList[[i]]
        for (cName in colnames(truthMat)) {
            lookup <- which(names(cv) == cName)
            if (length(lookup)) {
                truthMat[i,cName] <- cv[lookup]
            }
        }
    }
    return(truthMat)
}



