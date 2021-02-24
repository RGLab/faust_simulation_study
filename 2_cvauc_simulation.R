library(ncdfFlow)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(diffcyt)
library(citrus)
library(cvAUC)
library(Rclusterpp)
library(doParallel)
library(ClusterR)
library(flowMeans)
library(FlowSOM)
library(DepecheR)
library(Rphenograph)
library(lme4)
library(tidyr)
library(faust)
library(MCMCpack)
library(mvtnorm)
library(MASS)
library(flowCore)
library(flowWorkspace)
library(ggplot2)
library(reticulate)
parc <- import("parc",convert=FALSE)
source(file.path(normalizePath("."),"functionsForBenchmarking.R"))
source(file.path(normalizePath("."),"functionsForCvaucSim.R"))
########################################################################
#
#set simulation parameters
#
########################################################################
dimension <- 10 #number of columns
gaussianSwitch <- FALSE #if true, multivariate gaussian samples simulated.
sampleSize <- 25000
ssLB <- 25000 #forcing every sample to have the same size.
ssUB <- 25000 #forcing every sample to have the same size.
subSampleSize <- 25000 #turn off subsampling so all methods are applied to entire dataset.
opGridNum <- 15 #used to over-partition the flowSOM grid. 

########################################################################
#
# Fix mixture weights in advance. 
# This vector is explicity designed to have block exponential decay
# in mixture weights up to the 110th component.
# Then, mixture weights are deliberately modifeid so that the 
# [weight of component 110]/[weight of component 120] = 2.0
# In this simulation, the 120th component has its mass spiked to 
# components 110 in 10 of the samples. These spiked
# samples are the responders, the other 10 (where component 120 is untouched)
# are the non-responders.
#
########################################################################
clusterProbVec <- c(0.1425, 0.0712, 0.0712, 0.0356, 0.0356, 0.0356, 0.0356, 0.0178, 
                    0.0178, 0.0178, 0.0178, 0.0178, 0.0178, 0.0178, 0.0178, 0.0089, 
                    0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 
                    0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0089, 0.0045, 
                    0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 
                    0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 
                    0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 
                    0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0045, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 0.0022, 
                    0.0022, 0.0022, 0.0022, 0.0022, 0.0022)
#
#modify the weights to induce the desired foldchange, and put residual mass in first weight
#
clusterProbVec[110:120] <- seq(from=0.0022,to=0.0011,by=-0.00011)
clusterProbVec[121:125] <- clusterProbVec[120]
clusterProbVec[1] <- (clusterProbVec[1] + (1-sum(clusterProbVec)))
numberOfClusters <- length(clusterProbVec) 
targetSpikeRankList <- lapply(c(110,111,112,113,114,115),function(x){c(x,120)}) 
########################################################################
#
# specify the mean vectors. the dimension of specification 
# matrix determines the true annotation in the experiment.
# a mean vector for the multivariate normal is sampled from each column
# if the 0 is selected, the annotation is low for that variable.
# if the non-0 is selected, the annotation is high for that variable.
#
########################################################################
specMat <- matrix(0,nrow=2,ncol=dimension)
specMat[2,] <- c(8,8,8,8,8,8,8,8,8,8)
                                        
########################################################################
#
#Initialize data structures/folders to store simulation results
#
########################################################################
jobPath <- paste0("cvauc_sim_cytof_",dimension,"_dim")
dataTransformations <- list(function(x){return(gamma(1+abs(x/4)))},
                            function(x){return(gamma(1+abs(x/4)))},
                            function(x){return(gamma(1+abs(x/4)))})

if (gaussianSwitch) {
    jobPath <- paste0("cvauc_sim_gaussian_",dimension,"_dim")
    dataTransformations <- list(function(x){return((x))},
                                function(x){return((x))},
                                function(x){return((x))})
}

if (!dir.exists(file.path(normalizePath("."),jobPath))) {
    dir.create(file.path(normalizePath("."),jobPath))
}

if (!dir.exists(file.path(normalizePath("."),jobPath,"results"))) {
    dir.create(file.path(normalizePath("."),jobPath,"results"))
}

###########################################################################
#
#Assign responder status to the samples at random, 15 responders, 15 non
#
###########################################################################
print(paste0("Processing ",numberOfClusters," number of clusters."))        
    
################################################################################################
#
#Initialze containers for cvAUC/ARI of best cluster simulation
#
################################################################################################
rclusterppSimResults <- list()
phenographSimResults <- list()
PARCSimResults <- list()
depecheSimResults <- list()
kmeansSimResults <- list()
faustSimResults <- list()
oracleFlowsomSimResults <- list()
opFlowSOMSimResults <- list()
citrusSimResults <- list()
dsCitrusSimResults <- list()
diffcytSimResults <- list()
currentSeedValue <- 1234

for (tsprNum in seq(from=1,to=1)) {
for (iterNum in seq(25)) {
print(paste0("currentSeedValue: ",currentSeedValue))
set.seed(currentSeedValue)
responderStatus <- rep(0,20)
responderStatus[sort(sample(seq(20),10))] <- 1
#
#determine how much to increase the targeted phenotype for each iteration.
#    
targetSpikeRanks <- targetSpikeRankList[[tsprNum]] 
dir.create(file.path(normalizePath("."),jobPath,"citrusData"))
dir.create(file.path(normalizePath("."),jobPath,"citrusData","inputData"))
dir.create(file.path(normalizePath("."),jobPath,"citrusData","outputData"))

print(paste0("Current iteration number is ",iterNum))
################################################################################################
#
#Simulate an experimental dataset of 20 samples. 
#
################################################################################################
simmedExp <- simulateExperiment(
    meanVectorBoundsMatrix=specMat,
    numSamples=20,
    transformationList=dataTransformations,
    noiseDimension=0,
    probVecIn=clusterProbVec,
    minSampleSize=ssLB,
    maxSampleSize=ssUB,
    randomSeed=currentSeedValue,
    tncp=sampleSize,
    knockoutNumber=0,
    spikeInFlag=TRUE,
    targetRank=targetSpikeRanks,
    useBatchFlag=FALSE,
    batchEffectShift=1,
    fixedResFlag=TRUE,
    responderStatusVec=responderStatus
)

##############################################################
#
#apply faust to the simulated data
#
##############################################################
fs <- as(simmedExp[["flowFrameList"]], "flowSet")
gs <- flowWorkspace::GatingSet(fs)
print("starting faust")
faust(
    gatingSet=gs,
    startingCellPop="root",
    projectPath=file.path(normalizePath("."),jobPath),
    drawAnnotationHistograms = FALSE,
    threadNum = 10,
    annotationsApproved = TRUE,
    debugFlag = FALSE
)

##############################################################
#
#faust orders the selected variables in depth score order.
#this is not necessarily column order,
#so remap the true spiked population label into faust order
#    
##############################################################

faustCountMatrix <- readRDS(file.path(normalizePath("."),jobPath,"faustData","faustCountMatrix.rds"))
fo <- getStrOrder(colnames(faustCountMatrix)[1])
trueSpikedInPop <- truth2faust(simmedExp[["spikedPop"]],fo)
#
#record model responder type for later use
#
allSubjects <- rownames(faustCountMatrix)
resDF <- data.frame(
    modelRT=as.factor(as.numeric(responderStatus)),
    subjectID=allSubjects,
    stringsAsFactors=FALSE
)
#
#collect faust's per-sample clusterings into a single vector to score
#
faustSamplePaths <- list.files(file.path(normalizePath("."),
                                         jobPath,
                                         "faustData",
                                         "sampleData"),
                               full.names=TRUE)
firstSample <- TRUE
for (samplePath in faustSamplePaths) {
    fa <- read.table(paste0(samplePath,"/faustAnnotation.csv"),
                     header=F,
                     sep="`")[,1]
    if (firstSample) {
        faustClustering <- as.character(fa)
        firstSample <- FALSE
    }
    else {
        faustClustering <- append(faustClustering,as.character(fa))
    }
}
    
################################################################################################
#
#Concatenate data into single expression matrix in order to apply other methods.
#Create sub-sampled count matrix for long-running methods.
#
################################################################################################
totalL <- sum(unlist(lapply(simmedExp[["flowFrameList"]],nrow)))
trueLabelVec <- sampleLookups  <- rep("NoSample",totalL)
exprsMat <- matrix(0,nrow=totalL,ncol=dimension)
colnames(exprsMat) <- paste0("V",seq(dimension))
startL <- 1
for (sampleName in sampleNames(gs)) {
    newData <- exprs(gh_pop_get_data(gs[[sampleName]],"root"))
    #
    #transcribe data for a citrus analysis
    #
    fnOut <- file.path(normalizePath("."),jobPath,"citrusData","inputData",paste0(sampleName,".fcs"))
    ffOut <- flowFrame(exprs=newData)
    write.FCS(x=ffOut,filename=fnOut)
    #
    #extract all data
    #
    endL <- (startL + nrow(newData) - 1)
    exprsMat[startL:endL,] <- newData
    sampleLookups[startL:endL] <- sampleName
    trueLabelVec[startL:endL] <- simmedExp[["labelsList"]][[which(names(simmedExp[["labelsList"]])==sampleName)]]
    startL <- (endL + 1)
}
#
#record which elements of the true label vector correspond to the spiked population
#
spikedInClustering <- rep(0,length(trueLabelVec))
spikedInClustering[which(trueLabelVec==simmedExp[["spikedPop"]])] <- 1
################################################################################################
#
#score faust and map to interpretable labels from encoded labels.
#
################################################################################################
faustDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=faustClustering,
    truthSource=spikedInClustering
)
names(faustDerivedClusterARIs) <- gsub("~1~2~","-",names(faustDerivedClusterARIs))
names(faustDerivedClusterARIs) <- gsub("~2~2~","+",names(faustDerivedClusterARIs))

faustDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=faustClustering,
    expressionSource=exprsMat
)
names(faustDerivedClusterLabels) <- gsub("~1~2~","-",names(faustDerivedClusterLabels))
names(faustDerivedClusterLabels) <- gsub("~2~2~","+",names(faustDerivedClusterLabels))
faustResults <- computeAUCforBestCluster(
    countMatrix=faustCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=faustDerivedClusterARIs,
    derivedPhenotypes=faustDerivedClusterLabels
)
faustResults$numberOfClusters <- numberOfClusters
faustResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
faustResults$iterNum <- iterNum
faustResults$tsprNum <- tsprNum  
faustResults$spikedPop <- simmedExp[["spikedPop"]]
faustSimResults <- append(faustSimResults,list(faustResults))
saveRDS(faustSimResults,
        file.path(normalizePath("."),jobPath,"results","faust_results.rds"))

############################################################################################################
#
#collect relevant data about spiked population for discovery methods used in sequel.
#
############################################################################################################
sampleCountSummary <- table(sampleLookups)
daPopLookup <- which(trueLabelVec==simmedExp[["spikedPop"]])
daPopTrueCounts <- sampleCountSummary 
for (cSampleName in names(sampleCountSummary)) {
    daPopTrueCounts[cSampleName] <- length(intersect(which(sampleLookups==cSampleName),daPopLookup))
}
modelResDF <- data.frame(
    Responder=as.factor(as.numeric(responderStatus)),
    subjectID=names(sampleCountSummary),
    fileName=paste0(names(sampleCountSummary),".fcs"),
    parentCount=as.numeric(sampleCountSummary),
    childCount=as.numeric(daPopTrueCounts),
    stringsAsFactors=FALSE
)
modelResDF$prop <- modelResDF$childCount/modelResDF$parentCount
    
############################################################################################################
#
#run citrus analysis with default settings 
#
############################################################################################################
dataDirectory <- file.path(normalizePath("."),jobPath,"citrusData","inputData")
outputDirectory <- file.path(normalizePath("."),jobPath,"citrusData","outputData")
fileList <- modelResDF[,c("fileName"),drop=FALSE]
labels <- as.factor(modelResDF[,"Responder"])

#
#Run full citrusAnalysis
#
dsCitrusResults <- citrus.full(
    fileList = fileList, 
    labels = labels, 
    clusteringColumns = paste0("V",seq(10)),
    transformColumns = NULL,
    transformCofactor = NULL,
    scaleColumns = NULL,
    dataDirectory = dataDirectory, 
    outputDirectory = outputDirectory, 
    family = "classification",
    modelTypes = "glmnet",
    featureType = "abundances",
    nFolds = 1,
    fileSampleSize = subSampleSize
)

dsCitrusData <- dsCitrusResults$citrus.combinedFCSSet$data
#
#the key assumption here is that the rows of citrusPropDF
#map to sample001, sample002,...
#then bind the vector of proportions in this order for downstream glm modeling
#
dsCitrusPropDF <- modelResDF[,c("Responder"),drop=FALSE]
dsCitrusClusters <- dsCitrusResults$conditionRegressionResults$fileName$glmnet$differentialFeatures[["cv.min"]][["clusters"]]
print(paste0("Citrus cluster length: ",length(dsCitrusClusters)))    
for (dsCitrusClusterID in dsCitrusClusters) {
    dsCitrusClusterLookup <- sort(unique(dsCitrusResults$citrus.foldClustering$allClustering$clusterMembership[[as.numeric(dsCitrusClusterID)]]))
    dsCitrusDataSubset <- dsCitrusData[dsCitrusClusterLookup,,drop=FALSE]
    dsCitrusSamplesInCluster <- sort(unique(dsCitrusDataSubset[,"fileId",drop=FALSE]))
    dsCitrusClusterProps <- rep(0,20)
    for (cSampNum in dsCitrusSamplesInCluster) {
        dsCitrusClusterProps[cSampNum] <- (length(which(dsCitrusDataSubset[,"fileId"]==cSampNum))/subSampleSize)
    }
    dsCitrusPropDF <- cbind(dsCitrusPropDF,dsCitrusClusterProps)
    colnames(dsCitrusPropDF)[ncol(dsCitrusPropDF)] <- as.character(dsCitrusClusterID)
}
dsCitrusDerivedResults <- glmCVAUC(train = dsCitrusPropDF, y = "Responder", V = 5)
#
#derive phenotype and ARI for all clusters used by CITRUS    
#
allDsCitrusClusters <- rep(0,nrow(dsCitrusData))
for (dsCitrusClusterID in dsCitrusClusters) {
    dsCitrusClusterLookup <- sort(unique(dsCitrusResults$citrus.foldClustering$allClustering$clusterMembership[[as.numeric(dsCitrusClusterID)]]))
    allDsCitrusClusters[dsCitrusClusterLookup] <- 1
}
dsCitrusARIScore <- ClusterR::external_validation(
                                true_labels=spikedInClustering,
                                clusters=allDsCitrusClusters,
                                method="adjusted_rand_index",
                                summary_stats=FALSE
                            )
dsCitrusDerivedLabels <- getLabelsForClustering(
    clusterSource=allDsCitrusClusters,
    expressionSource=exprsMat
)
dsCitrusDerivedLabel <- dsCitrusDerivedLabels["1"]
#
#record result for output    
#
dsCitrusDerivedResults$bestCluster <- list(dsCitrusClusters)
dsCitrusDerivedResults$bestClusterARI <- dsCitrusARIScore
dsCitrusDerivedResults$bestClusterPhenotype <- dsCitrusDerivedLabel
dsCitrusDerivedResults$numberOfClusters <- numberOfClusters
dsCitrusDerivedResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
dsCitrusDerivedResults$iterNum <- iterNum
dsCitrusDerivedResults$tsprNum <- tsprNum
dsCitrusDerivedResults$spikedPop <- simmedExp[["spikedPop"]]
dsCitrusSimResults <- append(dsCitrusSimResults,list(dsCitrusDerivedResults))
saveRDS(dsCitrusSimResults,
        file.path(normalizePath("."),jobPath,"results","dsCitrus_results.rds"))
    
############################################################################################################
#
#run diffcyt analysis
#
############################################################################################################
marker_info <- data.frame(
    channel_name=paste0("V",seq(10)),
    marker_name=paste0("V",seq(10)),
    marker_class="type",
    stringsAsFactors=FALSE
)
experiment_info <- modelResDF[,c("subjectID"),drop=FALSE]
colnames(experiment_info) <- c("patient_id")
experiment_info$`sample_id` <- experiment_info$`patient_id`
experiment_info$`group_id` <- as.character(as.numeric(modelResDF$Responder)-1)
design <- createDesignMatrix(
  experiment_info, cols_design = c("group_id")
)
contrast <- createContrast(c(0, 1))
out_DA <- diffcyt(
  d_input = fs,
  experiment_info = experiment_info, 
  marker_info = marker_info, 
  design = design, 
  contrast = contrast, 
  analysis_type = "DA", 
  seed_clustering = currentSeedValue,
  xdim=11,
  ydim=12,
  transform=FALSE
)
bestDiffcytCluster <- topTable(out_DA,top_n=1)@rownames
diffcytDF <- as.data.frame(assay(out_DA$d_counts))
diffcytDF <- t(diffcytDF)
diffCytOutputDF <- data.frame(
    sampleID=names(diffcytDF[,bestDiffcytCluster,drop=TRUE]),
    childCount=diffcytDF[,bestDiffcytCluster,drop=TRUE],
    parentCount=apply(diffcytDF,1,sum),
    stringsAsFactors=FALSE
)
compareDF <- modelResDF[,c("subjectID","childCount","parentCount","Responder")]
colnames(compareDF) <- c("sampleID","trueChildCount","trueParentCount","Responder")
diffCytModelDF <- inner_join(diffCytOutputDF,compareDF,by=c("sampleID"))
diffCytModelDF$pct <- diffCytModelDF$childCount/diffCytModelDF$parentCount
diffCytAUCDF <- diffCytModelDF[,c("Responder","pct")]
diffCytResults <- glmCVAUC(train = diffCytAUCDF, y = "Responder", V = 5)

diffcytClustering <- as.data.frame(rowData(out_DA$d_se))[,"cluster_id"]
bestDiffcytClusterIndex <- rep(0,nrow(exprsMat))
bestDiffcytClusterIndex[which(diffcytClustering==bestDiffcytCluster)] <- 1
diffcytARIScore <- ClusterR::external_validation(
                                 true_labels=spikedInClustering,
                                 clusters=bestDiffcytClusterIndex,
                                 method="adjusted_rand_index",
                                 summary_stats=FALSE
                             )
diffcytDerivedLabels <- getLabelsForClustering(
    clusterSource=diffcytClustering,
    expressionSource=exprsMat
)
diffCytResults$bestCluster <- bestDiffcytCluster
diffCytResults$bestClusterARI <-diffcytARIScore
diffCytResults$bestClusterPhenotype <- diffcytDerivedLabels[bestDiffcytCluster]
diffCytResults$numberOfClusters <- numberOfClusters
diffCytResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
diffCytResults$iterNum <- iterNum
diffCytResults$tsprNum <- tsprNum
diffCytResults$spikedPop <- simmedExp[["spikedPop"]]
diffcytSimResults <- append(diffcytSimResults,list(diffCytResults))
saveRDS(diffcytSimResults,
        file.path(normalizePath("."),jobPath,"results","diffcyt_results.rds"))

################################################################################################
#
#First, we apply FlowSOM, and set the grid to 1 by number of true clusters (mixture components).
#The choice of the 1 by number of true clusters comes from
#Handbook of Cluster Analysis, Edited by Christian Hennig, Marina Meila
#Fionn Murtagh, and Roberto Rocci, c 2016 Taylor & Francis Group LLC 
#Version Date: 20151012, ISBN-13: 978-1-4665-5188-6
#SOM discussed pages 421-423. Clustering choice page 421.
#
################################################################################################
ffFlowSOM <- flowCore::flowFrame(exprsMat)
fSOM <- FlowSOM::ReadInput(ffFlowSOM, transform = FALSE, scale = FALSE)
fSOM <- FlowSOM::BuildSOM(fSOM,
                          colsToUse = seq(length(colnames(exprsMat))),
                          xdim = 1,
                          ydim = length(clusterProbVec))
fSOM <- FlowSOM::BuildMST(fSOM)
oracleClusterLabels <- fSOM$map$mapping[, 1]
flowsomCountMatrix <- getCountMatrixFromClustering(
    clusterSource=oracleClusterLabels,
    numberOfCols=length(clusterProbVec),
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
oracleDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=oracleClusterLabels,
    truthSource=spikedInClustering
)
oracleDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=oracleClusterLabels,
    expressionSource=exprsMat
)
oracleFlowSOMResults <- computeAUCforBestCluster(
    countMatrix=flowsomCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=oracleDerivedClusterARIs,
    derivedPhenotypes=oracleDerivedClusterLabels
)
oracleFlowSOMResults$numberOfClusters <- numberOfClusters
oracleFlowSOMResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
oracleFlowSOMResults$iterNum <- iterNum
oracleFlowSOMResults$tsprNum <- tsprNum    
oracleFlowSOMResults$spikedPop <- simmedExp[["spikedPop"]]
oracleFlowsomSimResults <- append(oracleFlowsomSimResults,list(oracleFlowSOMResults))
saveRDS(oracleFlowsomSimResults,
        file.path(normalizePath("."),jobPath,"results","oracle_FlowSOM_results.rds"))

################################################################################################
#
#We next set a grid that is deliberately over-partitioned, and cluster the data
#with FlowSOM again.
#
################################################################################################
opfSOM <- FlowSOM::BuildSOM(fSOM,
                            colsToUse = seq(length(colnames(exprsMat))),
                            xdim = opGridNum,
                            ydim = opGridNum)
opfSOM <- FlowSOM::BuildMST(opfSOM)
opClusterLabels <- opfSOM$map$mapping[, 1]
overpartitionedCountMatrix <- getCountMatrixFromClustering(
    clusterSource=opClusterLabels,
    numberOfCols=(opGridNum*opGridNum),
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
overpartitionedDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=opClusterLabels,
    truthSource=spikedInClustering
)
overpartitionedDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=opClusterLabels,
    expressionSource=exprsMat
)
overpartFlowSOMResults <- computeAUCforBestCluster(
    countMatrix=overpartitionedCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=overpartitionedDerivedClusterARIs,
    derivedPhenotypes=overpartitionedDerivedClusterLabels
)
overpartFlowSOMResults$numberOfClusters <- numberOfClusters
overpartFlowSOMResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
overpartFlowSOMResults$iterNum <- iterNum
overpartFlowSOMResults$tsprNum <- tsprNum
overpartFlowSOMResults$spikedPop <- simmedExp[["spikedPop"]]
opFlowSOMSimResults <- append(opFlowSOMSimResults,list(overpartFlowSOMResults))
saveRDS(opFlowSOMSimResults,
        file.path(normalizePath("."),jobPath,"results","overpart_FlowSOM_results.rds"))
    

################################################################################################
#
#apply kmeans to the simulated data using the true number of clusters
#
################################################################################################
kmeansResult <- kmeans(x=exprsMat,centers=numberOfClusters,iter.max=10000,algorithm="Lloyd")
kmeansClusterLabels <- as.numeric(kmeansResult$cluster)
kmeansNC <- length(table(kmeansClusterLabels))
kmeansCountMatrix <- getCountMatrixFromClustering(
    clusterSource=kmeansClusterLabels,
    numberOfCols=kmeansNC,
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
kmeansDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=kmeansClusterLabels,
    truthSource=spikedInClustering
)
kmeansDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=kmeansClusterLabels,
    expressionSource=exprsMat
)
kmeansResults <- computeAUCforBestCluster(
    countMatrix=kmeansCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=kmeansDerivedClusterARIs,
    derivedPhenotypes=kmeansDerivedClusterLabels
)
kmeansResults$numberOfClusters <- numberOfClusters
kmeansResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
kmeansResults$iterNum <- iterNum
kmeansResults$tsprNum <- tsprNum
kmeansResults$spikedPop <- simmedExp[["spikedPop"]]
kmeansSimResults <- append(kmeansSimResults,list(kmeansResults))
saveRDS(kmeansSimResults,
        file.path(normalizePath("."),jobPath,"results","kmeans_results.rds"))

################################################################################
#
#apply phenograph to the simulated data
#
################################################################################
phenographResult <- Rphenograph(data=exprsMat)
phenographClusterLabels <- as.numeric(membership(phenographResult[[2]]))
phenographNC <- length(table(phenographClusterLabels))
phenographCountMatrix <- getCountMatrixFromClustering(
    clusterSource=phenographClusterLabels,
    numberOfCols=phenographNC,
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
phenographDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=phenographClusterLabels,
    truthSource=spikedInClustering
)
phenographDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=phenographClusterLabels,
    expressionSource=exprsMat
)
phenographResults <- computeAUCforBestCluster(
    countMatrix=phenographCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=phenographDerivedClusterARIs,
    derivedPhenotypes=phenographDerivedClusterLabels
)
phenographResults$numberOfClusters <- numberOfClusters
phenographResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
phenographResults$iterNum <- iterNum
phenographResults$tsprNum <- tsprNum
phenographResults$spikedPop <- simmedExp[["spikedPop"]]
phenographSimResults <- append(phenographSimResults,list(phenographResults))
saveRDS(phenographSimResults,
        file.path(normalizePath("."),jobPath,"results","phenograph_results.rds"))


################################################################################
#
#apply phenograph to the simulated data
#
################################################################################
depecheResult <- depeche(exprsMat,nCores=10,k=(2*numberOfClusters))
depecheClusterLabels <- as.numeric(depecheResult$clusterVector)
depecheNC <- length(table(depecheClusterLabels))
depecheCountMatrix <- getCountMatrixFromClustering(
    clusterSource=depecheClusterLabels,
    numberOfCols=depecheNC,
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
depecheDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=depecheClusterLabels,
    truthSource=spikedInClustering
)
depecheDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=depecheClusterLabels,
    expressionSource=exprsMat
)
depecheResults <- computeAUCforBestCluster(
    countMatrix=depecheCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=depecheDerivedClusterARIs,
    derivedPhenotypes=depecheDerivedClusterLabels
)
depecheResults$numberOfClusters <- numberOfClusters
depecheResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
depecheResults$iterNum <- iterNum
depecheResults$tsprNum <- tsprNum
depecheResults$spikedPop <- simmedExp[["spikedPop"]]
depecheSimResults <- append(depecheSimResults,list(depecheResults))
saveRDS(depecheSimResults,
        file.path(normalizePath("."),jobPath,"results","depeche_results.rds"))

################################################################################
#
#apply parc to the simulated data
#
################################################################################

PARCResult <- parc$PARC(r_to_py(exprsMat))
PARCResult$run_PARC()
PARCClusterLabels <- unlist(py_to_r(PARCResult$labels))
PARCNC <- length(table(PARCClusterLabels))
PARCCountMatrix <- getCountMatrixFromClustering(
    clusterSource=PARCClusterLabels,
    numberOfCols=PARCNC,
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
PARCDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=PARCClusterLabels,
    truthSource=spikedInClustering
)
PARCDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=PARCClusterLabels,
    expressionSource=exprsMat
)
PARCResults <- computeAUCforBestCluster(
    countMatrix=PARCCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=PARCDerivedClusterARIs,
    derivedPhenotypes=PARCDerivedClusterLabels
)
PARCResults$numberOfClusters <- numberOfClusters
PARCResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
PARCResults$iterNum <- iterNum
PARCResults$tsprNum <- tsprNum
PARCResults$spikedPop <- simmedExp[["spikedPop"]]
PARCSimResults <- append(PARCSimResults,list(PARCResults))
saveRDS(PARCSimResults,
        file.path(normalizePath("."),jobPath,"results","PARC_results.rds"))


################################################################################################
#
#apply Rclusterpp to the simulated data
#
################################################################################################
rclusterppResult <- Rclusterpp.hclust(x=exprsMat)
#
#per the implementation, the k parameter in cutree cannot exceed the number of rows of the merge entry plus 1.
#
rclusterppClusterLabels <- cutree(rclusterppResult, k = min(numberOfClusters,(nrow(rclusterppResult$merge)+1)))
rclusterppNC <- length(table(rclusterppClusterLabels))
rclusterppCountMatrix <- getCountMatrixFromClustering(
    clusterSource=rclusterppClusterLabels,
    numberOfCols=rclusterppNC,
    sampleNames=names(table(sampleLookups)),
    sampleLookups=sampleLookups
)
rclusterppDerivedClusterARIs <- getCluteringARIsRelativeSpike(
    clusterSource=rclusterppClusterLabels,
    truthSource=spikedInClustering
)
rclusterppDerivedClusterLabels <- getLabelsForClustering(
    clusterSource=rclusterppClusterLabels,
    expressionSource=exprsMat
)
rclusterppResults <- computeAUCforBestCluster(
    countMatrix=rclusterppCountMatrix,
    responderStatusDF=resDF,
    derivedARIs=rclusterppDerivedClusterARIs,
    derivedPhenotypes=rclusterppDerivedClusterLabels
)
rclusterppResults$numberOfClusters <- numberOfClusters
rclusterppResults$expectedFoldChange <- ((clusterProbVec[targetSpikeRanks][1])/(clusterProbVec[targetSpikeRanks][2]))
rclusterppResults$iterNum <- iterNum
rclusterppResults$tsprNum <- tsprNum
rclusterppResults$spikedPop <- simmedExp[["spikedPop"]]
rclusterppSimResults <- append(rclusterppSimResults,list(rclusterppResults))
saveRDS(rclusterppSimResults,
        file.path(normalizePath("."),jobPath,"results","rclusterpp_results.rds"))

###########################################################################
#
#clean up and increment seed value
#
###########################################################################
unlink(file.path(normalizePath("."),jobPath,"faustData"),recursive=TRUE)
unlink(file.path(normalizePath("."),jobPath,"citrusData"),recursive=TRUE)
currentSeedValue <- (currentSeedValue + 1)
}
}
