library(ClusterR)
library(DepecheR)
library(flowMeans)
library(FlowSOM)
library(Rphenograph)
library(lme4)
library(tidyr)
library(dplyr)
library(faust)
library(MCMCpack)
library(mvtnorm)
library(MASS)
library(flowCore)
library(flowWorkspace)
library(ggplot2)
library(Rclusterpp)
library(reticulate)
parc <- import("parc",convert=FALSE)
source(file.path(normalizePath("."),"functionsForBenchmarking.R"))
set.seed(12345)
startingPath <- normalizePath(".")
#
#containers for results
#
parcNC <- rcppNC <- FlowORCNC <- kmeansNC <-  faustNC <- FlowSOMNC <- phenographNC <- depecheNC <- flowMeansNC <- c()
parcAll <- rcppAll <- FlowORCAll <- kmeansAll <-  faustAll <- FlowSOMAll <- phenographAll <- depecheAll <- flowMeansAll <- c()
parcPheno <- rcppPheno <- FlowORCPheno <- kmeansPheno <-  faustPheno <- FlowSOMPheno<- phenographPheno <- depechePheno <- flowMeansPheno <- c()
faustPctAnn <- c()
#
#simulation parameters
#
dimension <- 10 
sampleSize <- 25000 
ssLB <- 25000 #fixed so all samples are 25,000 obs
ssUB <- 25000 #fixed so all samples are 25,000 obs

#
#flip this switch to simulated well separated "flow == TRUE" or "cytof == FALSE"
#
gaussianSwitch <- FALSE

#
#specify the mean vectors. entry of 0 == "-" for the marker, entry of 8 == "+"
#
specMat <- matrix(0,nrow=2,ncol=dimension)
specMat[2,] <- c(8,8,8,8,8,8,8,8,8,8)

#
#directory for intermediate results and simulation output
#
jobPath <- paste0("num_of_clus_cytof_",dimension,"_",sampleSize)
dataTransformations <- list(function(x){return(gamma(1+abs(x/4)))},
                            function(x){return(gamma(1+abs(x/4)))},
                            function(x){return(gamma(1+abs(x/4)))})
if (gaussianSwitch) {
    jobPath <- paste0("num_of_clus_gaussian_",dimension,"_",sampleSize)
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

#
#fix a block exponential decay for the cluster proportions
#
startProb <- 0.1425
refClusterProbVec <- c(
    startProb,rep(startProb/2,2),rep(startProb/4,4),rep(startProb/8,8),rep(startProb/16,16),rep(startProb/32,32),rep(startProb/64,64)
)
#
#iterate over possible number of clusters
#
for (numberOfClusters in seq(5,125,by=10)) {
print(paste0("Processing ",numberOfClusters," number of clusters."))        
#
#5 simulations per each true numberOfClusters
#
for (iterNum in seq(5)) {
print(paste0("iter ",iterNum))
#
#from the refClusterProbVec, select weights according to the true numberOfClusters
#then uniformly add mass to make mixture weights that sum to 1
#
clusterProbVec <- sort(refClusterProbVec[seq(numberOfClusters)],decreasing=TRUE)
if (sum(clusterProbVec) < 1) {
    residualPV <- 1-sum(clusterProbVec)
    clusterProbVec <- (clusterProbVec + (residualPV/length(clusterProbVec)))
}
if (sum(clusterProbVec) != 1) {
    stop("Cluster weights don't sum to 1")
}
#
#Simulate an experiment. 
#
simmedExp <- simulateExperiment(
    meanVectorBoundsMatrix=specMat,
    numSamples=10,
    transformationList=dataTransformations,
    noiseDimension=0,
    probVecIn=clusterProbVec,
    minSampleSize=ssLB,
    maxSampleSize=ssUB,
    randomSeed=(iterNum*numberOfClusters),
    tncp=sampleSize
)

#
#process the experiment with faust
#
fs <- as(simmedExp[["flowFrameList"]], "flowSet")
gs <- flowWorkspace::GatingSet(fs)
print("starting faust")
faust(
    gatingSet=gs,
    startingCellPop="root",
    projectPath=file.path(normalizePath("."),jobPath),
    drawAnnotationHistograms = FALSE,
    threadNum = 10,
    annotationsApproved = TRUE
)
print("faust complete.")
#
#subtract 1 from faust to remove un-annotated cells 0_0_0_0_0 from the comparison
#
faustNC <- append(faustNC,
                  (ncol(readRDS(file.path(normalizePath("."),
                                          jobPath,
                                          "faustData",
                                          "faustCountMatrix.rds")))-1))
#
#Concatenate the simulated data in order to apply other methods
#
totalL <- sum(unlist(lapply(simmedExp[["flowFrameList"]],nrow)))
trueLabelVec <- sampleLookups  <- rep("NoSample",totalL)
exprsMat <- matrix(0,nrow=totalL,ncol=dimension)
colnames(exprsMat) <- paste0("V",seq(dimension))
startL <- 1
for (sampleName in sampleNames(gs)) {
    print(sampleName)
    newData <- exprs(gh_pop_get_data(gs[[sampleName]],"root"))
    endL <- (startL + nrow(newData) - 1)
    exprsMat[startL:endL,] <- newData
    sampleLookups[startL:endL] <- sampleName
    #
    #record the true labels for each simulated dataset
    #
    trueLabelVec[startL:endL] <- simmedExp[["labelsList"]][[which(names(simmedExp[["labelsList"]])==sampleName)]]
    startL <- (endL + 1)
}
#
#apply flowMeans to the simulated dataset
#
print("starting flowMeans")
flowMeansResult <- flowMeans(x=exprsMat,MaxN=(2*numberOfClusters))
print("flowMeans complete")
flowMeansClusterLabels <- as.numeric(flowMeansResult@Labels[[1]])
flowMeansNC <- append(flowMeansNC,length(table(flowMeansClusterLabels)))

#
#apply FlowSOM to the simulated dataset
#maxMeta uses ConsensusClusterPlus. This consturcts an object called 'this_dist'
#that bounds the number of elements returned by stats::hclust at 90
#this then throws an error in stats::cutree, since maxMeta > 90 leads to k > 90
#in that function.
#
print("Starting FlowSOM")
ffFlowSOM <- flowCore::flowFrame(exprsMat)
FlowSOMResult <- FlowSOM(input=ffFlowSOM,colsToUse=seq(ncol(exprsMat)),maxMeta=90)
print("FlowSOM complete")
FlowSOMClusterLabels <- as.numeric(GetMetaclusters(FlowSOMResult,FlowSOMResult$metaclustering))
FlowSOMNC <- append(FlowSOMNC,length(table(FlowSOMClusterLabels)))

#
#apply FlowSOM to the simulated dataset, using the true number of clusters
#need to use the member functions in place of the FlowSOM wrappers to
#avoid the ConsensusClusterPlut issue.
#
print("Starting FlowORC")
ffFlowORC <- flowCore::flowFrame(exprsMat)
fSOM <- FlowSOM::ReadInput(ffFlowORC, transform = FALSE, scale = FALSE)
fSOM <- FlowSOM::BuildSOM(
                     fSOM,
                     colsToUse = seq(length(colnames(exprsMat))),
                     xdim = 1,
                     ydim = numberOfClusters
                 )
fSOM <- FlowSOM::BuildMST(fSOM)
print("FlowORC complete")
FlowORCClusterLabels <- fSOM$map$mapping[, 1]
FlowORCNC <- append(FlowORCNC,length(table(FlowORCClusterLabels)))

#
#apply phenograph to the simulated dataset
#
print("Starting Phenograph")
phenographResult <- Rphenograph(data=exprsMat)
print("Phenograph complete")
phenographClusterLabels <- as.numeric(membership(phenographResult[[2]]))
phenographNC <- append(phenographNC,length(table(phenographClusterLabels)))

#
#apply Depeche to the simulated dataset, with starting k set to twice the true number of clusters
#
print("Starting Depeche")
depecheResult <- depeche(exprsMat,nCores=10,k=(2*numberOfClusters))
print("Depeche complete")
depecheClusterLabels <- as.numeric(depecheResult$clusterVector)
depecheNC <- append(depecheNC,length(table(depecheClusterLabels)))

#
#apply kmeans to the simulated datset, with k set to the true number clusters
#
print("Starting kmeans")
kmeansResult <- kmeans(x=exprsMat,centers=numberOfClusters,iter.max=10000,algorithm="Lloyd")
print("Kmeans complete")
kmeansClusterLabels <- as.numeric(kmeansResult$cluster)
kmeansNC <- append(kmeansNC,length(table(kmeansClusterLabels)))

#
#apply rcpp to the simulated datatset, with the tree cut at the true number of clusters
#
print("Starting rcpp")
rcppResult <- Rclusterpp.hclust(x=exprsMat)
print("Rcpp complete")
#
#per the implementation, the k parameter in cutree cannot exceed the number of rows of the merge entry plus 1.
#
cutreekVal <-  min(numberOfClusters,(nrow(rcppResult$merge)+1))
print(paste0("Using cutreekVal: ", cutreekVal))
rcppClusterLabels <- cutree(rcppResult, k = cutreekVal)
rcppNC <- append(rcppNC,length(table(rcppClusterLabels)))

#
#apply PARC to the simulated dataset
#
print("Starting parc")
parcResult <- parc$PARC(r_to_py(exprsMat))
parcResult$run_PARC()
print("parc complete")
parcClusterLabels <- unlist(py_to_r(parcResult$labels))
parcNC <- append(parcNC,length(table(parcClusterLabels)))

#
#compute ARI for the different methods relative to the simulated truth
#
trueLabels <- Reduce(append,simmedExp[["labelsList"]])
trueNumericLabels <- as.numeric(as.factor(trueLabels))

print("Start results all phenograph")
phenographScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                                  clusters=as.numeric(phenographClusterLabels),
                                                  method="adjusted_rand_index",
                                                  summary_stats=TRUE)
print("End results all phenograph")


print("Start results all flowMeans")
flowMeansScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                                 clusters=as.numeric(flowMeansClusterLabels),
                                                 method="adjusted_rand_index",
                                                 summary_stats=TRUE)
print("End results all flowMeans")


print("Start results all FlowSOM")
FlowSOMScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                               clusters=as.numeric(FlowSOMClusterLabels),
                                               method="adjusted_rand_index",
                                               summary_stats=TRUE)
print("End results all FlowSOM")

print("Start results all FlowORC")
FlowORCScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                               clusters=as.numeric(FlowORCClusterLabels),
                                               method="adjusted_rand_index",
                                               summary_stats=TRUE)
print("End results all FlowORC")


print("Start results all depeche")
depecheScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                               clusters=as.numeric(depecheClusterLabels),
                                               method="adjusted_rand_index",
                                               summary_stats=TRUE)
print("End results all depeche")

print("Start results all kmeans")
kmeansScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                               clusters=as.numeric(kmeansClusterLabels),
                                               method="adjusted_rand_index",
                                               summary_stats=TRUE)
print("End results all kmeans")

print("Start results all rcpp")
rcppScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                               clusters=as.numeric(rcppClusterLabels),
                                               method="adjusted_rand_index",
                                               summary_stats=TRUE)
print("End results all rcpp")


print("Start results all parc")
parcScores <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                            clusters=as.numeric(parcClusterLabels),
                                            method="adjusted_rand_index",
                                            summary_stats=TRUE)
print("End results all parc")

#
#collect faust per-sample clusterings into single vector to compute ARI
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

print("Start results all faust")
faustScoresAll <- ClusterR::external_validation(true_labels=trueNumericLabels,
                                               clusters=as.numeric(as.factor(faustClustering)),
                                               method="adjusted_rand_index",
                                               summary_stats=TRUE)
print("End results all faust")

#
#save results
#
faustAll <- append(faustAll,faustScoresAll)
FlowSOMAll <- append(FlowSOMAll,FlowSOMScores)
FlowORCAll <- append(FlowORCAll,FlowORCScores)
phenographAll <- append(phenographAll,phenographScores)
flowMeansAll <- append(flowMeansAll,flowMeansScores)
depecheAll <- append(depecheAll,depecheScores)
kmeansAll <- append(kmeansAll,kmeansScores)
rcppAll <- append(rcppAll,rcppScores)
parcAll <- append(parcAll,parcScores)

#
#look up the subset of events that FAUST annotates
#
faustPhenoLookup <- which(faustClustering != "0_0_0_0_0")
faustPctAnnotated <- length(faustPhenoLookup)/length(faustClustering)
faustPctAnn <- append(faustPctAnn,faustPctAnnotated)


#
#compute the ARI for all methods relative to the simulated truth on the annotated subset of events
#
print("Start results pheno faust")
faustScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                        clusters=as.numeric(as.factor(faustClustering[faustPhenoLookup])),
                                        method="adjusted_rand_index",
                                        summary_stats=TRUE)
print("End results all faust")

print("Start results pheno phenograph")
phenographScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                            clusters=phenographClusterLabels[faustPhenoLookup],
                                            method="adjusted_rand_index",
                                            summary_stats=TRUE)
print("End results pheno phenograph")


print("Start results pheno depeche")
depecheScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                          clusters=depecheClusterLabels[faustPhenoLookup],
                                          method="adjusted_rand_index",
                                          summary_stats=TRUE)
print("End results pheno depeche")

print("Start results pheno FlowSOM")
FlowSOMScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                          clusters=FlowSOMClusterLabels[faustPhenoLookup],
                                          method="adjusted_rand_index",
                                          summary_stats=TRUE)
print("End results pheno FlowSOM")


print("Start results pheno FlowORC")
FlowORCScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                          clusters=FlowORCClusterLabels[faustPhenoLookup],
                                          method="adjusted_rand_index",
                                          summary_stats=TRUE)
print("End results pheno FlowORC")

print("Start results pheno flowMeans")
flowMeansScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                            clusters=flowMeansClusterLabels[faustPhenoLookup],
                                            method="adjusted_rand_index",
                                            summary_stats=TRUE)
print("End results pheno flowMeans")

print("Start results pheno kmeans")
kmeansScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                            clusters=kmeansClusterLabels[faustPhenoLookup],
                                            method="adjusted_rand_index",
                                            summary_stats=TRUE)
print("End results pheno kmeans")

print("Start results pheno rcpp")
rcppScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                            clusters=rcppClusterLabels[faustPhenoLookup],
                                            method="adjusted_rand_index",
                                            summary_stats=TRUE)
print("End results pheno rcpp")


print("Start results pheno parc")
parcScoresPheno <- external_validation(true_labels=trueNumericLabels[faustPhenoLookup],
                                       clusters=parcClusterLabels[faustPhenoLookup],
                                       method="adjusted_rand_index",
                                       summary_stats=TRUE)
print("End results pheno parc")

#
#save results for annnotated subset
#
faustPheno <- append(faustPheno,faustScoresPheno)
FlowSOMPheno <- append(FlowSOMPheno,FlowSOMScoresPheno)
FlowORCPheno <- append(FlowORCPheno,FlowORCScoresPheno)
phenographPheno <- append(phenographPheno,phenographScoresPheno)
depechePheno <- append(depechePheno,depecheScoresPheno)
flowMeansPheno <- append(flowMeansPheno,flowMeansScoresPheno)
kmeansPheno <- append(kmeansPheno,kmeansScoresPheno)
rcppPheno <- append(rcppPheno,rcppScoresPheno)
parcPheno <- append(parcPheno,parcScoresPheno)

#
#clean up faust data directory for next iteration
#
unlink(file.path(normalizePath("."),jobPath,"faustData"),recursive=TRUE)
}
#
#checkpoint intermediate results
#
saveRDS(faustPctAnn,file.path(normalizePath("."),jobPath,"results","intermediate_faustPctAnn.rds"))

saveRDS(faustNC,file.path(normalizePath("."),jobPath,"results","intermediate_faustNC.rds"))
saveRDS(FlowSOMNC,file.path(normalizePath("."),jobPath,"results","intermediate_FlowSOMNC.rds"))
saveRDS(FlowORCNC,file.path(normalizePath("."),jobPath,"results","intermediate_FlowORCNC.rds"))
saveRDS(phenographNC,file.path(normalizePath("."),jobPath,"results","intermediate_phenographNC.rds"))
saveRDS(depecheNC,file.path(normalizePath("."),jobPath,"results","intermediate_depecheNC.rds"))
saveRDS(flowMeansNC,file.path(normalizePath("."),jobPath,"results","intermediate_flowMeansNC.rds"))
saveRDS(kmeansNC,file.path(normalizePath("."),jobPath,"results","intermediate_kmeansNC.rds"))
saveRDS(rcppNC,file.path(normalizePath("."),jobPath,"results","intermediate_rcppNC.rds"))
saveRDS(parcNC,file.path(normalizePath("."),jobPath,"results","intermediate_parcNC.rds"))

saveRDS(faustAll,file.path(normalizePath("."),jobPath,"results","intermediate_faustAll.rds"))
saveRDS(FlowSOMAll,file.path(normalizePath("."),jobPath,"results","intermediate_FlowSOMAll.rds"))
saveRDS(FlowORCAll,file.path(normalizePath("."),jobPath,"results","intermediate_FlowORCAll.rds"))
saveRDS(phenographAll,file.path(normalizePath("."),jobPath,"results","intermediate_phenographAll.rds"))
saveRDS(depecheAll,file.path(normalizePath("."),jobPath,"results","intermediate_depecheAll.rds"))
saveRDS(flowMeansAll,file.path(normalizePath("."),jobPath,"results","intermediate_flowMeansAll.rds"))
saveRDS(kmeansAll,file.path(normalizePath("."),jobPath,"results","intermediate_kmeansAll.rds"))
saveRDS(rcppAll,file.path(normalizePath("."),jobPath,"results","intermediate_rcppAll.rds"))
saveRDS(parcAll,file.path(normalizePath("."),jobPath,"results","intermediate_parcAll.rds"))

saveRDS(faustPheno,file.path(normalizePath("."),jobPath,"results","intermediate_faustPheno.rds"))
saveRDS(FlowSOMPheno,file.path(normalizePath("."),jobPath,"results","intermediate_FlowSOMPheno.rds"))
saveRDS(FlowORCPheno,file.path(normalizePath("."),jobPath,"results","intermediate_FlowORCPheno.rds"))
saveRDS(phenographPheno,file.path(normalizePath("."),jobPath,"results","intermediate_phenographPheno.rds"))
saveRDS(depechePheno,file.path(normalizePath("."),jobPath,"results","intermediate_depechePheno.rds"))
saveRDS(flowMeansPheno,file.path(normalizePath("."),jobPath,"results","intermediate_flowMeansPheno.rds"))
saveRDS(kmeansPheno,file.path(normalizePath("."),jobPath,"results","intermediate_kmeansPheno.rds"))
saveRDS(rcppPheno,file.path(normalizePath("."),jobPath,"results","intermediate_rcppPheno.rds"))
saveRDS(parcPheno,file.path(normalizePath("."),jobPath,"results","intermediate_parcPheno.rds"))
}
#
#save the final simulate results
#
saveRDS(faustPctAnn,file.path(normalizePath("."),jobPath,"results","final_faustPctAnn.rds"))

saveRDS(faustNC,file.path(normalizePath("."),jobPath,"results","final_faustNC.rds"))
saveRDS(FlowSOMNC,file.path(normalizePath("."),jobPath,"results","final_FlowSOMNC.rds"))
saveRDS(FlowORCNC,file.path(normalizePath("."),jobPath,"results","final_FlowORCNC.rds"))
saveRDS(phenographNC,file.path(normalizePath("."),jobPath,"results","final_phenographNC.rds"))
saveRDS(depecheNC,file.path(normalizePath("."),jobPath,"results","final_depecheNC.rds"))
saveRDS(flowMeansNC,file.path(normalizePath("."),jobPath,"results","final_flowMeansNC.rds"))
saveRDS(kmeansNC,file.path(normalizePath("."),jobPath,"results","final_kmeansNC.rds"))
saveRDS(rcppNC,file.path(normalizePath("."),jobPath,"results","final_rcppNC.rds"))
saveRDS(parcNC,file.path(normalizePath("."),jobPath,"results","final_parcNC.rds"))

saveRDS(faustAll,file.path(normalizePath("."),jobPath,"results","final_faustAll.rds"))
saveRDS(FlowSOMAll,file.path(normalizePath("."),jobPath,"results","final_FlowSOMAll.rds"))
saveRDS(FlowORCAll,file.path(normalizePath("."),jobPath,"results","final_FlowORCAll.rds"))
saveRDS(phenographAll,file.path(normalizePath("."),jobPath,"results","final_phenographAll.rds"))
saveRDS(depecheAll,file.path(normalizePath("."),jobPath,"results","final_depecheAll.rds"))
saveRDS(flowMeansAll,file.path(normalizePath("."),jobPath,"results","final_flowMeansAll.rds"))
saveRDS(kmeansAll,file.path(normalizePath("."),jobPath,"results","final_kmeansAll.rds"))
saveRDS(rcppAll,file.path(normalizePath("."),jobPath,"results","final_rcppAll.rds"))
saveRDS(parcAll,file.path(normalizePath("."),jobPath,"results","final_parcAll.rds"))

saveRDS(faustPheno,file.path(normalizePath("."),jobPath,"results","final_faustPheno.rds"))
saveRDS(FlowSOMPheno,file.path(normalizePath("."),jobPath,"results","final_FlowSOMPheno.rds"))
saveRDS(FlowORCPheno,file.path(normalizePath("."),jobPath,"results","final_FlowORCPheno.rds"))
saveRDS(phenographPheno,file.path(normalizePath("."),jobPath,"results","final_phenographPheno.rds"))
saveRDS(depechePheno,file.path(normalizePath("."),jobPath,"results","final_depechePheno.rds"))
saveRDS(flowMeansPheno,file.path(normalizePath("."),jobPath,"results","final_flowMeansPheno.rds"))
saveRDS(kmeansPheno,file.path(normalizePath("."),jobPath,"results","final_kmeansPheno.rds"))
saveRDS(rcppPheno,file.path(normalizePath("."),jobPath,"results","final_rcppPheno.rds"))
saveRDS(parcPheno,file.path(normalizePath("."),jobPath,"results","final_parcPheno.rds"))
