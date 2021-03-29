library(Polychrome)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)

#
#modify this path so it points to the location of the simulation results
#
resultsDir <- file.path(normalizePath("."),"simulation")

#
#get color palette, and define a common legend for displaying results.
#
set.seed(77)
myColors <- createPalette(N=13,seedcolors=brewer.pal(7, "Set1"),M=10000)
myColors <- myColors[sample(seq(13),13)]
names(myColors) <- c("k-means","FlowSOM","DEPECHE","Phenograph","flowMeans",
                     "Rclusterpp","FlowSOM\n(Oracle)","FAUST","FAUST\n(Annotated subset)",
                     "Citrus","FlowSOM\n(2*Oracle)","diffCyt","PARC")
myColors[["FlowSOM\n(2*Oracle)"]] <- "#B56340"
myColors[["FAUST\n(Annotated subset)"]] <- "#E0DED6"
myColors[["FlowSOM\n(Oracle)"]] <- "#E816F0"
myColors[["PARC"]] <- "#45FF00"
myColors[["k-means"]] <- "#FC00B1"
myColScale <- scale_colour_manual(name = "Method",values = myColors)
myFillScale <- scale_fill_manual(name = "Method",values = myColors)
legendDF <- data.frame(
              x=rep(c(0,1),13),
              y=seq(2*13),
              Method=c("k-means","FlowSOM","DEPECHE","Phenograph","flowMeans",
                       "Rclusterpp","FlowSOM\n(Oracle)","FAUST","FAUST\n(Annotated subset)",
                       "Citrus","FlowSOM\n(2*Oracle)","diffCyt","PARC")
)
legendDF$Method <- as.factor(legendDF$Method)
legendDF$Method <- relevel(legendDF$Method,ref="FAUST")
pLegend <- ggplot(legendDF,aes(x=x,y=y,color=Method,fill=Method))+
    geom_line(position=position_dodge(width=0.4))+
    theme_classic(base_size=8)+
    theme(
        legend.position="bottom",
        legend.text=element_text(size=8)
    )+
    scale_colour_manual(
        values = myColors,
        guide=guide_legend(
            title="Method",
            override.aes=list(
                size=8,
                shape=15
            )
        )
    )

       
#
#collect results from number of clusters simulation
#
numIter <- 5
baseClusterNum <- 5
computeNCDF <- function(dataSource,methodName){
    ub <- seq(numIter,length(dataSource),by=numIter)
    lb <- ub + 1
    lb <- c(1,lb[-length(ub)])
    methodClusterSD <- methodClusterAvg <- c()
    currentTruthNum <- 5
    for (i in seq(length(lb))) {
        methodClusterAvg <- append(methodClusterAvg,median(dataSource[seq(lb[i],ub[i])]))
        methodClusterSD <- append(methodClusterSD,sd(dataSource[seq(lb[i],ub[i])]))
        names(methodClusterAvg)[length(methodClusterAvg)] <- currentTruthNum 
        names(methodClusterSD)[length(methodClusterSD)] <- currentTruthNum 
        currentTruthNum <- (currentTruthNum + 10) #step size 10 for the sim
    }
    outputDF <- data.frame(
        truth=as.numeric(names(methodClusterAvg)),
        point.median=as.numeric(methodClusterAvg),
        point.sd=as.numeric(methodClusterSD),
        Method=methodName,
        stringsAsFactors=FALSE
    )
    return(outputDF)
}
computeMetricDF <- function(dataSource,methodName){
    ub <- seq(numIter,length(dataSource),by=numIter)
    lb <- ub + 1
    lb <- c(1,lb[-length(ub)])
    methodMetricSD <- methodMetricMed <- c()
    currentTruthNum <- 5
    for (i in seq(length(lb))) {
        methodMetricMed <- append(methodMetricMed,median(dataSource[seq(lb[i],ub[i])]))
        methodMetricSD <- append(methodMetricSD,sd(dataSource[seq(lb[i],ub[i])]))
        names(methodMetricMed)[length(methodMetricMed)] <- currentTruthNum 
        names(methodMetricSD)[length(methodMetricSD)] <- currentTruthNum 
        currentTruthNum <- (currentTruthNum + 10) #step size 10
    }
    outputDF <- data.frame(
        truth=as.numeric(names(methodMetricMed)),
        point.median=as.numeric(methodMetricMed),
        point.sd=as.numeric(methodMetricSD),
        Method=methodName,
        stringsAsFactors=FALSE
    )
    return(outputDF)
}

get_noc_plot <- function(simType) {
    faustNC <- readRDS(file.path(resultsDir,simType,"results","final_faustNC.rds"))
    FlowSOMNC <- readRDS(file.path(resultsDir,simType,"results","final_FlowSOMNC.rds"))
    FlowORCNC <- readRDS(file.path(resultsDir,simType,"results","final_FlowORCNC.rds"))
    flowMeansNC <- readRDS(file.path(resultsDir,simType,"results","final_flowMeansNC.rds"))
    depecheNC <- readRDS(file.path(resultsDir,simType,"results","final_depecheNC.rds"))
    phenographNC <- readRDS(file.path(resultsDir,simType,"results","final_phenographNC.rds"))
    kmeansNC <- readRDS(file.path(resultsDir,simType,"results","final_kmeansNC.rds"))
    rcppNC <- readRDS(file.path(resultsDir,simType,"results","final_rcppNC.rds"))
    parcNC <- readRDS(file.path(resultsDir,simType,"results","final_parcNC.rds"))
    plotDF <- Reduce(rbind,list(
                               faust=computeNCDF(faustNC,"FAUST"),
                               depeche=computeNCDF(depecheNC,"DEPECHE"),
                               phenograph=computeNCDF(phenographNC,"Phenograph"),
                               flowMeansDF=computeNCDF(flowMeansNC,"flowMeans"),
                               FlowSOMDF=computeNCDF(FlowSOMNC,"FlowSOM"),
                               parcDF=computeNCDF(parcNC,"PARC")
                           ))
    plotDF$Method <- as.factor(plotDF$Method)
    plotDF$Method <- relevel(plotDF$Method,ref="FAUST")
    p1 <- ggplot(plotDF,aes(x=truth,y=point.median,color=Method,fill=Method))+
        geom_point()+
        geom_line(position=position_dodge(width=0.4))+
        theme_classic(base_size=6)+
        xlab("True number of clusters")+
        ylab("Median estimated number of clusters")+
        geom_abline(intercept=0, slope=1, linetype="dashed", color="black", size=1)+
        ggtitle("Estimated number of clusters by method")+
        theme(plot.title = element_text(size = 6, face = "bold"))+
        coord_cartesian(ylim=c(0,150))+
        myColScale

    #
    #visualize ARI for all cells
    #
    FlowSOMMetric <- readRDS(file.path(resultsDir,simType,"results","final_FlowSOMAll.rds"))
    FlowORCMetric <- readRDS(file.path(resultsDir,simType,"results","final_FlowORCAll.rds"))
    flowMeansMetric <- readRDS(file.path(resultsDir,simType,"results","final_flowMeansAll.rds"))
    faustMetric <- readRDS(file.path(resultsDir,simType,"results","final_faustAll.rds"))
    faustMetricpheno <- readRDS(file.path(resultsDir,simType,"results","final_faustPheno.rds"))
    phenoMetric <- readRDS(file.path(resultsDir,simType,"results","final_phenographAll.rds"))
    depecheMetric <- readRDS(file.path(resultsDir,simType,"results","final_depecheAll.rds"))
    kmeansMetric <- readRDS(file.path(resultsDir,simType,"results","final_kmeansAll.rds"))
    rcppMetric <- readRDS(file.path(resultsDir,simType,"results","final_rcppAll.rds"))
    parcMetric <- readRDS(file.path(resultsDir,simType,"results","final_parcAll.rds"))

    metricPlotDF <- Reduce(rbind,list(
                                     faust=computeMetricDF(faustMetric,"FAUST"),
                                     faustp=computeMetricDF(faustMetricpheno,"FAUST\n(Annotated subset)"),
                                     depeche=computeMetricDF(depecheMetric,"DEPECHE"),
                                     phenograph=computeMetricDF(phenoMetric,"Phenograph"),
                                     flowMeansDF=computeMetricDF(flowMeansMetric,"flowMeans"),
                                     FlowSOMDF=computeMetricDF(FlowSOMMetric,"FlowSOM"),
                                     FlowORCDF=computeMetricDF(FlowORCMetric,"FlowSOM\n(Oracle)"),
                                     kmeansDF=computeMetricDF(kmeansMetric,"k-means"),
                                     rcppDF=computeMetricDF(rcppMetric,"Rclusterpp"),
                                     parcDF=computeMetricDF(parcMetric,"PARC")
                                 ))
    metricPlotDF$Method <- as.factor(metricPlotDF$Method)
    metricPlotDF$Method <- relevel(metricPlotDF$Method,ref="FAUST")

    p2 <- ggplot(metricPlotDF,aes(x=truth,y=point.median,color=Method))+
        geom_hline(yintercept=0.9,linetype="dashed",size=0.5,color="red")+
        geom_point(position=position_dodge(width=0.4))+
        geom_line(position=position_dodge(width=0.4))+
        theme_classic(base_size=6)+
        xlab("True number of clusters")+
        ylab("Median ARI")+
        coord_cartesian(ylim=c(0,1))+
        ggtitle("Median observed Adjusted Rand Index")+
        theme(plot.title = element_text(size = 6, face = "bold"))+
        myColScale
    plotList <- list()
    plotList <- append(plotList,list(p1))
    plotList <- append(plotList,list(p2))
    return(plotList)
}
pl_noc_gaussian <- get_noc_plot("num_of_clus_gaussian_10_25000")
pl_noc_cytof <- get_noc_plot("num_of_clus_cytof_10_25000")

#
#collection deterministic sim results
#
gatherResults <- function(methodType) {
    simResults <- readRDS(file.path(resPath,paste0(methodType,"_results.rds")))
    resultsMatrix <- matrix(nrow=0,ncol=8)
    colnames(resultsMatrix) <- c("tsprNum","simIter","expectedFoldChange","cvAUC","ciLow","ciHigh","ari","matchedPhenotype")
    for (simNum in seq(length(simResults))) {
        simRes <- simResults[[simNum]]
        iterRes <- c(
            simRes$tsprNum,
            simRes$iterNum,
            simRes$expectedFoldChange,
            simRes$cvAUC,
            simRes$ci[1],
            simRes$ci[2],
            simRes$bestClusterARI,
            as.logical(simRes$spikedPop==simRes$bestClusterPhenotype)
        )
        resultsMatrix <- rbind(resultsMatrix,iterRes)
    }
    resultsDF <- as.data.frame(resultsMatrix)
    resultsDF$Method <- methodType
    return(resultsDF)
}

########################################################################################################
#
#
#generate boxplots for the cvAUC simulation.
#points are jittered with 0.0125 height so there isn't crossing of predicted performance categories
#since we simulate a perfect predictor, warnings are due to a method finding it, recovering perfect
#performance, and then jittering above 1 for cvAUC; we constrain this from happening, and note so
#in the methods section of the manuscript.
#
#
########################################################################################################
#
#generate plots for the gaussian simulation
#
resPath <- file.path(resultsDir,"cvauc_sim_gaussian_10_dim","results")
methodsToCompare <- list.files(resPath)
plotDF <- Reduce(rbind,lapply(gsub("_results.rds","",methodsToCompare),gatherResults))
plotDF$Method <- gsub("faust","FAUST",plotDF$Method)
plotDF$Method <- gsub("phenograph","Phenograph",plotDF$Method)
plotDF$Method <- gsub("cellcnn","CellCNN",plotDF$Method)
plotDF$Method <- gsub("oracle_FlowSOM","FlowSOM\n(Oracle)",plotDF$Method)
plotDF$Method <- gsub("overpart_FlowSOM","FlowSOM\n(2*Oracle)",plotDF$Method)
plotDF$Method <- gsub("dsCitrus","Citrus",plotDF$Method)
plotDF$Method <- gsub("depeche","DEPECHE",plotDF$Method)
plotDF$Method <- gsub("kmeans","k-means",plotDF$Method)
plotDF$Method <- gsub("diffcyt","diffCyt",plotDF$Method)
plotDF$Method <- gsub("rclusterpp","Rclusterpp",plotDF$Method)
plotDF$Method <- as.factor(plotDF$Method)
plotDF$Method <- relevel(plotDF$Method,ref="FAUST")

p_gaussian_AUC <- ggplot(plotDF,aes(x=Method,y=cvAUC,fill=Method))+
    geom_hline(yintercept=0.9,linetype="dashed",size=0.5,color="red")+
    geom_boxplot(outlier.color=NA)+
    geom_jitter(size=0.5, height=0.0125)+
    theme_classic(base_size=6)+
    theme(
        legend.position="bottom",
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
    )+
    ylab("Top cluster 5-fold\ncross validated AUC")+
    xlab("")+
    ggtitle("")+
    myFillScale+
    ylim(c(0,1.0))

p_gaussian_ARI <- ggplot(plotDF,aes(x=Method,y=ari,fill=Method))+
    geom_hline(yintercept=0.9,linetype="dashed",size=0.5,color="red")+
    geom_boxplot(outlier.color=NA)+
    geom_jitter(size=0.5, height=0.0125)+
    theme_classic(base_size=6)+
    theme(
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
    )+
    ylab("Top cluster ARI relative to\nsimulated perfect predictor")+
    xlab("")+
    ggtitle("")+
    myFillScale

#
#generate plots for the cytof simulation
#
resPath <- file.path(resultsDir,"cvauc_sim_cytof_10_dim","results")
methodsToCompare <- list.files(resPath)
plotDF <- Reduce(rbind,lapply(gsub("_results.rds","",methodsToCompare),gatherResults))
plotDF$Method <- gsub("faust","FAUST",plotDF$Method)
plotDF$Method <- gsub("phenograph","Phenograph",plotDF$Method)
plotDF$Method <- gsub("cellcnn","CellCNN",plotDF$Method)
plotDF$Method <- gsub("oracle_FlowSOM","FlowSOM\n(Oracle)",plotDF$Method)
plotDF$Method <- gsub("overpart_FlowSOM","FlowSOM\n(2*Oracle)",plotDF$Method)
plotDF$Method <- gsub("dsCitrus","Citrus",plotDF$Method)
plotDF$Method <- gsub("depeche","DEPECHE",plotDF$Method)
plotDF$Method <- gsub("kmeans","k-means",plotDF$Method)
plotDF$Method <- gsub("diffcyt","diffCyt",plotDF$Method)
plotDF$Method <- gsub("rclusterpp","Rclusterpp",plotDF$Method)
plotDF$Method <- as.factor(plotDF$Method)
plotDF$Method <- relevel(plotDF$Method,ref="FAUST")

p_cytof_AUC <- ggplot(plotDF,aes(x=Method,y=cvAUC,fill=Method))+
    geom_hline(yintercept=0.9,linetype="dashed",size=0.5,color="red")+
    geom_boxplot(outlier.color=NA)+
    geom_jitter(size=0.5, height=0.0125)+
    theme_classic(base_size=6)+
    theme(
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        plot.title=element_text(hjust=0.5)
    )+
    ylab("Top cluster 5-fold\ncross validated AUC")+
    xlab("")+
    ggtitle("")+
    myFillScale+
    ylim(c(0,1))

p_cytof_ARI <- ggplot(plotDF,aes(x=Method,y=ari,fill=Method))+
    geom_hline(yintercept=0.9,linetype="dashed",size=0.5,color="red")+
    geom_boxplot(outlier.color=NA)+
    geom_jitter(size=0.5, height=0.0125)+
    theme_classic(base_size=6)+
    theme(
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
    )+
    ylab("Top cluster ARI relative to\nsimulated perfect predictor")+
    xlab("")+
    ggtitle("")+
    myFillScale


#
#collect and display results
#
gridOut <- plot_grid(
    (pl_noc_gaussian[[1]]+theme(legend.position="none")+ggtitle("")),
    (pl_noc_gaussian[[2]]+theme(legend.position="none")+ggtitle("")),
    (p_gaussian_AUC+theme(legend.position="none")+ggtitle("")),
    (p_gaussian_ARI+theme(legend.position="none")+ggtitle("")),
    (pl_noc_cytof[[1]]+theme(legend.position="none")+ggtitle("")),
    (pl_noc_cytof[[2]]+theme(legend.position="none")+ggtitle("")),
    (p_cytof_AUC+ggtitle("")+theme(legend.position="none")),
    (p_cytof_ARI+ggtitle("")+theme(legend.position="none")),
    labels=c("A","","B","","C","","D",""),
    ncol=2,
    nrow=4,
    label_size=9
)

shared_legend <- get_legend(pLegend)
pOut <- plot_grid(
    gridOut,
    shared_legend,
    ncol=1,
    rel_heights=c(1,(1/6))
)

cowplot::save_plot(
             filename="./figure_simulation.png",
             plot=pOut,
             base_width=6,
             base_height=9
         )
