# Plotting functions for choice history bias analysis

#' Plot psychometric curves
#'
#' Plot the proportion of right responses.
#' Also fit Weibull and show the fit
#'
#' @export
PlotPsychometricCurves <- function(inputData,
                                   figureWidth=21.37078, 	# Width of the plot
                                   figureHeight=9.034483, 	# Height of the plot
                                   plotByCondition=T,
                                   plotErrorBars=T,		# To plot error bars
                                   findAsymptotes=F,
                                   collapseByRun=T,
                                   lineColor="#c6dbef",	# Color of psychometric curve lines
                                   dotColor="#09519c",		# Color of dots for geom_point
                                   nrow=3) {
  ## Plot proportion correct of responding right to stimuli presented either to left or to right
  pcRight <- droplevels(inputData)
  ## Label left responses to right gratings with negative contrast. Right responses to right gratings will remain with positive sign
  pcRight[pcRight$VisualField == 1,]$Contrast <- pcRight[pcRight$VisualField == 1,]$Contrast * -1
  ## Summary of responses to gratings presented in the right
  pcRightSummary <- ddply(pcRight, .(SubjectID, SessionID, Condition, VisualField, Contrast), summarise,
                          nRightResp = sum(Response == 2),
                          nLeftResp = sum(Response == 1),
                          nRStim = sum(VisualField == 2),
                          nLStim = sum(VisualField == 1))
  pcRightSummary <- ddply(pcRightSummary, .(SubjectID, SessionID, Condition, VisualField, Contrast), summarise,
                          pRightCorrect = nRightResp / (nRStim + nLStim),
                          nYesR=nRightResp,
                          nNoR=nLeftResp)
  ## Convert contrast into %
  pcRightSummary$Contrast <- pcRightSummary$Contrast * 100

  ## Plot collapsed by SessionID
  if (findAsymptotes==TRUE) {
    models <- dlply(pcRightSummary, c("SubjectID", "SessionID", "Condition"), .fun=MakeModelWithAsymptote)  ## NEED TO FIX THIS Function. It is currently disabled
  } else {
    print("Fitting Probit function to data")
    models <- dlply(pcRightSummary, c("SubjectID", "SessionID", "Condition"), .fun=FitProbit)
  }
  predvals <- ldply(models, .fun=PredictvalsProbit, xvar="Contrast", yvar="pRightCorrect", type="response")
  #browser()
  g <- ggplot(data = pcRightSummary, aes(x = Contrast, y = pRightCorrect))
  g <- g + theme_few()
  g <- g + geom_vline(xintercept = 0.0, size = 0.2, colour = "grey30", linetype = "dashed")
  g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())
  g <- g + theme(axis.line = element_line(colour = "#a9a9a9", size = 0.3))
  g <- g + theme(axis.ticks.x = element_line(colour = "#a9a9a9", size = 0.3))

  g <- g + geom_line(data=predvals, aes(x=Contrast, y=pRightCorrect, group=SessionID), color=lineColor, size=0.3)
  if (plotErrorBars) g <- g + stat_summary(aes(group=SubjectID), fun.data = "mean_cl_boot", B=1000, conf.int = 0.68, geom = "errorbar", size = 0.5, width = 0.0, color='#9ecae1')
  g <- g + geom_point(aes(group=SubjectID), stat="summary", fun.y="mean", size=2.5, color=dotColor)

  if (plotByCondition) {
    nRunsLabel <- ddply(inputData, .(SubjectID), summarise, label=paste("r=", length(unique(SessionID)), sep=""),
                        nTrialsLab=paste("t=", length(SubjectID), sep=""))
  } else {
    nRunsLabel <- ddply(inputData, .(SubjectID), summarise, label=paste("r=", length(unique(SessionID)), sep=""),
                        nTrialsLab=paste("t=", length(SubjectID), sep=""))
  }

  ## x and y coordinates where labels will be shown
  ## y position of the label that shows the number of runs
  nRunsLabel$TextX <- min(pcRightSummary$Contrast) + 1.6
  nRunsLabel$TextY <- 0.90
  nRunsLabel$TextYTrls <- 0.97

  g <- g + geom_text(data=nRunsLabel, aes(x=TextX, y=TextY, label=label, group=SubjectID), color="#c8c8c8", size=3)
  g <- g + geom_text(data=nRunsLabel, aes(x=TextX, y=TextYTrls, label=nTrialsLab, group=SubjectID), color="#c8c8c8", size=3)

  g <- g + scale_color_brewer()
  #g <- g + geom_hline(yintercept = 0.5, size = 0.2, colour = "grey50", linetype = "dashed")
  g <- g + xlab('Contrast (%)') + ylab('Proportion Right Responses')

  #browser()
  if (plotByCondition) {
    g <- g + facet_grid(Condition~SubjectID, labeller=ConditionLabels)
  }
  else {
    g <- g + facet_wrap(~SubjectID, nrow=nrow)
  }
  dev.new(width=figureWidth, height=figureHeight)
  return(g)
}


#' Plot whole run as Tuft's sparklines of left/right choice and responses
#'
#' @export
PlotSparklineOfRun <- function(runData,   	 	# Input data in the format of "glmData"
                               plotResults=F, 	# If False, returns list of ggplot objects. If True, plots those objects, each subject in separate window
                               plotInColor=T)		# Success and failure are marked by two different colors. Otherwise, black and white
{
  ## Add column with trial numbers
  #runData <- droplevels(ddply(runData, mutate, TrialID=1:length(Response)))
  runData <- cbind(runData, TrialID=1:length(runData$Response))
  ## Define colors for Left and Right choices
  colChooseRight <- '#2b83ba'
  colChooseLeft <- '#d7191c'

  ## Generate ggplot graphs as binary sparklines
  g <- ggplot(data=runData, aes(x=TrialID, y=Response)) + theme_few()
  g <- g + geom_linerange(subset=.(Response==1), aes(ymin=1.5, ymax=Response, color="Left"), size=0.3) +
    geom_linerange(subset=.(Response==2), aes(ymin=1.5, ymax=Response, color="Right"), size=0.3) +
    geom_linerange(subset=.(CorrIncorr==0), aes(ymin=1.5, ymax=Response, color="Fail"), size=0.3) +
    scale_color_manual(values=c("black", colChooseLeft, colChooseRight)) +
    scale_x_continuous(expand=c(0.01,0.01)) +
    scale_y_continuous(breaks=c(1,2), labels=c("L","R")) +
    theme(axis.ticks.x=element_blank()) +
    theme(legend.title=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
    ylab("Response") +
    xlab("")

  if (plotResults) print(g)
  return(g)

}


#' Plot both contrast and choice history weights
#'
#' A beautiful plot of all weights (contrast and history). Contrast weights are
#' connected together, but history weights are not.
#' @export
PlotContrastAndHistoryWeights <- function(regularizedWeights,
                                          plotByCondition=T,		# Plot by condition
                                          legendPos=c(0.15,0.75),
                                          plotInColor=FALSE,
                                          showAxisLabels=TRUE,
                                          plotFigure=TRUE,      # Plot figure or return ggplot object
                                          figureWidth=6.382023,
                                          figureHeight=8.678161) {

  if (plotByCondition)
    idealSbjParams <- ddply(regularizedWeights, .(SubjectID, Parameter, Condition), summarise,
                            MeanWeight=mean(Weight),
                            std=sd(Weight, na.rm=TRUE),
                            n=sum(!is.na(Weight)),
                            se=std/sqrt(n))
  else idealSbjParams <- ddply(regularizedWeights, .(SubjectID, Parameter), summarise,
                               MeanWeight=mean(Weight),
                               std=sd(Weight, na.rm=TRUE),
                               n=sum(!is.na(Weight)),
                               se=std/sqrt(n))

  ## Sort contrast and history weights in the order that will make it intuitive to read the plot
#   paramNames <- levels(idealSbjParams$Parameter)
#   contrastNames <- sort(paramNames[grep("c0",paramNames)])
#   biasNames <- paramNames[!paramNames %in% contrastNames]
  idealSbjParams$Parameter <- factor(idealSbjParams$Parameter, levels=SortParams(idealSbjParams$Parameter))

  ## Add grouping parameter that will be used to plot different weights in different colors
  idealSbjParams$plotColor <-"Contrast"
  idealSbjParams$plotColor[idealSbjParams$Parameter=="(Intercept)"] <- "Intercept"
  idealSbjParams$plotColor[grep(successColName, idealSbjParams$Parameter)] <- "PrevSuccess"
  idealSbjParams$plotColor[grep(failColName, idealSbjParams$Parameter)] <- "PrevFail"
  #colorValues <- c("#a4a4a4",prevSuccessColor, prevFailColor, "#6f6f6f")
  colorValues <- c("#3A79B6","#3A79B6", "#3A79B6", "#3A79B6")

  names(colorValues) <- c("Intercept", "PrevSuccess", "PrevFail", "Contrast")
  idealSbjParams$plotOrder <- 4
  idealSbjParams$plotOrder[idealSbjParams$Parameter=="(Intercept)"] <- 1
  idealSbjParams$plotOrder[idealSbjParams$Parameter=="PrevSuccess"] <- 2
  idealSbjParams$plotOrder[idealSbjParams$Parameter=="PrevFail"] <- 3

  p <- ggplot(idealSbjParams, aes(x=Parameter, y=MeanWeight)) +
    theme_publish1() +
    theme(legend.position = legendPos) +
    geom_hline(yintercept = 0, size = 0.3, colour = "grey70", linetype = "dashed") +
    stat_summary(aes(group=Parameter), fun.data = "mean_cl_boot", fun.args=list(conf.int=0.68, B=1000), geom = "errorbar", size = 0.5, width = 0.0, color="#3A79B6") +
    geom_blank() +
    geom_line(data=subset(idealSbjParams, plotColor="Contrast"), aes(group=1), stat="summary", fun.y="mean", color='#3A79B6', size=1) +
    geom_point(aes(group=Parameter), size=7, color="white", stat="summary", fun.y="mean") +
    geom_point(aes(group=Parameter, color=plotColor), size=4, stat="summary", fun.y="mean") +
    scale_colour_manual(values=colorValues) +
    coord_cartesian(ylim = c(-1.2, 5.5)) +
    scale_y_continuous(breaks=seq(-1, 5.5, 1))

  if (showAxisLabels) {
    p <- p + theme(axis.text.x=element_text(angle = 90, vjust=0.5))
  } else {
    p <- p + theme(axis.text=element_blank(), axis.title=element_blank())
  }

  if (plotByCondition) p <- p + facet_grid(~Condition, labeller=ConditionLabels, scales="free")

  if (plotFigure) {
    dev.new(width=figureWidth, height=figureHeight)
    print(p)
  } else {
    return(p)
  }
}


#' A scatterplot of choice history biases
#'
#' Plot 'prevFail' vs 'prevSuccess' weights in forms of scatterplot.
#' @export
HistoryWeightsScatterplot <- function(weights,  				# Regularized weights
                                      figureWidth=15.516854, 	# Width of the plot
                                      figureHeight=5.689655, 	# Height of the plot
                                      xlims=c(-3,3),
                                      ylims=c(-3,3),
                                      conditionsToPlot, 		# Condition numbers to be plotted
                                      plotByCondition=TRUE,		# Plot by Condition or collapsed across Condition
                                      plotSubjectMeans=FALSE,	# Compute means across Sessions for each subject and plots those means
                                      plotSubjectNames=FALSE, 	# Plot names of subjects next to each dot. Better to use when plotSubjectMeans=TRUE
                                      plotMeans=FALSE         # Mean weights for each Condition
) {


  if (!missing(conditionsToPlot)) weights <- droplevels(subset(weights, Condition %in% conditionsToPlot))
  ## If asked for, only subset of conditions to plot

  historyWeights <- subset(weights, Parameter=="PrevCorr1" | Parameter=="PrevFail1")
  biases <- cast(historyWeights, SubjectID+SessionID+Condition ~ Parameter, value=.(Weight))

  ## Get means for each subject across all sessions
  if (plotSubjectMeans) {
    biases <- ddply(biases, .(SubjectID, Condition), summarise,
                    SePrevCorr1=sd(PrevCorr1)/sqrt(length(PrevCorr1)), 	# Standard error
                    SePrevFail1=sd(PrevFail1)/sqrt(length(PrevFail1)),	# Standard error
                    PrevCorr1=mean(PrevCorr1),
                    PrevFail1=mean(PrevFail1)
    )
  }

  plottingColor <- "#a9a9a9"
  p <- ggplot(data = biases, aes(x=PrevCorr1, y=PrevFail1))
  ## If only means for each subject are plotted, draw error bars first
  if (plotSubjectMeans) {
    p <- p + geom_errorbar(aes(ymin=PrevFail1-SePrevCorr1, ymax=PrevFail1+SePrevCorr1), color=plottingColor)
    p <- p + geom_errorbarh(aes(xmin=PrevCorr1-SePrevFail1, xmax=PrevCorr1+SePrevFail1), color=plottingColor)
    p <- p + scale_fill_manual(values=colorRampPalette(brewer.pal(12, "Set1"))(15))
    p <- p + geom_point(aes(fill=SubjectID), size=3.5, shape=21, color="white")
  } else {
    p <- p + geom_point(size=3.5, shape=21, color="white", fill="#6a6a6a")
  }

  p <- p + geom_hline(yintercept = 0.0, size = 0.3, colour = plottingColor, linetype = "solid") +
    geom_vline(xintercept = 0.0, size = 0.3, colour = plottingColor, linetype = "solid") +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    xlab("Weight success (Bs)") + ylab("Weight failure (Bf)") +
    theme_publish1()

  if (plotSubjectNames) {
    p <- p + geom_text(aes(label=SubjectID), size=3, color=plottingColor, vjust=2.5)
  }

  if (plotMeans) {
    biasMeans <- ddply(biases, .(Condition), summarise,
                       SePrevCorr1=sd(PrevCorr1)/sqrt(length(PrevCorr1)), 	# Standard error
                       SePrevFail1=sd(PrevFail1)/sqrt(length(PrevFail1)),	# Standard error
                       PrevCorr1=mean(PrevCorr1),
                       PrevFail1=mean(PrevFail1)
    )
    p <- p + geom_point(data=biasMeans, size=7, alpha=0.5)
  }
  if (plotByCondition) {
    p <- p + facet_grid(~Condition, labeller=ConditionLabels)
  }
  p <- p + coord_fixed(xlim=xlims, ylim=ylims)
  dev.new(width=figureWidth, height=figureHeight)
  return(p)

}



#' Scatterplot of choice history weights grouped by Education
#'
#' This is a plot for the manuscript
#' @export
HistoryScatterplotByEducation <- function(weights,    			# Regularized weights
                                          figureWidth=15.516854, 	  # Width of the plot
                                          figureHeight=5.689655, 	  # Height of the plot
                                          conditionsToPlot, 		    # Condition numbers to be plotted
                                          legendPos=c(0.15,0.75),   # Legend position
                                          plotMeansByEducation = F, # Plot mean weights grouped by Education
                                          plotByCondition=TRUE,		  # Plot by Condition or collapsed across Condition
                                          plotSubjectMeans=FALSE,	  # Compute means across Sessions for each subject and plots those means
                                          plotSubjectNames=FALSE 	  # Plot names of subjects next to each dot. Better to use when plotSubjectMeans=TRUE
) {

  # Get Education
  sbjEducation <- SubjectDemographics(getEducation = TRUE, getAge = FALSE)
  # Join demographics and history weights
  weights <- droplevels(join(weights, sbjEducation))
  if (!missing(conditionsToPlot)) weights <- droplevels(subset(weights, Condition %in% conditionsToPlot))
  ## If asked for, only subset of conditions to plot

  historyWeights <- subset(weights, Parameter=="PrevCorr1" | Parameter=="PrevFail1")
  biases <- cast(historyWeights, SubjectID+SessionID+Condition+Education ~ Parameter, value=.(Weight))

  ## Get means for each subject across all sessions
  if (plotSubjectMeans) {
    biases <- ddply(biases, .(SubjectID, Condition, Education), summarise,
                    SePrevCorr1=sd(PrevCorr1)/sqrt(length(PrevCorr1)), 	# Standard error
                    SePrevFail1=sd(PrevFail1)/sqrt(length(PrevFail1)),	# Standard error
                    PrevCorr1=mean(PrevCorr1),
                    PrevFail1=mean(PrevFail1))
  }
  plottingColor <- "#a9a9a9"
  errorbarColor <- "grey50"

  if (plotMeansByEducation) {
    geomPointSize <- 5
    alpha <- 1
  } else {
    geomPointSize <- 5
    alpha <- 1
  }
  p <- ggplot(data = biases, aes(x=PrevCorr1, y=PrevFail1))
  p <- p + geom_hline(yintercept = 0.0, size = 0.3, colour = plottingColor, linetype = "dashed") +
    geom_vline(xintercept = 0.0, size = 0.3, colour = plottingColor, linetype = "dashed")
  ## If only means for each subject are plotted, draw error bars first
  if (plotSubjectMeans) {
    #     if (!plotMeansByEducation)  {
    p <- p + geom_errorbar(aes(ymin=PrevFail1-SePrevCorr1, ymax=PrevFail1+SePrevCorr1), color=errorbarColor, width=0, size=0.3)
    p <- p + geom_errorbarh(aes(xmin=PrevCorr1-SePrevFail1, xmax=PrevCorr1+SePrevFail1), color=errorbarColor, height=0, size=0.3)
    #     }
    p <- p + scale_fill_manual(values=c("#fdae61", "#d7191c"))  # blue 2c7bb6
    p <- p + geom_point(aes(fill=Education), size=geomPointSize, shape=21, color="white", alpha=alpha)
  } else {
    p <- p + geom_point(aes(fill=Education), size=geomPointSize, shape=21, color="white", alpha=alpha)
  }

  p <- p + coord_fixed(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
    #p <- p + coord_fixed(xlim = c(-2.2, 2.2), ylim = c(-2.5, 1)) +
    xlab("Success bias") + ylab("Failure bias") +
    theme_publish1() +
    theme(legend.position = legendPos)
    #geom_rangeframe()
  #theme(panel.grid=element_blank(), panel.border=element_blank()) +
  #theme(axis.line = element_line(colour="grey", size=0.5)) +
  #theme(axis.text=element_text(size=12, face="bold"))

  if (plotSubjectNames) {
    p <- p + geom_text(aes(label=SubjectID), size=3, color=plottingColor, vjust=2.5)
  }
  # Show means by education
  #if (!plotMeansByEducation) {
  meansByEducation <- ddply(biases, .(Education), summarise,
                            SePrevCorr1=sd(PrevCorr1)/sqrt(length(PrevCorr1)),   # Standard error
                            SePrevFail1=sd(PrevFail1)/sqrt(length(PrevFail1)),	# Standard error
                            PrevCorr1=mean(PrevCorr1),
                            PrevFail1=mean(PrevFail1))
  #p <- p + geom_errorbar(data=meansByEducation, aes(ymin=PrevFail1-SePrevCorr1, ymax=PrevFail1+SePrevCorr1), color=errorbarColor, width=0)
  #p <- p + geom_errorbarh(data=meansByEducation, aes(xmin=PrevCorr1-SePrevFail1, xmax=PrevCorr1+SePrevFail1), color=errorbarColor, height=0)
  p <- p + geom_point(data=meansByEducation, aes(group=Education, fill=Education), size=9, shape=21, color="white", alpha=0.5)
  #}

  # Separate by Condition
  if (plotByCondition) {
    p <- p + facet_grid(~Condition, labeller=ConditionLabels)
  }
  return(p)
}

#' Function to plot Variance Inflation Factor (VIF)
#' @export
PlotVIFs <- function(allWeights) {
  vifsOnly <- allWeights[!is.na(allWeights$Vif),] # Pick rows that only contain VIFs
  vifData <- vifsOnly[c("SubjectID", "SessionID", "Condition", "Parameter", "Vif")]
  ggplot(data=vifData, aes(x=Parameter, y=Vif)) +
    theme_few() +
    ggtitle("Variance Inflation Factor (VIF). \n Non-shaded panels have fail/switch probability set at 80% (Condition 2). \n Shaded panels have success/stay probability set at 80% (Condition 3)") +
    geom_bar(stat="identity", position="dodge") +
    geom_rect(data = subset(vifData, Condition == 3), aes(fill=Condition), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.02) +
    geom_hline(yintercept = 1, size = 0.5, colour = "orange", linetype = "dashed") +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    theme(legend.position="none") +
    facet_wrap(SubjectID~Condition~SessionID, ncol=6)
}

#' Function to plot stimulus intensities during switching or staying behaviours
#' @export
PlotContrastPairs <- function(contrastPairs, plotTitle) {
  ggplot(data=contrastPairs, aes(x=factor(failIntensity*100), y= factor(nextIntensity*100))) +
    geom_tile(aes(fill= 100*prcntSwitch)) +
    theme_few() +
    facet_wrap(~SubjectID, scale="free") +
    #scale_fill_gradient(low="#f7f500", high="#d80570", na.value="black") +
    #scale_fill_gradient2(low="#FEE8C8", mid="#FDBB84", high="#E34A33") +
    #scale_fill_gradient(low="#FEE8C8", high="#E34A33", na.value="black") +
    scale_fill_gradientn(colours=c("#FEF0D9", "#FDCC8A", "#E34A33", "#B30000")) +
    labs(fill="% switches") +
    scale_x_discrete("Contrast when failed (%)", expand=c(0,0)) +
    scale_y_discrete("Contrast when switched (%)", expand=c(0,0)) +
    geom_text(aes(label=nTotalSwitches, color=nTotalSwitches), size=4) +
    scale_color_gradient(low="white", high="black", guide="none") +
    ggtitle(plotTitle)
}


#' History scatterplot that shows how history biases adapt
#'
#' Formerly "HistoryMigrationScatterplot"
#' Scatterplot of choice history weights grouped by Education
#' This is a plot for the manuscript
#' @export
HistoryAdaptationScatterplot <- function(weights,  				# Regularized weights
                                         figureWidth=15.516854, 	# Width of the plot
                                         figureHeight=5.689655, 	# Height of the plot
                                         conditionsToPlot, 		# Condition numbers to be plotted
                                         sourceCondition=NA,		# Condition from which arrow will extend to show transition
                                         plotErrorBars=FALSE,		#
                                         plotByCondition=FALSE,	# Plot by Condition or collapsed across Condition
                                         plotSubjectMeans=FALSE,	# Compute means across Sessions for each subject and plots those means
                                         plotSubjectNames=FALSE,  	# Plot names of subjects next to each dot. Better to use when plotSubjectMeans=TRUE
                                         plotLegend=TRUE)
{
  library(grid)
  plottingColor <- "#a9a9a9"
  errorbarColor <- "grey90"

  ## Add subject degree
  sbjDemographs <- SubjectDemographics()
  weights <- merge(weights, sbjDemographs, by='SubjectID')

  ## If asked for, only subset of conditions to plot
  if (!missing(conditionsToPlot)) weights <- droplevels(subset(weights, Condition %in% conditionsToPlot))
  ## If more than one condition, find common subjects
  ## Get subject labels for both conditions
  sbjInFirstCondition <- levels(droplevels(weights[weights$Condition==conditionsToPlot[1],]$SubjectID))
  sbjInSecondCondition <- levels(droplevels(weights[weights$Condition==conditionsToPlot[2],]$SubjectID))
  ## Find common subjects for both conditions. Only they will be plotted
  subjectsToPlot <- intersect(sbjInFirstCondition, sbjInSecondCondition)
  ## Select those subjects for further data processing and plotting
  weights <- droplevels(subset(weights, SubjectID %in% subjectsToPlot))

  historyWeights <- subset(weights, Parameter=="PrevCorr1" | Parameter=="PrevFail1")
  biases <- cast(historyWeights, SubjectID+SessionID+Condition+Education ~ Parameter, value=.(Weight))

  ## Get means for each subject across all sessions
  if (plotSubjectMeans) {
    biases <- ddply(biases, .(SubjectID, Condition, Education), summarise,
                    SePrevCorr1=sd(PrevCorr1)/sqrt(length(PrevCorr1)), 	# Standard error
                    SePrevFail1=sd(PrevFail1)/sqrt(length(PrevFail1)),	# Standard error
                    PrevCorr1=mean(PrevCorr1),
                    PrevFail1=mean(PrevFail1)
    )
    ## Means across subjects for each condition. This will be shown as a big vector
    overallMeans <- ddply(biases, .(Condition), summarise, MeanPrevCorr=mean(PrevCorr1), MeanPrevFail=mean(PrevFail1))
    ## To have the direction of the arrow set correctly, sort data to have "Condition" in proper order
    ## This only needs to be done if the Source condition has number greater than Target condition
    if (conditionsToPlot[1]>conditionsToPlot[2]) {
      biases <- biases[with(biases, order(SubjectID, -Condition)),]
      overallMeans <- overallMeans[with(overallMeans, order(-Condition)),]
    }
  }

  ## If it is requested to
  if (!missing(sourceCondition)) {
    sourceTargetOrder <- c(sourceCondition, conditionsToPlot[conditionsToPlot!=sourceCondition])
    ## Order factors according to sourceTargetOrder
    biases$Condition <- factor(biases$Condition, levels=sourceTargetOrder)
  } else {
    ## If source condition is not specified, the smallest number will be
    ## the source and the biggest number destination (where arrow will go)
    sourceTargetOrder <- sort(conditionsToPlot)
  }

  #browser()

  p <- ggplot(data = biases, aes(x=PrevCorr1, y=PrevFail1))
  p <- p + geom_hline(yintercept = 0.0, size = 0.3, colour = plottingColor, linetype = "dashed") +
    geom_vline(xintercept = 0.0, size = 0.3, colour = plottingColor, linetype = "dashed")
  ## If only means for each subject are plotted, draw error bars first
  if (plotSubjectMeans) {
    if (plotErrorBars) {
      p <- p + geom_errorbar(aes(ymin=PrevFail1-SePrevCorr1, ymax=PrevFail1+SePrevCorr1), color=errorbarColor)
      p <- p + geom_errorbarh(aes(xmin=PrevCorr1-SePrevFail1, xmax=PrevCorr1+SePrevFail1), color=errorbarColor)
    }
    #p <- p + scale_fill_brewer(palette="Set1")
    #p <- p + scale_fill_few()
    p <- p + scale_fill_manual(values=c("#fdae61", "#d7191c", "#2c7bb6"))
    p <- p + scale_colour_manual(values=c("#fdae61", "#d7191c", "#2c7bb6"))
    conditionLabels <- ConditionLabels("Condition", sourceTargetOrder)
    #ifelse (sourceTargetOrder[1]<sourceTargetOrder[2], shapeValues <- c(1,19), shapeValues <- c(19,1))
    #browser()
    p <- p + scale_shape_manual(values=c(1,19), labels= conditionLabels, name="Condition")
    #p <- p + scale_size_manual(values=c(5,5))
    #p <- p + geom_point(aes(fill=Education), size=3.5, shape=21, color="white", alpha=0.8)
    #p <- p + geom_point(aes(color=Education, shape=as.factor(Condition)), size=5, alpha=0.8)
    p <- p + geom_point(aes(color=Education, shape=as.factor(Condition)), size=5)
    p <- p + geom_path(aes(group=SubjectID), arrow = arrow(length = unit(0.2,"cm"), type="closed"), colour="grey50")
    #p <- p + geom_point(subset=.(Condition==sourceTargetOrder[1]), size=4.1, color="white")
    ## Show average change by computing means of failure and success weights for each condition
    p <- p + geom_path(data=overallMeans, aes(group=1, x= MeanPrevCorr, y= MeanPrevFail), arrow = arrow(length = unit(0.5,"cm"), type="closed"), size=7, alpha=0.3)
    #p <- p + geom_point(aes(fill=Education, shape=as.factor(Condition)), size=3.5)
  } else {
    p <- p + geom_point(aes(fill=Education), size=3.5, shape=21, color="white")
  }

  #geom_density2d(color="#6a6a6a", size=0.3) +
  #geom_point(size=5, color="white") +
  #geom_point(size=3, color="#6a6a6a") +
  #coord_cartesian(xlim = c(-2.2, 2.2), ylim = c(-1.9, 1)) +
  p <- p + coord_fixed(xlim = c(-2.2, 2.2), ylim = c(-1.9, 1)) +
    xlab("Success bias") + ylab("Failure bias") +
    theme_publish1()

  if (plotSubjectNames) {
    p <- p + geom_text(aes(label=SubjectID), size=3, color=plottingColor, vjust=2.5)
  }

  if (plotByCondition) {
    p <- p + facet_grid(~Condition, labeller=ConditionLabels)
  }

  if (!plotLegend) {
    p <- p + theme(legend.position="none")
  }
  #p <- p + ggtitle(paste(conditionLabels[1], "VS", conditionLabels[2]))
  # Compute paired t.tests for the manuscript
  print(t.test(data=biases, PrevFail1~Condition, paired=T))
  print(t.test(data=biases, PrevCorr1~Condition, paired=T))
  browser()
  return(p)
}

#' Show decline in visual sensitivity from choice history biases
#'
#' Displays decline in visual sensitivity as a result of choice history biases.
#' It can show reduction in contrast sensitivity (75% contrast threshold) or
#' change in slope of psychometric function (in case of probit, smaller slope
#' means less sensitivity).
#'
#' @param simThAndSlope simulated threshold and slope computed using function NoBiasVsSubjectBias (see example below).
#' @param whatToPlot either decline of slope or threshold can be plotted
#' @param whatToReturn return either 'plot' as ggplot object or 'data' prepared for plotting, which can be used for other operations, or 'stats'. 'stats'
#' parameter shows p value of one-sampled test of median change (whether sensitivity decline is significantly different from 0)
#'
#' @examples
#' load('20150919_allWeights_RIKEN_UCL_Stanford_cond1n13.RData',verbose=TRUE)
#' regWeights <- droplevels(subset(allWeights, Regularized=='yes'))
#' oneSbjWeights <- droplevels(subset(regWeights, SubjectID %in% c('s001', 's002', 's003', 's009', 's008', 's007')))
#' thAndSlope <- NoBiasVsSubjectBias(regWeights, nTrialsPerContrast=10, B=5) # B=5 is for fast simulation, for meaningful results use B>500.
#' PlotSensitivityDecline(thAndSlope) #
#' @export
PlotSensitivityDecline <- function(simThAndSlope,
                                   whatToPlot='slope', # can be either 'slope' or 'threshold'
                                   whatToReturn='plot') # can be plot, data or stats
{
  if (whatToPlot == 'slope') colName <- 'slope' else colName <- 'th75'
  # Make subject's (biased) and unbiased sensitivity column into two separate columns
  sensitivity <- cast(simThAndSlope, SubjectID+Degree+SimulationID~IsSubject, value=colName)
  # Clearer names for columns
  names(sensitivity)[names(sensitivity) == "FALSE"] <- "UnbiasedSensitivity"
  names(sensitivity)[names(sensitivity) == "TRUE"] <- "BiasedSensitivity"
  # Compute change in sensitivity
  sensitivity$Change <- 1-(sensitivity$UnbiasedSensitivity/sensitivity$BiasedSensitivity)
  # In case of probit slope, smaller slope means worse performance (and vice versa for th)
  # So need to invert the sign
  if (whatToPlot == 'slope') sensitivity$Change <- -1 * sensitivity$Change

  # Prepare for plotting
  sensitivityForPlot <- sensitivity
  # Median sensitivity change for each subject
  sensitivityForPlot <- ddply(sensitivityForPlot, .(SubjectID), summarise,
                              UnbiasedSensitivity=median(UnbiasedSensitivity),
                              BiasedSensitivity=median(BiasedSensitivity),
                              MedianChange=median(Change))
  # Order subjects by sensitivity change
  orderedSubjects <- with(sensitivityForPlot, SubjectID[order(MedianChange)])
  sensitivityForPlot$SubjectID <- factor(sensitivityForPlot$SubjectID, levels=orderedSubjects)
  sensitivity$SubjectID <- factor(sensitivity$SubjectID, levels=orderedSubjects)

  if (whatToReturn=='plot') {
    require(scales)
    gg <- ggplot(sensitivityForPlot, aes(x=SubjectID, y=MedianChange)) +
      theme_publish1() +
      theme(axis.text.y = element_text(size=8)) +
      xlab('Subject') + ylab('Decline in visual sensitivity (%)') +
      geom_hline(yintercept=0, linetype='dashed', color='grey60', size=0.1) +
      coord_cartesian() +
      scale_y_continuous(labels=percent, breaks=pretty_breaks(n=10)) +
      coord_flip() +
      geom_blank() +
      stat_summary(data=sensitivity, aes(x=SubjectID, y=Change, group=SubjectID),
                   fun.data="median_cl_boot", geom='errorbar', width=0.0, color = 'grey80') +
      geom_point(size=2, color='grey20')
    return(gg)
  }

  if (whatToReturn=='data') {
    return(sensitivityForPlot)
  }

  if (whatToReturn=='stats') {
    stat <- function(dat) {
      # Wilcoxon one-sample median test
      wilcoxonTest <- wilcox.test(dat$Change)
      #browser()
      data.frame(SubjectID=unique(dat$SubjectID), pVal=wilcoxonTest$p.value)
    }
    statDat <- ddply(sensitivity, .(SubjectID), .fun=stat)
    return(statDat)
  }

}


######################################################################################################
## Plot psychometric curves separated into preceding choice being left or right
## Blue color on the plot means preceding choice was right (that is, subject stayed on the same side,
## cause we are plotting proportion right responses)
## Red color means that preceding choice was left (subject switched)
PlotByPrecedingChoice <- function(inputData, 				      # This data should have glmData structure
                                  geomPointSize=3.5,		  # Size of geom_point
                                  showMidPoints=T, 		    # Shows vertical and horizontal lines at mid points on x and y axes
                                  confIntData,            # Simulated data in glmData format that will be used to show confidence intervals
                                  figureWidth=21.37078, 	# Width of the plot
                                  figureHeight=9.034483,  # Height of the plot
                                  plotFigure=TRUE) 	        # Plot figure or just return ggplot object otherwise
{
  # Function to mark trials by previous choice (L or R)
  # Adds new column called PrecedingChoice which is 1 when
  # preceding choice was L or 2 when R
  MarkTrialsByChoice <- function(dat) {
    markedTrials <- data.frame()	# Trials marked as preceded by L or R choice will be placed here
    for (ixSubject in levels(droplevels(dat$SubjectID))) {
      oneSubjectData <- droplevels(subset(dat, dat$SubjectID == ixSubject))
      #subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))
      #for (ixCondition in unique(oneSubjectData$Condition)) {
      for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
        oneSessionData <- subset(oneSubjectData, (oneSubjectData$SessionID==ixSession))
        ## Find trials that were preceded by Left choice
        ixLeftResponse <- which(oneSessionData$Response==1)
        ## If last trial was included, exclude it cause we cannot add 1 more trial to the session
        #browser()
        if (ixLeftResponse[length(ixLeftResponse)] == nrow(oneSessionData)) ixLeftResponse <- ixLeftResponse[-length(ixLeftResponse)]
        trialsPrecededByLChoice <- oneSessionData[ixLeftResponse+1, ]
        ## Find trials that were preceded by Right choice
        ixRightResponse <- which(oneSessionData$Response==2)
        if (ixRightResponse[length(ixRightResponse)] == nrow(oneSessionData)) ixRightResponse <- ixRightResponse[-length(ixRightResponse)]
        trialsPrecededByRChoice <- oneSessionData[ixRightResponse+1, ]
        trialsPrecededByLChoice$PrecedingChoice <- 1
        trialsPrecededByRChoice$PrecedingChoice <- 2
        rbind(trialsPrecededByLChoice, trialsPrecededByRChoice)
        ## Bind together trials marked by L and R preceding choices
        markedTrials <- rbind(markedTrials, trialsPrecededByLChoice, trialsPrecededByRChoice)
      }
      #}
    }
    markedTrials$PrecedingChoice <- as.factor(markedTrials$PrecedingChoice)
    return(markedTrials)
  }


  PercentRightChoices <- function(markedTrials)
  {
    ## Plot proportion correct of responding right to stimuli presented either to left or to right
    pcRight <- droplevels(markedTrials)
    ## Label left responses to right gratings with negative contrast. Right responses to right gratings will remain with positive sign
    pcRight[pcRight$VisualField == 1,]$Contrast <- pcRight[pcRight$VisualField == 1,]$Contrast * -1
    ## Summary of responses to gratings presented in the right
    pcRightSummary <- ddply(pcRight, .(SubjectID, SessionID, Condition, PrecedingChoice, VisualField, Contrast), summarise,
                            nRightResp = sum(Response == 2),
                            nLeftResp = sum(Response == 1),
                            nRStim = sum(VisualField == 2),
                            nLStim = sum(VisualField == 1))
    pcRightSummary <- ddply(pcRightSummary, .(SubjectID, SessionID, Condition, PrecedingChoice, VisualField, Contrast), summarise,
                            pRightCorrect = nRightResp / (nRStim + nLStim),
                            nYesR=nRightResp,
                            nNoR=nLeftResp)
    ## Convert contrast into %
    pcRightSummary$Contrast <- pcRightSummary$Contrast * 100
    return(pcRightSummary)

  }

  # Mark trials preceded by L and R choices from the data
  markedTrials <- MarkTrialsByChoice(inputData)
  markedTrialsPrcRChoice <- PercentRightChoices(markedTrials)

  # Mark trials preceded by L/R choices using simulated data to estimate confidence intervals
  if (!missing(confIntData)) {
    confIntMarkedTrials <- MarkTrialsByChoice(confIntData)
    # Compute proportion right choices using simulated data to estimate confidence intervals
    confIntPrcRChoice <- PercentRightChoices(confIntMarkedTrials)
    # Compute confidence intervals
    confInt <- ddply(confIntPrcRChoice, .
                     (SubjectID, Condition, PrecedingChoice, Contrast),
                     function(x) {c(quantile(x$pRightCorrect, c(0.16, 0.84), names=FALSE), # 68% confidence intervals.
                                    median(x$pRightCorrect),
                                    mean(x$pRightCorrect),
                                    mean(x$pRightCorrect)-sqrt(var(x$pRightCorrect)/length(x$pRightCorrect)),
                                    mean(x$pRightCorrect)+sqrt(var(x$pRightCorrect)/length(x$pRightCorrect)),
                                    mean_cl_boot(x$pRightCorrect, conf.int=0.95)$ymin,
                                    mean_cl_boot(x$pRightCorrect, conf.int=0.95)$ymax)
                     })
    colnames(confInt)[which(names(confInt) %in% c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"))] <- c("ciMin", "ciMax", "Median", "Mean", "seMin", "seMax", "ciBMin", "ciBMax")
    # This is a dummy variable needed for ggplot to be present, but ggplot will
    # only use ciMin and ciMax to plot the error bars
    confInt$pRightCorrect <- 0.5
  }

  g <- ggplot(data = markedTrialsPrcRChoice, aes(x = Contrast, y = pRightCorrect, color=PrecedingChoice))
  #g <- g + geom_line(aes(group=PrecedingChoice), stat="summary", fun.y="mean", size=0.3)
  g <- g + theme_publish2() + geom_rangeframe(color="grey30")
  if (showMidPoints) {
    g <- g + geom_vline(xintercept = 0.0, size = 0.2, colour = "grey30", linetype = "dashed")
    g <- g + geom_hline(yintercept = 0.5, size = 0.2, colour = "grey30", linetype = "dashed")
  }
  #g <- g + geom_errorbar(data=confInt, aes(ymin=ciMin, ymax=ciMax), width=0.1, alpha=0.3)
  if (!missing(confIntData)) g <- g + geom_smooth(data=confInt, aes(ymin=ciMin, ymax=ciMax, fill=PrecedingChoice), stat="identity", linetype=0, alpha=0.2)
  #g <- g + geom_smooth(data=confInt, aes(ymin=ciBMin, ymax=ciBMax, fill=PrecedingChoice), stat="identity", linetype=0, alpha=0.2)
  g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())
  #g <- g + theme(axis.line = element_line(colour = "#a9a9a9", size = 0.3))
  g <- g + theme(axis.ticks.x = element_line(colour = "#a9a9a9", size = 0.3))
  g <- g + stat_summary(aes(group=PrecedingChoice), fun.data = "mean_cl_boot", B=500, conf.int = 0.68, geom = "errorbar", size = 0.2, width = 0.0)
  if (!missing(confIntData)) g <- g + geom_line(data=confInt, aes(x=Contrast, y=Mean), alpha=0.5, size=0.2)
  #g <- g + geom_line(data=confInt, aes(x=Contrast, y=Median), alpha=0.4)
  g <- g + geom_point(aes(group=PrecedingChoice), stat="summary", fun.y="mean", size=geomPointSize)
  g <- g + facet_wrap(~SubjectID, scale="free_x")
  g <- g + scale_colour_manual(values=c(prevFailColor, prevSuccessColor))
  g <- g + ylab("Proportion rightward choices")
  g <- g + xlab("Contrast (%)")
  g <- g + theme(legend.position="none")

  if (plotFigure) {
    dev.new(width=figureWidth, height=figureHeight)
    g
  } else {
    return(g)
  }
}


######################################################################################################
## Plot how biases change as a function of run. Do subjects change their biases particularly during
## bias induction conditions?
PlotWeightChangeByRun <- function (weightsRegularized		# Bias weights
) {
  #weights <- droplevels(subset(weightsRegularized, Parameter!="(Intercept)"))
  weights <- droplevels(subset(weightsRegularized, Parameter %in% c("PrevCorr1", "PrevFail1")))
  ggplot(data= weights, aes(x=SessionID, y=Weight, color=Parameter)) +
    theme_few() +
    geom_path(aes(group=Parameter), stat="smooth", method="lm", color="#102d95", size=1.5, alpha=1.0, lineend="round") +
    geom_line(aes(group=Parameter)) +
    #geom_point(size=3, color="white", fill="white", shape=21) +
    geom_point() +
    #scale_colour_manual(values=c(prevSuccessColor, prevFailColor), labels=c("Prev success", "Prev failure"))  +
    xlab("Run") +
    ylab("Bias") +
    #coord_cartesian(ylim=c(-3.2, 3.2)) +
    #theme(axis.text.x = element_text(angle=-45)) +
    facet_grid(Condition~SubjectID, labeller=ConditionLabels, scales="free")
  #facet_grid(SubjectID~Condition, labeller=ConditionLabels, scales="free") +
  #guides(color=FALSE)
}

##########################################################
##
PlotBiasesAndTheirDiffs <- function(weights, 						## Weight to be plotted. Need to be from same condition and only one type, such as fail or success
                                    weightToPlot = 'PrevFail1',		# Name of the weigh to be plotted
                                    conditionsToPlot = c(1,2),		# Pair of conditions to be plotted
                                    sorted=T,						# Subjects will be sorted by the weight of first condition
                                    #subjectsToPlot = 'all', 		# Default is to plot all subjects
                                    geomPointColor='black',			# Color of geom_point
                                    figureWidth=7.095745, 			# Width of the plot
                                    figureHeight= 3.989362			# Height of the plot
) {

  dataToPlot <- droplevels(subset(weights, (Condition %in% conditionsToPlot) & (Parameter %in% weightToPlot)))
  ## Remove two unnecessary columns (Vif and Regularized)
  dataToPlot$Vif <- NULL
  dataToPlot$Regularized <- NULL
  ## Get subject labels for both conditions
  sbjInFirstCondition <- levels(droplevels(dataToPlot[dataToPlot$Condition==conditionsToPlot[1],]$SubjectID))
  sbjInSecondCondition <- levels(droplevels(dataToPlot[dataToPlot$Condition==conditionsToPlot[2],]$SubjectID))
  ## Find common subjects for both conditions. Only they will be plotted
  subjectsToPlot <- intersect(sbjInFirstCondition, sbjInSecondCondition)
  ## Select those subjects for further data processing and plotting
  dataToPlot <- droplevels(subset(dataToPlot, SubjectID %in% subjectsToPlot))

  ## Show failure rate for each condition (success rate is 100%-failRate)
  failRateFirstCondition <- ComputeFailRate(glmData, conditionsToPlot[1], subjectsToPlot)
  print(paste("Failure rate for first condition: ", sprintf("%.0f", failRateFirstCondition*100), "%", sep=""))
  failRateSecondCondition <- ComputeFailRate(glmData, conditionsToPlot[2], subjectsToPlot)
  print(paste("Failure rate for second condition: ", sprintf("%.0f", failRateSecondCondition*100), "%", sep=""))

  ## Compute mean weights
  meanDataToPlot <- ddply(dataToPlot, .(SubjectID, Parameter, Condition), summarise, MeanWeight=mean(Weight))
  ## Compute difference by subtracting weights from first condition from second condition
  weightsDifference <- ddply(meanDataToPlot, .(SubjectID, Parameter), summarise, Difference=diff(MeanWeight))
  ## Rename last column name to have the same name as meanDataToPlot data frame. This is to join them later.
  names(weightsDifference)[3] <- "MeanWeight"
  ## Add condition column and assign it 1000 which will indicate that it is the data containing difference of conditions
  weightsDifference$Condition <- 1000
  ## Bind mean weights of two conditions and their differences into one data frame for plotting
  meanDataToPlot <- rbind(meanDataToPlot, weightsDifference)

  ## If requested, sort subjects by first condition weights
  if (sorted) {
    ## Sort SubjectID so that subjects with lowest weight are plotted first and
    ## subjects with largest weight are plotted the last. Sorting is done based on
    ## on first condition passed to the function
    firstConditionData <- subset(meanDataToPlot, Condition== conditionsToPlot[1])
    weightsOrder <- order(-firstConditionData$MeanWeight)
    subjectsOrdered <- firstConditionData$SubjectID[weightsOrder]
    ## Change order of subjects in the data that will be plotted.
    meanDataToPlot$SubjectID <- factor(meanDataToPlot$SubjectID, levels=subjectsOrdered)
  }
  #2
  p <- ggplot(meanDataToPlot, aes(x=SubjectID, y=MeanWeight))
  p <- p + theme_few() +
    theme(strip.background=element_rect(colour="white", fill="white")) +
    theme(panel.background=element_rect(colour="white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
    theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
    scale_y_continuous(limits=c(-2.5, 2.5)) +
    geom_hline(yintercept=0.0, size=0.5, colour="#a9a9a9", linetype = "solid") +
    geom_segment(aes(xend=SubjectID), yend=0, colour=geomPointColor, size=1) +
    #geom_errorbar(aes(ymin=MeanWeight-se, ymax=MeanWeight+se), width=0.01, alpha=0.2) +
    #geom_point(size=1.5, )
    geom_point(size=5, color=geomPointColor) +
    geom_point(size=2.7, color='white') +
    xlab("")  +
    ylab("") +
    theme(axis.line = element_line(colour = "#a9a9a9", size = 0.3),axis.line.y = element_blank()) +
    coord_flip() +
    facet_grid(~Condition, labeller=ConditionLabels)

  dev.new(width=figureWidth, height=figureHeight)
  print(p)

}

