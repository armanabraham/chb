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
  g <- g + geom_vline(yintercept = 0.0, size = 0.2, colour = "grey30", linetype = "dashed")
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
                               plotInColor=T		# Success and failure are marked by two different colors. Otherwise, black and white
){
  ## Add column with trial numbers
  #runData <- droplevels(ddply(runData, mutate, TrialID=1:length(Response)))
  runData <- cbind(runData, TrialID=1:length(runData$Response))

  ## Generate ggplot graphs as binary sparklines
  g <- ggplot(data=runData, aes(x=TrialID, y=Response)) + theme_few()
  g <- g + geom_linerange(subset=.(Response==1), aes(ymin=1.5, ymax=Response, color="Left"), size=0.3, ) +
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

