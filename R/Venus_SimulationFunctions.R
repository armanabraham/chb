#' Estimate slope and threshold with and without history biases
#'
#' Run simulation using either subjects' natural history biases or without
#' history biases (success, failure and L/R biases are set to 0).
#' It sets success, failure and L/R biases of each subject
#' to 0 when simulating no-bias responses and leaves them intact when simulating
#' using subject's natural biases. The "IsSubject" field indicates if simulates
#' has subject's natural biases (TRUE) or has biases removed (FALSE).
#'
#' Note: Simulations run using subjects' mean weights
#'
#' @param sbjWeights regularized subject weights
#' @param nTrialsPerContrast number of trials for each contrast
#' @param B number of simulations
#' @param setLRBiasToZero removes L/R bias (intercept) from all simulations, including from subject's weights
#'
#' @examples
#' SimulateNoBiasVsSubjectBias(sbjWeights, nTrialsPerContrast = 100, B = 30)
#'
#' @export
SimulateNoBiasVsSubjectBias <- function(sbjWeights,               # Regularized weights
                                       nTrialsPerContrast = 100,  # Number of trials per contrast
                                       B = 30,                    # Number of simulations
                                       setLRBiasToZero=TRUE)      # Sets L/R bias to zero for all simulations, including for subject weights
{
  # Set L/R bias (intercept) to 0, if requested
  if (setLRBiasToZero) sbjWeights[sbjWeights$Parameter=='(Intercept)',]$Weight  <- 0

  # For each subject, simulate responses with
  # and without biases
  thsAndSlopes <- data.frame()
  for (ixSubject in levels(droplevels(sbjWeights$SubjectID))) {
    oneSbjWeights <- droplevels(subset(sbjWeights, SubjectID == ixSubject))
    execTime <- system.time({
      for (ixSim in 1:B) {
        # Simulate trials using subject history bias
        # and another set of trials after removing all biases (including left/right bias)
        simTrialsNoBias <- SimulateSubjectResponses(oneSbjWeights,
                                                    failWeightsToSim=0,
                                                    successWeightsToSim=0,
                                                    nTrialsPerContrast=nTrialsPerContrast, # Number of trials per contrast
                                                    setLRBiasToZero=setLRBiasToZero)
        # Compute 75% thresholds for biased and unbiased simulations
        oneSimThAndSlope <- ThAndSlopeForSimData(simTrialsNoBias)
        thsAndSlopes <- rbind(thsAndSlopes, data.frame(oneSimThAndSlope,
                                                       SimulationID=rep(ixSim, nrow(oneSimThAndSlope))))
      }
    }) # system.time
    print(paste("Subject #", ixSubject, "simulation #", ixSim, "lasted for", round(as.numeric(execTime[3]), 2), "sec"))
  }
  ## Add subject degree column
  sbjDemographs <- SubjectDemographics()
  thsAndSlopes <- merge(thsAndSlopes, sbjDemographs, by='SubjectID')
  colnames(thsAndSlopes)[which(names(thsAndSlopes) == "Education")] <- "Degree"

  return(thsAndSlopes)
}



#' Simulate subject responses for different set of failure weights
#'
#' Also controls setting L/R bias and success bias to 0
#' This function can be generalized because in current form
#' it only assumes that only failure biases change but success
#' biases remain either constant or absent (same for L/R bias)
#' @export
SimulateSubjectResponses <- function(subjectWeights,           # Weights of the model for each subject. Only regularized weights will be selected.
                                     failWeightsToSim,         # Failure weights for which simulations will be computed
                                     successWeightsToSim,      # Success weights to simulate
                                     nTrialsPerContrast=1000,                   # Number of trials one contrast intensity should be simulates (min should be 4 and preferably divisible by 4)
                                     setLRBiasToZero=TRUE,     # Set Left/Right subject bias to 0. Alternatively, subject's L/R bias will used
                                     conditionsToSimulate)     # Experimental conditions that will be used for simulation
{

  if (missing(failWeightsToSim)) stop("Please specify failWeightsToSim")
  if (missing(successWeightsToSim)) stop("Please specify successWeightsToSim")

  # Select conditions
  if (!missing(conditionsToSimulate)) subjectWeights <- subset(subjectWeights, Condition %in% conditionsToSimulate)
  # Ensure only regularized weights are selected
  regularizedWeights <- droplevels(subset(subjectWeights, Regularized=="yes"))

  # Compute mean subject weights for each condition
  meanSbjWeights <- ddply(regularizedWeights, .(SubjectID, Condition, Parameter), summarise,
                          MeanWeight=mean(Weight),
                          std=sd(Weight, na.rm=TRUE),
                          n=sum(!is.na(Weight)),
                          se=std/sqrt(n))
  meanSbjWeights <- droplevels(meanSbjWeights)
  ## Report number of simulations
  if (nTrialsPerContrast<4) nTrialsPerContrast <- 4
  print(paste('Requested', nTrialsPerContrast, 'trials, adjusted to ', as.integer(nTrialsPerContrast/4)*4))

  # Set L/R bias (intercept) to 0, if requested
  if (setLRBiasToZero) meanSbjWeights[meanSbjWeights$Parameter=="(Intercept)",]$MeanWeight <- 0

  simData <- data.frame() # data frame to store simulation results
  ## Loop through each "PrevFail" weights
  # Main loop where simulation runs
  for (ixSubject in levels(droplevels(meanSbjWeights$SubjectID))) {
    oneSbjConditions <- unique(meanSbjWeights$Condition[meanSbjWeights$SubjectID==ixSubject])
    for (ixCondition in oneSbjConditions) {
      subjectModel <- subset(meanSbjWeights, SubjectID==ixSubject & Condition==ixCondition)
      #browser()
      historyWeights <- expand.grid(failWeight=failWeightsToSim, successWeight=successWeightsToSim, IsSubject=FALSE)
      sbjFailWeight <- with(subjectModel, MeanWeight[Parameter=='PrevFail1'])
      sbjSuccessWeight <- with(subjectModel, MeanWeight[Parameter=='PrevCorr1'])

      historyWeights <- rbind(historyWeights,
                              data.frame(failWeight=sbjFailWeight, successWeight=sbjSuccessWeight, IsSubject=TRUE))
      for (ixWeight in 1:nrow(historyWeights)) {
        #browser()
        newModel <- subjectModel
        newModel$MeanWeight[newModel$Parameter=="PrevFail1"] <- historyWeights$failWeight[ixWeight]
        newModel$MeanWeight[newModel$Parameter=="PrevCorr1"] <- historyWeights$successWeight[ixWeight]
        ## Get contrast levels
        ixContrasts <- grep("c0", newModel$Parameter)
        contrastLevelsAsChar <- as.character(newModel$Parameter[ixContrasts])
        contrastLevelsAsNum <- sort(as.numeric(substring(contrastLevelsAsChar, 2)))
        ## Add new column with contrast levels
        #browser()
        simulationTrials <- do.call(rbind, lapply(contrastLevelsAsNum, PrepareTrials, subjectID=ixSubject, Condition=ixCondition))
        simulationTrials <- simulationTrials[rep(1:nrow(simulationTrials), nTrialsPerContrast/4), ] # nTrialsPerContrast/4 because left and right sides are already in simulationTrials as well as left and right drifts
        ## Randomly shuffle trials
        simulationTrials <- simulationTrials[sample(nrow(simulationTrials)), ]
        simTrialsForGlm <- BuilDataForGLM(simulationTrials, nHistoryBack=nBack, successColName=successColName, failColName=failColName)
        simTrialsForGlm <- droplevels(simTrialsForGlm)
        ## Add "Intercept" column
        simTrialsForGlm$"(Intercept)" <- 1
        ## Calculate the response to the first trial cause it doesn't have history term
        logOdds <- sum(simTrialsForGlm[1, as.character(newModel$Parameter)] * newModel$MeanWeight)
        p <- 1 / (1 + 1/exp(logOdds))
        simTrialsForGlm[1,]$Response <- c(1,2)[rbinom(1, size=1, prob=p)+1]
        simTrialsForGlm[1,]$CorrIncorr <- ifelse(simTrialsForGlm[1,]$VisualField==simTrialsForGlm[1,]$Response, 1, 0)
        simTrialsForGlm[1,]$y <- ifelse(simTrialsForGlm[1,]$Response==1, -1, 1)
        ## For each trial, get model responses
        for (ixTrial in 2:nrow(simTrialsForGlm)) {
          # If previous trial was failure
          if (simTrialsForGlm[ixTrial-1,]$CorrIncorr==0) {
            simTrialsForGlm[ixTrial,]$PrevFail1 <- ifelse(simTrialsForGlm[ixTrial-1,]$Response==1, -1, 1)
          } else {
            simTrialsForGlm[ixTrial,]$PrevCorr1 <- ifelse(simTrialsForGlm[ixTrial-1,]$Response==1, -1, 1)
          }
          currentTrial <- simTrialsForGlm[ixTrial, as.character(newModel$Parameter)]
          #print(currentTrial)
          logOdds <- sum(currentTrial * newModel$MeanWeight)
          ## Back transform logOdds into probability
          p <- 1 / (1 + 1/exp(logOdds))
          #browser()
          ## Toss a coin with that probability to identify left or right response to the stimulus
          simTrialsForGlm[ixTrial,]$Response <- c(1,2)[rbinom(1, size=1, prob=p)+1]
          #simTrialsForGlm[ixTrial,]$Response <- c(1,2)[rbinom(1, size=1, prob=p)+1]
          simTrialsForGlm[ixTrial,]$CorrIncorr <- ifelse(simTrialsForGlm[ixTrial,]$VisualField==simTrialsForGlm[ixTrial,]$Response, 1, 0)
          simTrialsForGlm[ixTrial,]$y <- ifelse(simTrialsForGlm[ixTrial,]$Response==1, -1, 1)
        }
        simData <- rbind.fill(simData, cbind(simTrialsForGlm,
                                             FailWeight=historyWeights$failWeight[ixWeight],
                                             SuccessWeight=historyWeights$successWeight[ixWeight],
                                             IsSubject=historyWeights$IsSubject[ixWeight]))
        print(paste(ixSubject, ', failWeight=',
                    as.character(historyWeights$failWeight[ixWeight]),
                    ', successWeight=',
                    historyWeights$successWeight[ixWeight],
                    ', Condition=', ixCondition, sep=''))
        #browser()
      }
    }
  }
  return(simData)
}


#' Simulates responses of one subject
#'
#' Simulate one subject using contrast weights, L/R bias and history weights
#' The simulation is computed ONLY ONCE using certain number of trials specified by
#' trialsPerContrast. You can call this function multiple times to estimates
#' confidence intervals of slope or threshold, for example
#'
#' @param subjectWeights contrast and history bias weights including L/R bias
#' @param SubjectID subject code, such as 's007'
#' @param Condition a number that specifiec condition of simulation
#' @param trialsPerContrast number of trials that each contrast intensity will be simulated
#' @param nBack history depth. At the moment this is not functional???
#'
#' @examples
#' SimulateOneSubject(oneSbjWeights, nTrialsPerContrast=50)
#' @export
SimulateOneSubject <- function(subjectWeights,           # This should include Intercept, history bias weights and contrast weights
                               SubjectID="temp",         # Subject ID
                               Condition= -777,          # Condition
                               trialsPerContrast=100,    # Number of trials per contrast
                               nBack=nBack)              # History depth
{
  sbjModel <- subjectWeights
  ## Get contrast levels
  ixContrasts <- grep("c0", sbjModel$Parameter)
  contrastLevelsAsChar <- as.character(sbjModel$Parameter[ixContrasts])
  contrastLevelsAsNum <- sort(as.numeric(substring(contrastLevelsAsChar, 2)))
  ## Add new column with contrast levels
  simulationTrials <- do.call(rbind, lapply(contrastLevelsAsNum, PrepareTrials, subjectID=SubjectID, Condition=Condition))
  simulationTrials <- simulationTrials[rep(1:nrow(simulationTrials), trialsPerContrast/4), ] # divided by 4 because left and right sides are already in simulationTrials as well as left and right drifts
  ## Randomly shuffle trials
  simulationTrials <- simulationTrials[sample(nrow(simulationTrials)), ]
  simTrialsForGlm <- BuilDataForGLM(simulationTrials, nBack, successColName=successColName, failColName=failColName)
  simTrialsForGlm <- droplevels(simTrialsForGlm)
  ## Add "Intercept" column
  simTrialsForGlm$"(Intercept)" <- 1
  ## Calculate the response to the first trial cause it doesn't have history term
  logOdds <- sum(simTrialsForGlm[1, as.character(sbjModel$Parameter)] * sbjModel$Weight)
  p <- 1 / (1 + 1/exp(logOdds))
  simTrialsForGlm[1,]$Response <- c(1,2)[rbinom(1, size=1, prob=p)+1]
  simTrialsForGlm[1,]$CorrIncorr <- ifelse(simTrialsForGlm[1,]$VisualField==simTrialsForGlm[1,]$Response, 1, 0)
  simTrialsForGlm[1,]$y <- ifelse(simTrialsForGlm[1,]$Response==1, -1, 1)
  ## For each trial, get model responses
  for (ixTrial in 2:nrow(simTrialsForGlm)) {
    # If previous trial was failure
    if (simTrialsForGlm[ixTrial-1,]$CorrIncorr==0) {
      simTrialsForGlm[ixTrial,]$PrevFail1 <- simTrialsForGlm[ixTrial-1,]$y
    } else {
      simTrialsForGlm[ixTrial,]$PrevCorr1 <- simTrialsForGlm[ixTrial-1,]$y
    }
    currentTrial <- simTrialsForGlm[ixTrial, as.character(sbjModel$Parameter)]
    logOdds <- sum(currentTrial * sbjModel$Weight)
    ## Back transform logOdds into probability
    p <- 1 / (1 + 1/exp(logOdds))
    ## Toss a coin with that probability to identify left or right response to the stimulus
    simTrialsForGlm[ixTrial,]$Response <- c(1,2)[rbinom(1, size=1, prob=p)+1]
    simTrialsForGlm[ixTrial,]$CorrIncorr <- ifelse(simTrialsForGlm[ixTrial,]$VisualField==simTrialsForGlm[ixTrial,]$Response, 1, 0)
    simTrialsForGlm[ixTrial,]$y <- ifelse(simTrialsForGlm[ixTrial,]$Response==1, -1, 1)
  }
  simTrialsForGlm$Condition <- as.factor(simTrialsForGlm$Condition)
  return(simTrialsForGlm)
  # Data frame to store the results
  #simulatedData <- rbind.fill(simulatedData, cbind(simTrialsForGlm, PrevFailWeight=prevFailValues$Weight[ixWeight], IsSubject=prevFailValues$IsSubject[ixWeight]))
  #print(paste(ixSubject, ', prevFailWeight=', as.character(prevFailValues$Weight[ixWeight]), ' Condition=', ixCondition, sep=''))
}


#' Simulate responses using given model weights
#'
#' Compared to SimulateSubjectResponses, this function is simpler
#' because it only deals with given set of weights and does not
#' allow running different model parameters. This function
#' can actually be called from within SimulateSubjectResponses
#' but for some reason it is not done. A TODO for future.
#'
#' @param subjectWeights
#' @param nSimulations
#' @param trialsPerContrast
#' @param Conditions
#'
#' @export
SimulateResponses <- function(subjectWeights,           # Model weights by subject and condition.
                              nSimulations=10,
                              trialsPerContrast=50,
                              Conditions)              # Conditions to simulate
{

  if (!missing(Conditions)) subjectWeights <- droplevels(subset(subjectWeights, Condition %in% Conditions))
  if (nrow(subjectWeights)==0) {
    cat("(SimulateResponses) subjectWeights in empty")
    return()
  }
  subjectWeights <- droplevels(subjectWeights)
  # Simple check to ensure non-regularized weights are excluded
  if ("Regularized" %in% colnames(subjectWeights))
    if (length(unique(subjectWeights$Regularized))>1) {
      print("Might have found non-regularized weights. They will be excluded.")
      subjectWeights <- droplevels(subset(subjectWeights, Regularized=="yes"))
    }
  # Progress bar indicator
  nSubjects <- length(levels(subjectWeights$SubjectID))
  if (nSubjects==1) nSubjects = 2
  progBar <- txtProgressBar(style=3, min=1, max=nSubjects, label="simulating")
  # Data storage
  simData <- data.frame()
  # Compute simulations
  for (ixSubject in levels(subjectWeights$SubjectID)) {
    sbjConditions <- unique(subjectWeights$Condition[subjectWeights$SubjectID==ixSubject])
    for (ixCondition in sbjConditions) {
      subjectModel <- subset(subjectWeights, SubjectID==ixSubject & Condition==ixCondition)
      for (ixSimulation in 1:nSimulations) {
        tmpData <- SimulateOneSubject(subjectModel, SubjectID=ixSubject, Condition=ixCondition, trialsPerContrast=trialsPerContrast)
        tmpData$SessionID <- ixSimulation
        simData <- rbind.fill(simData, tmpData)
        #print(paste(ixSubject, ', Condition=', ixCondition, sep=''))
      }
    }
    setTxtProgressBar(progBar, which(levels(subjectWeights$SubjectID)==ixSubject))
  }
  close(progBar)
  simData$SessionID <- as.factor(simData$SessionID)
  return(simData)
}

#' Compute threshold and slope from simulated trials
#'
#' There is another function called ThAndSlope which works
#' with subjects' data rather than simulated data, because simulated
#' data of glmData type includes additional column called IsSubject
#'
#' @export
ThAndSlopeForSimData <- function(simTrials) {
  pcRight <- simTrials
  ## Label left responses to right gratings with negative contrast. Right responses to right gratings will remain with positive sign
  pcRight[pcRight$VisualField == 1,]$Contrast <- pcRight[pcRight$VisualField == 1,]$Contrast * -1
  ## Summary of responses to gratings presented in the right
  pcRightSummary <- ddply(pcRight, .(SubjectID, SessionID, VisualField, Contrast, IsSubject, FailWeight, SuccessWeight), summarise,
                          nRightResp = sum(Response == 2),
                          nLeftResp = sum(Response == 1),
                          nRStim = sum(VisualField == 2),
                          nLStim = sum(VisualField == 1))
  pcRightSummary <- ddply(pcRightSummary, .(SubjectID, SessionID, VisualField, Contrast, IsSubject, FailWeight, SuccessWeight), summarise,
                          pRightCorrect = nRightResp / (nRStim + nLStim),
                          nYesR=nRightResp,
                          nNoR=nLeftResp)
  ## Convert contrast into %
  pcRightSummary$Contrast <- pcRightSummary$Contrast * 100

  models <- dlply(pcRightSummary, c("SubjectID", "IsSubject", "FailWeight", "SuccessWeight"), .fun=chb:::FitProbit)
  predvals <- ldply(models, .fun=chb:::PredictvalsProbit, xvar="Contrast", yvar="pRightCorrect", type="response")
  #browser()
  ## Further summarise results cause we need slopes only
  predvals1 <- ddply(predvals, .(SubjectID, IsSubject, FailWeight, SuccessWeight), summarise, slope=mean(slope), th50=mean(Th50), th75=mean(Th75))
  ## Get percent change of slope relative to slope when error weight is 0
  #predvals1 <- ddply(predvals1, .(SubjectID), transform, SlopeChange = slope[which(PrevFailWeight==0)] / slope)
  return(predvals1)
}


#' Prepare a dataframe with a given contrast intensity
#'
#' Makes all necessary columns that typically correspond to a trial.
#' However, this function returns four trials, by splitting
#' it into 2 VisualFields (left/right) and two drifts (left/right)
#' It is then possible to replicate rows of this dataframe
#' typically used to simulations
PrepareTrials <- function(contrast, subjectID="ERR", Condition=-1) {
  numReps <- 4
  trials <- data.frame(SubjectID=subjectID, SessionID=as.factor(1), Contrast=contrast, VisualField=as.numeric(1), Drift=as.numeric(1), Response=NaN, RT=NaN, CorrIncorr=NaN,  y=NaN, Condition=as.integer(1))
  trials <- trials[rep(1:nrow(trials), numReps),]
  trials$Drift[c(2, 4)] <- 2
  trials$VisualField[3:4] <- 2
  trials$Condition <- Condition
  return(trials)
}


#' Compute decline in sensitivity as a function of choice history biases
#'
#' Takes subject model weights, computes mean subject weights and then take the average
#' across all those subject means which results in one set of weights, which is average across
#' all subjects. Basically, it is like running a simulation with one "average" subject, who is the
#' average of all subjects. This approach saves time because this simulation take a while to run.
#'
#' This function stores thresholds and slopes of all simulations. To plot the results use
#' PlotDeclineMatrix function.
#'
#' @param regWeights subject weights
#' @param failWeightsToSim a matrix of failure weights
#' @param successWeightsToSim and matrix of success weights
#' @param nTrialPerContrast number of trials per contrast intensity
#' @param B number of simulations. A good number to start with is 30. Larger is better. T
#' @param setLRBiasToZero either L/R bias of the model as it is or set subjects' L/R bias to zero
#'
#' @return Returns a data frame containing threshold and slope values of each simulation done with all combinations of success and fail weights
#'
#' @examples
#' # This will run a massive simulation of subject responses
#' by using 13x13=169 matrix of fail and success weights and each
#' matrix (i.e., 169 simulations) will be run 50 times.
#' SensitivtyDeclineMatrix(regWeights, failWeightsToSim = seq(-2,2, length.out = 13),
#'                         successWeightsToSim = seq(-2,2, length.out = 13), nTrialsPerContrast=50, B=50)
#' # A quicker simulation
#' SensitivtyDeclineMatrix(regWeights, failWeightsToSim = c(-2, -1,0,1, 2),
#'                        successWeightsToSim = c(-2, -1,0,1, 2), nTrialsPerContrast=30, B=30)
#'
#' @export
SensitivtyDeclineMatrix <- function(regWeights,
                                    failWeightsToSim,
                                    successWeightsToSim,
                                    nTrialPerContrast,
                                    B=30,
                                    setLRBiasToZero = TRUE)
{
  regWeights <- droplevels(subset(regWeights, Regularized=='yes'))
  meanSbjWeights <- ddply(regWeights, .(SubjectID, Condition, Parameter), summarise,
                          Weight=mean(Weight))
  # Mean weight of all subjects. This will be used for simulation
  allSbjMeanWeights <- ddply(meanSbjWeights, .(Parameter), summarise, Weight=mean(Weight))
  # Adding some dummy columns which are needed by the SimulateSubjectResponses function
  allSbjMeanWeights$SubjectID <- as.factor('sAll')
  allSbjMeanWeights$Condition <- 1
  allSbjMeanWeights$Regularized <- 'yes'

  # To store all simulations
  allThsAndSlopes <- data.frame()
  # Simulation loop
  for (ixSim in 1:B) {
    print(paste('--------- Simulation #', ixSim))
    simDat <- SimulateSubjectResponses(allSbjMeanWeights,
                                       failWeightsToSim = failWeightsToSim,
                                       successWeightsToSim = failWeightsToSim,
                                       nTrialsPerContrast = nTrialPerContrast,
                                       setLRBiasToZero = setLRBiasToZero)
    thAndSlope <- ThAndSlopeForSimData(simDat)
    thsAndSlopeNoSbj <- droplevels(subset(thAndSlope, IsSubject==FALSE))
    #browser()
    allThsAndSlopes <- rbind(allThsAndSlopes, data.frame(thsAndSlopeNoSbj,
                                                         SimulationID=rep(ixSim, nrow(thsAndSlopeNoSbj))))
  }
  return(allThsAndSlopes)
}


