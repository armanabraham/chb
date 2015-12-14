## Functions for analysing choice history biases
#' @import ggplot2 ggthemes Hmisc aapack
#library(car)
#library(psyphy) ## To use fitting of curves with asymptotes
#library(RColorBrewer)
#library(plyr)
#library(grid)

# Themes for publication quality figures
#source('~/Desktop/Dropbox/R/Lib/PublishingThemes.R')
#source('~/Desktop/Dropbox/R/Lib/StandardErrorsByWinstonChan.R')

## Defining variables to use throughout
## When alpha is 0, run "ridge" regularization. When alpha is set to 1, run "lasso" regularization
alpha=0
# How many trials go back in history for prev faiure and success parameters
nBack <- 1
## Names of history columns
successColName <- "PrevCorr"
failColName <- "PrevFail"

## Colors to be used for plotting history weights
prevFailColor <- '#d7191c'
prevSuccessColor <- '#2b83ba'

## Old color values
# prevFailColor <- 'b62813'
# prevSuccessColor <- 'a7b112'

# Set of regularization parameters (lambdas) to test
lambdas <- c(exp(-8), exp(-7), exp(-6), exp(-5), exp(-4), exp(-3), 0.05, 0.1, 0.15, 0.25, 0.3, 0.5, 0.8, 1, 3, 7, 11, 15, 20)

#' Subject demographics such as age and education
#'
#' @export
SubjectDemographics <- function(getAge=TRUE,
                                getGender=TRUE,
                                getEducation=TRUE)
{
  subjectID <- c("s001", "s002", "s003", "s004", "s005", "s006", "s007", "s008",
                 "s009", "s010", "s011", "s012", "s013", "s014", "s015", "s016",
                 "s017", "s018", "s019", "s020", "s021", "s022", "s023", "s024",
                 "s025", "s026", "s027",
                 "s030", "s031", "s032", "s033", "s034", "s035", "s036", "s037",
                 "s038", "s039", "s040", "s041")
  demographics <- data.frame(SubjectID=subjectID)

  if (getAge) {
    age <- c(38, 36, 32, 21, 33, 31, 28, 29, 37, 28, 35, 35, 32, 22, 25,
             NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
             21, 36, 25, 19, 32, 24, 22, 21,
             20, 19, 27, 27)
    demographics <- cbind(demographics, Age=age)
  }
  # Subject degrees
  # s001 - Arman Abrahamyan					PhD
  # s002 - Georgios Keliris					PhD
  # s003 - Andrew Meso						PhD
  # s004 - Mackenzie Dolginow				Undergrad
  # s005 - Kenji Haruhana					Undergrad
  # s006 - Ilias Rentzeperis				PhD
  # s007 - Steeve Laquitaine				PhD
  # s008 - Katharina Dobs					Grad
  # s009 - Tancy Kao					PhD
  # s010 - Li Feng Yi					Grad
  # s011 - Erin Munro					PhD
  # s012 - Yuka Okazaki					PhD
  # s013 - Siva Kaveri					Grad
  # s014 - Darren Seibert (Summer School student)		Grad
  # s015 - Stephen Bruggemann (Summer School student)	Undergrad
  if (getEducation) {
    demographics <- cbind(demographics, Education=rep(NA, length(subjectID)))
    phds <- c("s001", "s002", "s003", "s006", "s007", "s008",
              "s009", "s010", "s011", "s012", "s013", "s014",
              "s016", "s017", "s018", "s019", "s020", "s031",
              "s032", "s036", "s041")
    demographics$Education[demographics$SubjectID %in% phds] <- "PhD"

    noPhds <- c("s004", "s005", "s015", "s021",
                "s022", "s023", "s024", "s025",
                "s026", "s027",
                "s030", "s031", "s033", "s034",
                "s035", "s037", "s038", "s039", "s040")
    demographics$Education[demographics$SubjectID %in% noPhds] <- "No PhD"
    demographics$Education <- as.factor(demographics$Education)

  }

  if (getGender) {
    demographics <- cbind(demographics, Gender=rep(NA, length(subjectID)))
    females <- c("s004", "s008", "s009", "s011", "s012", "s030", "s033", "s035", "s036", "s037", "s038", "s039", "s041")
    demographics$Gender[demographics$SubjectID %in% females] <- "F"
    males <- c("s001", "s002", "s003", "s005", "s006", "s007", "s010", "s013", "s014", "s015", "s031", "s032", "s034", "s040")
    demographics$Gender[demographics$SubjectID %in% males] <- "M"
  }

  return(demographics)
}

#' Compute correlation between age and bias
#'
#' @export
CorrelateAgeAndBias <- function(sbjInfo, biases) {
  meanBias <- subset(biases, Parameter %in% c('PrevFail1', 'PrevCorr1') & Regularized=='yes')
  meanBias <- ddply(meanBias, .(SubjectID, Condition, Parameter), numcolwise(mean))
  meanBiasWide <- cast(meanBias, SubjectID+Condition~Parameter, value=.(Weight))
  biasAndAge <- merge(meanBiasWide, sbjInfo)
  biasAndAge <- ddply(biasAndAge, .(SubjectID, Condition), transform, AbsBias=abs(PrevFail1) + abs(PrevCorr1))
  print(biasAndAge)
  # Plot correlations
  g1 <- qplot(biasAndAge$Age, biasAndAge$PrevFail1) + theme_publish1() + geom_smooth(method='lm', se=FALSE, color='grey50') + xlab("Age") + ylab("Fail bias")
  g2 <- qplot(biasAndAge$Age, biasAndAge$PrevCorr1) + theme_publish1() + geom_smooth(method='lm', se=FALSE, color='grey50') + xlab("Age") + ylab("Success bias")
  g3 <- qplot(biasAndAge$Age, biasAndAge$AbsBias) + theme_publish1() + geom_smooth(method='lm', se=FALSE, color='grey50') + xlab("Age") + ylab("Absolute bias")
  library(gridExtra)
  grid.arrange(g3, g1, g2, ncol=1)
  # Select only subjects whose age is registered
  biasAndAge <- biasAndAge[!is.na(biasAndAge$Age), ]
  print(cor.test(biasAndAge$Age, biasAndAge$AbsBias))
  print(cor.test( biasAndAge$Age, biasAndAge$PrevCorr1))
  print(cor.test(biasAndAge$Age, biasAndAge$PrevFail1))
  #   qplot(biasAndAge$Age, biasAndAge$PrevFail1) +
  #     stat_smooth(stat='lm', se=FALSE)
}


#' Stub for glmnet that accepts a formula
#'
#' Adapted from a book (NEED TO CHECK THE REF)
#' @import glmnet
#' @export
Glmnet <- function(formula, data, subset, na.action, ...) {
	call <- match.call() # returns the function call
	mf <- match.call(expand.dots = FALSE) # the function call w/o ...
	args <- match(c("formula", "data", "subset", "na.action"), names(mf), 0) # which arguments are present?
	mf <- mf[c(1, args)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval.parent(mf) # create a model frame
	terms <- attr(mf, "terms") # terms object for the model
	y <- model.response(mf) # response variable
	X <- model.matrix(terms, mf, contrasts) # model matrix
	## Arman: Remove Intercept from the model because glmnet includes it automatically
	##X <- X[,-1]
	glmnet::glmnet(X, y, ...)
	#browser()
}

#' Stub for glmnet cross-validation that accepts a formula
#'
#' @import glmnet
#' @export
CV.glmnet <- function(formula, data, subset, na.action, ...) {
	call <- match.call() # returns the function call
	mf <- match.call(expand.dots = FALSE) # the function call w/o ...
	args <- match(c("formula", "data", "subset", "na.action"), names(mf), 0) # which arguments are present?
	mf <- mf[c(1, args)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval.parent(mf) # create a model frame
	terms <- attr(mf, "terms") # terms object for the model
	y <- model.response(mf) # response variable
	X <- model.matrix(terms, mf, contrasts) # model matrix
	## Arman: Remove Intercept from the model because glmnet includes it automatically
	##X <- X[,-1]
	glmnet::cv.glmnet(X, y, ...)
	#browser()
}


#' Return root means square error
# rmse <- function(y, h) {
	# return(sqrt(mean((y - h)^2)))
# }


#' Log-likelihood of choosing right
#'
#' It is negative log-likelihood
Likelihood <- function(y,     # Participant responses
                       pRight)  # Probability of responding right
{
  -2*sum(log(ifelse(y==1, pRight, 1 - pRight)))
}

#' Prepare data for GLM analysis
#'
#' Build all necessary columns and fill in data for the GLM analysis
#' Previous success and previous failure (history terms) columns are also constructed here
#' @export
BuilDataForGLM <- function(rawData, 		# Input data frame
                           nHistoryBack, 	# Depth of constructing previous success and previous failure variables. 1 means one trial back, and 2 means 2 trials back
                           nDummyParams,	# Number of dummy parameters. Currently unused
                           successColName, 	# Name of previous success columns for GLM
                           failColName)		# Name of previous failure columns for GLM
{
  ## Prepare data for GLM fitting
  ## Convert contrast intensities into columns for glm fitting
  contrastLevels <- as.character(sort(unique(rawData$Contrast)))
  glmData <- rawData[, c("SubjectID", "SessionID", "Contrast", "VisualField", "Drift", "Response", "RT", "CorrIncorr", "y", "Condition")]
  ## Add contrast levels as columns
  contrastColumns <- data.frame(t(rep(0, length(contrastLevels))))
  glmData <- cbind(glmData, contrastColumns)

  ## Column names better not be numeric but start with a character
  contrastColumnNames <- paste("c", contrastLevels, sep = "")
  columnIndex <- (ncol(glmData) - length(contrastLevels) + 1):ncol(glmData)
  names(glmData)[columnIndex] <- contrastColumnNames
  ## Store information about columns where contrasts are stored in a separate data frame. Will use in GLM
  contrastColumnInfo <- data.frame(columnIndex, contrastLevels, contrastColumnNames)
  ## Find each contrast level across all sessions and subjects and assign them either -1 (for left VF) or +1 (for right VF)
  for (ixContrast in 1:length(contrastLevels)) {
    ## indices of contrast being processed
    ixThisContrastLeftVF <- which((as.character(glmData$Contrast) == contrastLevels[ixContrast]) & (glmData$VisualField == 1))
    ixThisContrastRightVF <- which((as.character(glmData$Contrast) == contrastLevels[ixContrast]) & (glmData$VisualField == 2))
    ## Column index in the data frame that corresponds to this contrast level
    thisContrastLevelColumn <- ncol(glmData) - length(contrastLevels) + ixContrast
    glmData[ixThisContrastLeftVF, thisContrastLevelColumn] <- -1
    glmData[ixThisContrastRightVF, thisContrastLevelColumn] <- 1
  }

  finalGLMData <- data.frame()
  ## Create history columns referring to past success and failure trials
  if (nHistoryBack > 0) {
    glmData <- cbind(glmData, data.frame(t(rep(0, 2*nHistoryBack))))
    histColumnNames <- paste(c(successColName, failColName), sort(rep(seq(1:nHistoryBack), 2)), sep="")
    names(glmData)[((ncol(glmData)-2*nHistoryBack)+1):ncol(glmData)] <- histColumnNames

    ## for each subject and each session,
    for (ixSubject in levels(droplevels(glmData$SubjectID))) {
      oneSubjectData <- subset(glmData, glmData$SubjectID == ixSubject)
      subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))
      for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
        oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)

        ## Past success/failure trials are processed one by one based
        ## on the number of trials that's required to go back in history
        for (ixGoBack in 1:nHistoryBack) {
          ixErrorTrials <- which(oneSessionData$CorrIncorr[1:(nrow(oneSessionData)-ixGoBack)] == 0)
          ixWhenPrevFailed <- ixErrorTrials + ixGoBack

          nBackFailColName <- paste(failColName, as.character(ixGoBack), sep="")
          #oneSessionData$PrevFail1[ixWhenPrevFailed] <- ifelse(oneSessionData$Response[ixErrorTrials] == 1, -1, 1)
          #browser()
          oneSessionData[, grep(nBackFailColName, colnames(oneSessionData))][ixWhenPrevFailed] <- ifelse(oneSessionData$Response[ixErrorTrials] == 1, -1, 1)

          ixSuccessTrials <- which(oneSessionData$CorrIncorr[1:(nrow(oneSessionData)-ixGoBack)] == 1)
          ixWhenPrevSuccess <- ixSuccessTrials + ixGoBack

          nBackSuccessColName <- paste(successColName, as.character(ixGoBack), sep="")
          oneSessionData[, grep(nBackSuccessColName, colnames(oneSessionData))][ixWhenPrevSuccess] <- ifelse(oneSessionData$Response[ixSuccessTrials] == 1, -1, 1)

          ## Define WinStayLoseSwitch (WSLS) variable as follows
          ## It will be -1 if subject meant to follow WSLS and respond left
          ## and it will be 1 if subject meant to follow WSLS and respond right
          ## Find trials on which subject was success on left or fail on right.
          ## After those trials, subject has to chose Left, according to WSLS
          oneSessionData$WSLS <- 0
          ixWSLSChooseLeft <- which(with(oneSessionData,
                                         (Response[1:(nrow(oneSessionData)-1)]==1 & CorrIncorr[1:(nrow(oneSessionData)-1)]==1) |
                                           (Response[1:(nrow(oneSessionData)-1)]==2 & CorrIncorr[1:(nrow(oneSessionData)-1)]==0))  )
          ixWSLSChooseLeft <- ixWSLSChooseLeft + 1	## Adding 1 to choose trials affected by WSLS
          ## Mark trials affected by WLSL rule after which subject meant to choose left
          oneSessionData$WSLS[ixWSLSChooseLeft] <- -1

          ## Now let's do the same for right side
          ixWSLSChooseRight <- which(with(oneSessionData,
                                          (Response[1:(nrow(oneSessionData)-1)]==2 & CorrIncorr[1:(nrow(oneSessionData)-1)]==1) |
                                            (Response[1:(nrow(oneSessionData)-1)]==1 & CorrIncorr[1:(nrow(oneSessionData)-1)]==0))  )

          ixWSLSChooseRight <- ixWSLSChooseRight + 1
          oneSessionData$WSLS[ixWSLSChooseRight] <- 1
          #browser()
          ### This is the old way to define WSLS, in which 1 represented that subjected followed WSLS rule
          ### and -1 indicated that subjected followed the opposite of WSLS (WinSwitch,LoseStay)
          ### This way of coding the variable is not correct, because our response variable (y) is -1 or 1,
          ### meaning Left or Right, respectively. Thus, we need to code WSLS using this notation. The updated
          ### code above does just that.
          # oneSessionData$WSLS <- 0
          # ixSuccessStay <- which(with(oneSessionData, (CorrIncorr[1:(nrow(oneSessionData)-1)]==1) & (Response[1:(nrow(oneSessionData)-1)] == Response[2:nrow(oneSessionData)])))
          # ixSuccessStay <- ixSuccessStay + 1 ## +1 because those indexes refer to previous trial which was success/fail and caused switch stay. We want current trial that was a results of success/stay or fail/switch
          # ixFailSwitch <- which(with(oneSessionData, (CorrIncorr[1:(nrow(oneSessionData)-1)]==0) & (Response[1:(nrow(oneSessionData)-1)] != Response[2:nrow(oneSessionData)])) )
          # ixFailSwitch <- ixFailSwitch + 1 ## +1 because those indexes refer to previous trial which was success/fail and caused switch stay. We want current trial that was a results of success/stay or fail/switch
          # oneSessionData$WSLS[c(ixSuccessStay,ixFailSwitch)] <- 1
          # ixSuccessSwitch <- which(with(oneSessionData, (CorrIncorr[1:(nrow(oneSessionData)-1)]==1) & (Response[1:(nrow(oneSessionData)-1)] != Response[2:nrow(oneSessionData)])) )
          # ixSuccessSwitch <- ixSuccessSwitch + 1
          # ixFailStay <- which(with(oneSessionData, (CorrIncorr[1:(nrow(oneSessionData)-1)]==0) & (Response[1:(nrow(oneSessionData)-1)] == Response[2:nrow(oneSessionData)])) )
          # ixFailStay <- ixFailStay + 1
          # oneSessionData$WSLS[c(ixSuccessSwitch, ixFailStay)] <- -1
        }
        finalGLMData <- rbind(finalGLMData, oneSessionData)
      }
    }

  }
  #browser()


  return(finalGLMData)
}

#' Enrich rawData with helper columns
#'
#' Previous success and previous failure (history terms) columns are also constructed here
#' @export
PrepareRawData <- function(rawData 		# Input data frame
) {
  ## Here to order SessionID starting from 1
  ## However, store the original SessionIDs just in case
  rawData$OriginalSessionID <- rawData$SessionID
  ## Change SessionIDs to start from 1 for each subject
  for (ixSubject in levels(rawData$SubjectID)) {
    oneSubjectData <- subset(rawData, rawData$SubjectID==ixSubject)
    existingSessionIDs <- unique(oneSubjectData$SessionID)
    existingSessionIDs <- sort(existingSessionIDs)
    newSessionIDs <- 1:length(existingSessionIDs)
    for (ixSessionID in newSessionIDs) {
      rawData[((rawData$SubjectID==ixSubject) & (rawData$SessionID==existingSessionIDs[ixSessionID])),]$SessionID <- newSessionIDs[ixSessionID]
    }
  }
  ## Convert SessionID from number into factor
  rawData$SessionID <- as.factor(rawData$SessionID)
  ## Prepare GLM y response
  rawData$y <- rawData$Response
  rawData$y[rawData$y == 1] <- -1
  rawData$y[rawData$y == 2] <- 1

  ## Make a new column for correct/incorrect and fill it out with 0(incorrect) or 1(correct)
  rawData$CorrIncorr <- 1
  ## Indices of correct responses
  ixCorrect <- which(rawData$VisualField == rawData$Response)
  ## This is redundant but for clarity
  rawData$CorrIncorr[ixCorrect] = 1
  ixIncorrect <- which(rawData$VisualField != rawData$Response)
  rawData$CorrIncorr[ixIncorrect] = 0
  ## When responses are NaN or other than left or right assign NaN to CorrIncorr
  ixIlligalResponses <- which(is.nan(rawData$Response) | ((rawData$Response!=1) & (rawData$Response!=2)))
  rawData$CorrIncorr[ixIlligalResponses] = NaN
  rawData$y[ixIlligalResponses] = NaN

  return(rawData)
}

#' Fit model parameters using regularized regression
#'
#' Uses regulirized logistic regression (glmnet)
#' TODO: describe fully how training and validation done to find the best lambda
#' @export
RegularizedRegression <- function (glmData, 			          # subject responses
                                   fitWithLapseRate=FALSE,  # If TRUE, probabilities are adjusted (made smaller) using lapse rate
                                   lambdas, 			          #
                                   alpha=0,			            # When alpha is 0, run "ridge" regularization. When alpha is set to 1, run "lasso" regularization
                                   B,                       # number of simulation to estimate lambda that provides lowest log-likelihood
                                   fitHistoryBias=TRUE) 	  # Include history parameters. When FALSE, only contrast and general bias parameters are computed
{
  maximumSessions <- max(as.numeric(levels(glmData$SessionID)))
  nParticipants <- length(levels(glmData$SubjectID))
  ## Data frame to store each fit and related info
  fittedParams <- data.frame(ixSubject = character(), ixSession = character, TermsLogistic = character(), Terms = character(), CoefLogistic = double(), CoefGaussian = double())
  modelsPredictions <- data.frame()
  ## Create an array to store GLMs. We will use them later on for simulations
  maximumSessions <- max(as.numeric(levels(glmData$SessionID)))
  nParticipants <- length(levels(glmData$SubjectID))
  glmsFullModel <- array(list(), c(nParticipants, maximumSessions)) # In this two-dimensional list we store all GLM for each session
  ## Fit GLM to each session data for each participant
  performance <- data.frame()

  #glmnetFits <- array(list(), c(nParticipants, maximumSessions))  # In this two-dimensional list we store all GLM for each session
  for (ixSubject in levels(droplevels(glmData$SubjectID))) {
    oneSubjectData <- subset(glmData, glmData$SubjectID == ixSubject)
    subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))

    for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
      oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)
      ifelse(!is.null(oneSessionData$Condition), numCondition <- unique(oneSessionData$Condition), numCondition <- -1)

      # Adjustment for lapse rate, if necessary
      if (fitWithLapseRate) {
        lapseRate <- ThSlopeAndLapse(oneSessionData)$LapseRate
      } else {
        lapseRate <- 0
      }

      n <- nrow(oneSessionData)
      for (ixSimulation in 1:B) {
        indices <- sort(sample(1:n, round(0.8 * n)))
        training.data <- oneSessionData[indices, ]
        test.data <- oneSessionData[-indices, ]

        contrasts <- levels(as.factor(oneSessionData$Contrast))
        contrasts <- paste("c", contrasts, sep = "")
        ## Find out number of prev success and failure parameters
        nHistoryBack <- length(grep(chb:::successColName, colnames(glmData)))
        if (nHistoryBack > 0 & fitHistoryBias) {
          histColumnNames <- paste(c(chb:::successColName, chb:::failColName), sort(rep(seq(1:nHistoryBack), 2)), sep="")
          coefs <- paste(c(histColumnNames, contrasts), sep = "")
        } else {
          coefs <- contrasts
        }

        formulaLogistic <- as.formula(paste("as.factor(y) ~", paste(coefs, collapse = "+")))
        glmnetCurrentFit <- Glmnet(formula=formulaLogistic, family='binomial', data=training.data, lambda=lambdas, alpha=alpha)

        ## Prepare X predictor values to run on test dataset
        #browser()
        newx <- cbind(`(Intercept)` = rep(1, nrow(test.data)), data.matrix(test.data[, coefs]))
        likelihoods <- sapply(lambdas, function(lambda) {
          pRight <- predict(glmnetCurrentFit, newx, s = lambda, type = "response")
          responses <- test.data$y
          return(Likelihood(responses, pRight, lapseRate))
        })
        #browser()
        nLambdas <- length(lambdas)
        performance <- rbind(performance, data.frame(SubjectID=rep(ixSubject, nLambdas),
                                                     SessionID=rep(ixSession, nLambdas),
                                                     Condition=rep(numCondition, nLambdas),
                                                     SimulationID=rep(ixSimulation, nLambdas),
                                                     Lambda = lambdas,
                                                     LapseRate = rep(lapseRate, nLambdas),
                                                     Likelihood = likelihoods))
      }
    }
  }
  return(performance)
}


#' Find best weights
#'
#' Chooses best weights out of all weights generated for each lambda
#' of regularized regression.
#' Then, compute glm coefficients on full dataset using the best lambda
#'
#' @param performance all results from regularized regression after running RegularizedRegression
#' @param lambdas list of lambda parameters that was used to run RegularizedRegression
#' @param glmData trial-by-trial subject responses
#' @param nBack number of trials to go back in history. Default is 1.
#'
#' @export
BestWeightsByLambda <- function (performance, 		# GLM weights with different lambda values
                                 lambdas, 		    # Lambda values used in regularized regression
                                 alpha=0,         # Regularization type
                                 glmData,         # trial by trial subject responses
                                 nBack) 			    # n history back
{
  performanceSummary <- ddply(performance, .(SubjectID, SessionID, Lambda), summarise, MedianLhood=median(Likelihood), LapseRate=mean(LapseRate))
  bestLambdas <- ddply(performanceSummary, .(SubjectID, SessionID), summarise, Lambda=Lambda[which.min(MedianLhood)], MinLhood=min(MedianLhood), LapseRate=mean(LapseRate))
  print(bestLambdas)

  allWeights <- data.frame() # For storing final results
  maximumSessions <- max(as.numeric(levels(glmData$SessionID)))
  nParticipants <- length(levels(glmData$SubjectID))
  glmnetFits <- array(list(), c(nParticipants, maximumSessions))

  for (ixSubject in levels(droplevels(glmData$SubjectID))) {
    oneSubjectData <- subset(glmData, glmData$SubjectID == ixSubject)
    subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))
    for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
      oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)
      ifelse(!is.null(oneSessionData$Condition), numCondition <- unique(oneSessionData$Condition), numCondition <- -1)
      contrasts <- levels(as.factor(oneSessionData$Contrast))
      contrasts <- paste("c", contrasts, sep = "")
      ## Find out number of prev success and failure parameters
      if (missing(nBack)) {
        nHistoryBack <- length(grep(chb:::successColName, colnames(glmData)))
      } else {
        nHistoryBack <- nBack
      }
      if (nHistoryBack > 0) {
        histColumnNames <- paste(c(chb:::successColName, chb:::failColName), sort(rep(seq(1:nHistoryBack), 2)), sep="")
        coefs <- paste(c(histColumnNames, contrasts), sep = "")
      } else {
        coefs <- contrasts
      }
      #browser()
      formulaLogistic <- as.formula(paste("as.factor(y) ~", paste(coefs, collapse = "+")))
      ## Fit the model
      glmnetFit <- Glmnet(formula=formulaLogistic, family='binomial', data = oneSessionData, lambda=lambdas, alpha=alpha)
      bestLambda <- with(bestLambdas, bestLambdas[(SubjectID==ixSubject & SessionID==ixSession), ])$Lambda
      lapseRate <- with(bestLambdas, bestLambdas[(SubjectID==ixSubject & SessionID==ixSession), ])$LapseRate
      # Get likelihood computed with best lambda
      glmnetBestWeights <- coef(glmnetFit, s=bestLambda)
      ## Let's store glmnet values for later use
      glmnetFits[[subjectNumber, as.numeric(ixSession)]] <- glmnetFit
      ## Remove redundant "intercept". There is one there already
      glmnetBestWeights <- glmnetBestWeights[-2,]
      nWeights <- length(glmnetBestWeights)

      # Compute AIC. logLhoodBestLambda is already computed as -2*log(L)
      newx <- cbind(`(Intercept)` = rep(1, nrow(oneSessionData)), data.matrix(oneSessionData[, coefs]))
      pRight <- predict(glmnetFit, newx, s = bestLambda, type='response')
      responses <- oneSessionData$y
      #browser()
      logLhoodBestLambda <- Likelihood(responses, pRight, lapseRate)
      aicReg <- logLhoodBestLambda + 2*nWeights # AIC computed for the regularized regression
      # AICc
      aiccReg <- aicReg + ((2*nWeights*(nWeights+1)) / (nrow(oneSessionData) - nWeights - 1))
      # Get logLik for regularized regression
      logLikReg <- logLhoodBestLambda/(-2)

      ## Store best weights for each subject and each session
      allWeights <- rbind(allWeights, data.frame(SubjectID=rep(ixSubject, nWeights),
                                                 SessionID=rep(ixSession, nWeights),
                                                 Condition=rep(numCondition, nWeights),
                                                 Regularized=rep("yes", nWeights),
                                                 Parameter=rownames(data.frame(glmnetBestWeights)),
                                                 Weight=as.numeric(glmnetBestWeights),
                                                 LapseRate=rep(lapseRate, nWeights),
                                                 Vif=NA,
                                                 AIC=aicReg,
                                                 AICc=aiccReg,
                                                 logLhood=logLikReg,
                                                 df=nWeights-1))

      glmFit <- glm(formula=formulaLogistic, family='binomial', data = oneSessionData)
      glmWeights <- coef(glmFit)
      nWeights <- length(glmWeights)
      ## Also calculate Variance Inflation Factor (VIF) to estimate multicolinearity present in data
      getVif <- car::vif(glmFit)
      vifDf <- data.frame(Parameter=names(getVif), Vif=as.numeric(getVif))
      # AIC for unregularized GLM model
      aicUnreg <- summary(glmFit)$aic
      aiccUnreg <- aicUnreg + ((2*nWeights*(nWeights+1)) / (nrow(oneSessionData) - nWeights - 1))
      # Get log likelihood of the model fit
      logLikUnreg <- logLik(glmFit)
      temp <- data.frame(SubjectID=rep(ixSubject, nWeights),
                         SessionID=rep(ixSession, nWeights),
                         Condition=rep(numCondition, nWeights),
                         Regularized=rep("no", nWeights),
                         Parameter=rownames(data.frame(glmWeights)),
                         Weight=as.numeric(glmWeights),
                         LapseRate=rep(lapseRate, nWeights),
                         AIC=aicUnreg,
                         AICc=aiccUnreg,
                         logLhood=logLikUnreg,
                         df=nWeights-1)
      glmFitWeightsAndVif <- merge(temp, vifDf, by="Parameter", all.x=T)
      ## VIF is only available for non-regularized GLM (that is glmFit)
      ## When regularized is used, VIF is set to NA. Also, it is NA for "(Intercept)"
      allWeights <- rbind.fill(allWeights, glmFitWeightsAndVif)
    }
  }
  return(allWeights)
}





#' Model weight for each lambda
#'
#' Function that computes weights for all lambda regularization
#' parameters
#' @export
WeightsWithAllLambdas <- function(performance,
                                  lambdas, 		# Lambda values used in regularized
                                  glmData) 			# subject responses
{

  # bestLambdas <- ddply(performanceSummary, .(SubjectID, SessionID), summarise, Lambda=Lambda[which.min(MedianLhood)], MinLhood=min(MedianLhood))
  # allWeights <- data.frame(ixSubject =  character(), ixSession = character, Lambda=double(), Regularized=character(), Parameter = character(), Weight = double())
  weightsByLambda <- data.frame(ixSubject =  character(), ixSession = character(), Condition=integer(), Lambda=double(), Parameter = character(), Weight = double())
  # glmnetFits is a two-dimensional list to store computed GLMs for each subject and each session
  maximumSessions <- max(as.numeric(levels(glmData$SessionID)))
  nParticipants <- length(levels(glmData$SubjectID))
  glmnetFits <- array(list(), c(nParticipants, maximumSessions))
  for (ixSubject in levels(droplevels(glmData$SubjectID))) {
    oneSubjectData <- subset(glmData, glmData$SubjectID == ixSubject)
    subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))
    for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
      oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)
      numCondition <- unique(oneSessionData$Condition) # The Condition number
      contrasts <- levels(as.factor(oneSessionData$Contrast))
      contrasts <- paste("c", contrasts, sep = "")
      ## Find out number of prev success and failure parameters
      nHistoryBack <- length(grep(successColName, colnames(glmData)))
      if (nHistoryBack > 0) {
        histColumnNames <- paste(c(successColName, failColName), sort(rep(seq(1:nHistoryBack), 2)), sep="")
        coefs <- paste(c(histColumnNames, contrasts), sep = "")
      } else {
        coefs <- contrasts
      }

      formulaLogistic <- as.formula(paste("as.factor(y) ~", paste(coefs, collapse = "+")))
      ## Fit the model
      glmnetFit <- Glmnet(formula=formulaLogistic, family='binomial', data = oneSessionData, lambda=lambdas, alpha=alpha)
      # bestLambda <- with(bestLambdas, bestLambdas[(SubjectID==ixSubject & SessionID==ixSession), ])$Lambda
      # glmnetBestWeights <- coef(glmnetFit, s=bestLambda)
      ## Let's store glmnet values for later use
      glmnetFits[[subjectNumber, as.numeric(ixSession)]] <- glmnetFit
      ## Remove redundant "intercept". There is one there already
      #glmnetBestWeights <- glmnetBestWeights[-2,]
      #nWeights <- length(glmnetBestWeights)
      ## Store best weights for each subject and each session
      # allWeights <- rbind(allWeights, data.frame(SubjectID=rep(ixSubject, nWeights),
      # 											SessionID=rep(ixSession, nWeights),
      # 											Regularized=rep("yes", nWeights),
      # 											Parameter=rownames(data.frame(glmnetBestWeights)),
      # 											Weight=as.numeric(glmnetBestWeights)))
      ## Also, store all other weights (not just the best one) calculated for each lambda
      coefsByLambdas <- coef(glmnetFit, s=lambdas)[-2,]
      nCoefsByLambdas <- length(rownames(coefsByLambdas))
      for (ixLambda in lambdas) {
        weightsByLambda <- rbind(weightsByLambda, data.frame(SubjectID=rep(ixSubject, nCoefsByLambdas),
                                                             SessionID=rep(ixSession, nCoefsByLambdas),
                                                             Condition=rep(numCondition, nCoefsByLambdas),
                                                             Lambda=rep(lambdas[lambdas==ixLambda], nCoefsByLambdas),
                                                             Parameter=rownames(coefsByLambdas),
                                                             Weight=as.numeric(coefsByLambdas[,lambdas==ixLambda])))
      }
    }
  }
  return(weightsByLambda)
}



#' Compute correlation matrices for each parameter of the GLM model
#'
#' Helps to check if there are any strong effects of multicolinearity
#' @export
CorrelationMatrixOfParameters <- function(glmData) {
	allCorrelations <- data.frame(SubjectID=as.character(), SessionID=as.character(), Condition=as.integer())
	for (ixSubject in levels(droplevels(glmData$SubjectID))) {
		oneSubjectData <- subset(glmData, SubjectID == ixSubject)
		subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))
		for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
			oneSessionData <- subset(oneSubjectData, SessionID==ixSession)
			numCondition <- unique(oneSessionData$Condition) # The Condition number
			contrastColumns <- grep('c0', colnames(oneSessionData))
			ixModelParams <- c(contrastColumns, grep(successColName, colnames(oneSessionData)), grep(failColName, colnames(oneSessionData)))
			modelParams <- oneSessionData[, ixModelParams]
			# Exclude contrast intensities which come from other runs and not relevant to this run
			modelParams <- modelParams[, sapply(abs(modelParams), sum) != 0]
			c <- cor(modelParams)
			c.m <- melt(c)
			# Remove correlations with itself
			c.m[which(c.m$X1==c.m$X2),]$value <- NA
			oneSessionCorrelations <- cbind(data.frame(SubjectID=ixSubject, SessionID=ixSession, Condition=numCondition), c.m)
			allCorrelations <- rbind(allCorrelations, oneSessionCorrelations)
		}
	}
	return(allCorrelations)
}


#' Compute length of sequences on left or right sides
#'
#' This is to figure out if there are long sequences, particularly when introducing biases
#' such as succeed/stay bias which can generate long sequences on one side
#' @export
StimulusSideSequences <- function(glmData) {
  allSequences <- data.frame(SubjectID=as.character(), SessionID=as.character(), Condition=integer())
  for (ixSubject in levels(droplevels(glmData$SubjectID))) {
    oneSubjectData <- subset(glmData, glmData$SubjectID == ixSubject)
    #subjectNumber <- which(ixSubject == levels(droplevels(glmData$SubjectID)))
    for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
      oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)
      numCondition <- unique(oneSessionData$Condition)
      visualField <- oneSessionData$VisualField
      ## calculate length of sequences of 1 (L) or 2 (R)
      stimSideAndSequence <- rle(visualField)
      sequenceLength <- stimSideAndSequence$lengths
      stimSide <- stimSideAndSequence$values
      sequenceNumber <- 1:length(sequenceLength)
      oneSessionSequences <- cbind(data.frame(SubjectID=ixSubject, SessionID=ixSession, Condition=numCondition),
                                   SequenceNumber=sequenceNumber, StimSide=stimSide, SequenceLength=sequenceLength)
      allSequences <- rbind(allSequences, oneSessionSequences)
    }
  }
  return(allSequences)
}


#' Compute weights after randomizing trials or responses
#'
#' This is checker function to ensure that history weights do indeed
#' disapper after trial order is made random or subjects responses are made
#' random. This randomization effectively destroys history effects effectively
#' setting history weights to 0.
#' !!! CHECK IF YOU HAVE ANOTHER FUNCTION THAT DOES THIS!!!! It is a suspect, because
#' uses rawData as input and not glmData, which is much more clean data.
#' Again, compute history weights after either shuffling trial order or shuffling
#' responses to each trial, by either scrambling the order of
#' trials, or by generating responses randomly. This procedure should
#' decrease the history weights to 0. If we run this function many times we
#' will get confidence intervals for history biases
#' @export
BiasAfterRandomization <- function(rawData, 		# rawData as input
                                   randomizationType=1, # Use 1 to scramble trials; use 2 to randomly generate responses to each trial (this is depricated); use 3 to scramble responses but keep trial order
                                   alpha=0, 		# when alpha is 0 use "ridge" regularization. When 1 run with "Lasso" regularization
                                   nBack=1, 		# how many trials go back in history
                                   B=10 			  # Number of simulations to find best lambda for the regularized regression
) {
  ## Names of history columns
  successColName <- "PrevCorr"
  failColName <- "PrevFail"

  ## Scramble rawData. This destroys history information
  if (randomizationType == 1) {
    scrambledRawData <- rawData[sample(nrow(rawData)),]
  }

  ## Generate random responses. This also destroys history information. THIS IS DEPRICATED
  if (randomizationType == 2) {
    print('Randomization type 2 is depricated')
    return()
    scrambledRawData <- rawData
    scrambledRawData$Response <- c(1,2)[rbinom(nrow(scrambledRawData), size=1, 0.5)+1]
    #browser()
  }

  ## Scramble responses
  if (randomizationType == 3) {
    scrambledRawData <- rawData
    scrambledRawData$Response <- sample(scrambledRawData$Response)
  }

  # Change CorrIncorr based on newly generated responses
  scrambledRawData$CorrIncorr <- 1
  ## Indices of correct responses
  ixCorrect <- which(scrambledRawData$VisualField == scrambledRawData$Response)
  ## This is redundant but for clarity
  scrambledRawData$CorrIncorr[ixCorrect] = 1
  ixIncorrect <- which(scrambledRawData$VisualField != scrambledRawData$Response)
  scrambledRawData$CorrIncorr[ixIncorrect] = 0
  ## When responses are NaN or other than left or right assign NaN to CorrIncorr
  ixIlligalResponses <- which(is.nan(scrambledRawData$Response) | ((scrambledRawData$Response!=1) & (scrambledRawData$Response!=2)))
  scrambledRawData$CorrIncorr[ixIlligalResponses] = NaN
  scrambledRawData$y[ixIlligalResponses] = NaN

  # Assign y based on new set of simulated responses
  # Prepare GLM y response
  scrambledRawData$y <- scrambledRawData$Response
  scrambledRawData$y[scrambledRawData$y == 1] <- -1
  scrambledRawData$y[scrambledRawData$y == 2] <- 1

  ## Prepare data to run logistic regression
  scrambledGlmData <- BuilDataForGLM(scrambledRawData, nHistoryBack=nBack, nDummyParams=0, successColName=successColName, failColName=failColName)
  ## Remove trials when participants didn't respond or pressed buttons other than left or right
  scrambledGlmData <- scrambledGlmData[!is.nan(scrambledGlmData$CorrIncorr), ]
  ## Fit regularized logistic model to the scrambled trials
  performance <- FitRegulirizedLogisticRegression(scrambledGlmData, lambdas=lambdas, alpha=alpha, B=B)
  ## Compute logistic regression weights using "best lambdas"
  weights <- ComputeWeightsWithBestLambda(performance, lambdas, scrambledGlmData)
  #browser()
  weights <- subset(weights, weights$Regularized=="yes")

  return(weights)
}

#' Labels assigned to each experimental condition
#'
#' On 09 Jan 2014, swapped labels for fail/stay and fail/success conditions.
#' Initially, Condition 2 was fail/switch, but then realized that switching trial
#' position after a failure actually encourages to stay on the same side where failure
#' happened. Conversly, fail/stay (stim presented on same side) actuallly encourages
#' to switch choice on next trial cause failure happened on the opposite side.
#' @export
ConditionLabels <- function(variable, value) {
  if (variable=='Condition' | variable=='Condition.y') {
    newValue <- c()
    newValue[which(value==1)] <- 'Natural bias\n'
    newValue[which(value==2)] <- 'Fail-stay bias'
    newValue[which(value==3)] <- 'Success-stay bias'
    newValue[which(value==4)] <- 'Natural bias: stim diam\n 6deg, eccnt 8deg'
    newValue[which(value==5)] <- 'Natural bias: stim diam\n 12deg, eccnt 12deg'
    newValue[which(value==6)] <- 'Natural bias: stim diam\n 6deg, eccnt 10deg'
    newValue[which(value==7)] <- 'Fail-switch bias'
    newValue[which(value==8)] <- 'Success-switch bias'
    newValue[which(value==9)] <- 'Induced bias corr~incorr: \n p stay after failure 80%'
    newValue[which(value==10)] <- 'Induced bias corr~incorr: \n p stay after success 80%'
    newValue[which(value==11)] <- 'Induced bias corr~incorr: \n p switch after failure 80%'
    newValue[which(value==12)] <- 'Induced bias corr~incorr: \n p switch after success 80%'
    newValue[which(value==13)] <- 'Natural bias, UCL'
    newValue[which(value==1000)] <- 'Difference'
    value <- newValue
  }
  return(as.character(value))
}


#' Plot model parameters as a distance from vertical
#'
#' Forest plot provides convenient way of plotting model parameters to see how far
#' they deviate from 0
#' @export
ForestPlot <- function(plotData, 						# Data to plot
                       xlims=c(-3, 3),					# Range of values along x axis
                       xbreaks=c(-2, -1, 0, 1, 2), 	# Break to be displayed on x axis
                       figureWidth=4.820225, 			# Width of the plot
                       figureHeight=7.896552, 			# Height of the plot
                       plotMeans=TRUE,					# Mean values of weights
                       plotMeanPointsHollow=TRUE, 		# Plot mean points as hollow dots which is easier to visually read
                       plotEachWeight=TRUE,			# Plot weights from each run
                       plot95ConfInt=FALSE,			# Plot error bars
                       plotBackground=FALSE,			# Show or not a grey background for easy reading of visual info
                       plotByCondition=TRUE, 			# Plot by each condition
                       plotByParameter=FALSE) {		# Plot each parameter as a separate column
  dev.new(width=figureWidth, height=figureHeight)
  thePlot <- ggplot(data=plotData, aes(x=Parameter, y=Weight, group=SessionID, colour=Parameter)) +
    geom_hline(yintercept = 0, size = 0.5, colour = "grey30", linetype = "dashed") +
    #geom_segment(aes(xend=Parameter), yend=0, colour="grey50") +
    theme_few() +
    theme(strip.background=element_rect(colour="white", fill="white")) +
    theme(panel.background=element_rect(colour="white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
    scale_colour_manual(values=c("#a4a4a4",prevSuccessColor, prevFailColor), breaks=c("B0", "pS", "pF"), labels=c("Bias", "Prev success", "Prev failure")) +
    ylab("Weight") +
    xlab("") +
    scale_x_discrete(limits=c("PrevFail1", "PrevCorr1", "(Intercept)"), breaks=NULL) +
    scale_y_continuous(limits=xlims, breaks=xbreaks) +
    theme(axis.line = element_line(colour = "black", size = 0.3),axis.line.y = element_blank()) +
    theme(legend.position = c(0.15,0.5)) +
    coord_flip()

  if (plotEachWeight) thePlot <- thePlot + geom_point(size=1.0, alpha=0.8)
  if (plot95ConfInt) thePlot <- thePlot + stat_summary(aes(group=Parameter),
                                                       fun.data = "mean_cl_boot",
                                                       B=500, conf.int = 0.95,
                                                       geom = "errorbar", size = 0.7, width = 0.0)
  if (plotMeans)	{
    if (plotMeanPointsHollow) {
      thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=5)
      if (plotBackground)
        thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=3, colour="#f0f0f0")
      else thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=3, colour="white")
    }
    else thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=3)
  }
  if (plotBackground) thePlot <- thePlot + theme(panel.background=element_rect(fill="#f0f0f0"))
  if (plotByCondition) {
    thePlot <- thePlot + facet_grid(SubjectID~Condition, labeller=ConditionLabels)
  }
  if (plotByParameter) thePlot <- thePlot + facet_grid(SubjectID~Parameter)
  if ((!plotBackground) & (!plotByCondition)) thePlot <- thePlot + facet_grid(SubjectID~.)

  thePlot
}




#' Fit probit function
#'
#' Fit probit (cummmulative Gaussian). The best fit after Weibull (with
#' log transformed contrast values)
FitProbit <- function(dat) {
	glm(data=dat, cbind(nYesR, nNoR)~Contrast, binomial(probit))
}

#' Fit probit function
#'
#' Fit logistic function
FitLogit <- function(dat) {
	glm(data=dat, cbind(nYesR, nNoR)~Contrast, binomial(logit))
}

#' Fit Weibull function
#'
#' Fit Weibull without log transforming contrast values
#' AVOID using this because the fit is poor
FitWeibull <- function(dat) {
  ## AVOID using this because the fit is poor
  ## Instead use FitWeibullLogContrast
  glm(data=dat, cbind(nYesR, nNoR)~Contrast, binomial(cloglog))
}


#' Fit Weibull after log transforming contrast values
#'
#' This produces the best fit out of all psychometric curves
#' NOTE: when contrast is negative, an easier approach is to fit probit
#' Tried Gumbel, but it didn't provide a good fit.
FitWeibullLogContrast <- function(data) {
	## Log transform Contrast before fitting
	data$Contrast <- log10(data$Contrast)  ## Applied log transform here cause when done in formula, the name of the column appears clumsy.
	glm(data=data, cbind(nYesR, nNoR)~Contrast, binomial(cloglog))
}


#' Compute contrast threshold using fitted probit model
#'
#' Estimates contrast intensities which generate response
#' at probablity p (can be an array). The fitted psychometric
#' parameters are in model glm.
#'
#' @param p either a single value or an array of probabilities for which contrast thresholds will be computed
#' @param model a glm object of fitted probit psychometric curve
ProbitThreshold <- function(p,        # probability
                            model){   # glm object
  coefs <- coef(model)
  nhu <- -coefs[1]/coefs[2]   # mean of Gaussian
  sigma <- 1/coefs[2]         # standard deviation of Gaussian
  thP <- qnorm(p, nhu, sigma)
  th50 <- qnorm(0.5, nhu, sigma)
  # Threshold is contrast increment that will raise the threshold from 50% to p percent
  th <- thP - th50
  # If p=50% was requested, return 50% threshold as is
  # which will make it 0
  th[p==0.5] <- th50
  return(th)
}

#' Estimate threshold and slope using Weibull
#'
#' DON"T USE THIS. The most accurate threshold estimation is probit.
#' Compute threshold and slope of psychometric function which was fitted with a
#' Weilbull. takes in a glm object and probability at which to estimate the threshold
#' Adapted from "Modelling Psychophysical Data in R", p. 155
WeibullThAndSlope <- function(p,        # probability at which to compute the contrast
                              model)		  # glm object
{  # This was added to GLM model to aid log transform of negative contrast values
  if (length(p) > 1) {
    print("*WeibullThAndSlope* Please provide only one p value")
    return()
  }
  ## Extracting fitted parameters
  ln10 <- log(10) ## This is to backtransform Contrast from log10 to original values
  coefs <- coef(model)
  th <- qweibull(p, shape=coefs[2]/ln10, scale=exp(-ln10 * coefs[1]/coefs[2])) - offset
  browser()
  weibParams <- c(th, coefs[2]/ln10)
  names(weibParams) <- c("th", "slope")
  weibParams
}

############################################
## Given a model, predict values of yvar from xvar
## This supports one predictor and one predicted variable
## xrange: If NULL, determine the x range from the model object. If a vector with # two numbers, use those as the min and max of the prediction range.
## samples: Number of samples across the x range.
## ...: Further arguments to be passed to predict()
PredictvalsProbit <- function(model, xvar, yvar, xrange=NULL, samples=100, ...) {
  ## If xrange isn't passed in, determine xrange from the models.
  ## Different ways of extracting the x range, depending on model type
  if (is.null(xrange)) {
    if (any(class(model) %in% c("lm", "glm"))) xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess")) xrange <- range(model$x)
  }
  newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)

  newdata[["slope"]] <- coef(model)[2]
  ## Get 50% and 75% thresholds
  ths <- ProbitThreshold(c(0.5, 0.75), model)

  newdata[["Th50"]] <- ths[1]
  newdata[["Th75"]] <- ths[2]
  ##newdata[["slope"]] <- exp(coef(model)[2])  ## This might need to be used with Weibull function but not with probit (see Modelling Psychophysical Data in R, p. 152)
  #browser()
  newdata
}

############################################
## Given a model, predict values of yvar from xvar
## The slope and threshold parameters for Weibull are
## calculated differently to logistic or porbit
## See "Modelling Psychophysical Data in R", p. 154
PredictvalsWeibull <- function(model, xvar, yvar, xrange=NULL, samples=100, ...) {
  ## If xrange isn't passed in, determine xrange from the models.
  ## Different ways of extracting the x range, depending on model type
  if (is.null(xrange)) {
    if (any(class(model) %in% c("lm", "glm"))) xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess")) xrange <- range(model$x)
  }
  newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  ## Let's now backtransform the data
  newdata$Contrast <- 10^newdata$Contrast
  ## Get 75% threshold and slope
  thAndSlope <- WeibullThAndSlope(0.75, model)
  newdata[["Th75"]] <- thAndSlope[1]
  newdata[["slope"]] <- thAndSlope[2]
  ##newdata[["slope"]] <- exp(coef(model)[2])  ## This might need to be used with Weibull function but not with probit (see Modelling Psychophysical Data in R, p. 152)
  newdata
}




####################################################################################
## Compute and display run and subject statistics
DataStatistics <- function(weights,  				# Regularized weights
							figureWidth=5.68, 		# Width of the plot
							figureHeight=15.51, 	# Height of the plot
							plotByCondition=TRUE,	# Plot by Condition. Otherwise, collapse across Condition
							plotSubjectMeans		# Compute means across Sessions for each subject and plots those means
									   ) {
# TODO


}


####################################################################################
## Sorted plot of history weights
## THIS IS NOT WORKING CURRENTLY. NEEDS FIXING
SortedPlotOfHistoryWeights <- function(weights,  				# Regularized weights
									  figureWidth=5.68, 	# Width of the plot
									  figureHeight=15.51, 	# Height of the plot
									  plotByCondition=TRUE,		# Plot by Condition or collapsed across Condition
									  plotSubjectMeans=TRUE,			# Compute means across Sessions for each subject and plots those means
									  plotEachWeight=T

									   ) {

thePlot <- ggplot(data= weights, aes(x=Parameter, y=Weight, group=SessionID, colour=Parameter)) +
	geom_hline(yintercept = 0, size = 0.5, colour = "grey30", linetype = "dashed") +
	theme_few() +
	theme(strip.background=element_rect(colour="white", fill="white")) +
	theme(panel.background=element_rect(colour="white")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
	scale_colour_manual(values=c("#a4a4a4",prevSuccessColor, prevFailColor), breaks=c("B0", "pS", "pF"), labels=c("Bias", "Prev success", "Prev failure")) +
	ylab("Weight") +
	xlab("") +
	scale_x_discrete(limits=c("PrevFail1", "PrevCorr1", "(Intercept)"), breaks=NULL) +
	#scale_y_continuous(limits=xlims, breaks=xbreaks) +
	theme(axis.line = element_line(colour = "black", size = 0.3),axis.line.y = element_blank()) +
	theme(legend.position = c(0.15,0.5)) +
	geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=5) +
	coord_flip()

	# if (plotEachWeight) thePlot <- thePlot + geom_point(size=1.0, alpha=0.8)
 	# if (plot95ConfInt) thePlot <- thePlot + stat_summary(aes(group=Parameter),
 														# fun.data = "mean_cl_boot",
 														# B=500, conf.int = 0.95,
 														# geom = "errorbar", size = 0.7, width = 0.0)
	# if (plotMeans)	{
		# if (plotMeanPointsHollow) {
	        # thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=5)
	        # if (plotBackground)
	        	# thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=3, colour="#f0f0f0")
			# else thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=3, colour="white")
		# }
		# else thePlot <- thePlot + geom_point(aes(group= Parameter), stat="summary", fun.y="mean", alpha=1.0, size=3)
	# }
    # if (plotBackground) thePlot <- thePlot + theme(panel.background=element_rect(fill="#f0f0f0"))
	# if (plotByCondition) thePlot <- thePlot + facet_grid(SubjectID~Condition, labeller=ConditionLabels)
	# if (plotByParameter) thePlot <- thePlot + facet_grid(SubjectID~Parameter)
	# if ((!plotBackground) & (!plotByCondition)) thePlot <- thePlot + facet_grid(SubjectID~.)

    thePlot
}


#######################################
## WRITE A DESCRIPTION FOR THIS ONE
PlotOrderedWeights <- function(weights, 	## Weight to be plotted. Need to be from same condition and only one type, such as fail or success
								weightToPlot = 'PrevFail1',		# Name of the weigh to be plotted
								conditionToPlot = 1,			# Single condition to be plotted
								subjectsToPlot = 'all', 		# Default is to plot all subjects
								geomPointColor='black',			# Color of geom_point
								figureWidth=3.851064, 			# Width of the plot
								figureHeight=5.755319 			# Height of the plot
								) {

	## Subset data based on weight to be plotted and condition that needs to be plotted
	#weightsToPlot <- subset(weights, Parameter==weightToPlot & Condition==conditionToPlot)
	weightsToPlot <- subset(weights, Condition==conditionToPlot)

	## Select subjects to be plotted
	if (subjectsToPlot != 'all') weightsToPlot <- subset(weightsToPlot, SubjectID %in% subjectsToPlot)

	## Compute means weighs and standard error
	meanWeights <- ddply(weightsToPlot, .(SubjectID, Parameter), summarise, MeanWeight=mean(Weight), StDev=sd(Weight), se=StDev/sqrt(length(Weight)))
	meanWeights <- meanWeights[, c(1,2,3)]
	colnames(meanWeights) <- c("SubjectID", "variable", "value")
	meanWeights <- droplevels(cast(meanWeights))

	#if (weightToPlot=='PrevFail1') meanWeights$SubjectID <- reorder(meanWeights$SubjectID, -meanWeights$MeanWeight)
	#if (weightToPlot=='PrevCorr1') meanWeights$SubjectID <- reorder(meanWeights$SubjectID, -as.numeric(c(15, 12,  5,  6,  4,  7,  8, 13,  9, 14,  1,  2,  3, 10, 11)))
	#if (weightToPlot=='(Intercept)') meanWeights$SubjectID <- reorder(meanWeights$SubjectID, -c(15, 12,  5,  6,  4,  7,  8, 13,  9, 14,  1,  2,  3, 10, 11))

	#browser()

	dev.new(width=figureWidth, height=figureHeight)
	#if (weightToPlot=='PrevFail1') p <- ggplot(data=meanWeights, aes(y=PrevFail1, x=(reorder(SubjectID, -PrevFail1))))
	#if (weightToPlot=='PrevCorr1') p <- ggplot(data=meanWeights, aes(y=PrevCorr1, x=(reorder(SubjectID, -PrevFail1))))
	#if (weightToPlot=='(Intercept)') p <- ggplot(data=meanWeights, aes(y=(Intercept), x=(reorder(SubjectID, -PrevFail1))))


	if (weightToPlot=='PrevFail1') p <- ggplot(data=meanWeights, aes(y=PrevFail1, x=SubjectID))
	if (weightToPlot=='PrevCorr1') p <- ggplot(data=meanWeights, aes(y=PrevCorr1, x=SubjectID))

	#p <- ggplot(data=meanWeights, aes(y=PrevFail1, x=SubjectID))
	#p <- ggplot(data=meanWeights, aes(y=PrevCorr1, x=SubjectID))


	#ggplot(data=meanWeights, aes(y=MeanWeight, x=reorder(SubjectID, -MeanWeight))) +
	#ggplot(data=meanWeights, aes(y=MeanWeight, x=SubjectID)) +
	p <- p + theme_few() +
		theme(strip.background=element_rect(colour="white", fill="white")) +
		theme(panel.background=element_rect(colour="white")) +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
		theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
		scale_y_continuous(limits=c(-2.2, 2.2)) +
		geom_hline(yintercept=0.0, size=0.5, colour="#a9a9a9", linetype = "solid") +
		geom_segment(aes(xend=SubjectID), yend=0, colour=geomPointColor, size=1.5) +
	  	#geom_errorbar(aes(ymin=MeanWeight-se, ymax=MeanWeight+se), width=0.01, alpha=0.2) +
	  	#geom_point(size=1.5, )
	  	geom_point(size=6, color=geomPointColor) +
	  	geom_point(size=3, color='white') +
	  	xlab("")  +
	  	ylab("") +
		theme(axis.line = element_line(colour = "#a9a9a9", size = 0.3),axis.line.y = element_blank()) +
	  	coord_flip()

	print(p)
}




#######################################
## REDESIGNING THE ABOVE
PlotOrderedWeights_1 <- function(weights, 	## Weight to be plotted. Need to be from same condition and only one type, such as fail or success
                                 weightToPlot = 'PrevFail1',		# Name of the weigh to be plotted
                                 conditionToPlot = 1,			# Single condition to be plotted
                                 subjectsToPlot = 'all', 		# Default is to plot all subjects
                                 geomPointColor='black',			# Color of geom_point
                                 figureWidth=3.851064, 			# Width of the plot
                                 figureHeight=5.755319 			# Height of the plot
) {

  ## Subset data based on weight to be plotted and condition that needs to be plotted
  weightsToPlot <- subset(weights, Parameter==weightToPlot & Condition==conditionToPlot)

  ## Select subjects to be plotted
  if (subjectsToPlot != 'all') weightsToPlot <- subset(weightsToPlot, SubjectID %in% subjectsToPlot)

  ## Compute means weighs and standard error
  meanWeights <- ddply(weightsToPlot, .(SubjectID, Parameter), summarise, MeanWeight=mean(Weight), StDev=sd(Weight), se=StDev/sqrt(length(Weight)))

  dev.new(width=figureWidth, height=figureHeight)
  ggplot(data=meanWeights, aes(y=MeanWeight, x=reorder(SubjectID, -MeanWeight))) +
    theme_few() +
    theme(strip.background=element_rect(colour="white", fill="white")) +
    theme(panel.background=element_rect(colour="white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
    theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
    scale_y_continuous(limits=c(-2.2, 2.2)) +
    geom_hline(yintercept=0.0, size=0.5, colour="#a9a9a9", linetype = "solid") +
    geom_segment(aes(xend=SubjectID), yend=0, colour=geomPointColor, size=1.5) +
    #geom_errorbar(aes(ymin=MeanWeight-se, ymax=MeanWeight+se), width=0.01, alpha=0.2) +
    #geom_point(size=1.5, )
    geom_point(size=6, color=geomPointColor) +
    geom_point(size=3, color='white') +
    xlab("")  +
    ylab("") +
    theme(axis.line = element_line(colour = "#a9a9a9", size = 0.3),axis.line.y = element_blank()) +
    coord_flip()
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


##########################################################
## This is a code to check if last 10 to 20 trials, which where were sometimes only presented on one side
## particularly during the success/stay bias induction condition, didn't cause any trouble when computing weights
RemoveLastTrialsFromRun <- function(dat) {
	if (unique(dat$Condition) %in% c(1,2,3,7,8)) {
		dat <- head(dat, -20) ## Return data minus last 20 trials
	}
	return(dat)
}

## PlotSuccessFailOnEachSide -----------------------------------------------
## Intuitive plot of stimulus presentation side over time.
## Helps to discover unusual patterns such as long streaks
## when inducing bias in success/stay condition.
PlotSuccessFailOnEachSide <- function(inputData,   	 	# Input data in the format of "glmData"
                                      plotResults=F, 	# If False, returns list of ggplot objects. If True, plots those objects, each subject in separate window
                                      plotInColor=T		# Success and failure are marked by two different colors. Otherwise, black and white
){
  ## Add column with trial numbers
  inputData <- droplevels(ddply(inputData, .(SubjectID, SessionID), mutate, TrialID=1:length(VisualField)))

  ## Generate ggplot graphs as binary sparklines
  ## Store those graphs for individual plotting because plotting them
  ## altogether has many empty spaces between sessions
  nSubjects <- length(levels(inputData$SubjectID))
  g <- vector(mode="list", length=nSubjects)
  for (ixSubject in levels(droplevels(inputData$SubjectID))) {
    oneSubjectData <- subset(inputData, SubjectID==ixSubject)
    subjectIndex <- which(levels(inputData$SubjectID)==ixSubject)

    g[[subjectIndex]] <- ggplot(data=oneSubjectData, aes(x=TrialID, y=VisualField)) +
						      theme_few()

	if (plotInColor) {
#		g[[subjectIndex]] <- g[[subjectIndex]] + geom_linerange(subset=.(CorrIncorr==0), aes(ymin=1.5, ymax=VisualField), size=0.3, color='#b62813') +
#		      									 geom_linerange(subset=.(CorrIncorr==1), aes(ymin=1.5, ymax=VisualField), size=0.3, color='#a7b112')

		g[[subjectIndex]] <- g[[subjectIndex]] + geom_linerange(subset=.(CorrIncorr==0), aes(ymin=1.5, ymax=VisualField), size=0.3, color='#d7191c') +
		      									 geom_linerange(subset=.(CorrIncorr==1), aes(ymin=1.5, ymax=VisualField), size=0.3, color='#2b83ba')


	}
	else {
		g[[subjectIndex]] <- g[[subjectIndex]] + geom_linerange(aes(ymin=1.5, ymax=VisualField), size=0.3)
	}

	g[[subjectIndex]] <- g[[subjectIndex]] + scale_x_continuous(expand=c(0.01,0.01)) +
						      scale_y_continuous(breaks=c(1,2), labels=c("L","R"), limits=c(0.8, 2.2)) +
						      #facet_grid(SessionID+Condition~., labeller=ConditionLabels) +
						      #facet_wrap(SubjectID~Condition+SessionID) +
						      facet_grid(SubjectID~Condition+SessionID, labeller=ConditionLabels) +
						      theme(strip.text=element_text(size=8, angle=90)) +
						      theme(axis.text=element_text(size=8)) +
						      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) +
						      ylab("Stimulus presentation side") +
						      xlab("Trial number") +
						      coord_flip()
  }
	## Plot each subject in a separate window.
	## Each window width depends on number of runs for that subject
	if (plotResults){
		## Results are plotted for each subject
		for (subjectIndex in 1:nSubjects) {
			nRuns <- length(unique(g[[subjectIndex]]$data$SessionID)) ## a shortcut to get number of runs from ggplot object
			plotWidth <- nRuns * 0.4
			dev.new(width=plotWidth, height=7.22)
			print(g[subjectIndex])
		}
	}
} ## -- PlotSuccessFailOneachSide --



################################## ComputeFailRate #########################################
## Compute failure rate.
ComputeFailRate <- function(inputData, 	# trail by trial data
							condition,  # One condition for which stats will be computed
							subjects	# Subject for whom stats will be computed
							) {

	selectData <- subset(inputData, Condition %in% condition)
	selectData <- subset(selectData, SubjectID %in% subjects)
	selectData <- droplevels(selectData)
	res <- table(selectData$CorrIncorr) / nrow(selectData)
	failRate <- res[1]
	return(failRate)
}


################################## PlotSlopeVsBias #########################################
# Plot change of slope as a function if failure bias
PlotSlopeVsBias <- function(inputData,          # Data of type slopeSimData
                            weightsToPlot) {    # Weights that will be plotted. Default all weights will be plotted. Subjects will always be plotted.
  ## Plot proportion correct of responding right to stimuli presented either to left or to right
  pcRight <- inputData
  ## Label left responses to right gratings with negative contrast. Right responses to right gratings will remain with positive sign
  pcRight[pcRight$VisualField == 1,]$Contrast <- pcRight[pcRight$VisualField == 1,]$Contrast * -1
  ## Summary of responses to gratings presented in the right
  pcRightSummary <- ddply(pcRight, .(SubjectID, SessionID, VisualField, Contrast, IsSubject, PrevFailWeight), summarise,
                          nRightResp = sum(Response == 2),
                          nLeftResp = sum(Response == 1),
                          nRStim = sum(VisualField == 2),
                          nLStim = sum(VisualField == 1))
  pcRightSummary <- ddply(pcRightSummary, .(SubjectID, SessionID, VisualField, Contrast, IsSubject, PrevFailWeight), summarise,
                          pRightCorrect = nRightResp / (nRStim + nLStim),
                          nYesR=nRightResp,
                          nNoR=nLeftResp)
  ## Convert contrast into %
  pcRightSummary$Contrast <- pcRightSummary$Contrast * 100
  models <- dlply(pcRightSummary, c("SubjectID", "IsSubject", "PrevFailWeight"), .fun=FitProbit)  ## Probit analysis

  predvals <- ldply(models, .fun=PredictvalsProbit, xvar="Contrast", yvar="pRightCorrect", type="response")
  ## Further summarise results cause we need slopes only
  predvals1 <- ddply(predvals, .(SubjectID, IsSubject, PrevFailWeight), summarise, slope=mean(slope), th75=mean(Th75))
  ## Get percent change of slope relative to slope when error weight is 0
  predvals1 <- ddply(predvals1, .(SubjectID), transform, SlopeChange = 1- (slope[which(PrevFailWeight==0)] / slope))

  ## Plot change of slope as a function of failure bias
  if (!missing(weightsToPlot)) {
    datForPlot <- subset(predvals1, (PrevFailWeight %in% weightsToPlot) | (IsSubject==T))
  } else datForPlot <- predvals1
  #browser()
  errbarDat <- subset(droplevels(datForPlot), PrevFailWeight==0 | IsSubject==T)
  # Compute within-subjects error bars using Winston Chan's functions
  datForPlot$PrevFailWeightPlot <- datForPlot$PrevFailWeight
  datForPlot$PrevFailWeightPlot[datForPlot$IsSubject] <- mean(datForPlot$PrevFailWeight[datForPlot$IsSubject])
  datForPlot <- droplevels(datForPlot)

  source('~/Desktop/Dropbox/R/Lib/StandardErrorsByWinstonChan.R')
  datForPlotErrBars <- summarySEwithin(datForPlot, measurevar='slope', withinvars=c('PrevFailWeightPlot'), idvar='SubjectID')
  datForPlotErrBars$PrevFailWeightPlot <- as.numeric(as.character(datForPlotErrBars$PrevFailWeightPlot))

  ggplot(data=datForPlot, aes(x=PrevFailWeightPlot, y=slope)) +
    theme_publish2() +
    xlab("Failure bias (weight)") + ylab("Slope") +
    coord_cartesian(ylim=c(0.468, 0.52), xlim=c(-2.5, 2.5)) +
    geom_errorbar(data=datForPlotErrBars, aes(ymin = slope - ci, ymax = slope + ci), width = 0.00, size = 0.2) +
    geom_line(subset=.(!IsSubject), stat="summary", fun.y=mean, size=1, color="grey80") +
    geom_point(subset=.(!IsSubject), stat="summary", fun.y=mean, size=7, color="white") +
    geom_point(subset=.(IsSubject), aes(x=mean(PrevFailWeight), y=slope), stat="summary", fun.y=mean, size=7, color="white") +
    geom_point(subset=.(!IsSubject), stat="summary", fun.y=mean, size=4, shape=1) +
    geom_point(subset=.(IsSubject), aes(x=mean(PrevFailWeight), y=slope), stat="summary", fun.y=mean, size=4)

  browser()

  # Compute some stats for the paper
  sbjAndIdealSbj <- droplevels(subset(datForPlot, IsSubject==TRUE | PrevFailWeight==0))
  sbjAndIdealSbj$DataSet <- NA
  RIKENSbjCodes <- paste('s', sprintf('%03d',seq(1,15)), sep='')
  UCLSbjCodes <- paste('s', sprintf('%03d',seq(16,27)), sep='')
  colnames(castedSlope)  <- c('SubjectID', 'NoBias', 'Bias')

  sbjAndIdealSbj$DataSet[which(sbjAndIdealSbj$SubjectID %in% RIKENSbjCodes)] <- 'RIKEN'
  sbjAndIdealSbj$DataSet[which(sbjAndIdealSbj$SubjectID %in% UCLSbjCodes)] <- 'UCL'

  castedSlope <- cast(data=sbjAndIdealSbj, SubjectID+DataSet~IsSubject, value=.(slope))
  t.test(castedSlope$NoBias, castedSlope$Bias, paired=TRUE)

  castedTh <- droplevels(cast(data=sbjAndIdealSbj, SubjectID+DataSet~IsSubject, value=.(th75)))
  colnames(castedTh)  <- c('SubjectID', 'DataSet', 'NoBias', 'Bias')
  t.test(castedTh$Bias, castedTh$NoBias, paired=TRUE)

  # Analyse the rest of biases, not just subject and no bias
  datForPlot$FailWeightSbjUnited <-datForPlot$PrevFailWeight
  datForPlot$FailWeightSbjUnited[datForPlot$IsSubject==TRUE] <- -500
  datForPlot$DataSet <- NA
  datForPlot$DataSet[datForPlot$SubjectID %in% RIKENSbjCodes] <- 'RIKEN'
  datForPlot$DataSet[datForPlot$SubjectID %in% UCLSbjCodes] <- 'UCL'

  ddply(datForPlot, .(FailWeightSbjUnited, DataSet), numcolwise(mean))
  datForPlotWideSlope <- cast(datForPlot, SubjectID+DataSet~FailWeightSbjUnited, value=.(slope))
}


################################## PlotSlopeVsBias #########################################
# Plot threshold vs subject bias to see if the magnitude of the bias
# depends on subject sensitivity. Subjects with higher visual sensitivity might
# have smaller biases because they have to rely less on priors and more on sensory evidence
PlotThVsBias <- function(inputData,          # Data of type glmData
                         modelWeights,       # Weights that will be plotted. Default all weights will be plotted. Subjects will always be plotted.
                         conditionsToPlot=1)   # Which conditions to be plotted
{
  # Select only bias weights
  ## Plot proportion correct of responding right to stimuli presented either to left or to right
  pcRight <- droplevels(subset(inputData, Condition %in% conditionsToPlot))
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
  print("Fitting Probit function to data")
  models <- dlply(pcRightSummary, c("SubjectID", "SessionID", "Condition"), .fun=FitProbit)
  predvals <- ldply(models, .fun=PredictvalsProbit, xvar="Contrast", yvar="pRightCorrect", type="response")
  ## Compute mean
  meanPredvals <- ddply(predvals, .(SubjectID, Condition), summarise, slope=mean(slope), th75=mean(Th75), th50=mean(Th50))

  # Select only bias weights from the set of all weights
  biasWeights <- droplevels(subset(modelWeights, Parameter %in% c("PrevFail1", "PrevCorr1") & Regularized=="yes"))
  biasWeights <- droplevels(subset(biasWeights, Condition %in% conditionsToPlot))
  meanBiasWeights <- ddply(biasWeights, .(SubjectID, Condition, Regularized, Parameter), summarise, Weight=mean(Weight)) # Compute mean subject weight
  meanBiasWeightsWide <- cast(meanBiasWeights, SubjectID+Condition~Parameter, value=.(Weight))

  # Get reaction time for each subject
  meanRTs <- ddply(pcRight, .(SubjectID), summarise, MeanRT=mean(RT))

  # Merge biases and thresholds into one data frame. Slopes will also be here.
  thAndBiases <- merge(meanBiasWeightsWide, meanPredvals)
  thAndBiases <- merge(thAndBiases, meanRTs)
  thAndBiases <- ddply(thAndBiases, .(SubjectID), transform,
                        HistoryBias=abs(PrevCorr1) + abs(PrevFail1))

  browser()
  #   ## Get percent change of slope relative to slope when error weight is 0
  #   predvals1 <- ddply(predvals1, .(SubjectID), transform, SlopeChange = 1- (slope[which(PrevFailWeight==0)] / slope))
  #
  ## Plot change of slope as a function of failure bias
  #   if (!missing(weightsToPlot)) {
  #     datForPlot <- subset(predvals1, (PrevFailWeight %in% weightsToPlot) | (IsSubject==T))
  #   } else datForPlot <- predvals1

  ggplot(data=thAndBiases, aes(x=HistoryBias, y=MeanRT)) +
    theme_few() +
    #scale_y_log10() +
    geom_point(size=3) +
    geom_text(aes(label=SubjectID), size=3, color='grey50', vjust=2.5) +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~Condition, scales="free")
ddply(thAndBiases, .(Condition), summarise, cor(HistoryBias, th75, method='pearson'))
ddply(thAndBiases, .(Condition), summarise, cor(PrevCorr1, th75, method='pearson'))

# Partial correlation between RT and bias while controlling for th
cor.test(thAndBiases$HistoryBias, thAndBiases$MeanRT)


ggplot(data=thAndBiases, aes(x=MeanRT, y=th75)) +
  theme_few() +
  #scale_y_log10() +
  geom_point(size=3) +
  geom_text(aes(label=SubjectID), size=3, color='grey50', vjust=2.5) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~Condition, scales="free")

  ggplot(data=thAndBiases, aes(x=PrevCorr1, y=PrevFail1)) +
    theme_few() +
    xlim(-1.3,1.3) + ylim(-2.0,.5) +
    geom_hline(vintercept=0, color="grey70", linetype='dashed') +
    geom_vline(hintercept=0, color="grey70", linetype='dashed') +
    #scale_y_log10() +
    geom_point(aes(size=th75)) +
    #geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~Condition, scales="free")

  # 3D plot
#   library(Rcmdr)
#   scatter3d(data=thAndBiases, th75~PrevFail1 + PrevCorr1)
#   scatter3d(data=thAndBiases, slope~PrevFail1 + PrevCorr1)
#

#   errbarDat <- subset(datForPlot, PrevFailWeight==0 | IsSubject==T)
#   ggplot(data=datForPlot, aes(x=PrevFailWeight, y=slope)) +
#     theme_publish2() +
#     xlab("Failure bias (weight)") + ylab("Slope") +
#     coord_cartesian(ylim=c(0.468, 0.52), xlim=c(-2.5, 2.5)) +
#     geom_line(subset=.(!IsSubject), stat="summary", fun.y=mean, size=1, color="grey80") +
#     geom_point(subset=.(!IsSubject), stat="summary", fun.y=mean, size=8, color="white") +
#     geom_point(subset=.(IsSubject), aes(x=mean(PrevFailWeight), y=slope), stat="summary", fun.y=mean, size=8, color="white") +
#     geom_point(subset=.(!IsSubject), stat="summary", fun.y=mean, size=5, shape=1) +
#     geom_point(subset=.(IsSubject), aes(x=mean(PrevFailWeight), y=slope), stat="summary", fun.y=mean, size=5)
}


#' Get threshold, slope and lapse rate from glmData
#'
#' Works better with subjects' real data rather than simulated trials
#' For simulated trials,use ThAndSlopeForSimData
#'
#' @param inputData data structured as glmData
#' @param returnMeans if TRUE, computes average for each subject. Otherwise, returns by-run results
#' @param whichConditions conditions to analyse. For natural bias, always include conditions 1 and 13
#'
#' @return data frame containing threshold, slope and lapse rate of each subject. If requested, can return by-run results, otherwise means of each subject
#'
#' @examples
#' ThSlopeAndLapse(glmData, returnMeans=TRUE, whichCondition=c(1,13))
#'
#' @export
ThSlopeAndLapse <- function(inputData,                  # Data of type glmData
                            returnMeans=TRUE,           # When FALSE, returns per-run data, otherwise, averages across runs per subject
                            whichConditions=c(1,13))    # condition or conditions to process
{
  # Select only bias weights
  ## Plot proportion correct of responding right to stimuli presented either to left or to right
  pcRight <- droplevels(subset(inputData, Condition %in% whichConditions))
  #pcRight <- droplevels(subset(inputData, Condition %in% conditionsToPlot))
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
  # Get lapse rate
  lapses <- ddply(pcRightSummary, .(SubjectID, SessionID, Condition), .fun=LapseRateFromHighestContrast)
  ## Convert contrast into %
  pcRightSummary$Contrast <- pcRightSummary$Contrast * 100
  print("Fitting Probit function to data")
  models <- dlply(pcRightSummary, c("SubjectID", "SessionID", "Condition"), .fun=chb:::FitProbit)
  predvals <- ldply(models, .fun=chb:::PredictvalsProbit, xvar="Contrast", yvar="pRightCorrect", type="response")

  # Return either means, or per run threshold and slope
  if (returnMeans) {
    cols <- c("SubjectID", "Condition")
  } else {
    cols <- c("SubjectID", "Condition", "SessionID")
  }
  thAndSlope <- ddply(predvals, cols, summarise, slope=mean(slope), th75=mean(Th75), th50=mean(Th50))
  lapses <- ddply(lapses, cols, summarise,  LapseRate=mean(LapseRate))
  # Combine threshold, slope and lapses
  thSlopeAndLapse <- merge(thAndSlope, lapses)
  return(thSlopeAndLapse)
}


#' Compute lapse rate from highest contrast intensity
#'
#' Lapses occur due to inattention and are not related to decision making processes
#' in the brain. Here, the lapse rate is estimated by finding the strongest contrast intensity
#'  - when missing stimulus is likely to be due to inattention - and computing the error rate
#'  at that stimulus intensity
#'
#' This lapse rate estimation approach is based on common sense, experience and supported by:
#' Prins, N. (2012). The psychometric function: the lapse rate revisited. Journal of Vision, 12(6). http://doi.org/10.1167/12.6.25
#'
#' NOTE: If your data does not strong contrast intensities, using this method of lapse rate estimation is not recommended.
#' In such a case, it is better to use
#'
#' @param rightwardResponses data in the format of proportion of rightward responses (also called pcRightSummary) for one subject. It can either be average across many runs or one run. The function doesn't care about it.
#'
#' @return lapse rate as proportion of errors at the highest contrast intensity
LapseRateFromHighestContrast <- function(rightwardResponses) # proportion rightward responses (also called pcRightSummary in other places)
{
  # Find highest contrast intensity
  maxContrast <- max(abs(rightwardResponses$Contrast))
  stimOnLeft <- with(rightwardResponses, rightwardResponses[Contrast == -maxContrast,]) # Highest contrast stim presented on the left
  stimOnRight <- with(rightwardResponses, rightwardResponses[Contrast == maxContrast,]) # Highest contrast stim presented on the right
  # Lapse rate as average error rate in response to the strongest stimulus
  lapseRate <- mean(c(stimOnLeft$pRightCorrect, 1-stimOnRight$pRightCorrect))
  return(data.frame(LapseRate=lapseRate))
}


#' Sort parameters of the model in intuitive way
#'
#' Model parameter names can get scrambled and when plotting them
#' it can be difficult to understand the graph. This function
#' convniently groups model parameters by name such that contrast
#' intensities are followed by bias weights.
#'
#' @param paramNames model parameter names as factor (see usage)
#'
#' @examples
#' modelWeights$Parameters <- sortParams(modelWeights$Parameters)
SortParams <- function(paramNames) {
  ## Sort contrast and history weights in the order that will make it intuitive to read the plot
  paramNamesSorted <- levels(paramNames)
  contrastNames <- sort(paramNamesSorted[grep("c0",paramNamesSorted)])
  biasNames <- paramNamesSorted[!paramNamesSorted %in% contrastNames]
  paramNamesSorted <- unique(c(biasNames, contrastNames))
  return(paramNamesSorted)
}

