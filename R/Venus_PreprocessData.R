#' Function to load data. It needs a vector containing 
#' full names of files. 
#'
#' LoadData(c('file1', 'file', 'file3'))
#'
#' @export

###########################################
## Read data from all subject and all sessions into memory
LoadData <- function(fileList) {
  ## Set the working directory
  #setwd('~/Desktop/Dropbox/Work/RIKEN_BSI_2012/Venus/')
  ## List of files to load

  ## Read the content of the first file
  loadedData <- read.csv(fileList[1], header = TRUE)
  ## Then the rest of files
  for (fileName in fileList[2:length(fileList)]) {
    tmp <- read.csv(fileName, header = TRUE)
    loadedData <- rbind(loadedData, tmp)
  }
  return(loadedData)
}

#' Filter out unnecessary or pilot data
#' Then assign "conditions"
#' Condition 1: No bias condition. fail/switch and succeed/stay probabilities are set to 50% 
#' Condition 2: fail/switch probability set at 80%, succeed/stay at 50%
#' Condition 3: fail/switch probability set at 50%, succeed/stay at 80%
#' Function to load data. It only loads what's here
#'
#' FilterRawData()
#'

#' @export
FilterRawData <- function(rawData)   	# Input data set that will be filtered
{	
  # ## Exclude these session because wanted to concentrate on certain range of contrast values
  rawData <- subset(rawData, !((rawData$SessionID == 2) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 3) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 4) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 5) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 19) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 20) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 17) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 18) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 19) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 20) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 21) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 22) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 23) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 24) & (rawData$SubjectID == 's001')))
  # ## These three sessions were collected using different kind of stimulus - gratings with smooth edges (28 & 29) and shorter presentation time (30)
  rawData <- subset(rawData, !((rawData$SessionID == 28) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 29) & (rawData$SubjectID == 's001')))
  rawData <- subset(rawData, !((rawData$SessionID == 30) & (rawData$SubjectID == 's001')))
  ## For subject s002, 4% contrast was too high as even using 3% responses saturated. Exclude 4%
  rawData <- subset(rawData, !((rawData$Contrast == 0.04) & (rawData$SubjectID=='s002')))
  
  # ## Select a particular subject if needed
  # #rawData <- subset(rawData, rawData$SubjectID == "s003")
  # ## Exclude very low contrast intensities for now
  # rawData <- subset(rawData, rawData$Contrast != 0.0005)
  # rawData <- subset(rawData, rawData$Contrast != 0.001)
  # rawData <- subset(rawData, rawData$Contrast != 0.002)
  # #rawData <- subset(rawData, rawData$Contrast != 0.13)
  # #rawData <- subset(rawData, rawData$Contrast != 0.003)
  # #rawData <- subset(rawData, rawData$Contrast != 0.005)
  # #rawData <- subset(rawData, rawData$Contrast != 0.008)
  # ## Exclude first session for s001 because it was his training sesssion 
  # rawData <- subset(rawData, !((rawData$SessionID == 1) & (rawData$SubjectID == 's001')))
  
  # ## Uncommenting will collapse all sessions into one
  # #rawData$SessionID <- 1
  
  return(rawData)
}

######################### --- ClassifyRawData --- #############################
#' Assign "Conditions"
#' Condition 1: No bias condition. fail/switch and succeed/stay probabilities are set to 50% 
#' Condition 2: fail/switch probability set at 80%, succeed/stay at 50%. !!!IMPORTANT. This is to induce fail/stay bias, because when we switch after fail, we end up on same side where fail was. 
#' Condition 3: fail/switch probability set at 50%, succeed/stay at 80%. This is to induce success/stay bias.
#' Condition 7: succeed/stay probability set to 50%, fail/switch is set to 20%. This is to induce fail/switch bias, because after fail stimulus will likely to remain on the opposite side of fail, inviting the subbject to switch side. 
#' Condition 8: succeed/stay probability=20%, fail/switch=50%. Lots of switching after success to induce success/switch bias. 
#' Condition 9: This is same as Condition 2, but only low contrasts used to have similar number of failure and success. Succeed/stay probability=50%, fail/switch=80%.
#' Condition 10: This is same as Condition 3, but only low contrasts used to have similar number of failure and success. Succeed/stay probability=80%, fail/switch=50%.
#' Condition 11: This is same as Condition 7, but only low contrasts used to have similar number of failure and success. Succeed/stay probability=50%, fail/switch=20%.
#' Condition 12: This is same as Condition 8, but only low contrasts used to have similar number of failure and success. Succeed/stay probability=20%, fail/switch=50%.
#' Condition 4: Stimulus diameter is 6??, eccentricity is 8??
#' Condition 5: Stimulus diameter is 12??, eccentricity is 12??
#' Condition 6: Stimulus diameter is 6??, eccentricity is 10??
#' @export
ClassifyRawData <- function(rawData)   	# Input data set that will be filtered
{	
  rawData$Condition <- -1 ## -1 means no condition was assigned
  #### Subject 1
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(43, 44)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(45, 55, 56)))] <- 3
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(48, 52, 54)))] <- 4
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(49, 51)))] <- 5
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(50, 53)))] <- 6
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(59, 60, 61)))] <- 7
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(57, 58)))] <- 8
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(62, 63, 64)))] <- 9
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(65, 66, 67)))] <- 10
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(68)))] <- 11
  rawData$Condition[with(rawData, (SubjectID=='s001') & (OriginalSessionID %in% c(69, 70)))] <- 12
  
  #### Subject 2	
  rawData$Condition[with(rawData, (SubjectID=='s002') & (OriginalSessionID %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10)))] <- 1
  #### Subject 3
  rawData$Condition[with(rawData, (SubjectID=='s003') & (OriginalSessionID %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)))] <- 1
  #### Subject 4	
  rawData$Condition[with(rawData, (SubjectID=='s004') & (OriginalSessionID %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)))] <- 1
  #### Subject 5	
  rawData$Condition[with(rawData, (SubjectID=='s005') & (OriginalSessionID %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s005') & (OriginalSessionID %in% c(15,16,17)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s005') & (OriginalSessionID %in% c(18,19,20,21,22,23)))] <- 3
  rawData$Condition[with(rawData, (SubjectID=='s005') & (OriginalSessionID %in% c(24,25,26)))] <- 7
  rawData$Condition[with(rawData, (SubjectID=='s005') & (OriginalSessionID %in% c(27,28,29,30,31,32)))] <- 8
  #### Subject 6	
  rawData$Condition[with(rawData, (SubjectID=='s006') & (OriginalSessionID %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)))] <- 1
  #### Subject 7	
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(16,17,18)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(13,14,15, 19, 20, 21)))] <- 3
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(22)))] <- 5
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(23)))] <- 4
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(27, 28, 29)))] <- 7
  rawData$Condition[with(rawData, (SubjectID=='s007') & (OriginalSessionID %in% c(24, 25, 26)))] <- 8
  #### Subject 8	
  rawData$Condition[with(rawData, (SubjectID=='s008') & (OriginalSessionID %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)))] <- 1
  #### Subject 9	
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(14,15,16)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(11,12,13)))] <- 3
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(17)))] <- 4
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(19)))] <- 5
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(18)))] <- 6
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(23, 24, 25)))] <- 7
  rawData$Condition[with(rawData, (SubjectID=='s009') & (OriginalSessionID %in% c(20, 21, 22)))] <- 8
  #### Subject 10
  rawData$Condition[with(rawData, (SubjectID=='s010') & (OriginalSessionID %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s010') & (OriginalSessionID %in% c(10,11,12)))] <- 3
  rawData$Condition[with(rawData, (SubjectID=='s010') & (OriginalSessionID %in% c(13)))] <- 5
  rawData$Condition[with(rawData, (SubjectID=='s010') & (OriginalSessionID %in% c(14)))] <- 4
  rawData$Condition[with(rawData, (SubjectID=='s010') & (OriginalSessionID %in% c(18,19,20)))] <- 7	
  rawData$Condition[with(rawData, (SubjectID=='s010') & (OriginalSessionID %in% c(15,16,17)))] <- 8	
  #### Subject 11
  rawData$Condition[with(rawData, (SubjectID=='s011') & (OriginalSessionID %in% c(1, 2, 3, 4, 5)))] <- 1
  #### Subject 12
  rawData$Condition[with(rawData, (SubjectID=='s012') & (OriginalSessionID %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s012') & (OriginalSessionID %in% c(10,11,12)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s012') & (OriginalSessionID %in% c(19,20,21)))] <- 3
  rawData$Condition[with(rawData, (SubjectID=='s012') & (OriginalSessionID %in% c(16, 17, 18)))] <- 7
  rawData$Condition[with(rawData, (SubjectID=='s012') & (OriginalSessionID %in% c(13, 14, 15)))] <- 8
  #### Subject 13
  rawData$Condition[with(rawData, (SubjectID=='s013') & (OriginalSessionID %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9)))] <- 1
  #### Subject 14
  rawData$Condition[with(rawData, (SubjectID=='s014') & (OriginalSessionID %in% c(1, 2, 3, 4, 5, 6)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s014') & (OriginalSessionID %in% c(10, 11, 12)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s014') & (OriginalSessionID %in% c(7, 8, 9)))] <- 3
  #### Subject 15
  rawData$Condition[with(rawData, (SubjectID=='s015') & (OriginalSessionID %in% c(1, 2, 3, 7, 8, 9)))] <- 1
  rawData$Condition[with(rawData, (SubjectID=='s015') & (OriginalSessionID %in% c(10, 11, 12)))] <- 2
  rawData$Condition[with(rawData, (SubjectID=='s015') & (OriginalSessionID %in% c(4, 5, 6)))] <- 3
  
  #### Subject 16 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s016') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 17 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s017') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 18 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s018') & (OriginalSessionID %in% c(4, 5, 3)))] <- 13
  #### Subject 19 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s019') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 20 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s020') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 21 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s021') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 22 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s022') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 23 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s023') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 24 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s024') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 25 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s025') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 26 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s026') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  #### Subject 27 (UCL)
  rawData$Condition[with(rawData, (SubjectID=='s027') & (OriginalSessionID %in% c(1, 2, 3)))] <- 13
  
  return(rawData)
}

#' Assign condition to sessions  
#' 
#' Important: SessionID should be provided based on OriginalSessionID
#' @export
AssignCondition <- function(dat,     # Data in rawData or glmData format
                            SubjectID,   # Subject ID
                            SessionID,   # Run numbers based on OriginalSessionID (not SessionID, which changes at some point)
                            Condition) {   # Conditions to assign to those run numbers
# Find the subject, get those sessions and assign the condition
  dat$Condition[(dat$SubjectID == SubjectID) & 
                           (dat$OriginalSessionID %in% SessionID)] <- Condition
    return(dat)
}

