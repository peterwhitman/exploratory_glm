library(stats)
library(car)

ExploratoryGLM <- function(data, response, family, linkArg, offsetVariable, offsetFormula, numCombinations, vifThresh)
{
  
  #set up function to deal with missing arguments
  `%||%` <- function(lhs, rhs) {
    if (!missing(lhs)) {
      lhs
    } else {
      rhs
    }
  }
  
  #use missing argument function to set default parameters
  family$link <- linkArg %||% family$link
  removeList <- response
  removeList[2] <- offsetVariable %||% response
  offsetFormula <- offsetFormula %||% ""
  if(nchar(offsetFormula) > 0){offsetFormula <- paste(offsetFormula,"+")}
  if(ncol(data)-length(unique(removeList)) > 5){altNum <- 5} else{altNum <- ncol(data)-length(unique(removeList))}
  numCombinations <- numCombinations %||% altNum
  vifThresh <- vifThresh %||% 2
  
  #create null model and combinations of explanatory variables
  selectedData <- data[, !names(data) %in% removeList]
  myVariables <- names(selectedData)
  null <- paste(response, " ~ ", offsetFormula, "1") #Specify null model
  out <- unlist(lapply(1:numCombinations, function(n) combn(myVariables, n, 
                                            FUN=function(row) paste0(substr(null, 1, nchar(null)-1), 
                                                                     paste0(row, collapse = "+")))))

  null[2:(length(out)+1)] <- out
  
  #Create models for all combinations of explanatory variables
  myFormulas <- as.list(null)
  myModels <- list()
  for(i in 1:length(myFormulas))
  {
    myModels[[i]] <- glm(as.formula(myFormulas[[i]]), data = data, family = family)
  }
  names(myModels) <- null[1:length(myModels)] #name outputed models with its formula

  #Remove models with multicollinearity
  RefinedList <- myModels[1:(length(myVariables)+1)]
  NameList <- myFormulas[1:(length(myVariables)+1)]
  for(i in (length(myVariables)+2):length(myModels))
  {
    t <- vif(myModels[[i]])
    if(max(t) <= vifThresh) {
      RefinedList <- append(RefinedList, myModels[i])
      NameList <- append(NameList, myFormulas[i])
    }
  }
  names(RefinedList) <- NameList

  #Compute some statistics for each model
  dev <- vector()
  df_resid <- vector()
  aic <- vector()
  loglike <- vector()
  for(i in 1:length(RefinedList))
  {
    dev[i] <- deviance(RefinedList[[i]])
    df_resid[i] <- ncol(RefinedList[[i]][["model"]])-length(unique(removeList))
    aic[i] <- AIC(RefinedList[[i]])
    loglike[i] <- logLik(RefinedList[[i]])
  }

  #Combine statistics into a table with model formula
  ModelComparison <- data.frame(cbind(NameList[1:length(NameList)], aic, loglike, dev, df_resid))
  ModelComparison$aic <- as.numeric(as.character(ModelComparison$aic))
  ModelComparison$df_resid <- as.numeric(as.character(ModelComparison$df_resid))
  ModelComparison$ModelID <- row.names(ModelComparison)

  #Select models with smallest AIC at each unique number of degrees of freedom
  SmallestAIC <- data.frame(tapply(ModelComparison$aic, ModelComparison$df_resid, min))
  names(SmallestAIC)[1] <- "aic"
  SmallestAIC$df_resid <- as.numeric(row.names(SmallestAIC))
  SmallestAIC <- merge(SmallestAIC, ModelComparison, by = c("df_resid","aic"), all = FALSE) #extract remaining stats/info for the selected models from original table
  SmallestAIC <- SmallestAIC[,c(6,3,1:2,4:5)] #reorder table
  names(SmallestAIC)[2] <- "Formula"
  
  #output list of filtered models and table of the "best" models to global environment
  RefinedModelList <<- RefinedList
  SelectedModels <<- SmallestAIC
}

library(ResourceSelection)
df <- goats
ExploratoryGLM(df[,c(1,3:8)], response = names(df)[1], family = binomial())
