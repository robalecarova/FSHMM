# Feature Selection for multiple condition gene expression 
# Script with a Hidden Markov Model (R-package RcppHMM)


#  Dataset documentation
#' Gene expression of a two condition time series experiment
#'
#' A dataset containing the gene expression of a time series experiment
#' with 6 time points (0,2,6,12,18,24), 2 conditions and 3 biological replicates
#'
#'
#' @format A data frame with 12762 rows and 36 samples
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75417}
"data.ann"
#> [1] "data.ann"


# Constants to use
LOG2PI <- log(2.0 * pi) 

#' @title  Probability density function for multivariate normal distribution
#'
#' @description 
#' This function computes the probability density function for a multivariate normal distribution
#' @param x Observation vector
#' @param mean Multivariate normal distribution mean vector
#' @param sigma Multivariate normal distribution variance/covariance matrix
#' @param log The probability value returns as log-value. Defaults to TRUE.
#' @return The probability value of 'x'
#' @keywords probabilityDensityFunction
#' @export
#' @examples
#' x <- c(1.5, 4.5)
#' mean <- c(1,4)
#' sigma <- matrix(c(1,0.7,0.7,1), byrow=TRUE, ncol = 2)
#' dmvn(x, mean , sigma)
dmvn <- function(x, mean, sigma, log = TRUE)
{
  xdim <- length(x)
  out <- 0.0
  
  rooti <- chol(sigma)
  rooti[lower.tri(rooti)] <- 0.0
  rooti <- t(solve(rooti))
  
  rootisum <- sum(log(diag(rooti)))
  constants  <- -((xdim)/2.0) * LOG2PI
  
  temp <- rooti %*% (x - mean) ;
  out <- constants - 0.5 * sum(temp * temp) + rootisum;
  
  if (log == F) {
    out = exp(out);
  }
  return(out);
}

#' @title One condition feature selection hidden Markov model
#'
#' @description 
#' This function computes a feature selection filtering via a hidden Markov model
#' 
#' @details 
#' This function uses a hidden Markov model with an observation vector 
#' that can be modeled with a multivariate normal distribution for one condition
#' to select the best features of the observation matrix. 
#' 
#' @param data          numeric matrix - Expression data (Log2 data)
#' @param times             numeric vector - All the time points when the data was measured
#' @param reps              numeric - How many biological replicates do we have per time point?
#' @param LOG               logical - Do you wanna use the data in its Log2 transformation or set them in its natural range? Defaults to TRUE.
#' @param t.size            numeric - What percentage of data will be used to get estimate the parameters? Range - (0,1]. Defaults to .
#' @param iters             numeric - How many iterations will be used to estimate the model parameters? Defaults to 50.
#' @param USE.REP           logical - Want to use all the replicates as input? Defaults to TRUE.
#' @param func              function - If USE.REP == TRUE, then how are they going to be summarized? (mean, median, customFunction). Defaults to median.
#' @param factor            numeric - Factor to multiply the data. Usually used recommended to be used with the LOG = T. Defaults to 10.
#' @param prior             logical - Do you want to use a prior for the model parameters? Defaults to FALSE.
#' @param initial.vector    vector - stochastic vector with the model initial probability vector prior. Defaults to NULL.
#' @param transition.matrix matrix - stochastic matrix with the model transition probability matrix prior. Defaults to NULL.
#' 
#' @return final.results a list consisting of: 
#' \item{data}{matrix - relevant and filtered data} 
#' \item{names}{vector - relevant features names} 
#' \item{model}{list - HMM used to select the most relevant data} 
#' \item{scoreReplicates}{vector - The score of biological replicates per feature} 
#' \item{scoreChanges}{vector - the score of the feature given the number of times the it traverse a 'change state'} 
#' \item{scoreDelta}{vector - the score given the proportion of the changes in time per feature} 
#' \item{scoreMaxDelta}{ - vector - the score of the feature given its maximum change in time } 
#' \item{ranking}{vector - the rank of the feature given the replicate, changes and delta scores} 
#' \item{decodedStates}{matrix - each row has the decoded states of each change per feature }
#' 
#' @seealso \code{\link{featureSelectionHMM}} for multiple condition feature selection.
#' @keywords FeatureSelection
#' @export
#' @examples
#' 
#' 
#' data("ratTimeSeries")
#' 
#' # All the necessary variables to run the script
#' times <- c(0,2,6,12,18,24)  # study time-points vector
#' conds <- c("CTRL","IKAROS") # Condition vector
#' reps <- 3                   # number of biological replicates
#' \dontrun{
#' # It has some random sampling inside, so a set.seed is recommended
#' # Case: Biological replicates summarized by median
#' set.seed(100)
#' summ.results <- featureSelectionOneConditionHMM(data.ann[,19:36],      # Gene expression data matrix
#'                                                 times,         # time-point vector
#'                                                reps,          # Number of replicates
#'                                                USE.REP  = F,  # Do not use all replicate data
#'                                                func = median  # Summarize the data
#' )
#'
#' # Case: Biological replicates not summarized
#' set.seed(100) 
#' mv.results <- featureSelectionHMM(data.ann[,19:36],      # Gene expression data matrix
#'                                  times,         # time-point vector
#'                                  reps           # Number of replicates
#' )  
#'
#' plotMultivariate(mv.results, times, conds[2], reps, ".", "Fat.Rat.multi.all.pdf")
#' plotUnivariate(summ.results, times, conds[2], reps, ".", "Fat.Rat.uni.res.pdf")
#'}
#' 
featureSelectionOneConditionHMM <- function(data,                    
                                            times,                   
                                            reps,                    
                                            LOG = TRUE,                 
                                            t.size = 1,              
                                            iters = 50,              
                                            USE.REP = T,             
                                            func = median,           
                                            factor = 10,             
                                            prior = FALSE,               
                                            initial.vector = NULL,   
                                            transition.matrix = NULL 
) 
{
  # Install and load necessary libraries
  requireNamespace("RcppHMM")
  
  size <- length(times)*reps
  noTrans <- nrow(data)
  conds.data <- data
  
  ###########################################################################
  ###########################################################################
  ###########################################################################
  # Data preprocessing
  
  print(" Starting Data preprocessing")
  
  # Data reordering
  # Set the data from 2D matrix to 3D matrix
  
  data.ord <- array(0, dim = c(reps, length(times), noTrans) )
  # Each row represents a replicate. 
  # If the data is going to be summarized, then we only need one row.
  data.ord.sort <- array(0, dim = c( ifelse(USE.REP, reps, 1) , 
                                     length(times), 
                                     noTrans) )
  
  for(i in 1:noTrans)
  {
    data.ord[,,i] <-  matrix(as.numeric(conds.data[i,]), nrow = reps)
    # Once the data is in a 3D matrix. We need other preprocessing step
    for(j in 1:length(times))
    {
      if(USE.REP)
        # we can summarize the biological replicates 
        data.ord.sort[,j,i] <- sort(data.ord[,j,i])
      else
        # Or sort them to have less variablity
        data.ord.sort[,j,i] <- func(data.ord[,j,i])
    }
  }
  
  
  # Then, we get the value of the gene expression difference between 2 consecutive time points
  # This is used to remove the offset value of each gene, and make them comparable 
  data.deltas.ord <- array(0, dim = c(ifelse(USE.REP, reps, 1) ,
                                      length(times)-1,
                                      noTrans) )
  
  # Are they going to be analized as LOG2 values?
  if( LOG == T)
  {
    for(i in 1:noTrans)
    {
      for(j in 2:length(times))
      {
        data.deltas.ord[,j-1,i] <-  data.ord.sort[,j,i] - data.ord.sort[,j-1,i]
      }
    }
  } else
  {
    for(i in 1:noTrans)
    {
      for(j in 2:length(times))
      {
        data.deltas.ord[,j-1,i] <-  2^data.ord.sort[,j,i] - 2^data.ord.sort[,j-1,i]
      }
    }
  }
  
  print(" Finished Data preprocessing")
  
  # We select how many genes are going to be used as the training set
  # All the conditions are going to be compared versus the Control set
  # With this, we can compare multiple conditions with only one model
  m.size <- floor(noTrans * t.size)
  indx <- sample(1:noTrans, m.size)
  matrix.training <- array(0, dim = c(ifelse(USE.REP, reps, 1) ,
                                      length(times) - 1,
                                      m.size) )
  
  
  temp.data.deltas.ord <- data.deltas.ord
  for(i in 1:m.size)
  {
    matrix.training[,,i] <-  temp.data.deltas.ord[,,indx[i]]
  }  
  
  factor <- ifelse(LOG == T, factor, 1)
  matrix <- matrix.training * factor
  M <- nrow(matrix)
  
  i <- 2
  
  print("Estimating Model parameters")
  # Estimate model parameters
  model <- initGHMM(i, M)
  
  # If the user inputs a prior, then it is necessary to replace the paramaters
  if(prior)
  {
    if(!is.null(initial.vector))
    {
      model$Pi <- initial.vector
    }
    
    if(!is.null(transition.matrix))
    {
      model$A <- transition.matrix
    }
    
    model <- verifyModel(model)
  }
  
  
  model <- learnEM(model, matrix, 
                   iter = iters, delta = 1e-09, print = T) 
  
  # Select the state that represents the No Change Sate (N) - Low variance
  mu.1 <- mean(model$Sigma[,,1])
  mu.2 <- mean(model$Sigma[,,2])
  
  if(mu.1 < mu.2){
    names <- c("N","C")
  }else{
    names <- c("C","N")
  }
  
  model <- setNames(model, list( 'StateNames' = names))
  
  print("Decoding hidden states")
  # Decode hidden states
  decoded.states.collapsed.cond <- matrix("N", ncol = (length(times)-1), nrow = dim(data.deltas.ord)[3])
  for(i in 1:noTrans)
  {
    temp <- viterbi(model, matrix(data.deltas.ord[,,i] * factor, nrow= ifelse(USE.REP, reps, 1) )  )
    # decoded.states.collapsed.cond <- rbind(decoded.states.collapsed.cond, temp)
    decoded.states.collapsed.cond[i, ] <- temp
  }
  
  cond.decoded <- apply(decoded.states.collapsed.cond,
                        1,
                        function(x){paste(x, collapse="")})
  
  # decoded.states.collapsed[[z]] <- cond.decoded
  # decoded.states[[z]] <- decoded.states.collapsed.cond
  command <- paste(rep("N", length(times) - 1), collapse = "")
  
  # All the Control genes that are not chaning in all the experiment
  # And all the Condition genes that change even one in all the experiment
  important.indx <- cond.decoded != command
  important.features <- data[important.indx,]
  
  noTransImp <- sum(important.indx)
  ########################################
  # Let's put a score and then rank them
  
  # The Score will be based in the Porbability of observing each gene sequence
  # We get a score per condition, and then we summarize them
  
  print("Score and ranking...")
  
  score.states.cond <- vector(mode = "numeric", length = noTransImp )
  delta.score.cond <- vector(mode = "numeric", length = noTransImp )
  max.delta.score.cond <- vector(mode = "numeric", length = noTransImp )
  decoded.important <- matrix("N", ncol = (length(times)-1), noTransImp)
  counter <- 1
  
  # For each gene
  for(i in 1:noTrans)
  {
    # Each genes that is "important"
    if(important.indx[i] == TRUE)
    {
      x <- matrix(data.deltas.ord[,,i] * factor, nrow= ifelse(USE.REP, reps, 1))
      hiddenPath <- decoded.states.collapsed.cond[i,]
      hiddenPath <- match(hiddenPath, model$StateNames)
      score <- 0
      for(t in 1:length(hiddenPath))
      {
        if(USE.REP)
        {
          score <- score + dmvn(x[,t], 
                                mean = model$Mu[,hiddenPath[t]],
                                sigma = model$Sigma[,,hiddenPath[t]],
                                log = T )
        }else{
          score <- score + dnorm(x[t], 
                                 mean = model$Mu[hiddenPath[t]],
                                 sd = model$Sigma[,,hiddenPath[t]],
                                 log = T )
        }
      }
      score.states.cond[counter] <-  score
      decoded.important[counter, ] <-  decoded.states.collapsed.cond[i,]
      
      temp <- abs(data.deltas.ord[,,i])
      delta.score.cond[counter] <- sum(temp)
      if(USE.REP == TRUE)
      {
        max.delta.score.cond[counter] <- max(apply(temp, 2, mean))  
      }else
      {
        max.delta.score.cond[counter] <- max(temp)  
      }
      
      counter <- counter + 1
    }
  }

  all.scores.changes <- apply(decoded.important, 1, function(x){
    sum(x == "C")
  })
  
  ranking <- order(all.scores.changes, delta.score.cond, score.states.cond, decreasing = T)
  
  #################################
  
  # Return all the important variables as a list
  print("Finished")
  final.results <- list()
  final.results[["data"]]    <- important.features
  final.results[["names"]]   <- rownames(important.features)
  final.results[["model"]]   <- model
  final.results[["scoreReplicates"]]   <- score.states.cond
  final.results[["scoreChanges"]]   <- all.scores.changes
  final.results[["ranking"]] <- ranking
  final.results[["decodedStates"]] <- decoded.important
  final.results[["scoreDelta"]] <- delta.score.cond
  final.results[["scoreMaxDelta"]] <- max.delta.score.cond
  
  return(final.results)
}

#' @title Multiple condition feature selection hidden Markov model
#'
#' @description 
#' This function computes a feature selection filtering via a hidden Markov model
#' 
#' @details 
#' This function uses a hidden Markov model with an observation vector 
#' that can be modeled with a multivariate normal distribution for multiples condition.
#' It compares the studied conditions versus a control condition. And it 
#' select the best features of the observation matrix. 
#' 
#' 
#' @param data              numeric matrix - Expression data (Log2 data)
#' @param times             numeric vector - All the time points when the data was measured
#' @param conds             character vector - All the condition names. The first one is the Control/Baseline condition.
#' @param reps              numeric - How many biological replicates do we have per time point?
#' @param LOG               logical - Do you wanna use the data in its Log2 transformation or set them in its natural range? Defaults to TRUE.
#' @param t.size            numeric - What percentage of data will be used to get estimate the parameters? Range - (0,1]. Defaults to .
#' @param iters             numeric - How many iterations will be used to estimate the model parameters? Defaults to 50.
#' @param USE.REP           logical - Want to use all the replicates as input? Defaults to TRUE.
#' @param func              function - If USE.REP == TRUE, then how are they going to be summarized? (mean, median, customFunction). Defaults to median.
#' @param factor            numeric - Factor to multiply the data. Usually used recommended to be used with the LOG = T. Defaults to 10.
#' @param prior             logical - Do you want to use a prior for the model parameters? Defaults to FALSE.
#' @param initial.vector    vector - stochastic vector with the model initial probability vector prior. Defaults to NULL.
#' @param transition.matrix matrix - stochastic matrix with the model transition probability matrix prior. Defaults to NULL.
#' 
#' @return final.results a list consisting of: 
#' \item{data}{matrix - relevant and filtered data} 
#' \item{names}{vector - relevant features names} 
#' \item{model}{list - HMM used to select the most relevant data} 
#' \item{scoreReplicates}{vector - The score of biological replicates per feature} 
#' \item{scoreChanges}{vector - the score of the feature given the number of times the it traverse a 'change state'} 
#' \item{scoreDelta}{vector - the score given the proportion of the changes in time per feature} 
#' \item{scoreMaxDelta}{ - vector - the score of the feature given its maximum change in time } 
#' \item{ranking}{vector - the rank of the feature given the replicate, changes and delta scores} 
#' \item{decodedStates}{list -  each element of the list represents a condition and has a matrix with each row has the decoded states of each change per feature }
#'  
#' @seealso \code{\link{featureSelectionHMM}} for multiple condition feature selection.
#' @keywords FeatureSelection
#' @export
#' @examples
#' 
#' data("ratTimeSeries")
#' 
#' # All the necessary variables to run the script
#' times <- c(0,2,6,12,18,24)  # study time-points vector
#' conds <- c("CTRL","IKAROS") # Condition vector
#' reps <- 3                   # number of biological replicates
#' \dontrun{
#' # It has some random sampling inside, so a set.seed is recommended
#' # Case: Biological replicates summarized by median
#' set.seed(100)
#' summ.results <- featureSelectionHMM(data.ann,      # Gene expression data matrix
#'                                    times,         # time-point vector
#'                                    conds,         # Condition vector
#'                                    reps,          # Number of replicates
#'                                    USE.REP  = F,  # Do not use all replicate data
#'                                    func = median  # Summarize the data 
#' )  
#'
#' # Case: Biological replicates not summarized
#' set.seed(100) 
#' mv.results <- featureSelectionHMM(data.ann,      # Gene expression data matrix
#'                                  times,         # time-point vector
#'                                  conds,         # Condition vector
#'                                  reps           # Number of replicates
#' )  
#'
#'
#' plotMultivariate(mv.results, times, conds[2], reps, ".", "Fat.Rat.multi.all.pdf")
#' plotUnivariate(summ.results, times, conds[2], reps, ".", "Fat.Rat.uni.res.pdf")
#' }
#' 
featureSelectionHMM <- function(data,                
                                times,                   
                                conds,                   
                                reps,                    
                                LOG = TRUE,                 
                                t.size = 1,              
                                iters = 50,              
                                USE.REP = T,             
                                func = median,           
                                factor = 10,             
                                prior = FALSE,               
                                initial.vector = NULL,   
                                transition.matrix = NULL 
) 
{
  # Install and load necessary libraries
  requireNamespace("RcppHMM")
  
  # It is necessary to be in interactive mode to see the plots
  if(interactive()) {
    print("In interactive mode")
  } else {
    print("Not in interactive mode")
  }
  
  size <- length(times)*reps
  noTrans <- nrow(data)
  
  ###########################################################################
  ###########################################################################
  ###########################################################################
  # Data preprocessing
  
  print(" Starting Data preprocessing")
  
  # Data reordering
  # Set the data from 2D matrix to 3D matrix
  conds.data <- list()
  data.ord.sort.conds <- list()
  data.deltas.ord.conds <- list()
  
  for(z in 1:length(conds))
  { 
    ind <- (1:size) + (z-1)*size
    conds.data[[z]] <- as.matrix(data[,ind])
    noTrans <- nrow(conds.data[[z]] )
    
    data.ord <- array(0, dim = c(reps, length(times), noTrans) )
    # Each row represents a replicate. 
    # If the data is going to be summarized, then we only need one row.
    data.ord.sort <- array(0, dim = c( ifelse(USE.REP, reps, 1) , 
                                       length(times), 
                                       noTrans) )
    
    
    for(i in 1:noTrans)
    {
      data.ord[,,i] <-  matrix(conds.data[[z]][i,], nrow = reps)
      # Once the data is in a 3D matrix. We need other preprocessing step
      for(j in 1:length(times))
      {
        if(USE.REP)
          # we can summarize the biological replicates 
          data.ord.sort[,j,i] <- sort(data.ord[,j,i])
        else
          # Or sort them to have less variablity
          data.ord.sort[,j,i] <- func(data.ord[,j,i])
      }
    }
    
    data.ord.sort.conds[[z]] <- data.ord.sort
    
    # Then, we get the value of the gene expression difference between 2 consecutive time points
    # This is used to remove the offset value of each gene, and make them comparable 
    data.deltas.ord <- array(0, dim = c(ifelse(USE.REP, reps, 1) ,
                                        length(times)-1,
                                        noTrans) )
    
    # Are they going to be analized as LOG2 values?
    if( LOG == T)
    {
      for(i in 1:noTrans)
      {
        for(j in 2:length(times))
        {
          data.deltas.ord[,j-1,i] <-  data.ord.sort[,j,i] - data.ord.sort[,j-1,i]
        }
      }
    } else
    {
      for(i in 1:noTrans)
      {
        for(j in 2:length(times))
        {
          data.deltas.ord[,j-1,i] <-  2^data.ord.sort[,j,i] - 2^data.ord.sort[,j-1,i]
        }
      }
    }
    
    data.deltas.ord.conds[[z]] <- data.deltas.ord
  }
  
  print(" Finished Data preprocessing")
  
  # We select how many genes are going to be used as the training set
  # All the conditions are going to be compared versus the Control set
  # With this, we can compare multiple conditions with only one model
  m.size <- floor(noTrans * t.size)
  indx <- sample(1:noTrans, m.size)
  matrix.training <- array(0, dim = c(ifelse(USE.REP, reps, 1) ,
                                      length(times) - 1,
                                      m.size * length(conds)) )
  
  # The first condition must be the Control
  # All the conditions must be taken to train the model. To not be biased
  for(z in 1:length(conds))
  {
    temp.data.deltas.ord <- data.deltas.ord.conds[[z]]
    offset <- (m.size*(z-1))
    for(i in 1:m.size)
    {
      matrix.training[,,i + offset] <-  temp.data.deltas.ord[,,indx[i]]
    }  
  }
  
  
  factor <- ifelse(LOG == T, factor, 1)
  matrix <- matrix.training * factor
  M <- nrow(matrix)
  
  i <- 2
  
  print("Estimating Model parameters")
  # Estimate model parameters
  model <- initGHMM(i, M)
  
  # If the user inputs a prior, then it is necessary to replace the paramaters
  if(prior)
  {
    if(!is.null(initial.vector))
    {
      model$Pi <- initial.vector
    }
    
    if(!is.null(transition.matrix))
    {
      model$A <- transition.matrix
    }
    
    model <- verifyModel(model)
  }
  
  model <- learnEM(model, matrix, 
                   iter = iters, delta = 1e-09, print = T) 
  
  # Select the state that represents the No Change Sate (N) - Low variance
  mu.1 <- mean(model$Sigma[,,1])
  mu.2 <- mean(model$Sigma[,,2])
  
  if(mu.1 < mu.2){
    names <- c("N","C")
  }else{
    names <- c("C","N")
  }
  
  model <- setNames(model, list( 'StateNames' = names))
  
  print("Decoding hidden states")
  # Decode hidden states
  decoded.states <- list()
  decoded.states.collapsed <- list()
  
  for(z in 1:length(conds)){
    data.deltas.ord <- data.deltas.ord.conds[[z]]
    decoded.states.collapsed.cond <- matrix("N", ncol = (length(times)-1), nrow = dim(data.deltas.ord)[3])
    for(i in 1:noTrans)
    {
      temp <- viterbi(model, matrix(data.deltas.ord[,,i] * factor, nrow= ifelse(USE.REP, reps, 1) )  )
      # decoded.states.collapsed.cond <- rbind(decoded.states.collapsed.cond, temp)
      decoded.states.collapsed.cond[i, ] <- temp
    }
    
    cond.decoded <- apply(decoded.states.collapsed.cond,
                          1,
                          function(x){paste(x, collapse="")})
    decoded.states.collapsed[[z]] <- cond.decoded
    decoded.states[[z]] <- decoded.states.collapsed.cond
  }
  
  command <- paste(rep("N", length(times) - 1), collapse = "")
  conds.inds <- matrix(F, nrow = noTrans, ncol = length(conds) )
  
  conds.inds[,1] <- decoded.states.collapsed[[1]] == command
  for(z in 2:length(conds))
  {
    conds.inds[,z] <- decoded.states.collapsed[[z]] != command
  }
  
  # All the Control genes that are not changing in all the experiment
  # And all the Condition genes that change even one in all the experiment
  important.indx <- apply(conds.inds, 1, 
                          function(x){
                            # At least one condition is changing (OR)
                            # And the Control gene is not Changing (AND)
                            Reduce('||', x[-1]) && x[1]
                          })
  
  ########################################
  # Let's put a score and then rank them
  
  # The Score will be based in the Porbability of observing each gene sequence
  # We get a score per condition, and then we summarize them
  
  print("Score and ranking...")
  score.cond <- list()
  decoded.cond <- list()
  delta.cond <- list()
  max.delta.cond <- list()
  
  # For each condition we get a score
  for(z in 1:length(conds)){
    data.deltas.ord <- data.deltas.ord.conds[[z]]
    
    score.states.cond <- vector(mode = "numeric", length = sum(important.indx) )
    delta.score.cond <- vector(mode = "numeric", length = sum(important.indx) )
    max.delta.score.cond <- vector(mode = "numeric", length = sum(important.indx) )
    decoded.important <- matrix("N", ncol = (length(times)-1), nrow = sum(important.indx))
    counter <- 1
    
    # For each gene
    for(i in 1:noTrans)
    {
      # Each genes that is "important"
      if(important.indx[i] == TRUE)
      {
        x <- matrix(data.deltas.ord[,,i] * factor, nrow= ifelse(USE.REP, reps, 1))
        hiddenPath <- decoded.states[[z]][i,]
        hiddenPath <- match(hiddenPath, model$StateNames)
        score <- 0
        for(t in 1:length(hiddenPath))
        {
          if(USE.REP)
          {
            score <- score + dmvn(x[,t], 
                                  mean = model$Mu[,hiddenPath[t]],
                                  sigma = model$Sigma[,,hiddenPath[t]],
                                  log = T )
          }else{
            score <- score + dnorm(x[t], 
                                   mean = model$Mu[hiddenPath[t]],
                                   sd = model$Sigma[,,hiddenPath[t]],
                                   log = T )
          }
        }
        score.states.cond[counter] <-  score
        decoded.important[counter, ] <-  decoded.states[[z]][i,]
        
        temp <- abs(data.deltas.ord[,,i])
        delta.score.cond[counter] <- sum(temp)
        if(USE.REP == TRUE)
        {
          max.delta.score.cond[counter] <- max(apply(temp, 2, mean))  
        }else
        {
          max.delta.score.cond[counter] <- max(temp)  
        }
        
        counter <- counter + 1
      }
    }
    decoded.cond[[conds[z]]] <- decoded.important
    score.cond[[z]] <- score.states.cond
    delta.cond[[z]] <- delta.score.cond
    max.delta.cond[[z]] <- max.delta.score.cond
  }
  
  
  # Sum all the scores into 1 per gene
  all.score.rep <- vector(mode = "numeric", length = sum(important.indx) )
  all.scores.changes <- vector(mode = "numeric", length = sum(important.indx) )
  
  all.score.delta <- vector(mode = "numeric", length = sum(important.indx) )
  all.score.delta.max <- vector(mode = "numeric", length = sum(important.indx) )
  
  for(i in 1:sum(important.indx))
  {
    score <- 0
    score.changes <- 0
    score.delta <- 0
    score.delta.max <- vector(mode = "numeric", length = length(conds))
    for(z in 1:length(conds))
    {
      score <- score + score.cond[[z]][i]
      score.changes <- score.changes + sum(decoded.cond[[conds[z]]][i, ] == "C")
      score.delta <- score.delta + delta.cond[[z]][i]
      score.delta.max[z] <- max.delta.cond[[z]][i]
    }
    all.score.rep[i] <-  score
    all.scores.changes[i] <- score.changes
    all.score.delta[i] <- score.delta
    all.score.delta.max[i] <- max(score.delta.max)
  }
  
  ranking <- order(all.scores.changes,all.score.delta, all.score.rep, decreasing = T)
  
  #################################
  
  important.data <- as.matrix(data[important.indx,])
  
  # Return all the important variables as a list
  print("Finished")
  final.results <- list()
  final.results[["data"]]    <- important.data
  final.results[["names"]]   <- rownames(important.data)
  final.results[["model"]]   <- model
  final.results[["scoreReplicates"]]   <- all.score.rep
  final.results[["scoreChanges"]]   <- all.scores.changes
  final.results[["ranking"]] <- ranking
  final.results[["decodedStates"]] <- decoded.cond
  final.results[["scoreDelta"]] <- all.score.delta
  final.results[["scoreMaxDelta"]] <- all.score.delta.max
  
  return(final.results)
}

#' @title Plot univariate feature selection object.
#'
#' @description 
#' This function plots in a pdf file the top ranked features stored in a
#' feature selection object that summarized the biological replicates 
#' via a statistical measure.
#' 
#' @param results  the output list from the featureSelectionHMM/featureSelectionOneConditionHMM function
#' @param times    vector with all the time-points
#' @param conds    vector with the condition names
#' @param rep      Number of biological replicates
#' @param path     Directory path to set the pdf. Defaults to '.'
#' @param pdf.name Character - pdf name. Defaults to 'topResults.pdf'
#' @param top      Number of genes to be plotted. Defaults to 100
#' 
#' @seealso \code{\link{plotMultivariate}} for not summarized biological replicates.
#' @keywords FS, 
#' @export
#' @examples
#' 
#' data("ratTimeSeries")
#' 
#' # All the necessary variables to run the script
#' times <- c(0,2,6,12,18,24)  # study time-points vector
#' conds <- c("CTRL","IKAROS") # Condition vector
#' reps <- 3                   # number of biological replicates
#' \dontrun{
#' # It has some random sampling inside, so a set.seed is recommended
#' # Case: Biological replicates summarized by median
#' set.seed(100)
#' summ.results <- featureSelectionHMM(data.ann,     # Gene expression data matrix
#'                                    times,         # time-point vector
#'                                    conds,         # Condition vector
#'                                    reps,          # Number of replicates
#'                                    USE.REP  = F,  # Do not use all replicate data
#'                                    func = median  # Summarize the data 
#' )  
#' 
#' plotUnivariate(summ.results, times, conds[2], reps, ".", "Fat.Rat.uni.res.pdf")
#' }
#' 
plotUnivariate <- function(results,                       
                           times,                         
                           conds,                         
                           rep,                           
                           path = ".",                    
                           pdf.name = "topResults.pdf",   
                           top = 100)                     
{
  data.graph <- results$data[results$ranking, ]
  top <- min(top, nrow(data.graph))
  data.graph <- data.graph[1:top, ]
  
  par(mfrow= c(1, 1))
  size <- rep*length(times)
  xlim <- c(min(times), max(times))
  
  relative.path <- getwd()
  setwd(path)
  pdf(file= pdf.name,height=6.5,width=20)
  
  for(i in 1:nrow(data.graph))
  {
    x <- data.graph[i, ]
    
    ylim <- c(min(x), max(x))  # this offset is for the legend
    plot(0,0,type = "n", 
         xlim = xlim,
         ylim = ylim,
         main = paste(rownames(data.graph)[i], "rank", i), xlab= "Time", ylab="Expression")
    
    for(z in 1:length(conds))
    {
      ind <- (1:size) + (z-1)*size
      x.x <- matrix(as.numeric(x[ind]), nrow=rep)
      x.x <- apply(x.x, 2, median)
      lines(times, x.x, col = z+1, lty = z)
      points(times, x.x, col = z+1, pch = z)
    }
    legend(xlim[1],ylim[2], conds, pch = 1:length(conds), lty = 1:length(conds), col = (1:length(conds)+1))
  }
  dev.off()
  setwd(relative.path)
}

#' @title Plot multivariate feature selection object.
#'
#' @description 
#' This function plots in a pdf file the top ranked features stored in a
#' feature selection object that did not summarized the biological replicates 
#' via a statistical measure.
#' 
#' @param results  the output list from the featureSelectionHMM/featureSelectionOneConditionHMM function
#' @param times    vector with all the time-points
#' @param conds    vector with the condition names
#' @param rep      Number of biological replicates
#' @param path     Directory path to set the pdf. Defaults to '.'
#' @param pdf.name Character - pdf name. Defaults to 'topResults.pdf'
#' @param top      Number of genes to be plotted. Defaults to 100
#' 
#' @seealso \code{\link{plotUnivariate}} for summarized biological replicates.
#' @keywords FS, 
#' @export
#' @examples
#' 
#' #' data("ratTimeSeries")
#' 
#' # All the necessary variables to run the script
#' times <- c(0,2,6,12,18,24)  # study time-points vector
#' conds <- c("CTRL","IKAROS") # Condition vector
#' reps <- 3                   # number of biological replicates
#' \dontrun{
#' # It has some random sampling inside, so a set.seed is recommended
#' # Case: Biological replicates not summarized
#' set.seed(100) 
#' mv.results <- featureSelectionHMM(data.ann,      # Gene expression data matrix
#'                                  times,         # time-point vector
#'                                  conds,         # Condition vector
#'                                  reps           # Number of replicates
#' )  
#'
#'
#' plotMultivariate(mv.results, times, conds[2], reps, ".", "Fat.Rat.multi.all.pdf")
#' }
#' 
plotMultivariate <- function(results,                       
                             times,                         
                             conds,                         
                             rep,                           
                             path = ".",                    
                             pdf.name = "topResults.pdf",   
                             top = 100)
{
  data.graph <- results$data[results$ranking, ]
  top <- min(top, nrow(data.graph))
  data.graph <- data.graph[1:top, ]
  
  par(mfrow= c(1, 1))
  size <- rep*length(times)
  xlim <- c(min(times), max(times))
  
  relative.path <- getwd()
  setwd(path)
  pdf(file= pdf.name,height=6.5,width=20)
  
  for(i in 1:nrow(data.graph))
  {
    x <- data.graph[i, ]
    
    ylim <- c(min(x), max(x))  # this offset is for the legend
    plot(0,0,type = "n", 
         xlim = xlim,
         ylim = ylim,
         main = paste(rownames(data.graph)[i], "rank", i), xlab= "Time", ylab="Expression")
    
    
    for(z in 1:length(conds))
    {
      ind <- (1:size) + (z-1)*size
      x.x <- matrix(as.numeric(x[ind]), nrow=rep)
      for(r in 1:rep)
      {
        lines(times, x.x[r,], col = z+1, lty = r)
        points(times, x.x[r,], col = z+1, pch = z)
      } 
    }
    legend(xlim[1],ylim[2], conds, pch = 1:length(conds), lty = 1, col = (1:length(conds)+1))
    
  }
  
  dev.off()
  setwd(relative.path)
}