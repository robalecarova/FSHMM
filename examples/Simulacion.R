##########
setwd("D:/GitHub/FSHMM/examples")
#############
# Part 1
# Data simulation
rm(list = ls())

library(MASS)

# Functions
generateVariables <- function(states, A, mean, sigma, Pi, times){
  sequence <- rep(0, times)
  state <- which.max(cumsum(Pi) > runif(1))
  sequence[1] <- states[state]
  for (i in 2:times){
    state <- which.max(cumsum(A[state, ]) > runif(1)) 
    sequence[i] <- states[state]
  }
  sequence
}

# model parameters
# sigma = c(1/100, 20/100) #sim 2
sigma = c(4/100, 1) #sim 1
sigma

mean = c(0,0)
A = matrix(c( 0.975, 0.025,
              0.0775,0.925),
           ncol = 2, byrow = T)
Pi = c(0.0, 1.0)
states = c("N", "C")
genes <- 1e4
times = 5
reps = 3

rand = 0.1 # sim1 ~ 0.1
changingGenes <- 1e3

replicatesCtrl = matrix(0, nrow = genes, ncol = (times+1)*reps)

replicatesCondition = matrix(0, nrow = genes, ncol = (times+1)*reps)
geneStatesCondition = matrix(0, nrow = genes, ncol = times)

# Temporals
replicaCond = rep(0, times)
replicaCont = rep(0, times)
replicasCond = rep(0,reps*(times+1))
replicasCont = rep(0,reps*(times+1))

set.seed(100)


# Changing genes
for(j in 1:changingGenes){
  condition <- generateVariables(states, A, mean, sigma, Pi, times)
  for (i in 1:times){
    estadoNumero = which.max(condition[i] == states)
    replicaCond[i] = rnorm(n = 1, mean[estadoNumero], sd = sqrt(sigma[estadoNumero])) # factor
    replicaCont[i] = rnorm(n = 1, mean[1], sd = sqrt(sigma[1]))  # factor
  }
  
  offset = rnorm(1,10,0.1)
  replicaCond <- cumsum(c(offset, replicaCond))
  replicaCont <- cumsum(c(offset, replicaCont))
  
  
  for(i in 1:(times+1)){
   replicasCond[(((i-1)*3)+1):(i*3)]  = rnorm(reps, replicaCond[i], rand) 
   replicasCont[(((i-1)*3)+1):(i*3)]  = rnorm(reps, replicaCont[i], rand) 
  }
  
  replicatesCondition[j,] = replicasCond
  replicatesCtrl[j,] = replicasCont
  geneStatesCondition[j,] = condition
  
}

# Not changing genes
 for(j in (changingGenes+1):genes){
  for (i in 1:times){
    replicaCond[i] = rnorm(n = 1, mean[1], sd = sqrt(sigma[1]))  # factor
    replicaCont[i] = rnorm(n = 1, mean[1], sd = sqrt(sigma[1]))  # factor
  }
  
  offset = rnorm(1,10,0.1)
  replicaCond <- cumsum(c(offset, replicaCond))
  replicaCont <- cumsum(c(offset, replicaCont))
  
  
  for(i in 1:(times+1)){
    replicasCond[(((i-1)*3)+1):(i*3)]  = rnorm(reps, replicaCond[i], rand) 
    replicasCont[(((i-1)*3)+1):(i*3)]  = rnorm(reps, replicaCont[i], rand) 
  }
  
  replicatesCondition[j,] = replicasCond
  replicatesCtrl[j,] = replicasCont
  geneStatesCondition[j,] = rep("N", times)
}


exp.matrix <- cbind(replicatesCtrl, replicatesCondition)
times <- seq(0,times)
conds <- c("Normal", "Condition") # Condition vector 
decoded <- apply(geneStatesCondition, 1, function(x) { paste(x, collapse = "")})
table(decoded)
hist(exp.matrix, breaks = 50)
save(exp.matrix, decoded, reps, times, conds, file = "Simulation.RData")


#############
# Part 2
# Model benchmark
rm(list = ls())
load("Simulation.RData")

# Load library
library(FSHMM)
rownames(exp.matrix) <- as.character(1:length(decoded))
hist(exp.matrix, breaks = 50)

# It has some random sampling inside, so a set.seed is recommended
set.seed(100)
summ.results <- featureSelectionHMM(exp.matrix,    # Gene expression data matrix
                                    times,         # time-point vector
                                    conds,         # Condition vector
                                    reps,          # Number of replicates
                                    iters = 100,
                                    USE.REP  = T,  # Do not use all replicate data
                                    func = median  # function to use to summarize the data 
)

# Lets see the model
{
# Lets compare
# Simulated data
simData <- decoded
tabla <- table(simData)
print(tabla)
totalSim <- length(simData) - tabla["NNNNN"]
totalSim

# Model data
totalModel <- nrow(summ.results$data)
percentage <- totalModel/totalSim   


modelData <- apply(summ.results$decodedStates$Condition, 1, function(x){paste(x, collapse = "")})
table(modelData)

decodedModel <- rownames(summ.results$data)
decodedSim <- decoded != "NNNNN"
decodedSim <- (1:length(decodedSim))[decodedSim]
decodedSim <- as.character(decodedSim)
}
summ.results$model
print(percentage)
sum(decodedModel %in% decodedSim)/length(decodedModel) 
# plotUnivariate(summ.results, times, conds, reps, ".", "prueba4.pdf")
plotMultivariate(summ.results, times, conds, reps, ".", "prueba4mv.pdf")

