
# clean all varibale envirnment
rm(list = ls())

# install packages
# point this to the source file location
install.packages("./FSHMM_1.0.3.tar.gz", repos = NULL, type = "source")

# Load library
library(FSHMM)

## DataMouse variables: GSE39549
# It loads the annotated data matrix "data.ann"
load("DataMouse.RData") 
times <- c(0,2,4,6,8,12,16,20,24) # study time-points vector
reps <- 3                         # number of biological replicates
conds <- c("Normal", "Fat")       # Condition vector


# DataRat variables: Japanese Toxicogenomics project
# It loads the annotated data matrix "data.ann"
load("DataRat.RData") 
times <- c(2,8,24)                 # study time-points vector
conds <- c("CTRL","LOW", "MED")    # Condition vector
reps <- 2                          # number of biological replicates

# DataRatRNAseq variables: GSE75417
# It loads the annotated data matrix "data.ann"
load("DataMouseRNAseq.RData") 
times <- c(0,2,6,12,18,24)                 # study time-points vector
conds <- c("Control", "IKAROS")    # Condition vector
reps <- 3                          # number of biological replicates


# It has some random sampling inside, so a set.seed is recommended
set.seed(100)
summ.results <- featureSelectionHMM(data.ann,      # Gene expression data matrix
                                    times,         # time-point vector
                                    conds,         # Condition vector
                                    reps,          # Number of replicates
                                    USE.REP  = F,  # Do not use all replicate data
                                    func = median  # function to use to summarize the data 
)

set.seed(100) 
mv.results <- featureSelectionHMM(data.ann,      # Gene expression data matrix
                                  times,         # time-point vector
                                  conds,         # Condition vector
                                  reps           # Number of replicates
)  



plotMultivariate(mv.results, times, conds, reps, ".", "multivariate.pdf")
plotUnivariate(summ.results, times,conds, reps, ".", "univariate.pdf")

