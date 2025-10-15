#source('experiment1.R') # regenerates files
library(ScITree)

dir <- "./spatial-experiment/kappa-0/ScITree/inputs/"
covariates     <- read.csv(paste0(dir,"covariates.csv"),     stringsAsFactors = FALSE, check.names = FALSE)
moves.inputs   <- read.csv(paste0(dir,"moves.inputs.csv"),   stringsAsFactors = FALSE, check.names = FALSE)
pars.aux       <- read.csv(paste0(dir,"pars.aux.csv"),       stringsAsFactors = FALSE, check.names = FALSE)
para.key.inits <- read.csv(paste0(dir,"para.key.inits.csv"), stringsAsFactors = FALSE, check.names = FALSE)
para.priors    <- read.csv(paste0(dir,"para.priors.csv"),    stringsAsFactors = FALSE, check.names = FALSE)
para.sf        <- read.csv(paste0(dir,"para.sf.csv"),        stringsAsFactors = FALSE, check.names = FALSE)

dnaPath <- "./spatial-experiment/kappa-0/ScITree"
sink(file = "./spatial-experiment/kappa-0/ScITree/inference_output.txt")
infer.out<-infer(covariates = covariates,
                 moves.inputs = moves.inputs,
                 parsAux = pars.aux,
                 keyInits = para.key.inits,
                 priors = para.priors,
                 scalingFactors = para.sf,
                 seed = 1,
                 accTable = covariates[,9],
                 t.sample = covariates[,6],
                 inputPath = "./spatial-experiment/kappa-0/ScITree/inputs",
                 outputPath = "./spatial-experiment/kappa-0/ScITree/outputs",
                 dnaPath = "./spatial-experiment/kappa-0/ScITree/gen_inputs",
                 dnaReference = F)
sink(file=NULL)