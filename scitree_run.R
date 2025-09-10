library(ScITree)

dir <- "./ScITree/inputs/"
covariates     <- read.csv(paste(dir,"covariates.csv",sep=''),     stringsAsFactors = FALSE, check.names = FALSE)
moves.inputs   <- read.csv(paste(dir,"moves.inputs.csv",sep=''),   stringsAsFactors = FALSE, check.names = FALSE)
pars.aux       <- read.csv(paste(dir,"pars.aux.csv",sep=''),       stringsAsFactors = FALSE, check.names = FALSE)
para.key.inits <- read.csv(paste(dir,"para.key.inits.csv",sep=''), stringsAsFactors = FALSE, check.names = FALSE)
para.priors    <- read.csv(paste(dir,"para.priors.csv",sep=''),    stringsAsFactors = FALSE, check.names = FALSE)
para.sf        <- read.csv(paste(dir,"para.sf.csv",sep=''),        stringsAsFactors = FALSE, check.names = FALSE)

dnaPath <- "./ScITree"
sink(file = "ScITree/inference_output.txt")
infer.out<-infer(covariates = covariates,
                 moves.inputs = moves.inputs,
                 parsAux = pars.aux,
                 keyInits = para.key.inits,
                 priors = para.priors,
                 scalingFactors = para.sf,
                 seed = 1,
                 accTable = covariates[9],
                 t.sample = covariates[,6],
                 inputPath = "./ScITree/inputs",
                 outputPath = "./ScITree/outputs",
                 dnaPath = dnaPath,
                 dnaReference = F)
sink(file=NULL)