library(ScITree)
dnaPath <- "C:/Users/tobyh/Documents/epidemic_simulator"
sink(file = "inference_output.txt")
infer.out<-infer(covariates = covariates,
                 moves.inputs = moves.inputs,
                 parsAux = pars.aux,
                 keyInits = para.key.inits,
                 priors = para.priors,
                 scalingFactors = para.sf,
                 seed = 1,
                 accTable = covariates[9],
                 t.sample = covariates[,6],
                 inputPath = "./inputs",
                 outputPath = "./outputs",
                 dnaPath = dnaPath,
                 dnaReference = F)
sink(file=NULL)