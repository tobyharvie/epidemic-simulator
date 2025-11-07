prepare_data <- function (epi,sequences,ps,dir,v){
  
  nexus_str<-sequences
  # PREPROCESSING
  # Example string
  # Split into lines
  lines <- unlist(strsplit(nexus_str, "\n"))
  # Keep only labelled rows (start with digits or letters)
  labelled <- grep("^[[:alnum:]_]+", lines, value = TRUE)
  # Extract numeric labels
  labels <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", labelled)))
  # Filter multiples of (n_pop*100)
  keep <- labelled[!is.na(labels) & labels %% (100) == 0]
  # Replace the numeric labels with divided-by-(n_pop*100) labels
  keep_div <- sapply(keep, function(line) {
    num <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    new_num <- num / (100)
    sub("^[0-9]+", new_num, line)
  })
  # Collapse back into a string
  cleaned_str <- paste(keep_div, collapse = "\n")
  #convert back to NEXUS
  matrix_lines <- unlist(strsplit(cleaned_str, "\n"))
  matrix_lines <- matrix_lines[nchar(matrix_lines) > 0]
  
  # BREATH requires NEXUS
  # switch unit to days
  time_range <- 1000 # in days
  max_time <- max(epi[is.finite(epi)], na.rm = TRUE, finite = TRUE)
  start_date <- as.Date("2020-01-01") # arbitrary starting date
  # include time data
  matrix_lines_new <- sapply(matrix_lines, function(line) {
    node <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    seq <- sub("^[0-9]+\\s+", "", line)
    # Calculate timestamp as a Date
    days <- (epi[epi[, 1] == node, 6] / max_time) * time_range
    ts <- start_date + round(days)
    # Build new label
    paste0(node, "_", ts, "    ", seq)
  })
  nexus_out <- c(
    "#NEXUS","","BEGIN DATA;",
    paste0("    DIMENSIONS NTAX=", length(matrix_lines_new), " NCHAR=", max(nchar(gsub(".*?\\s+", "", matrix_lines_new))), ";"),
    "    FORMAT DATATYPE=DNA MISSING=? GAP=-;",
    "MATRIX",matrix_lines_new,";","END;"
  )
  #print(matrix_lines_new)
  #writeLines(nexus_out, paste(dir,"BREATH/seqs.nex",sep=''))
  
  # calculate sampling shape and rate for sampling hazard
  stimes <- (epi[,6]-epi[,3])/max_time*1000/365
  stimes <- stimes[!is.na(stimes)]
  sfit <- fitdistr(stimes, densfun = "gamma")
  sa <-sfit$estimate["shape"]
  sb <-sfit$estimate["rate"]
  #sk <-mean(stimes,na.rm=TRUE)
  #sv <-var(stimes,na.rm=TRUE)
  #sa <-sk^2/sv #shape
  #sb <-sk/sv #rate
  sp <- length(epi[!is.na(epi[, 6]), 6])/length(epi[!is.na(epi[, 3]), 3])
  #print(paste('Sampling prob:', sp))
  #print(paste('Sampling rate:', sb))
  #print(paste('Sampling shape:', sa))
  
  epidf <- as.data.frame(epi)
  df_merged <- merge(epidf, epidf[, c("Node ID", "Etime")],
                     by.x = "Parent", by.y = "Node ID",
                     all.x = TRUE, suffixes = c("", "_parent"))
  # we need these in unit time years
  
  # Calculate transmission interval
  df_merged$transmission_interval <- df_merged$Etime - df_merged$Etime_parent
  #df_merged$transmission_interval<-df_merged$transmission_interval/max_time*1000
  #ttimes <- c(df_merged$transmission_interval,(epi[,6]-epi[,3]))/max_time*1000/365
  ttimes<-df_merged$transmission_interval/max_time*1000/365
  ttimes <- ttimes[!is.na(ttimes)]
  tfit <- fitdistr(ttimes, densfun = "gamma")
  ta <-tfit$estimate["shape"]
  tb <-tfit$estimate["rate"]
  
  # Average transmission time (exclude NAs)
  #tk <- mean(ttimes, na.rm = TRUE)
  #tk <-mean(ttimes,na.rm=TRUE)
  #tv <- var(ttimes, na.rm = TRUE)
  #tv <-var(ttimes,na.rm=TRUE)
  #ta <-tk^2/tv #shape
  #tb <-tk/tv #rate
  
  #print(paste('Transmission rate:', tb))
  #print(paste('Transmission shape:', ta))
  # average number of people that a single node infects
  tc <- mean(table(factor(epi[,2],levels=1:length(epi[,1]))))
  #print(paste("transmission constant: ", tc))
  ni <-  nrow((epi[!is.na(epi[,6]),]))
  breathparamnames <- c('sampling prob','sampling rate','sampling shape','transmission rate','transmission shape','transmission constant','initial pop size')
  breathparams <- c(sp,sb,sa,tb,ta,tc,ni)
  breathinput <- data.frame(breathparamnames,breathparams)
  #file_conn <- file(paste0(dir,"BREATH/params.txt"), open = "w")
  #writeLines(as.character(breathparamnames), con = file_conn)
  #writeLines(as.character(breathparams), con = file_conn)
  #close(file_conn)
  
  # GENERATE BREATH XML FILE
  chainLength <- 4000000
  storeEvery <- 1000
  logEvery <- 2000
  if (1) # WRITE XML FILE
  {
    fn=paste0(dir,'breath-',ps*100,'-',v,'.xml')
    if (file.exists(fn)) {
      file.remove(fn)
    }
    f <- file(fn, "a") 
    
    xml_content <- r"(<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate='Standard' beautistatus='' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.base.evolution.alignment:beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.evolution.tree.coalescent:beast.pkgmgmt:beast.base.core:beast.base.inference.util:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood" required="BEAST.base v2.7.8:BREATH v0.0.1" version="2.7">)"
    writeLines(xml_content, f)
    xml_content <- r"(<data
    id="seqs"
    spec="Alignment"
    name="alignment">)"
    writeLines(xml_content, f)
    for (i in 1:length(matrix_lines_new)){
      x <- strsplit(matrix_lines_new[[i]], " +")
      xml_content <- sprintf(
        '<sequence id="%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>',
        paste0('seq_',x[[1]][1]),x[[1]][1],x[[1]][2]
      )
      writeLines(xml_content, f)
    }
    # sequences
    xml_content <- r"(
      </data>
      <map name="Uniform" >beast.base.inference.distribution.Uniform</map>
      <map name="Exponential" >beast.base.inference.distribution.Exponential</map>
      <map name="LogNormal" >beast.base.inference.distribution.LogNormalDistributionModel</map>
      <map name="Normal" >beast.base.inference.distribution.Normal</map>
      <map name="Beta" >beast.base.inference.distribution.Beta</map>
      <map name="Gamma" >beast.base.inference.distribution.Gamma</map>
      <map name="LaplaceDistribution" >beast.base.inference.distribution.LaplaceDistribution</map>
      <map name="prior" >beast.base.inference.distribution.Prior</map>
      <map name="InverseGamma" >beast.base.inference.distribution.InverseGamma</map>
      <map name="OneOnX" >beast.base.inference.distribution.OneOnX</map>
  )"
    writeLines(xml_content, f)
    ###  ----------- MCMC ------------
    xml_content <- sprintf(
      '<run id="mcmc" spec="MCMC" chainLength="%s">',
      chainLength
    )
    writeLines(xml_content, f)
    ###  ----------- STORE EVERY ------------
    xml_content <- sprintf(
      '<state id="state" spec="State" storeEvery="%s">',
      storeEvery
    )
    writeLines(xml_content, f)
    xml_content <- '<tree id="Tree.t:seqs" spec="beast.base.evolution.tree.Tree" name="stateNode">'
    writeLines(xml_content, f)
    ### node names
    node_names <- epi[!is.na(epi[, 6]), 1]
    time_values <- epi[!is.na(epi[, 6]), 6]
    days <- (time_values / max_time) * time_range
    ts <- start_date + round(days)
    date_strings <- format(ts, "%Y-%m-%d")
    value_pairs <- paste0(node_names, "_", date_strings, "=", date_strings)
    value_string <- paste(value_pairs, collapse = ",")
    xml_trait <- sprintf(
      '<trait id="dateTrait.t:seqs" spec="beast.base.evolution.tree.TraitSet" dateFormat="yyyy-M-dd" traitname="date" value="%s">',
      value_string
    )
    writeLines(xml_trait, f)
    xml_content <- r'(<taxa id="TaxonSet.seqs" spec="TaxonSet">
      <alignment idref="seqs"/>
      </taxa>
      </trait>
      <taxonset idref="TaxonSet.seqs"/>
      </tree>
      <parameter id="clockRate.c:seqs" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>
      <parameter id="kappa.s:seqs" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>
      <parameter id="freqParameter.s:seqs" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>)'
    writeLines(xml_content, f)
    
    # BLOCK COUNTS
    block_start <- 2*nrow(epi[!is.na(epi[,6]),])-2
    blockcount <- block_start + 1
    xml_content <- sprintf(
      '<parameter id="blockstart.t:seqs" spec="parameter.RealParameter" dimension="%s" lower="0.0" name="stateNode" upper="1.0">0.5</parameter>
        <parameter id="blockend.t:seqs" spec="parameter.RealParameter" dimension="%s" lower="0.0" name="stateNode" upper="1.0">0.5</parameter>
        <stateNode id="blockcount.t:seqs" spec="parameter.IntegerParameter" dimension="%s" lower="-1" upper="100000">0</stateNode>
    ',
    block_start,block_start,blockcount
    )
    writeLines(xml_content, f)
    
    xml_content <- r'(
      <parameter id="transmissionPopSize.t:seqs" spec="parameter.RealParameter" lower="0.1" name="stateNode" upper="60.0">0.5</parameter>
      <parameter id="transmissionOrigin.t:seqs" spec="parameter.RealParameter" name="stateNode">10.0</parameter>
      </state>
      <init id="RandomTree.t:seqs" spec="RandomTree" estimate="false" initial="@Tree.t:seqs" taxa="@seqs">
      <populationModel id="ConstantPopulation0.t:seqs" spec="ConstantPopulation">
      <parameter id="randomPopSize.t:seqs" spec="parameter.RealParameter" name="popSize">1.0</parameter>
      </populationModel>
      </init>
      <distribution id="posterior" spec="CompoundDistribution">
      <distribution id="prior" spec="CompoundDistribution">
      <distribution id="transmissionLikelihood.t:seqs" spec="breath.distribution.TransmissionTreeLikelihood" blockcount="@blockcount.t:seqs" blockend="@blockend.t:seqs" blockstart="@blockstart.t:seqs" origin="@transmissionOrigin.t:seqs" tree="@Tree.t:seqs">
      <populationModel id="ConstantTransmissionPopulation.t:seqs" spec="ConstantPopulation" popSize="@transmissionPopSize.t:seqs"/>
      <parameter id="endTime.t:seqs" spec="parameter.RealParameter" name="endTime">-1.0</parameter>
    )'
    writeLines(xml_content, f)
    
    # SAMPLING AND TRANSMISSION
    xml_content <- sprintf('
      <samplingHazard id="samplingHazard.t:seqs" spec="breath.distribution.GammaHazardFunction">
      <shape id="Function$Constant" spec="Function$Constant" value="%s"/>
      <rate id="Function$Constant1" spec="Function$Constant" value="%s"/>
      <C id="Function$Constant2" spec="Function$Constant" value="%s"/>
      </samplingHazard>
      <transmissionHazard id="transmissionHazard.t:seqs" spec="breath.distribution.GammaHazardFunction">
      <shape id="Function$Constant3" spec="Function$Constant" value="%s"/>
      <rate id="Function$Constant4" spec="Function$Constant" value="%s"/>
      <C id="Function$Constant5" spec="Function$Constant" value="2.0"/>
      </transmissionHazard>',
      sa,sb,sp,ta,tb)
    writeLines(xml_content, f)
    
    xml_content <- r'(</distribution>
      <prior id="BlockCountPrior.t:seqs" name="distribution" x="@blockcount.t:seqs">
      <Uniform id="Uniform.3" lower="-1.0" name="distr" upper="4.0"/>
      </prior>
      <prior id="BlockEndPrior.t:seqs" name="distribution" x="@blockend.t:seqs">
      <Uniform id="Uniform.5" name="distr"/>
      </prior>
      <prior id="BlockStartPrior.t:seqs" name="distribution" x="@blockstart.t:seqs">
      <Uniform id="Uniform.4" name="distr"/>
      </prior>
      <prior id="ClockPrior.c:seqs" name="distribution" x="@clockRate.c:seqs">
      <Uniform id="Uniform.0" name="distr" upper="Infinity"/>
      </prior>
      <prior id="FrequenciesPrior.s:seqs" name="distribution" x="@freqParameter.s:seqs">
      <distr id="Dirichlet.0" spec="distribution.Dirichlet">
      <parameter id="RealParameter.3" spec="parameter.RealParameter" dimension="4" estimate="false" name="alpha">4.0 4.0 4.0 4.0</parameter>
      </distr>
      </prior>
      <prior id="KappaPrior.s:seqs" name="distribution" x="@kappa.s:seqs">
      <LogNormal id="LogNormalDistributionModel.0" name="distr">
      <parameter id="RealParameter.1" spec="parameter.RealParameter" estimate="false" name="M">1.0</parameter>
      <parameter id="RealParameter.2" spec="parameter.RealParameter" estimate="false" name="S">1.25</parameter>
      </LogNormal>
      </prior>
      <prior id="TransmissionOriginPrior.t:seqs" name="distribution" x="@transmissionOrigin.t:seqs">
      <Uniform id="Uniform.7" name="distr" upper="Infinity"/>
      </prior>
      <prior id="TransmissionPopSizePrior.t:seqs" name="distribution" x="@transmissionPopSize.t:seqs">
      <Uniform id="Uniform.6" lower="0.1" name="distr" upper="60.0"/>
      </prior>
      </distribution>
      <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">
      <distribution id="treeLikelihood.seqs" spec="ThreadedTreeLikelihood" data="@seqs" tree="@Tree.t:seqs">
      <siteModel id="SiteModel.s:seqs" spec="SiteModel">
      <parameter id="mutationRate.s:seqs" spec="parameter.RealParameter" estimate="false" lower="0.0" name="mutationRate">1.0</parameter>
      <parameter id="gammaShape.s:seqs" spec="parameter.RealParameter" estimate="false" lower="0.1" name="shape">1.0</parameter>
      <parameter id="proportionInvariant.s:seqs" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
      <substModel id="hky.s:seqs" spec="HKY" kappa="@kappa.s:seqs">
      <frequencies id="estimatedFreqs.s:seqs" spec="Frequencies" frequencies="@freqParameter.s:seqs"/>
      </substModel>
      </siteModel>
      <branchRateModel id="StrictClock.c:seqs" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:seqs"/>
      </distribution>
      </distribution>
      <distribution id="fossilCalibrations" spec="CompoundDistribution"/>
      </distribution>
      <operator id="StrictClockRateScaler.c:seqs" spec="AdaptableOperatorSampler" weight="1.5">
      <parameter idref="clockRate.c:seqs"/>
      <operator id="AVMNOperator.seqs" spec="kernel.AdaptableVarianceMultivariateNormalOperator" allowNonsense="true" beta="0.05" burnin="400" initial="800" weight="0.1">
      <transformations id="AVMNSumTransform.seqs" spec="operator.kernel.Transform$LogConstrainedSumTransform">
      <f idref="freqParameter.s:seqs"/></transformations>
      <transformations id="AVMNLogTransform.seqs" spec="operator.kernel.Transform$LogTransform">
      <f idref="clockRate.c:seqs"/>
      <f idref="kappa.s:seqs"/>
      </transformations>
      <transformations id="AVMNNoTransform.seqs" spec="operator.kernel.Transform$NoTransform">
      <f idref="Tree.t:seqs"/>
      </transformations>
      </operator>
      <operator id="StrictClockRateScalerX.c:seqs" spec="kernel.BactrianScaleOperator" parameter="@clockRate.c:seqs" upper="10.0" weight="3.0"/>
      </operator>
      <operator id="strictClockUpDownOperator.c:seqs" spec="AdaptableOperatorSampler" weight="1.5">
      <parameter idref="clockRate.c:seqs"/>
      <tree idref="Tree.t:seqs"/>
      <operator idref="AVMNOperator.seqs"/>
      <operator id="strictClockUpDownOperatorX.c:seqs" spec="operator.kernel.BactrianUpDownOperator" scaleFactor="0.75" weight="3.0">
      <up idref="clockRate.c:seqs"/>
      <down idref="Tree.t:seqs"/>
      </operator>
      </operator>
      <operator id="KappaScaler.s:seqs" spec="AdaptableOperatorSampler" weight="0.05">
      <parameter idref="kappa.s:seqs"/>
      <operator idref="AVMNOperator.seqs"/>
      <operator id="KappaScalerX.s:seqs" spec="kernel.BactrianScaleOperator" parameter="@kappa.s:seqs" scaleFactor="0.1" upper="10.0" weight="0.1"/>
      </operator>
      <operator id="FrequenciesExchanger.s:seqs" spec="AdaptableOperatorSampler" weight="0.05">
      <parameter idref="freqParameter.s:seqs"/>
      <operator idref="AVMNOperator.seqs"/>
      <operator id="FrequenciesExchangerX.s:seqs" spec="operator.kernel.BactrianDeltaExchangeOperator" delta="0.01" weight="0.1">
      <parameter idref="freqParameter.s:seqs"/>
      </operator>
      </operator>
      <operator id="transmissionLikelihoodBICEPSEpochTop.t:seqs" spec="EpochFlexOperator" scaleFactor="0.1" tree="@Tree.t:seqs" weight="2.0"/>
      <operator id="transmissionLikelihoodBICEPSEpochAll.t:seqs" spec="EpochFlexOperator" fromOldestTipOnly="false" scaleFactor="0.1" tree="@Tree.t:seqs" weight="2.0"/>
      <operator id="transmissionLikelihoodBICEPSTreeFlex.t:seqs" spec="TreeStretchOperator" scaleFactor="0.01" tree="@Tree.t:seqs" weight="2.0"/>
      <operator id="transmissionLikelihoodTreeRootScaler.t:seqs" spec="kernel.BactrianScaleOperator" rootOnly="true" scaleFactor="0.1" tree="@Tree.t:seqs" upper="10.0" weight="3.0"/>
      <operator id="transmissionLikelihoodUniformOperator.t:seqs" spec="kernel.BactrianNodeOperator" tree="@Tree.t:seqs" weight="30.0"/>
      <operator id="transmissionLikelihoodSubtreeSlide.t:seqs" spec="kernel.BactrianSubtreeSlide" tree="@Tree.t:seqs" weight="15.0"/>
      <operator id="transmissionLikelihoodNarrow.t:seqs" spec="Exchange" tree="@Tree.t:seqs" weight="15.0"/>
      <operator id="transmissionLikelihoodWide.t:seqs" spec="Exchange" isNarrow="false" tree="@Tree.t:seqs" weight="3.0"/>
      <operator id="transmissionLikelihoodWilsonBalding.t:seqs" spec="WilsonBalding" tree="@Tree.t:seqs" weight="3.0"/>
      <operator id="transmissionInfectionMover.t:seqs" spec="breath.operator.InfectionMover" blockcount="@blockcount.t:seqs" blockend="@blockend.t:seqs" blockstart="@blockstart.t:seqs" likelihood="@transmissionLikelihood.t:seqs" weight="50.0"/>
      <operator id="transmissionBlockOperator.t:seqs" spec="breath.operator.BlockOperator" blockcount="@blockcount.t:seqs" blockend="@blockend.t:seqs" blockstart="@blockstart.t:seqs" tree="@Tree.t:seqs" weight="50.0"/>
      <operator id="transmissionPopSizeScaler.t:seqs" spec="kernel.BactrianScaleOperator" parameter="@transmissionPopSize.t:seqs" scaleFactor="0.1" upper="10.0" weight="5.0"/>
      <operator id="transmissionOriginScaler.t:seqs" spec="kernel.BactrianScaleOperator" parameter="@transmissionOrigin.t:seqs" scaleFactor="0.1" upper="10.0" weight="0.5"/>)'
    writeLines(xml_content, f)
    
    # log and store
    xml_content <- sprintf('<logger id="tracelog" spec="Logger" fileName="breath-logs/$(filebase).log" logEvery="%s" model="@posterior" sanitiseHeaders="true" sort="smart">
                <log idref="posterior"/>
                <log idref="likelihood"/>
                <log idref="prior"/>
                <log idref="treeLikelihood.seqs"/>
                <log id="TreeHeight.t:seqs" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t:seqs"/>
                <log idref="clockRate.c:seqs"/>
                <log idref="kappa.s:seqs"/>
                <log idref="freqParameter.s:seqs"/>
                <log idref="transmissionLikelihood.t:seqs"/>
                <log idref="transmissionPopSize.t:seqs"/>
                <log idref="transmissionOrigin.t:seqs"/>
                </logger>
                <logger id="screenlog" spec="Logger" logEvery="%s">
                <log idref="posterior"/>
                <log idref="likelihood"/>
                <log idref="prior"/>
                </logger>
                <logger id="treelog.t:seqs" spec="Logger" fileName="breath-trees/$(filebase)-$(tree).trees" logEvery="%s" mode="tree">
                <log id="TreeWithMetaDataLogger.t:seqs" spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree.t:seqs">
                <metadata idref="blockstart.t:seqs"/>
                <metadata idref="blockend.t:seqs"/>
                <metadata idref="blockcount.t:seqs"/>
                </log>
                </logger>
                <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>
                </run>
                </beast>',
                logEvery, logEvery, logEvery
    )
    writeLines(xml_content, f)
    close(f)
  }
  
  # sanity checks
  #curve(dgamma(x, shape = sa, scale = 1/sb),from = 0, to = max_time)
  #curve(dgamma(x, shape = ta, scale = 1/tb),from = 0, to = max_time)
  
  # ---------------
  # SCOTTI Pipeline
  # ---------------
  # FASTA
  matrix_lines_new <- sapply(matrix_lines, function(line) {
    node <- as.numeric(sub("^([0-9]+).*", "\\1", line))
    seq <- sub("^[0-9]+\\s+", "", line)
    # Build new label
    if (!is.na(epi[epi[, 1] == node, 6])){
      paste0('>', node, "\n", seq)
    }
    else{ paste0('')}
  })
  fasta_out <- c(
    matrix_lines_new
  )
  writeLines(fasta_out, paste0(dir,"seqs.fasta"))
  # sampling dates
  epi[1,2]<-1 # this allows us to use complete.cases and is reverted after writing these fiels
  write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,6]), file = paste(dir,"dates.csv",sep=''),sep = ",",row.names = FALSE, col.names = FALSE)
  # hosts .. unneeded?
  write.table(cbind(epi[complete.cases(epi),][,1],epi[complete.cases(epi),][,1]), file = paste(dir,"hosts.csv",sep=''),sep = ",",row.names = FALSE, col.names = FALSE)
  # earliest and latest times host in infection
  # assume sampled during infectious period
  # thus earliest possible time infectious period + exposed period at 99th percentile
  lower <- qgamma(p = 0.99, shape=ke, scale = thetae) + qgamma(p = 0.99, shape=ki, scale = thetai)
  # for upper, we assume that we remove after sampling
  upper <- qgamma(p = 0.99, shape=ki, scale = thetai)
  write.table(cbind(epi[complete.cases(epi),][,1],pmax(0,epi[complete.cases(epi),][,6]-lower),epi[complete.cases(epi),][,6]+upper), file = paste(dir,"hostTimes.csv",sep=''),sep = ",",row.names = FALSE, col.names = FALSE)
  # use py script to generate xml
  epi[1,2]<-NA
  args <- c(
    "SCOTTI_generate_xml.py",
    "--fasta", paste0(dir,"seqs.fasta"),
    "--dates", paste0(dir,"dates.csv"),
    "--hosts", paste0(dir,"hosts.csv"),
    "--hostTimes", paste0(dir,"hostTimes.csv"),
    "--output", paste0(dir,"scotti-",ps*100,"-",v),
    "--maxHosts",  dim(epi[!is.na(epi[,5]),])[1],
    "--numIter", '5000000',
    "--tracelog", 1000,
    "--treelog", 1000,
    "--screenlog", 1000
  )
  
  if (file.exists(paste0(dir,'scotti-',ps*100,'.xml'))) {
    #Delete file if it exists
    file.remove(paste0(dir,'scotti-',ps*100,'.xml'))
  }
  
  # Call (modified) python script to generate xml file
  result<-system2("python", args = args, stdout = TRUE, stderr = TRUE)
  
  # delete other files
  #file.remove(paste0(dir,'dates.csv'))
  #file.remove(paste0(dir,'hosts.csv'))
  #file.remove(paste0(dir,'hostTimes.csv'))
  #file.remove(paste0(dir,'seqs.fasta'))
}
