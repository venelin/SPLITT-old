# the following code is executed by a simName.R script found in one of the subdirectories below
# the directory DATA. 
# It assumes that there are the following global variables:
# tree : a phylo object
# tMean : Numeric the mean root-tip distance in tree;
# N: number of tips in tree.

# used to generate parameters from the prior distribution; called at the beginning of 
# each replication.
if(simType == "POUMM") {
  rParPrior <- function() {
    c(alpha = rexp(1, rate = .1), 
      theta = runif(1, 2, 8), 
      sigma = rexp(1, .4), 
      sigmae = rexp(1, rate = 1),
      g0 = rnorm(1, mean = 5, sd = 5))
  }
} else if(simType == "PMM") {
  rParPrior <- function() {
    c(alpha = 0, #rexp(1, rate = .1), 
      theta = runif(1, 2, 8), 
      sigma = rexp(1, .4), 
      sigmae = rexp(1, rate = 1),
      g0 = rnorm(1, mean = 5, sd = 5))
  }
}

specPOUMM <- POUMM::specifyPOUMM_ATSSeG0(
  tMean = tMean,
  parLower = c(alpha = 0, theta = 0, sigma = 0, sigmae = 0, 
               g0 = qnorm(0.001, 5, 5)), 
  
  parUpper = c(alpha = 50, theta = 8, sigma = qexp(0.999, 0.4), sigmae = qexp(.999, 1), 
               g0 = qnorm(0.999, 5, 5)),
  
  parPriorMCMC = function(par) {
    dexp(par[1], rate = .1, log = TRUE) +
      dunif(par[2], min = 2, max = 8, log = TRUE) +
      dexp(par[3], rate = 0.4, log = TRUE) +
      dexp(par[4], rate = 1, log = TRUE) + 
      dnorm(par[5], mean = 5, sd = 5, log = TRUE)
  },
  
  parInitMCMC = function(chainNo, fitML = NULL) {
    if(!is.null(fitML)) {
      parML <- c(fitML$par[c('alpha', 'theta', 'sigma', 'sigmae', 'g0')])
      names(parML) <- c('alpha', 'theta', 'sigma', 'sigmae', 'g0')
    } else {
      parML <- NULL
    }
    
    init <- rbind(
      c(alpha = 0, theta = 0, sigma = qexp(.9, 0.4), sigmae = 0, g0 = 0),
      parML,
      c(alpha = 0, theta = 4, sigma = qexp(.1, 0.4), sigmae = 1, g0 = 6)
    )
    
    init[(chainNo - 1) %% nrow(init) + 1, ]
  }, 
  
  g0Prior = list(mean=5, var=25),
  nSamplesMCMC = 1e6,  thinMCMC = 10
  
)

specOU <- POUMM::specifyPOUMM_ATSG0(
  tMean = tMean,
  parLower = c(alpha = 0, theta = 0, sigma = 0, g0 = qnorm(0.001, 5, 5)), 
  
  parUpper = c(alpha = 50, theta = 8, sigma = qexp(0.999, 0.4), g0 = qnorm(0.999, 5, 5)),
  
  parPriorMCMC = function(par) {
    dexp(par[1], rate = .1, log = TRUE) +
      dunif(par[2], min = 2, max = 8, log = TRUE) +
      dexp(par[3], rate = 0.4, log = TRUE) +
      dnorm(par[4], mean = 5, sd = 5, log = TRUE)
  },
  
  parInitMCMC = function(chainNo, fitML = NULL) {
    if(!is.null(fitML)) {
      parML <- c(fitML$par[c('alpha', 'theta', 'sigma', 'g0')])
      names(parML) <- c('alpha', 'theta', 'sigma', 'g0')
    } else {
      parML <- NULL
    }
    
    init <- rbind(
      c(alpha = 0, theta = 0, sigma = qexp(.9, 0.4), g0 = 0),
      parML,
      c(alpha = 0, theta = 4, sigma = qexp(.1, 0.4), g0 = 6)
    )
    
    init[(chainNo - 1) %% nrow(init) + 1, ]
  }, 
  
  g0Prior = list(mean = 5, var = 25),
  nSamplesMCMC = 1e6, thinMCMC = 10,
  sigmaeFixed = 0
)

specPMM <- POUMM::specifyPMM_SSeG0(
  tMean = tMean,
  parLower = c(sigma = 0, sigmae = 0, g0 = qnorm(0.001, 5, 5)), 
  
  parUpper = c(sigma = qexp(.999, 0.4), sigmae = qexp(.999, 1), g0 = qnorm(0.999, 5, 5)),
  
  parPriorMCMC = function(par) {
    dexp(par[1], rate = 0.4, log = TRUE) +
      dexp(par[2], rate = 1, log = TRUE) + 
      dnorm(par[3], mean = 5, sd = 5, log = TRUE)
  },
  
  parInitMCMC = function(chainNo, fitML = NULL) {
    if(!is.null(fitML)) {
      parML <- fitML$par
    } else {
      parML <- NULL
    }
    
    init <- rbind(
      c(sigma = qexp(.9, rate = 0.4), sigmae = 0, g0 = 0),
      parML,
      c(sigma = qexp(.1, rate = 0.4), sigmae = 1, g0 = 6)
    )
    
    init[(chainNo - 1) %% nrow(init) + 1, ]
  }, 
  
  g0Prior = list(mean=5, var=25),
  
  nSamplesMCMC = 1e6, thinMCMC = 10
)


# the argument p is a list containign a single integer entry called id. This entry
# is used as a set.seed.
replication <- function(p) {
  # load tree from file found in the same directory
  p <- as.list(p)
  p$id <- as.integer(p$id)
  
  N <- length(tree$tip.label)

  set.seed(p$id)
  parTrue <- rParPrior()
  
  g <- rVNodesGivenTreePOUMM(
    tree, z0 = parTrue["g0"], parTrue["alpha"], parTrue["theta"], parTrue["sigma"])
  
  e <- rnorm(length(g), 0, parTrue["sigmae"])
  z <- g + e
  
  H2eTrue <- var(g[1:N])/var(z[1:N])
  
  
  parTrueFull <- parTrue
  
  parTrueFull["H2tMean"] <- 
    POUMM::H2(parTrue["alpha"], parTrue["sigma"],  parTrue["sigmae"], tMean)
  
  parTrueFull["H2e"] <- POUMM::H2e(z[1:N], parTrue["sigmae"])
  
  parTrueFull["H2tInf"] <- POUMM::H2(alpha = parTrueFull["alpha"], 
                                     sigma = parTrueFull["sigma"], 
                                     sigmae = parTrueFull["sigmae"], t = Inf)
  
  cat("True parameters (full):\n")
  print(parTrueFull)
  
  cat("True heritability:\n")
  print(H2eTrue)
  
  cat("Log-likelihood at true parameters:\n") 
  
  print(logLikAtTrueParam <- likPOUMMGivenTreeVTips(
    z, tree, parTrueFull["alpha"], parTrueFull["theta"], parTrueFull["sigma"], 
    parTrueFull["sigmae"], g0 = parTrueFull["g0"], 
    g0Prior = list(mean=0, var = 4)))
  
  pruneInfo <- pruneTree(tree, z)
  
  print(microbenchmark::microbenchmark(
    "POUMM: C++, Armadillo" = {
      likPOUMMGivenTreeVTipsC(integrator = pruneInfo$integrator,
                               parTrueFull["alpha"], parTrueFull["theta"], parTrueFull["sigma"], 
                               parTrueFull["sigmae"], g0 = parTrueFull["g0"], 
                               g0Prior = list(mean=0, var = 4))
      },
    "POUMM: C++, omp" = {
      likPOUMMGivenTreeVTipsC4(integrator = pruneInfo$integrator,
        parTrueFull["alpha"], parTrueFull["theta"], parTrueFull["sigma"], 
        parTrueFull["sigmae"], g0 = parTrueFull["g0"], 
        g0Prior = list(mean=0, var = 4))
    }))
  
  cat("Fitting POUMM on z ...\n")
  fitPOUMM <- POUMM(z[1:N], tree, spec = c(specPOUMM), verbose = TRUE)
  summaryShortPOUMM <- summary(fitPOUMM, mode = "short",
                               startMCMC = specPOUMM$nSamplesMCMC / 2)
  summaryLongPOUMM <- summary(fitPOUMM, mode = "long",
                              startMCMC = specPOUMM$nSamplesMCMC / 2,
                              thinMCMC = 10)
  summaryExpertPOUMM <- summary(fitPOUMM, mode = "expert",
                                startMCMC = specPOUMM$nSamplesMCMC / 2)
  rm(fitPOUMM)
  print(summaryShortPOUMM)
  
  cat("Fitting POUMM on g ...\n")
  fitPOUMMonG <- POUMM(g[1:N], tree, spec = c(specPOUMM), verbose = TRUE)
  summaryShortPOUMMonG <- summary(fitPOUMMonG, mode = "short",
                               startMCMC = specPOUMM$nSamplesMCMC / 2)
  summaryLongPOUMMonG <- summary(fitPOUMMonG, mode = "long",
                              startMCMC = specPOUMM$nSamplesMCMC / 2,
                              thinMCMC = 10)
  summaryExpertPOUMMonG <- summary(fitPOUMMonG, mode = "expert",
                                startMCMC = specPOUMM$nSamplesMCMC / 2)
  rm(fitPOUMMonG)
  print(summaryShortPOUMMonG)
  
  
  cat("Fitting OU on z ...\n")
  fitOU <- POUMM(z[1:N], tree, spec = c(specOU), verbose = TRUE)
  summaryShortOU <- summary(fitOU, mode = "short",
                               startMCMC = specOU$nSamplesMCMC / 2)
  summaryLongOU <- summary(fitOU, mode = "long",
                              startMCMC = specOU$nSamplesMCMC / 2,
                              thinMCMC = 10)
  summaryExpertOU <- summary(fitOU, mode = "expert",
                                startMCMC = specOU$nSamplesMCMC / 2)
  rm(fitOU)
  print(summaryShortOU)
  
  cat("Fitting OU on G ...\n")
  fitOUonG <- POUMM(g[1:N], tree, spec = c(specOU), verbose = TRUE)
  summaryShortOUonG <- summary(fitOUonG, mode = "short",
                            startMCMC = specOU$nSamplesMCMC / 2)
  summaryLongOUonG <- summary(fitOUonG, mode = "long",
                           startMCMC = specOU$nSamplesMCMC / 2,
                           thinMCMC = 10)
  summaryExpertOUonG <- summary(fitOUonG, mode = "expert",
                             startMCMC = specOU$nSamplesMCMC / 2)
  rm(fitOUonG)
  print(summaryShortOUonG)
  
  
  cat("Fitting PMM ...\n")
  fitPMM <- POUMM(z[1:N], tree, spec = specPMM, verbose = TRUE)
  
  summaryShortPMM <- summary(fitPMM, mode = "short",
                             startMCMC = specPMM$nSamplesMCMC / 2)
  summaryLongPMM <- summary(fitPMM, mode = "long",
                            startMCMC = specPMM$nSamplesMCMC / 2,
                            thinMCMC = 10)
  summaryExpertPMM <- summary(fitPMM, mode = "expert",
                              startMCMC = specPMM$nSamplesMCMC / 2)
  rm(fitPMM)
  
  print(summaryShortPMM)
  
  posteriorQuantilesPOUMM <- list()
  posteriorQuantilesPOUMMonG <- list()
  posteriorQuantilesOU <- list()
  posteriorQuantilesOUonG <- list()
  posteriorQuantilesPMM <- list()
  
  for(name in names(parTrueFull)) {
    sample <- as.vector(summaryLongPOUMM[stat == name, mcmc[[1]]])
    posteriorQuantilesPOUMM[[name]] <-
      sum(sample < parTrueFull[name]) / length(sample)
    
    sample <- as.vector(summaryLongPOUMMonG[stat == name, mcmc[[1]]])
    if(name == "sigmae") {
      posteriorQuantilesPOUMMonG[[name]] <-
        sum(sample < 0) / length(sample)
    } else if(name %in% c("H2tMean", "H2e", "H2tInf")) {
      # this was introduced after the execution, so the resulting quantiles
      # are set artificially to 1 in the simStats data.table (see Manuscript.Rmd,
      # chunk simStats)
      posteriorQuantilesPOUMMonG[[name]] <-
        sum(sample < 1) / length(sample)
    } else {
      posteriorQuantilesPOUMMonG[[name]] <-
        sum(sample < parTrueFull[name]) / length(sample)
    }
    
    sample <- as.vector(summaryLongOU[stat == name, mcmc[[1]]])
    posteriorQuantilesOU[[name]] <-
      sum(sample < parTrueFull[name]) / length(sample)
    
    sample <- as.vector(summaryLongOUonG[stat == name, mcmc[[1]]])
    if(name == "sigmae") {
      posteriorQuantilesOUonG[[name]] <-
        sum(sample < 0) / length(sample)
    } else {
      posteriorQuantilesOUonG[[name]] <-
        sum(sample < parTrueFull[name]) / length(sample)
    }
    
    sample <- as.vector(summaryLongPMM[stat == name, mcmc[[1]]])
    posteriorQuantilesPMM[[name]] <-
      sum(sample < parTrueFull[name]) / length(sample)
  }
  
  list(id = id, parTrueFull = parTrueFull, H2eTrue = H2eTrue,
       logLikAtTrueParam = logLikAtTrueParam, 
       g = g, z = z,
       summaryShortPOUMM = summaryShortPOUMM, summaryExpertPOUMM = summaryExpertPOUMM,
       summaryShortPOUMMonG = summaryShortPOUMMonG, summaryExpertPOUMMonG = summaryExpertPOUMMonG,
       summaryShortOU = summaryShortOU, summaryExpertOU = summaryExpertOU,
       summaryShortOUonG = summaryShortOUonG, summaryExpertOUonG = summaryExpertOUonG,
       summaryShortPMM = summaryShortPMM, summaryExpertPMM = summaryExpertPMM,
       posteriorQuantilesPOUMM = posteriorQuantilesPOUMM,
       posteriorQuantilesPOUMMonG = posteriorQuantilesPOUMMonG,
       posteriorQuantilesOU = posteriorQuantilesOU,
       posteriorQuantilesOUonG = posteriorQuantilesOUonG,
       posteriorQuantilesPMM = posteriorQuantilesPMM
  )
}
