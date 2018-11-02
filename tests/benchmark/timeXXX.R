
getModeStr <- function(mode) {
  if(mode == 11) "serial"
  else if(mode == 21) "parallel range"
  else if(mode == 24) "parallel queue"
  else if(mode == 25) "parallel range no except"
  else if(mode == 31) "hybrid range"
  else stop("mode not supported");
}

timeRphylopars <- function(X, tree) {
  cat("timeRphylopars... ")
  N <- length(tree$tip.label)
  Z <- X
  data <- cbind(data.table(species = tree$tip.label), as.data.table(t(Z)[1:N,]))
  time <- system.time(output <- capture.output(
    p_OU <- phylopars(trait_data = data, tree = tree, model="OU", optim_verbose = TRUE)
  )
  )[3]*1000
  #print(output)
  timeLik <- mean(sapply(output, function(o) {
    as.double(strsplit(o, " ")[[1]][3])
  }))
  cat(timeLik, "ms\n")

  timeExtra <- if(is.na(timeLik)) {
    as.double(NA)
  } else {
    (time - length(output)*timeLik) / length(output)
  }

  data.table(impl = paste0("Rphylopars: C++, 3-point"), time = timeLik,
             timeExtra = timeExtra, mode=11)
}

timePCMBaseR <- function(X, tree, model) {
  cat("timePCMBaseR... ")

  Z <- X[,1:PCMTreeNumTips(tree), drop=FALSE]
  timeExtra <- system.time(
    {
      metaI <- PCMInfo(X = Z, tree = tree, model = model)
    })[3]*1000

  timeLik <-
    system.time(
      for(i in 1:nTests)
        PCMBase::PCMLik(Z[, 1:metaI$N,drop=FALSE], tree, model, metaI= metaI)
    )[3]*1000/nTests

  cat(timeLik, "ms\n")

  data.table(impl = paste0("PCMBase: R, Poly"), time = timeLik, timeExtra = timeExtra, mode=11)
}

timePCMBaseCpp <- function(X, tree, model, mode = 0) {
  cat("timePCMBaseCpp... ")

  options("PCMBase.Lmr.mode" = mode)
  Z <- X[,1:PCMTreeNumTips(tree), drop=FALSE]
  timeExtra <- system.time(
    {
      metaI <- PCMInfoCpp(X = Z, tree = tree, model = model)
    })[3]*1000

  # warm-up
  for(i in 1:nTests)
    PCMBase::PCMLik(Z, tree, model, metaI= metaI)


  timeLik <-
    system.time(
      for(i in 1:nTests)
        PCMBase::PCMLik(Z, tree, model, metaI= metaI)
    )[3]*1000/nTests

  cat(timeLik, "ms\n")

  data.table(impl = paste0("PCMBase: C++, Poly, ", getModeStr(mode)),
             time = timeLik,
             timeExtra = timeExtra,
             mode=mode)
}

timePOUMMR <- function(X, tree, model) {
  cat("timePOUMMR... ")

  N <- length(tree$tip.label)
  Z <- X[,1:N]
  if(is.matrix(Z)) {
    list(impl = paste0("POUMM: R, Poly"), time = as.double(NA), timeExtra = as.double(NA))
    cat(NA, "ms\n")
  } else {
    timeExtra <- system.time(
      {
        pruneI <- POUMM:::pruneTree(tree, Z)
      }
    )[3]*1000

    timeLik <- system.time(
      for(i in 1:nTests) POUMM::likPOUMMGivenTreeVTips(Z, tree, alpha = model$H[1,1,1], theta=model$Theta[1,1], sigma = model$Sigma_x[1,1,1]^2, sigmae = model$Sigmae_x[1,1,1]^2, g0=model$X0[1],  pruneInfo = pruneI)
    )[3]*1000/nTests
    cat(timeLik, "ms\n")

    data.table(impl = paste0("POUMM: R, Poly"), time = timeLik, timeExtra = timeExtra, mode=11)
  }
}

timePOUMMCpp <- function(X, tree, model, mode) {
  options(SPLITT.postorder.mode = mode)
  #options(splittree.postorder.mode = mode)
  N <- length(tree$tip.label)
  Z <- X[,1:N]
  if(is.matrix(Z)) {
    data.table(impl = paste0("POUMM: C++, Poly, ", getModeStr(mode)), time = as.double(NA), timeExtra = as.double(NA))
  } else {
    timeExtra <- system.time(
      {
        pruneI <- POUMM:::pruneTree(tree, Z)
      }
    )[3]*1000

    alpha = model$H[1,1,1]
    theta=model$Theta[1,1]
    sigma = model$Sigma_x[1,1,1]^2
    sigmae = model$Sigmae_x[1,1,1]^2
    g0=model$X0[1]

    # warm-up
    for(i in 1:(10*nTests)) POUMM::likPOUMMGivenTreeVTipsC(integrator = pruneI$integrator,
                                                           alpha = alpha,
                                                           theta = theta,
                                                           sigma = sigma,
                                                           sigmae = sigmae,
                                                           g0 = g0)


    timeLik <- system.time(
      for(i in 1:(10*nTests)) POUMM::likPOUMMGivenTreeVTipsC(integrator = pruneI$integrator,
                                                             alpha = alpha,
                                                             theta = theta,
                                                             sigma = sigma,
                                                             sigmae = sigmae,
                                                             g0 = g0)
    )[3]*1000/(10*nTests)

    data.table(impl = paste0("POUMM: C++, Poly, ", getModeStr(mode)), time = timeLik, timeExtra = timeExtra, mode=mode)
  }
}
