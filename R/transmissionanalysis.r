

transmission_analysis <- function(epidata, distmat, seqIDs, resmat, iterations=100000, augmoves=10, feedb=1, 
                                  tag=1, noaug=1, model=2, path=NULL, sensprior=c(1,1), impprior=c(1,1), 
                                  betaprior=1E6, gammaprior=c(1,1), gammaGprior=c(1,1), tparprior=1E6, 
                                  sigma=c(0.004, 0.03, 0.005, 0.25)) {

  n <- nrow(epidata) # no. patients
  if (is.null(path)) {
    stop("Please specify a path to which output can be saved.")
  } else if (substring(path,first=nchar(path))!="/") {
    path <- paste(path, "/", sep="")
  }
  if (!file.exists(path)) {
    stop("Path does not exist")
  }
  
  pID <- epidata[,1] # patient ID
  adm <- epidata[,2] # day of admission
  dis <- epidata[,3] # day of discharge
  col_t <- numeric(n) # day of infection (0 if never infected)
  group <- numeric(n) # infection group (0 if never infected)
  psource <- numeric(n) # source of infection (ID of infector, 0 if never infected)
  
  GD <- distmat # genetic distance matrix
  K <- cbind(1:n, resmat) # Tack on pat IDs to resmat
  res <- as.numeric(t(K)) # convert results matrix to vector, to pass to C
  noseqs <- length(seqIDs) # number of sequences
  maxD <- max(dis) # final day of study
  
  # Tsim object contains true infection times, routes. We ignore these and set 
  # them to default values (the initial tree for the MCMC), allowing the algorithm 
  # to recover the tree based upon data which would be observed in reality.
  for (i in 1:n) {
    if (1%in%resmat[i,]) { # if positive result
      col_t[i] <- adm[i] # set time of infection to admission day (initial settings)
      group[i] <- i # set each host as its own group (initial settings)
      psource[i] <- i # set each host to be importation (their own source of infection - initial settings)
    } else { # if no positive result, assume not infected
      col_t[i] <- 0
      group[i] <- 0
      psource[i] <- 0
    }
  }
  
  parnames <- c("p", "z", "beta", "gamma", "gamma_G", "tpar") # Names of parameters
  # Prior distribution parameters
  priors <- c(impprior, sensprior, betaprior, gammaprior, gammaGprior, tparprior)
  
  ########################
  # LOAD & RUN MCMC CODE #
  ########################

  WGSmcmc <- .C("network_mcmc", as.integer(n), as.integer(pID), as.integer(adm), 
                as.integer(dis), as.integer(col_t), as.integer(psource), as.integer(group),
                as.integer(res), as.integer(model), as.integer(noseqs), as.integer(seqIDs), as.integer(GD), 
                as.integer(iterations), as.integer(augmoves), as.integer(maxD), as.double(sigma),
                as.double(priors), as.integer(tag), as.character(path), 
                as.double(feedb), as.integer(noaug))
  
  mcmcoutput <- read.table(paste(path, "WGS_mcmc", tag, ".txt", sep=""))
  colnames(mcmcoutput) <- c(parnames, "importations", "acquisitions", "n_groups", "loglikelihood",
                            paste("source", 1:n, sep=""), paste("group", 1:n, sep=""))
  
  return(mcmcoutput)
}
