
simulate_data <- function(n,D,LOS=7,p=0.05,z=0.8,b=0.005,gamma=0.3,gamma_gl=0.03, 
                           genpar=0.8, testdays=3, model=2) {
  
  print(paste("Simulating", n, "patients over", D, "days"))
  # set up
  day_adm <- numeric(n)
  day_dis <- numeric(n)
  col_t <- numeric(n)
  psource <- numeric(n)
  group <- numeric(n)
  pID <- 1:n
  
  day_adm <- sort(sample(1:D,n,replace=T))
  day_dis <- day_adm + rpois(n,LOS)
  col_t <- rbinom(n,1,p)
  psource <- -1*col_t
  col_t <- col_t*day_adm
  
  maxD <- max(day_dis)
  testdays <- seq(1,maxD,testdays)
  
  print("Generating importations...")
  
  imports <- which(col_t!=0)
  if (model==1) {
    for (i in imports) {
      if (runif(1,0,1)<genpar || sum(group!=0)==0) {
        group[i] <- pID[i]
      } else {
        group[i] <- group[sample(imports[which(col_t[imports]<=col_t[i])],1)]
      }
      if (group[i] == 0) {
        group[i] = pID[i]
      }
    }
  } else {
    for (i in imports) {
      group[i] <- pID[i]
    }
  }
  
  cx <- 0
  for (j in which(psource==-1)) { # which patients are pos on adm
    cx <- cx+1
    cat(cx, " ")
    t <- day_adm[j]   
  }
  cat("\n")
  
  if (sum(psource==-1)==0) {
    stop("No infections generated!")
  }
  
  # calculate initial colonised population vector
  colp <- numeric(maxD)
  pop <- numeric(maxD)
  for (t in 1:maxD) {
    for (j in 1:n) {
      if (day_adm[j]<=t && day_dis[j] >= t && col_t[j] != 0) {
        if (col_t[j] < t || (psource[j]==-1 && col_t[j] == t)) {
          colp[t] <- colp[t]+1
        }
      }
      if (day_adm[j]<=t && day_dis[j] >= t) {
        pop[t] <- pop[t]+1
      }
    }
  }
  colcan <- colp
  # Add acquisitions
  print("Generating acquisitions...")
  cx <- 0
  for (t in 1:maxD) {
    colp <- colcan
    for (j in 1:n) {
      if (day_adm[j] <= t && day_dis[j] >= t) {
        if (col_t[j] == 0) {
          if (runif(1,0,1)<(1-exp(-b*colp[t]))) { # acquisition
            # choose source
            cx <- cx+1
            cat(cx, " ")
            src <- sample(colp[t],1)
            k <- 0
            for (m in 1:n) {
              if (day_adm[m]<=t && day_dis[m] >= t && col_t[m] !=0) {
                if (col_t[m] < t || (psource[m]==-1 && col_t[m] == t)) {
                  k <- k+1
                }
              }
              if (k==src) { # mth patient is infector
                psource[j] <- pID[m]
                group[j] <- group[m]
                break
              }
            }
            col_t[j] <- t
            colcan[col_t[j]:day_dis[j]] <- colp[col_t[j]:day_dis[j]]+rep(1,day_dis[j]-col_t[j]+1)
          }
        }
      }
    }
  }
  cat("\n")
  #Add observations
  resultsmatrix <- matrix(-1,n,maxD)
  print("Generating observed results...")
  for (t in testdays) {
    for (j in 1:n) {
      if (day_adm[j] <= t && day_dis[j] >= t) {
        if (col_t[j] <= t && col_t[j]!=0 && runif(1,0,1)<z) {
          resultsmatrix[j,t] <- 1
        } else {
          resultsmatrix[j,t] <- 0
        }
      }
    }
  }
  
  positives <- which(psource!=0)
  npositives <- length(positives)
  positives <- which(psource!=0)
  npositives <-length(positives)
  nsequences <-sum(resultsmatrix==1)
  
  distmat <- matrix(0,nsequences,nsequences)
  patientseqIDs <- numeric(nsequences)
  mark<-0
  
  k<-1
  for (i in which(col_t>0)) {
    for (t in day_adm[i]:day_dis[i]) {
      if (resultsmatrix[i,t]==1) {
        patientseqIDs[k] <- pID[i]
        k <- k+1
      }
    }
  }
  # Pairwise distance matrix
  print("Generating genetic distance matrix...")
  
  for (i in 2:nsequences) {
    for (j in 1:(i-1)) {
      if (model==1) {
        if (group[patientseqIDs[i]]==group[patientseqIDs[j]]) {
          distmat[i,j] <- rgeom(1,gamma)
          distmat[j,i] <- distmat[i,j]
        } else {
          distmat[i,j] <- rgeom(1,gamma_gl)
          distmat[j,i] <- distmat[i,j]
        }
      } else {
        if (translinks(patientseqIDs[i], patientseqIDs[j], pID, psource)>=0) {
          distmat[i,j] <- rgeom(1,gamma*genpar^translinks(patientseqIDs[i], patientseqIDs[j], pID, psource))
          distmat[j,i] <- distmat[i,j]
        } else {
          distmat[i,j] <- rgeom(1,gamma_gl)
          distmat[j,i] <- distmat[i,j]
        }
      }
    }
  }
  
  epi <- cbind(pID,day_adm,day_dis,col_t, psource,group)
  gen <- cbind(patientseqIDs,distmat)
  
  return(invisible(list(epi=epi,
                        resmat=resultsmatrix, distmat=distmat, 
                        patientseqIDs=patientseqIDs)))
}


translinks <- function(p1, p2, patID, psource) {
  ch1 <- p1
  ch2 <- p2
  if (ch1!=ch2) {
    cur <- p1
    while (1) {
      cur <- psource[which(patID==cur)]
      if (cur==-1) {
        break
      } else {
        ch1 <- c(ch1, cur)
      }
    }
    
    cur <- p2
    while (1) {
      cur <- psource[which(patID==cur)]
      if (cur==-1) {
        break
      } else {
        ch2 <- c(ch2, cur)
      }
    }
    if (sum(ch1%in%ch2)==0) {
      gen <- -1
    } else {
      p1 <- which(ch1%in%ch2)[1]-1
      p2 <- which(ch2%in%ch1)[1]-1
      gen <- p1+p2
    }
  } else {
    gen <- 0
  }
  
  return(gen)
}

simulate_data_dates <- function(day_adm,day_dis,p=0.05,z=0.8,b=0.005,
                                gamma=0.3,gamma_gl=0.03, genpar=0.8, testdays=3, model=2) {

  n <- length(day_adm)
  col_t <- numeric(n)
  psource <- numeric(n)
  group <- numeric(n)
  pID <- 1:n
  col_t <- rbinom(n,1,p)
  psource <- -1*col_t
  col_t <- col_t*day_adm
  
  maxD <- max(day_dis)
  testdays <- seq(1,maxD,testdays)
  
  print("Generating importations...")
  
  imports <- which(col_t!=0)
  if (model==1) {
    for (i in imports) {
      if (runif(1,0,1)<genpar || sum(group!=0)==0) {
        group[i] <- pID[i]
      } else {
        group[i] <- group[sample(imports[which(col_t[imports]<=col_t[i])],1)]
      }
      if (group[i] == 0) {
        group[i] = pID[i]
      }
    }
  } else {
    for (i in imports) {
      group[i] <- pID[i]
    }
  }
  
  cx <- 0
  for (j in which(psource==-1)) { # which patients are pos on adm
    cx <- cx+1
    cat(cx, " ")
    t <- day_adm[j]   
  }
  cat("\n")
  
  # calculate initial colonised population vector
  colp <- numeric(maxD)
  pop <- numeric(maxD)
  for (t in 1:maxD) {
    for (j in 1:n) {
      if (day_adm[j]<=t && day_dis[j] >= t && col_t[j] != 0) {
        if (col_t[j] < t || (psource[j]==-1 && col_t[j] == t)) {
          colp[t] <- colp[t]+1
        }
      }
      if (day_adm[j]<=t && day_dis[j] >= t) {
        pop[t] <- pop[t]+1
      }
    }
  }
  colcan <- colp
  # Add acquisitions
  print("Generating acquisitions...")
  cx <- 0
  for (t in 1:maxD) {
    colp <- colcan
    for (j in 1:n) {
      if (day_adm[j] <= t && day_dis[j] >= t) {
        if (col_t[j] == 0) {
          if (runif(1,0,1)<(1-exp(-b*colp[t]))) { # acquisition
            # choose source
            cx <- cx+1
            cat(cx, " ")
            src <- sample(colp[t],1)
            k <- 0
            for (m in 1:n) {
              if (day_adm[m]<=t && day_dis[m] >= t && col_t[m] !=0) {
                if (col_t[m] < t || (psource[m]==-1 && col_t[m] == t)) {
                  k <- k+1
                }
              }
              if (k==src) { # mth patient is infector
                psource[j] <- pID[m]
                group[j] <- group[m]
                break
              }
            }
            col_t[j] <- t
            colcan[col_t[j]:day_dis[j]] <- colp[col_t[j]:day_dis[j]]+rep(1,day_dis[j]-col_t[j]+1)
          }
        }
      }
    }
  }
  cat("\n")
  #Add observations
  resultsmatrix <- matrix(-1,n,maxD)
  print("Generating observed results...")
  for (t in testdays) {
    for (j in 1:n) {
      if (day_adm[j] <= t && day_dis[j] >= t) {
        if (col_t[j] <= t && col_t[j]!=0 && runif(1,0,1)<z) {
          resultsmatrix[j,t] <- 1
        } else {
          resultsmatrix[j,t] <- 0
        }
      }
    }
  }
  
  positives <- which(psource!=0)
  npositives <- length(positives)
  positives <- which(psource!=0)
  npositives <-length(positives)
  nsequences <-sum(resultsmatrix==1)
  
  distmat <- matrix(0,nsequences,nsequences)
  patientseqIDs <- numeric(nsequences)
  mark<-0
  
  k<-1
  for (i in which(col_t>0)) {
    for (t in day_adm[i]:day_dis[i]) {
      if (resultsmatrix[i,t]==1) {
        patientseqIDs[k] <- pID[i]
        k <- k+1
      }
    }
  }
  # Pairwise distance matrix
  print("Generating genetic distance matrix...")
  
  for (i in 2:nsequences) {
    for (j in 1:(i-1)) {
      if (model==1) {
        if (group[patientseqIDs[i]]==group[patientseqIDs[j]]) {
          distmat[i,j] <- rgeom(1,gamma)
          distmat[j,i] <- distmat[i,j]
        } else {
          distmat[i,j] <- rgeom(1,gamma_gl)
          distmat[j,i] <- distmat[i,j]
        }
      } else {
        if (translinks(patientseqIDs[i], patientseqIDs[j], pID, psource)>=0) {
          distmat[i,j] <- rgeom(1,gamma*genpar^translinks(patientseqIDs[i], patientseqIDs[j], pID, psource))
          distmat[j,i] <- distmat[i,j]
        } else {
          distmat[i,j] <- rgeom(1,gamma_gl)
          distmat[j,i] <- distmat[i,j]
        }
      }
    }
  }
  
  epi <- cbind(pID,day_adm,day_dis,col_t, psource,group)
  gen <- cbind(patientseqIDs,distmat)
  
  return(invisible(list(epi=epi,
                        resmat=resultsmatrix, distmat=distmat, 
                        patientseqIDs=patientseqIDs)))
}