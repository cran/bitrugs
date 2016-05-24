
traceplots <- function(mcmcoutput, tracecol="lightgrey", labels=c("p","z","beta","gamma","gamma_G","genpar"), burn=1, filter=1) {
    
  layout(rbind(1:3,4:6))
  for (i in 1:length(labels)) {
    stream <- mcmcoutput[burn:nrow(mcmcoutput),i]
    stream <- stream[seq(from=1, to=length(stream), by=filter)]
    plot(NULL, xlim=c(0,length(stream)), ylim=range(stream), xlab="Iterations", ylab=labels[i])
    lines(1:length(stream), stream, col=tracecol)
    abline(h=quantile(mcmcoutput[,i], c(0.5,0.025,0.975)), lty=c(1,2,2), col="red")
  }

}


plot_transnetwork <- function(mcmcoutput, epidata=NULL, type=c(1,2,3), plotthresh=0.05, 
                              labels="ID", text.cex=1, adj=0.25, ID=NULL, n=NULL) {

  if (is.null(n)) {
    n <- (ncol(mcmcoutput)-10)/2
  }
  
  if (is.null(ID)) {
    ID <- sample(1:n,n)
  }
  
  namevec <- c("Inferred", "Naive", "True")
  if (!is.null(epidata)) {
    infsource <- epidata[,5]
  }
  
  if (is.null(infsource) && 3%in%type) {
    stop("To plot the true network (type=3), please specify the true infection routes (infsource)")
  }
  
  if (1%in%type) {
    iterations <- nrow(mcmcoutput)
    # After column 10, output represents source of infection of reach host
    colsource <- 10 # IN WHICH COLUMN IN THE OUTPUT DO THE INFERRED SOURCES BEGIN?
    colgroup <- colsource+n
    colprob <- numeric(n)
    for (i in 1:n) {
      colprob[i] <- sum(mcmcoutput[,colsource+i]!=0)/iterations # posterior prob of each patient being colonized
    }
    hicolprob <- colprob[which(colprob==1)] # Those probabilities which are high enough
    hiprobID <- which(colprob==1) # IDs for those patients
  
    hiprobno <- length(hiprobID) # No. patients to plot
    #pIDsamp <- sample(hiprobID,hiprobno) # Randomize the order of patients
    pIDsamp <- ID[which(ID%in%hiprobID)] # Randomize the order of patients
    n_inf <- length(pIDsamp)
    net_inf <- matrix(0,n,n)
    for (i in 1:n) {
      for (j in 1:n) {
        net_inf[j,i] <- sum(mcmcoutput[,colsource+i]==j)/iterations
      }
    }
    
  } else {
    pIDsamp <- sample(ID[which(infsource!=0)],sum(infsource!=0))
    pIDsamp <- ID[which(ID%in%which(infsource!=0))]
    n_inf <- length(pIDsamp)
  }

  if (2%in%type) {
    net_nv <- matrix(0,n,n)
    for (i in 1:n) {
      if (infsource[i]!=0) {
        colpos <- sum(epidata[,2]<=epidata[i,4] & epidata[,3]>=epidata[i,4] & epidata[,4]>0)+1
        net_nv[which(epidata[,2]<=epidata[i,4] & epidata[,3]>=epidata[i,4] & epidata[,4]>0),i] <- 1/colpos
      }
    }
  }
  
  if (3%in%type) {
    positives <- which(infsource!=0)
    npositives <-length(positives)
    net_true <- matrix(0,n,n)
    infsource[which(infsource==-1)] <- which(infsource==-1)
    for (i in 1:n) {
      if (infsource[i]>=1){
        net_true[infsource[i],i]=1
      } else if (infsource[i]==-1) {
        net_true[i,i] <- 1
      }
    }
  }

  
  layout(t(1:length(type)))
  
  for (k in 1:length(type)) {
    plot(NULL, xlim=c(-3.5,3.5), ylim=c(-3.5,3.5), yaxt="n", xaxt="n", bty="n",
         xlab="", ylab="", main=paste(namevec[type[k]], "transmission network"))
    if (type[k]==1) {
      net <- net_inf
    } else if (type[k]==2) {
      net <- net_nv
    } else {
      net <- net_true
    }
    for (i in 1:n_inf) {
      points(2.3*cos(2*pi*i/n_inf), 2.3*sin(2*pi*i/n_inf), 
             cex=sqrt(sum(net[,pIDsamp[i]])), pch=16, col=rgb(1-net[pIDsamp[i],pIDsamp[i]],0,0))
      if (labels=="ID") {
        text((2.5+adj)*cos(2*pi*i/n_inf), (2.5+adj)*sin(2*pi*i/n_inf), 
             pIDsamp[i], cex=text.cex, srt=(360*i/n_inf)+180*ifelse(i/n_inf>0.25&&i/n_inf<0.75,1,0))
      } else if (labels=="import") {
        text((2.5+adj)*cos(2*pi*i/n_inf), (2.5+adj)*sin(2*pi*i/n_inf), 
             format(round(net[pIDsamp[i],pIDsamp[i]],2), nsmall=2), cex=text.cex, 
             srt=(360*i/n_inf)+180*ifelse(i/n_inf>0.25&&i/n_inf<0.75,1,0))
      } else if (labels=="secondary") {
        text((2.5+adj)*cos(2*pi*i/n_inf), (2.5+adj)*sin(2*pi*i/n_inf), 
             format(round(sum(net[pIDsamp[i],-pIDsamp[i]]),2), nsmall=2), cex=text.cex, 
             srt=(360*i/n_inf)+180*ifelse(i/n_inf>0.25&&i/n_inf<0.75,1,0))        
      } else {
        stop("'labels' must be either 'ID', 'import' or 'secondary'")
      }
    }
    for (i in 1:n_inf) {
      for (j in 1:n_inf) {
        if (i!=j) {
          man <- net[pIDsamp[j],pIDsamp[i]]
          if (man > plotthresh) {
            arrows(2.3*cos(2*pi*j/n_inf),2.3*sin(2*pi*j/n_inf),
                   2.3*cos(2*pi*i/n_inf),2.3*sin(2*pi*i/n_inf), 
                   col=rgb(0,1-man,man,man), length=0.1, lwd=sqrt(man)*3)
          }
        }
      } 
    }
  }
}
