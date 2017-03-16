###########################
.indexing.markets <- function( cdidOld ) {
  unique.ids <- unique(cdidOld)
  nmkt <- length(unique.ids)
  cdid <- numeric(length(cdidOld))
  cdindex <- numeric(nmkt+1) #;% gives the last product position in the market
  cdindex[1]<- 0
  for(i in 1:nmkt ){
    relevantMarkets <- cdidOld == unique.ids[i]
    cdid[ relevantMarkets ] <- i
    cdindex[i+1]<- max(which(cdid==i))
  }
  return( list("cdid" = cdid,
               "cdindex" = cdindex,
               "nmkt" = nmkt))
}


###########################
.get.theta2.reshape <- function(theta2.in, totalRC, total.demogr.in,
                                indices.in, fill, printLevel) {

  theta2.matrix.out <- matrix(fill, nrow = totalRC, ncol = total.demogr.in +
                                1)

  for (i in 1:length(theta2.in)) {
    theta2.matrix.out[indices.in[i, 1], indices.in[i, 2]] <- theta2.in[i]
  }



  return(theta2.matrix.out)
}

############################

.MLHS <- function(D, N) {

  draws  <- numeric(N)
  shuffleddraws <- matrix(NA, nrow=N, ncol=D)

  for (i in 1:D) {
    draws <-  ((1:N)-1)/N + runif(1)/N;
    shuffle = sample(N)
    shuffleddraws[ , i] <- draws[shuffle]

  }
  return( shuffleddraws )
}
############################

.Halton<-function(D, N , randomized){

  if(randomized == F){
    out<- halton( N , D )} else {
      out<- (halton(N, D ) + matrix( rep(runif(D,0,1), N ),nrow=N,byrow = T))%%1
    }
  return(out)
}
