#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats pchisq
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats na.fail
#' @importFrom stats optim
#' @importFrom Formula as.Formula
#' @importFrom mvQuad createNIGrid
#' @importFrom mvQuad rescale
#' @importFrom mvQuad getWeights
#' @importFrom mvQuad getNodes
#' @importFrom numDeriv hessian
#' @importFrom randtoolbox halton


.get.theta2.reshape <- function(theta2.in, totalRC, total.demogr.in,
                                indices.in, fill, printLevel) {

  theta2.matrix.out <- matrix(fill, nrow = totalRC, ncol = total.demogr.in +
                                1)

  for (i in 1:length(theta2.in)) {
    theta2.matrix.out[indices.in[i, 1], indices.in[i, 2]] <- theta2.in[i]
  }



  return(theta2.matrix.out)
}



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


.Halton<-function(D, N , randomized){

  if(randomized == F){
    out<- halton( N , D )} else {
      out<- (halton(N, D ) + matrix( rep(runif(D,0,1), N ),nrow=N,byrow = T))%%1
    }
  return(out)
}






  get_integration_input <- function(dim, method,
                                    accuracy, nmkt, seed) {

    if (!missing(seed) && !is.na(seed) )
      set.seed(seed)

    if (!(method %in% c("MLHS","Halton","randHalton","MC","sgGH","sgNH"))){
      stop("Integration method not available! Choose *MLHS*, *Halton*, *randHalton*, *MC*, *sgGH*, *sgNH*.")
    }

    ## Calc. Nodes and Weights ----
    output <- list()
    attributes(output) <- list(method = method,
                               accuracy = accuracy)
    accuracy <- na.fail(accuracy)


    # MLHS
    if (method == "MLHS") {
      nodes <- replicate(nmkt, .MLHS(dim, accuracy))
      nodes_allMkt <- qnorm(t(sapply(1:nmkt,function(x){ c(nodes[,,x]) } )))
      weights <- matrix(1/accuracy, ncol = 1, nrow= accuracy)
    }


    # Halton
    if (method == "Halton") {
      nodes <- replicate(nmkt, .Halton(dim, accuracy, randomized = FALSE))
      nodes_allMkt <- qnorm(t(sapply(1:nmkt,function(x){ c(nodes[,,x]) } )))
      weights <- matrix(1/accuracy, ncol = 1, nrow= accuracy)
    }

    # randomized Halton
    if (method == "randHalton") {
      nodes <- replicate(nmkt, .Halton(dim, accuracy, randomized = TRUE))
      nodes_allMkt <- qnorm(t(sapply(1:nmkt,function(x){ c(nodes[,,x]) } )))
      weights <- matrix(1/accuracy, ncol = 1, nrow= accuracy)
    }

    # MC
    if (method == "MC") {
      nodes_allMkt <- matrix(rnorm(nmkt * accuracy * dim), nrow=nmkt)
      weights <- matrix(1/accuracy, ncol = 1, nrow= accuracy)
    }

    # Gauss Hermite
    if (method == "sgGH") {
      grid <- createNIGrid(dim, type = "GHN", level = accuracy, ndConstruction = "sparse")
      nodes <- c(getNodes(grid))
      nodes_allMkt <- t(replicate(nmkt, nodes))
      weights<-as.matrix(getWeights(grid), ncol = 1)
    }

    # Gauss Hermite nested
    if (method == "sgNH") {
      grid <- createNIGrid(dim, "nHN", accuracy, ndConstruction = "sparse")
      nodes <- c(getNodes(grid))
      nodes_allMkt <- t(replicate(nmkt, nodes))
      weights<-as.matrix(getWeights(grid), ncol = 1)
    }

    # Sparse Grids Trapezoidal
    if (method == "sgTr1") {
      grid <- createNIGrid(dim, 'cNC1', level = accuracy, ndConstruction = "sparse", level.trans = TRUE)
      rescale(grid,cbind(rep(-3.5,dim),rep(+3.5,dim)))
      nodes <- getNodes(grid)
      nodes_allMkt <- t(replicate(nmkt, c(nodes)))
      weights<-as.matrix(getWeights(grid), ncol = 1)
      tmp <- apply(nodes, MARGIN = 1,
                   FUN = function(x){
                     prod(dnorm(x))
                   })
      weights <- weights * tmp
    }


    if (method == "sgTr") {
      grid <- createNIGrid(dim, 'Trapez', level = accuracy, ndConstruction = "sparse", level.trans = TRUE)
      rescale(grid,cbind(rep(-3.5,dim),rep(+3.5,dim)))
      nodes <- getNodes(grid)
      nodes_allMkt <- t(replicate(nmkt, c(nodes)))
      weights<-as.matrix(getWeights(grid), ncol = 1)
      tmp <- apply(nodes, MARGIN = 1,
                   FUN = function(x){
                     prod(dnorm(x))
                   })
      weights <- weights * tmp
    }


    output$nodesMktShape <- nodes_allMkt
    output$weights <- weights

    return(output)
  }


  .draws_listToMatrix <- function( drawList, amountDraws,
                                 market_identifier_pd ,
                                 market_identifier_list_name,
                                 use ){
    list_names <- names( drawList )

    # checks
    drawsInList<- unlist(
      lapply( drawList,
              function(i){
                draws <- ncol(i) -1 # minus 1 because of market id col.
                return(draws)
              }))

    if( !all( drawsInList ==  amountDraws ) ){
      cat( paste(use ,":") )
      stop("Number of draws for at least one list entry is
           smaller as the provided integration accuracy.
           Include draws accordingly.")
    }

    # extract data
    tmp<- lapply( list_names,
                  function(i){
                    drawset_i <- na.fail( as.data.frame( get( i , drawList ) ) )
                    marketid_list_i <- get( market_identifier_list_name, drawset_i)

                    if( any( table(marketid_list_i) >1) ){
                      cat( paste(use ,":") )
                      stop("Draws in one list entry are not unique for at least one market.")
                    }

                    # reorder according to ordered markets
                    pD_order <- match( unique(market_identifier_pd),
                                       marketid_list_i)
                    if( any( is.na( pD_order ) ) ){
                      cat( paste(use ,":") )
                      stop("Draws in one list entry are not available for at least one market.")
                    }

                    drawset_i[[market_identifier_list_name]] <- NULL
                    tmpList <- as.matrix( drawset_i[ pD_order , 1:amountDraws ]) # now, tmpDemogr.data has roworder as pD
                    colnames(tmpList) <- paste0(use ,"_draw",1:amountDraws,"_", i )

                    if( ! all( is.finite( tmpList ) ) ){
                      cat( paste(use ,":") )
                      stop("Draws contain non-numeric values.")
                    }


                    return(tmpList)
                  } )



    # Demographic data is arranged in a matrix, with the draws of
    # different variables next to each other :
    D <- do.call(cbind,tmp)
    return(D)

  }





  .indexing.markets <- function( cdidOld ) {
    cdidOld <- as.character(cdidOld)
    unique.ids <-  unique(cdidOld)
    nmkt <- length(unique.ids)
    cdid <- numeric(length(cdidOld))

    for(i in 1:nmkt ){
      relevantMarkets <- cdidOld == unique.ids[i]
      cdid[ relevantMarkets ] <- i

    }
    return( cdid )
  }




  .prepare_theta2 <- function(par_theta2,
                              final_col_names_par, final_row_names_par ,
                              K , M){
    if( missing(par_theta2) ){
      par_theta2 <- matrix( 0, nrow = K, ncol = 1 + M )
    } else{


      if( !setequal( colnames(par_theta2) , final_col_names_par ))
        stop("Colnames of par_theta2 do not match with names of obs. and unobs. heterogeneity. Remember to name column of unobs. heterogeneity \"unobs_sd\" .\n")

      if( !setequal( rownames(par_theta2) , final_row_names_par ))
        stop("Rownames of par_theta2 do not match with random coefficients. Remember to name any constant \"(Intercept)\" .\n")

      #reorder accoring to demographic_names and random coefficients from formula
      col_order <- match( final_col_names_par, colnames(par_theta2) )
      row_order <- match( final_row_names_par, rownames(par_theta2) )

      par_theta2 <- par_theta2[ row_order , col_order ,drop=FALSE]
    }

    rownames(par_theta2) <- final_row_names_par
    colnames(par_theta2) <- final_col_names_par

    # indices for saving col and row structure
    indices <- which( !is.na( par_theta2 ),
                      arr.ind = TRUE ) # NA's are excluded from the estimation algorithm

    par_theta2 <- na.omit( c( par_theta2) )
    total_par <- length(par_theta2)
    if( any (!is.finite( par_theta2 ) ))
      stop( "Provide finite starting guesses for par_theta2.")


    out <- list(par_theta2=par_theta2,
                indices=indices,
                total_par =total_par,
                final_row_names_par= final_row_names_par,
                final_col_names_par=final_col_names_par
    )


    return( out )



  }


