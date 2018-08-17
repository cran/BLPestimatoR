#' @useDynLib BLPestimatoR
#' @importFrom Rcpp sourceCpp
NULL


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
#'
gmm_obj <- function(par_theta2 , indices, blp_results, blp_data,
                    printLevel){

  if(class(blp_data) != "blp_data")
    stop("Input has wrong class. Call BLP_data() first.")



  deltaOld <- blp_results$deltaOld # delta vector from previous run (without the blp.results environment)

  theta2Mat<- .get.theta2.reshape(theta2.in       = par_theta2,
                                  totalRC         = blp_data$parameters$K,
                                  total.demogr.in = blp_data$parameters$total_demogr,
                                  indices.in      = indices,
                                  fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu

  #Call C++ function:
  tmp <- getDelta(  theta2    = theta2Mat,
                    cdid      = blp_data$parameters$market_id,
                    cdindex   = blp_data$parameters$cdindex,
                    innerCrit = blp_data$parameters$inner_tol,
                    indices   = indices,
                    innerMaxit= blp_data$parameters$inner_maxit,
                    Xrandom          = blp_data$data$X_rand,
                    obsshare         = blp_data$data$shares,
                    deltaOld         = blp_results$deltaOld,
                    nodesDemMktShape = blp_data$integration$drawsDemMktShape ,
                    nodesRcMktShape  = blp_data$integration$drawsRcMktShape,
                    weights          = blp_data$integration$weights,
                    printLevel = printLevel)


  delta <- tmp$delta
  counter <- tmp$counter

  bet <- NA
  xi <- NA
  sij<- NA
  jacobian <- NA
  gradient <- NA
  delta_has_na <- any(is.na(delta))

  if (delta_has_na) {
    gradient <- rep(Inf, length(par_theta2))
    f <- Inf
  } else {
    #save delta as start solution for the next run,
    #only if contraction mapping converged in the previous step:
    blp_results$deltaOld <- delta
    ## Objective
    bet <- blp_data$data$invxzwzx %*% (blp_data$data$xzwz %*% delta)
    xi <- delta - blp_data$data$X_lin %*% bet
    tmp2 <- t(xi) %*% blp_data$data$Z
    f <- c(tmp2 %*% blp_data$data$W %*% t(tmp2))

    ## Gradient
    sij <- matrix(NA_real_, nrow = blp_data$parameters$nobs ,
                  ncol = blp_data$integration$amountDraws )
    sij[ , ] <- tmp$sij
    jacobian <- jacob_c(sij = sij,
                        indices   = indices,
                        blp_data = blp_data$data,
                        blp_parameters = blp_data$parameters,
                        blp_integration = blp_data$integration,
                        printLevel = printLevel)

    if (any(is.na(jacobian))) {
      if (printLevel > 0) { cat("\t gradient contains Na's --> objective value and gradient replaced by Inf \n")}
      gradient <- rep(Inf, length(par_theta2))
      f <- Inf
    }else{
      gradient <- 2 * t(jacobian) %*% blp_data$data$Z %*% blp_data$data$W %*%
        t(blp_data$data$Z) %*% xi
    } #is.na(jacobian)

  } #end is.na(delta)



  if (printLevel >= 1) {
    cat("gmm objective:", round(f, 4))
    if ( delta_has_na )  cat(" [delta contains NaN's] ")
    cat("\n")
  }
  if (printLevel >= 2) {
    cat("\t theta (RC): ")
    cat( round(theta2Mat[ ,1] , 2), "\n")
    if( length(par_theta2) >0 ){
      cat("\t theta (demogr.): ")
      cat( round(c(theta2Mat[ ,-1]) , 2), "\n")
    }
    cat("\t inner iterations: ")
    cat( counter, "\n")
    cat("\t gradient: " )
    cat( round(gradient,4), "\n")
  }

  # save to environment (not provided as return object, because optim only accepts a single number as output)
  blp_results$bet <- bet
  blp_results$gradient <- gradient
  blp_results$jacobian <- jacobian
  blp_results$xi <- xi
  blp_results$innerItAll <- c(blp_results$innerItAll, tmp$counter)
  blp_results$sij <- sij
  if(tmp$negShares==TRUE) blp_results$negShares <- TRUE

  return(f)
}


gmm_gr <- function(par_theta2 , indices, blp_results, blp_data,
                   printLevel) {

  return(blp_results$gradient)

}


