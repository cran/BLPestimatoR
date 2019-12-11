#' @useDynLib BLPestimatoR
#' @importFrom Rcpp sourceCpp
NULL

#' Calculates derivatives of all shares with respect to all non-linear parameters in a given market.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param par_theta2 matrix with column and rownames providing a starting value for the optimization routine (see details),
#' @param market character specifying the market in which derivatives are calculated,
#' @param printLevel level of output information (default = 1)
#'
#' @return Returns a numeric matrix with derivatives.
#' Cell in row i and col j is the derivative of share i with respect to parameter j.
#'
#' @details NA's in \code{par_theta2} entries indicate the exclusion from estimation, i.e. the coefficient is assumed to be zero.
#' If only unobserved heterogeneity is used (no demographics), the column name of \code{par_theta2} must be "unobs_sd".
#' With demographics the colnames must match the names of provided demographics (as in \code{demographic_draws}) and "unobs_sd".
#' Row names of \code{par_theta2} must match random coefficients as specified in \code{model}. Constants must be named "(Intercept)".
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                         Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                         Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                         Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                         true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                Xlin.true.price = -0.2,
#'                                                Xrandom.true = rep(2,K),
#'                                                instrument.effects = rep(2,10),
#'                                                instrument.Xexo.effects = rep(1,5)),
#'                         price.endogeneity = list( mean.xi = -2,
#'                                                   mean.eita = 0,
#'                                                   cov = cbind( c(1,0.7), c(0.7,1))),
#'                         printlevel = 0, seed = 234234 )
#'
#'
#' model <- as.formula("shares ~  price + x1 + x2 + x3 + x4 + x5 |
#'                     x1 + x2 + x3 + x4 + x5 |
#'                     0+ x1 + x2 |
#'                     iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )
#'
#' blp_data <- BLP_data(model = model, market_identifier="cdid",
#'                      product_id = "prod_id",
#'                      productData = data,
#'                      integration_method = "MLHS" ,
#'                      integration_accuracy = 40,
#'                      integration_seed = 1)
#'
#' theta2 <- matrix(c(0.5,2), nrow=2)
#' rownames(theta2) <- c("x1","x2")
#' colnames(theta2) <- "unobs_sd"
#'
#' derivatives1 <- dstdtheta_wrap(  blp_data=blp_data,
#'                                   par_theta2 = theta2,
#'                                  market = 2)
#'
#' @export
dstdtheta_wrap <- function(  blp_data, par_theta2, market, printLevel = 1 ){
  nobs <- blp_data$parameters$nobs
  K <- blp_data$parameters$K
  original_market_id <- blp_data$parameters$market_id_char_in
  original_product_id <- blp_data$parameters$product_id
  market <- as.character(market)

  ## check inouts
  if(class(blp_data) != "blp_data") stop("Input has wrong class. Call BLP_data() first.")
  if( printLevel > 0) cat("Mean utility (delta) is used as provided in the BLP_data() function.")
  if( missing(market)) stop("Specify a valid market.")
  if( (length(market) != 1)  ) stop("Only one market can be specified.")
  if( !any(market %in% original_market_id) ) stop("Market is not available in provided dataset.")

  ## calculate share evaluations (all markets)
  current_delta <- blp_data$data$delta
  start_theta2 <- .prepare_theta2(par_theta2,
                                  final_col_names_par = c( "unobs_sd" ,
                                                           blp_data$parameters$demographic_names),
                                  final_row_names_par = colnames(blp_data$data$X_rand),
                                  K = blp_data$parameters$K,
                                  M = blp_data$parameters$total_demogr)
  theta2Mat<- .get.theta2.reshape(theta2.in       = start_theta2$par_theta2,
                                  totalRC         = blp_data$parameters$K,
                                  total.demogr.in = blp_data$parameters$total_demogr,
                                  indices.in      = start_theta2$indices,
                                  fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu
  expmu <-  getExpMu( theta2Matrix = theta2Mat,
                      qv = blp_data$integration$drawsRcMktShape,
                      Xrandom = blp_data$data$X_rand,
                      cdid = blp_data$parameters$market_id,
                      demographics = blp_data$integration$drawsDemMktShape)
  sij <-  getSij(expmu = expmu,
                 expdelta = exp(current_delta),
                 cdindex = blp_data$parameters$cdindex )

  ## extract market information
  indicator_prod <- which( market == original_market_id)
  indicator_mkt <- which( market == unique(original_market_id))

  X_rand_mkt <- blp_data$data$X_rand[indicator_prod,,drop=FALSE]
  drawsRcMktShape_mkt <- blp_data$integration$drawsRcMktShape[indicator_mkt,,drop=FALSE]
  if(blp_data$parameters$total_demogr > 0){
    drawsDemMktShape_mkt <- blp_data$integration$drawsDemMktShape[indicator_mkt,,drop=FALSE]
  }else{
    drawsDemMktShape_mkt<- matrix(NA)
  }
  sij_mkt <- sij[indicator_prod,,drop=FALSE]

  ## calculate dstdtheta_c in the given market
  out <- dstdtheta_c( sijt_arma = sij_mkt,
                      indices = start_theta2$indices,
                      xt_arma = blp_data$data$X_rand[indicator_prod,,drop=FALSE],
                      qvt_arma = drawsRcMktShape_mkt,
                      dt_arma = drawsDemMktShape_mkt,
                      weights_arma = blp_data$integration$weights )

  ## preparing output
  rownames(out) <- paste0("share_" ,original_product_id[indicator_prod])
  names_par <- kronecker( start_theta2$final_col_names_par ,
                          start_theta2$final_row_names_par, paste, sep="*")
  relevantRcDem_index <- start_theta2$indices[,"row"] +
    max( start_theta2$indices[,"row"] ) * ( start_theta2$indices[,"col"] - 1 )
  colnames(out) <- names_par[relevantRcDem_index]

  return(out)

}



#' Calculates derivatives of all shares with respect to all mean utilities in a given market.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param par_theta2 matrix with column and rownames providing a starting value for the optimization routine (see details),
#' @param market character specifying the market in which derivatives are calculated,
#' @param printLevel level of output information (default = 1)
#'
#' @return Returns a numeric matrix with derivatives.
#' Cell in row i and col j is the derivative of share i with respect to mean utility j.
#'
#' @details NA's in \code{par_theta2} entries indicate the exclusion from estimation, i.e. the coefficient is assumed to be zero.
#' If only unobserved heterogeneity is used (no demographics), the column name of \code{par_theta2} must be "unobs_sd".
#' With demographics the colnames must match the names of provided demographics (as in \code{demographic_draws}) and "unobs_sd".
#' Row names of \code{par_theta2} must match random coefficients as specified in \code{model}. Constants must be named "(Intercept)".
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                         Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                         Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                         Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                         true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                Xlin.true.price = -0.2,
#'                                                Xrandom.true = rep(2,K),
#'                                                instrument.effects = rep(2,10),
#'                                                instrument.Xexo.effects = rep(1,5)),
#'                         price.endogeneity = list( mean.xi = -2,
#'                                                   mean.eita = 0,
#'                                                   cov = cbind( c(1,0.7), c(0.7,1))),
#'                         printlevel = 0, seed = 234234 )
#'
#'
#' model <- as.formula("shares ~  price + x1 + x2 + x3 + x4 + x5 |
#'                     x1 + x2 + x3 + x4 + x5 |
#'                     0+ x1 + x2 |
#'                     iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )
#'
#' blp_data <- BLP_data(model = model, market_identifier="cdid",
#'                      product_id = "prod_id",
#'                      productData = data,
#'                      integration_method = "MLHS" ,
#'                      integration_accuracy = 40,
#'                      integration_seed = 1)
#'
#' theta2 <- matrix(c(0.5,2), nrow=2)
#' rownames(theta2) <- c("x1","x2")
#' colnames(theta2) <- "unobs_sd"
#'
#' derivatives2 <- dstddelta_wrap(  blp_data=blp_data,
#'                                   par_theta2 = theta2,
#'                                  market = 2)
#' @export
dstddelta_wrap <- function(   blp_data, par_theta2, market, printLevel = 1 ){
  nobs <- blp_data$parameters$nobs
  K <- blp_data$parameters$K
  original_market_id <- blp_data$parameters$market_id_char_in
  original_product_id <- blp_data$parameters$product_id
  market <- as.character(market)

  ## check inouts
  if(class(blp_data) != "blp_data") stop("Input has wrong class. Call BLP_data() first.")
  if( printLevel > 0) cat("Mean utility (delta) is used as provided in the BLP_data() function.")
  if( missing(market)) stop("Specify a valid market.")
  if( (length(market) != 1)  ) stop("Only one market can be specified.")
  if( !any(market %in% original_market_id) ) stop("Market is not available in provided dataset.")

  ## calculate share evaluations (all markets)
  current_delta <- blp_data$data$delta
  start_theta2 <- .prepare_theta2(par_theta2,
                                  final_col_names_par = c( "unobs_sd" ,
                                                           blp_data$parameters$demographic_names),
                                  final_row_names_par = colnames(blp_data$data$X_rand),
                                  K = blp_data$parameters$K,
                                  M = blp_data$parameters$total_demogr)
  theta2Mat<- .get.theta2.reshape(theta2.in       = start_theta2$par_theta2,
                                  totalRC         = blp_data$parameters$K,
                                  total.demogr.in = blp_data$parameters$total_demogr,
                                  indices.in      = start_theta2$indices,
                                  fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu
  expmu <-  getExpMu( theta2Matrix = theta2Mat,
                      qv = blp_data$integration$drawsRcMktShape,
                      Xrandom = blp_data$data$X_rand,
                      cdid = blp_data$parameters$market_id,
                      demographics = blp_data$integration$drawsDemMktShape)
  sij <-  getSij(expmu = expmu,
                 expdelta = exp(current_delta),
                 cdindex = blp_data$parameters$cdindex )

  ## extract market information
  indicator_prod <- which( market == original_market_id)
  indicator_mkt <- which( market == unique(original_market_id))

  X_rand_mkt <- blp_data$data$X_rand[indicator_prod,,drop=FALSE]
  drawsRcMktShape_mkt <- blp_data$integration$drawsRcMktShape[indicator_mkt,,drop=FALSE]
  if(blp_data$parameters$total_demogr > 0){
    drawsDemMktShape_mkt <- blp_data$integration$drawsDemMktShape[indicator_mkt,,drop=FALSE]
  }else{
    drawsDemMktShape_mkt<- matrix(NA)
  }
  sij_mkt <- sij[indicator_prod,,drop=FALSE]


  ## calculate dstdtheta_c in the given market
  out <- dstddelta_c( sijt = sij_mkt,
                      weights= blp_data$integration$weights )
  colnames(out) <- paste0("meanUtility_" ,original_product_id[indicator_prod])
  rownames(out) <- paste0("share_" ,original_product_id[indicator_prod])
  return(out)

}

#' Calculating the GMM objective for a given set of non-linear parameters.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param par_theta2 matrix with column and rownames providing a starting value for the optimization routine (see details),
#' @param printLevel level of output information ranges from 1 (no GMM results) to 4 (every norm in the contraction mapping)
#'
#' @return Returns a list with results from the GMM evaluation.
#' \describe{
#' \item{\code{local_min}}{GMM point evaluation}
#' \item{\code{gradient}}{GMM derivative with respect to non-linear parameters}
#' \item{\code{delta}}{result of the contraction mapping}
#' \item{\code{xi}}{residuals of GMM evaluation} }
#'
#' @details NA's in \code{par_theta2} entries indicate the exclusion from estimation, i.e. the coefficient is assumed to be zero.
#' If only unobserved heterogeneity is used (no demographics), the column name of \code{par_theta2} must be "unobs_sd".
#' With demographics the colnames must match the names of provided demographics (as in \code{demographic_draws}) and "unobs_sd".
#' Row names of \code{par_theta2} must match random coefficients as specified in \code{model}. Constants must be named "(Intercept)".
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                         Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                         Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                         Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                         true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                Xlin.true.price = -0.2,
#'                                                Xrandom.true = rep(2,K),
#'                                                instrument.effects = rep(2,10),
#'                                                instrument.Xexo.effects = rep(1,5)),
#'                         price.endogeneity = list( mean.xi = -2,
#'                                                   mean.eita = 0,
#'                                                   cov = cbind( c(1,0.7), c(0.7,1))),
#'                         printlevel = 0, seed = 234234 )
#'
#'
#' model <- as.formula("shares ~  price + x1 + x2 + x3 + x4 + x5 |
#'                     x1 + x2 + x3 + x4 + x5 |
#'                     0+ x1 + x2 |
#'                     iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )
#'
#' blp_data <- BLP_data(model = model, market_identifier="cdid",
#'                      product_id = "prod_id",
#'                      productData = data,
#'                      integration_method = "MLHS" ,
#'                      integration_accuracy = 40,
#'                      integration_seed = 1)
#'
#' theta_guesses <- matrix(c(0.5,2), nrow=2)
#' rownames(theta_guesses) <- c("x1","x2")
#' colnames(theta_guesses) <- "unobs_sd"
#'
#' gmm <- gmm_obj_wrap(  blp_data=blp_data,
#'                       par_theta2 = theta_guesses,
#'                       printLevel = 2)
#' gmm$local_min
#'
#'
#' @export
gmm_obj_wrap <- function( blp_data, par_theta2, printLevel = 2){

  nobs <- blp_data$parameters$nobs
  K <- blp_data$parameters$K
  ## BLP_data class
  if(class(blp_data) != "blp_data")
    stop("Input has wrong class. Call BLP_data() first.")

  ## calc matrices
  Z <- blp_data$data$Z
  W <-  try( solve((t(Z) %*% Z)) )
  if (any(class(W) == "try-error"))
    stop("Problems with singular matrizes. This might be caused by (nearly) linear dependent regressors or weak instruments.")
  xzwz <- t(blp_data$data$X_lin) %*% Z %*% W %*% t(Z)
  xzwzx <- xzwz %*% blp_data$data$X_lin
  invxzwzx <- try( solve(xzwzx) )
  if (any(class(invxzwzx) == "try-error"))
    stop("Problems with singular matrices. This might be caused by (nearly) linear dependent regressors or weak instruments.")

  blp_data$data$W <- W
  blp_data$data$xzwz <- xzwz
  blp_data$data$invxzwzx <- invxzwzx

  ## check and prepare par_theta2

  start_theta2 <- .prepare_theta2(par_theta2,
                                  final_col_names_par = c( "unobs_sd" ,
                                                           blp_data$parameters$demographic_names),
                                  final_row_names_par = colnames(blp_data$data$X_rand),
                                  K = blp_data$parameters$K,
                                  M = blp_data$parameters$total_demogr)




  # global variables with blp_results:
  blp_results <- new.env( parent = emptyenv())
  blp_results$deltaOld <- blp_data$data$delta
  blp_results$innerItAll <- c()
  blp_results$negShares<- FALSE
  blp_results$gradient <- rep(NA_real_, start_theta2$total_par )
  innerItAll_out <- blp_results$innerItAll

  start <- Sys.time()
  finalTmp <- gmm_obj(par_theta2 = start_theta2$par_theta2,####
                      indices=start_theta2$indices,
                      blp_results=blp_results,
                      blp_data=blp_data,
                      printLevel=printLevel)
  end <- Sys.time()

  delta_out<- blp_results$deltaOld
  theta_rc_out <- start_theta2$par_theta2
  theta_lin_out <- blp_results$bet
  sij_out <- blp_results$sij
  local_min_out <- finalTmp
  gradient_out <- blp_results$gradient
  jacob_out <- blp_results$jacobian
  xi_out <- blp_results$xi

  if( printLevel > 0) print(local_min_out)
  out <- list("delta" = delta_out,
              "theta_rc" = theta_rc_out ,
              "theta_lin" = theta_lin_out ,
              "sij" = sij_out ,
              "local_min" = local_min_out ,
              "gradient" = gradient_out ,
              "jacob" = jacob_out,
              "xi" = xi_out,
              "time" = end - start)
  names(out$delta) <- paste0( blp_data$parameters$product_id , "_" ,
                              blp_data$parameters$market_id_char_in )

  rownames(out$sij) <- paste0( blp_data$parameters$product_id , "_" ,
                               blp_data$parameters$market_id_char_in )
  colnames(out$sij) <- paste0("individual_", 1: blp_data$integration$amountDraws)

  names_par <- kronecker( start_theta2$final_col_names_par ,
                          start_theta2$final_row_names_par , paste, sep="*")
  relevantRcDem_index <- start_theta2$indices[,"row"] +
    max( start_theta2$indices[,"row"] ) * ( start_theta2$indices[,"col"] - 1 )
  rownames(out$gradient) <- names_par[relevantRcDem_index]

  return( out )


}

#' Performs a contration mapping for a given set of non-linear parameters.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param par_theta2 matrix with column and rownames providing a starting value for the optimization routine (see details),
#' @param printLevel level of output information (default = 1)
#'
#' @return Returns an object of class "blp_cm" with results from the contraction mapping.
#' \describe{
#' \item{\code{delta}}{resulting vector of mean utilities after the contraction mapping}
#' \item{\code{counter}}{inner iterations needed to convergence}
#' \item{\code{sij}}{market share integral evaluations for each product (in rows) for the final mean utility} }
#'
#' @details NA's in \code{par_theta2} entries indicate the exclusion from estimation, i.e. the coefficient is assumed to be zero.
#' If only unobserved heterogeneity is used (no demographics), the column name of \code{par_theta2} must be "unobs_sd".
#' With demographics the colnames must match the names of provided demographics (as in \code{demographic_draws}) and "unobs_sd".
#' Row names of \code{par_theta2} must match random coefficients as specified in \code{model}. Constants must be named "(Intercept)".
#'
#' Starting guesses for the contraction mapping are provided with \code{BLP_data}.
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                         Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                         Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                         Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                         true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                Xlin.true.price = -0.2,
#'                                                Xrandom.true = rep(2,K),
#'                                                instrument.effects = rep(2,10),
#'                                                instrument.Xexo.effects = rep(1,5)),
#'                         price.endogeneity = list( mean.xi = -2,
#'                                                   mean.eita = 0,
#'                                                   cov = cbind( c(1,0.7), c(0.7,1))),
#'                         printlevel = 0, seed = 234234 )
#'
#'
#' model <- as.formula("shares ~  price + x1 + x2 + x3 + x4 + x5 |
#'                     x1 + x2 + x3 + x4 + x5 |
#'                     0+ x1 + x2 |
#'                     iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )
#'
#' blp_data <- BLP_data(model = model, market_identifier="cdid",
#'                      product_id = "prod_id",
#'                      productData = data,
#'                      integration_method = "MLHS" ,
#'                      integration_accuracy = 40,
#'                      integration_seed = 1)
#'
#' theta_guesses <- matrix(c(0.5,2), nrow=2)
#' rownames(theta_guesses) <- c("x1","x2")
#' colnames(theta_guesses) <- "unobs_sd"
#'
#' delta_eval <- getDelta_wrap(  blp_data=blp_data,
#'                               par_theta2 = theta_guesses,
#'                               printLevel = 4)
#'
#' @export
getDelta_wrap <- function(blp_data, par_theta2, printLevel = 1){
  nobs <- blp_data$parameters$nobs
  K <- blp_data$parameters$K
  ## BLP_data class
  if(class(blp_data) != "blp_data")
    stop("Input has wrong class. Call BLP_data() first.")


  ## check and prepare par_theta2

  start_theta2 <- .prepare_theta2(par_theta2,
                                  final_col_names_par = c( "unobs_sd" ,
                                                           blp_data$parameters$demographic_names),
                                  final_row_names_par = colnames(blp_data$data$X_rand),
                                  K = blp_data$parameters$K,
                                  M = blp_data$parameters$total_demogr)

  deltaOld <- blp_data$data$delta

  theta2Mat<- .get.theta2.reshape(theta2.in       = start_theta2$par_theta2,
                                  totalRC         = blp_data$parameters$K,
                                  total.demogr.in = blp_data$parameters$total_demogr,
                                  indices.in      = start_theta2$indices,
                                  fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu

  #Call C++ function:
  tmp <- getDelta(  theta2    = theta2Mat,
                    cdid      = blp_data$parameters$market_id,
                    cdindex   = blp_data$parameters$cdindex,
                    innerCrit = blp_data$parameters$inner_tol,
                    indices   = start_theta2$indices,
                    innerMaxit= blp_data$parameters$inner_maxit,
                    Xrandom          = blp_data$data$X_rand,
                    obsshare         = blp_data$data$shares,
                    deltaOld         = deltaOld,
                    nodesDemMktShape = blp_data$integration$drawsDemMktShape ,
                    nodesRcMktShape  = blp_data$integration$drawsRcMktShape,
                    weights          = blp_data$integration$weights,
                    printLevel = printLevel)


  names(tmp$delta) <-  paste0( blp_data$parameters$product_id , "_" ,
                               blp_data$parameters$market_id_char_in )
  rownames(tmp$sij) <- paste0("share_",
                              blp_data$parameters$product_id , "_" ,
                              blp_data$parameters$market_id_char_in )
  colnames(tmp$sij) <- paste0("individual_", 1: blp_data$integration$amountDraws)

  out <- list( delta = tmp$delta,
               counter = tmp$counter,
               sij = tmp$sij,
               theta2 = start_theta2 )
  class(out) <- "blp_cm"
  return(out)


}


#' Calculates information related to predicted shares for a given set of non-linear parameters and data.
#'
#' @param blp_data data object created by the function \code{BLP_data} (provides, among others, mean utilitys and integration draws),
#' @param par_theta2 matrix with column and rownames providing the evaluation point (see details),
#' @param printLevel level of output information (default = 1)
#'
#' @return Returns a list with information related to predicted shares.
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                         Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                         Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                         Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                         true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                Xlin.true.price = -0.2,
#'                                                Xrandom.true = rep(2,K),
#'                                                instrument.effects = rep(2,10),
#'                                                instrument.Xexo.effects = rep(1,5)),
#'                         price.endogeneity = list( mean.xi = -2,
#'                                                   mean.eita = 0,
#'                                                   cov = cbind( c(1,0.7), c(0.7,1))),
#'                         printlevel = 0, seed = 234234 )
#'
#' model <- as.formula("shares ~  price + x1 + x2 + x3 + x4 + x5 |
#'                     x1 + x2 + x3 + x4 + x5 |
#'                     0+ x1 + x2 |
#'                     iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )
#'
#' blp_data <- BLP_data(model = model, market_identifier="cdid",
#'                      product_id = "prod_id",
#'                      productData = data,
#'                      integration_method = "MLHS" ,
#'                      integration_accuracy = 40,
#'                      integration_seed = 1)
#'
#' theta_guesses <- matrix(c(0.5,2), nrow=2)
#' rownames(theta_guesses) <- c("x1","x2")
#' colnames(theta_guesses) <- "unobs_sd"
#'
#' shares <- getShareInfo(  blp_data=blp_data,
#'                            par_theta2 = theta_guesses,
#'                            printLevel = 4)
#'
#' @export
getShareInfo <- function(blp_data, par_theta2, printLevel = 1){
  nobs <- blp_data$parameters$nobs
  K <- blp_data$parameters$K
  ## BLP_data class
  if(class(blp_data) != "blp_data")
    stop("Input has wrong class. Call BLP_data() first.")

  ## mean utility
  if( printLevel > 0){
    cat("Mean utility (delta) is used as provided in the BLP_data() function.")
  }

  current_delta <- blp_data$data$delta


  ## check and prepare par_theta2

  start_theta2 <- .prepare_theta2(par_theta2,
                                  final_col_names_par = c( "unobs_sd" ,
                                                           blp_data$parameters$demographic_names),
                                  final_row_names_par = colnames(blp_data$data$X_rand),
                                  K = blp_data$parameters$K,
                                  M = blp_data$parameters$total_demogr)

  theta2Mat<- .get.theta2.reshape(theta2.in       = start_theta2$par_theta2,
                                  totalRC         = blp_data$parameters$K,
                                  total.demogr.in = blp_data$parameters$total_demogr,
                                  indices.in      = start_theta2$indices,
                                  fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu

  colnames(theta2Mat) <- c( "unobs_sd" ,
                            blp_data$parameters$demographic_names)
  rownames(theta2Mat) <- colnames(blp_data$data$X_rand)


  # get the exp of individual part of utility:
  expMu <- getExpMu( theta2Matrix = theta2Mat,
                     qv = blp_data$integration$drawsRcMktShape,
                     Xrandom = blp_data$data$X_rand,
                     cdid = blp_data$parameters$market_id,
                     demographics = blp_data$integration$drawsDemMktShape ) ;

  rownames(expMu) <- paste0("expMu_",
                          blp_data$parameters$product_id , "_" ,
                          blp_data$parameters$market_id_char_in )
  colnames(expMu) <- paste0("individual_", 1: blp_data$integration$amountDraws)

  # calculate individual choice probabilities
  Sij <- getSij(expmu = expMu,
                expdelta = exp(current_delta),
                cdindex = blp_data$parameters$cdindex )
  rownames(Sij) <- paste0("share_",
                              blp_data$parameters$product_id , "_" ,
                              blp_data$parameters$market_id_char_in )
  colnames(Sij) <- paste0("individual_", 1: blp_data$integration$amountDraws)

  # calc. aggregated choice probabilities, i.e. shares
  shares <- c( Sij %*% blp_data$integration$weights)


  names(shares) <- paste0( blp_data$parameters$product_id , "_" ,
                           blp_data$parameters$market_id_char_in )

  out <- list( "shares"=shares,
               "sij" = Sij,
               "expMu" =  expMu,
               "theta2" = theta2Mat)

  class(out) <- "shareInfo"

  return(out)



}



#' Calculating the Jacobian for a given set of non-linear parameters and mean utilities.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param par_theta2 matrix with column and rownames providing the evaluation point (see details),
#' @param printLevel level of output information (default = 1)
#'
#' @return Returns a matrix with the jacobian (products in rows, parameters in columns).
#'
#' @details NA's in \code{par_theta2} entries indicate the exclusion from estimation, i.e. the coefficient is assumed to be zero.
#' If only unobserved heterogeneity is used (no demographics), the column name of \code{par_theta2} must be "unobs_sd".
#' With demographics the colnames must match the names of provided demographics (as in \code{demographic_draws}) and "unobs_sd".
#' Row names of \code{par_theta2} must match random coefficients as specified in \code{model}. Constants must be named "(Intercept)".
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                         Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                         Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                         Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                         true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                Xlin.true.price = -0.2,
#'                                                Xrandom.true = rep(2,K),
#'                                                instrument.effects = rep(2,10),
#'                                                instrument.Xexo.effects = rep(1,5)),
#'                         price.endogeneity = list( mean.xi = -2,
#'                                                   mean.eita = 0,
#'                                                   cov = cbind( c(1,0.7), c(0.7,1))),
#'                         printlevel = 0, seed = 234234 )
#'
#'
#' model <- as.formula("shares ~  price + x1 + x2 + x3 + x4 + x5 |
#'                     x1 + x2 + x3 + x4 + x5 |
#'                     0+ x1 + x2 |
#'                     iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )
#'
#' blp_data <- BLP_data(model = model, market_identifier="cdid",
#'                      product_id = "prod_id",
#'                      productData = data,
#'                      integration_method = "MLHS" ,
#'                      integration_accuracy = 40,
#'                      integration_seed = 1)
#'
#' theta_guesses <- matrix(c(0.5,2), nrow=2)
#' rownames(theta_guesses) <- c("x1","x2")
#' colnames(theta_guesses) <- "unobs_sd"
#'
#' jacobian <- getJacobian_wrap(blp_data=blp_data,
#'                              par_theta2 = theta_guesses,
#'                              printLevel = 2)
#' head(jacobian)
#' @export
getJacobian_wrap <- function( blp_data, par_theta2, printLevel = 1){

  nobs <- blp_data$parameters$nobs
  K <- blp_data$parameters$K
  ## BLP_data class
  if(class(blp_data) != "blp_data")
    stop("Input has wrong class. Call BLP_data() first.")

  ## mean utility
  if( printLevel > 0){
    cat("Mean utility (delta) is used as provided in the BLP_data() function.")
  }

  current_delta <- blp_data$data$delta

  ## check and prepare par_theta2

  start_theta2 <- .prepare_theta2(par_theta2,
                                  final_col_names_par = c( "unobs_sd" ,
                                                           blp_data$parameters$demographic_names),
                                  final_row_names_par = colnames(blp_data$data$X_rand),
                                  K = blp_data$parameters$K,
                                  M = blp_data$parameters$total_demogr)

  theta2Mat<- .get.theta2.reshape(theta2.in       = start_theta2$par_theta2,
                                  totalRC         = blp_data$parameters$K,
                                  total.demogr.in = blp_data$parameters$total_demogr,
                                  indices.in      = start_theta2$indices,
                                  fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu


  expmu <-  getExpMu( theta2Matrix = theta2Mat,
                      qv = blp_data$integration$drawsRcMktShape,
                      Xrandom = blp_data$data$X_rand,
                      cdid = blp_data$parameters$market_id,
                      demographics = blp_data$integration$drawsDemMktShape)

  sij <-  getSij(expmu = expmu,
                 expdelta = exp(current_delta),
                 cdindex = blp_data$parameters$cdindex )

  jacobian <- jacob_c(sij = sij,
                      indices   = start_theta2$indices,
                      blp_data = blp_data$data,
                      blp_parameters = blp_data$parameters,
                      blp_integration = blp_data$integration,
                      printLevel = printLevel)

  rownames(jacobian) <-  paste0( blp_data$parameters$product_id , "_" ,
                                 blp_data$parameters$market_id_char_in )

  names_par <- kronecker( start_theta2$final_col_names_par ,
                          start_theta2$final_row_names_par , paste, sep="*")
  relevantRcDem_index <- start_theta2$indices[,"row"] +
    max( start_theta2$indices[,"row"] ) * ( start_theta2$indices[,"col"] - 1 )

  colnames(jacobian) <- names_par[relevantRcDem_index]

  return( jacobian )


}
