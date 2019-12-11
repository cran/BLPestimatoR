#' @useDynLib BLPestimatoR
#' @importFrom Rcpp sourceCpp
NULL

#' Performs a BLP demand estimation.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param par_theta2 matrix with column and rownames providing a starting value for the optimization routine (see details),
#' @param solver_method character specifying the solver method in \code{optim} (further arguments can be passed to \code{optim} by ...)
#' @param solver_maxit integer specifying maximum iterations for the optimization routine (default=10000),
#' @param solver_reltol integer specifying tolerance for the optimization routine (default= 1e-6),
#' @param standardError character specifying assumptions about the GMM residual (homoskedastic , heteroskedastic (default), or cluster)
#' @param extremumCheck if \code{TRUE}, second derivatives are checked for the existence of minimum at the point estimate (default = FALSE),
#' @param printLevel level of output information ranges from 0 (no GMM results) to 4 (every norm in the contraction mapping)
#' @param ... additional arguments for \code{optim}

#'
#' @return Returns an object of class "blp_est". This object contains, among others, all estimates for preference parameters and standard errors.
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
#' blp_est <- estimateBLP(blp_data =blp_data,
#'                        par_theta2 = theta_guesses,
#'                        extremumCheck = FALSE ,
#'                        printLevel = 1 )
#' summary(blp_est)
#'
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
#' @importFrom Matrix Diagonal
#'
#' @export
estimateBLP <- function( blp_data,
                         par_theta2,
                         solver_method = "BFGS", solver_maxit = 10000, solver_reltol = 1e-6,

                         standardError = "heteroskedastic", extremumCheck = FALSE,
                         printLevel = 2, ... ) {


  call_arguments <- match.call(expand.dots = TRUE) # capture the call used to create the model

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



  cat("blp_data were prepared with the following arguments:\n")

  print(blp_data$call_arguments)



  ## Initialising and optimisation

  if (printLevel > 0) {
    cat("Starting a BLP demand estimation with ", blp_data$parameters$nobs, " observations in ", blp_data$parameters$nmkt, " markets...\n")
    cat("[integration::method", blp_data$integration$integration_method, " integration::amountDraws", blp_data$integration$amountDraws, "]\n")
    cat("[blp::inner_tol", blp_data$parameters$inner_tol, " blp::inner_maxit", blp_data$parameters$inner_maxit, "]\n")
    cat("[solver::method", solver_method, " solver::maxit", solver_maxit," solver::reltol", solver_reltol, "]\n")
  }

  # making use of global variables (environments), because optim allows just a scalar as output of gmm_obj:
  blp_results <- new.env( parent = emptyenv())
  blp_results$deltaOld <- blp_data$data$delta
  blp_results$innerItAll <- c()
  blp_results$negShares<- FALSE
  blp_results$gradient <- rep(NA_real_, start_theta2$total_par )



  ## Estimation ----

  start_time <- Sys.time()

  # optimization

  res <- optim(par = start_theta2$par_theta2,
               fn = gmm_obj, gr = gmm_gr,

               method = solver_method,
               control = list( reltol = solver_reltol,
                               maxit = solver_maxit ),

               indices=start_theta2$indices,
               blp_results=blp_results,
               blp_data=blp_data,
               printLevel=printLevel,
               ... )

  solverMessage<- if( res$convergence==0 ) "Successful convergence" else paste("See error code (optim package)", res$convergence )
  outer_it_out <-  res$counts[1]


  end_time <- Sys.time()
  time <- end_time - start_time

  cat("------------------------------------------ \n")
  cat(paste("Solver message:", solverMessage ,"\n") )
  if( !( solverMessage=="Successful convergence" ) ) stop( "Cannot compute post estimation results due to failed minimization routine." )
  cat("------------------------------------------ \n")
  cat("Final GMM evaluation at optimal parameters: \n")


  # the next call ensures that values that are written to environments
  # and are used for post estimation analysis are based on the optimal
  # set of parameters and not just the last step of the solver:

  innerItAll_out <- blp_results$innerItAll
  blp_results$deltaOld <- rep(0,nobs) #reset environment for final evaluation

  finalTmp <- gmm_obj(par_theta2 = res$par,####
                      indices=start_theta2$indices,
                      blp_results=blp_results,
                      blp_data=blp_data,
                      printLevel=3)

  delta_out<- blp_results$deltaOld
  theta_rc_out <- res$par
  theta_lin_out <- blp_results$bet
  sij_out <- blp_results$sij
  local_min_out <- finalTmp
  gradient_out <- blp_results$gradient
  jacob_out <- blp_results$jacobian
  xi_out <- blp_results$xi

  # naming of rc
  names_rc <- kronecker( start_theta2$final_col_names_par ,
                         start_theta2$final_row_names_par, paste, sep="*")
  relevantRcDem_index <- start_theta2$indices[,"row"] +
    max( start_theta2$indices[,"row"] ) * ( start_theta2$indices[,"col"] - 1 )
  names(theta_rc_out) <- names_rc[relevantRcDem_index]

  ### Post estimation----

  ## standard errors
  X_lin <- blp_data$data$X_lin
  Z <- blp_data$data$Z
  W <- blp_data$data$W

  a <- t(cbind(X_lin, jacob_out)) %*% Z
  tmpSE <- try( solve(a %*% W %*% t(a)) )
  lin_len <- dim(X_lin)[2]

  valid_SE <- (standardError %in% c("heteroskedastic","homoskedastic","cluster")) && (length(standardError)==1)

  if(!valid_SE){
    message("Invalid standard error option is provided. Switching to heteroskedastic standard errors...")
    standardError <- "heteroskedastic"
  }

  if (any(class(tmpSE) == "try-error"))
    stop("Standard errors cannot be computed due to singular matrices.")

  if( standardError == "heteroskedastic") { # default
    cat("Using the heteroskedastic asymptotic variance-covariance matrix... \n")
    #omega <- diag( diag( xi_out %*% t(xi_out) ) )
    #b<- t(Z) %*% omega %*% Z
    omega <- xi_out^2
    b<- as.matrix(t(Z) %*% Diagonal(length(omega),omega) %*% Z)
    COV <- tmpSE %*% a %*% W %*% b %*% W %*% t(a) %*% tmpSE
  }
  if( standardError == "homoskedastic") {
    cat("Using the homoskedastic asymptotic variance-covariance matrix... \n")
    COV <- c( (t(xi_out) %*% xi_out)/nobs ) * tmpSE
  }

  if( standardError == "cluster") {
    group_structure <- blp_data$data$group_structure

    if( any(is.na(group_structure)) || is.null(group_structure) )
      stop("Valid group_structure is not availalbe in blp_data. Clustered standard errors require a variable that describes the cluster affiliation.")

    group_structure <- data.frame(group_structure=group_structure)
    group_matrix <- model.matrix( as.Formula("~0+group_structure"), group_structure )

    tmp <- c(xi_out) * group_matrix

    omega <- tmp %*% t(tmp)
    b<- t(Z) %*% omega %*% Z
    COV <- tmpSE %*% a %*% W %*% b %*% W %*% t(a) %*% tmpSE
  }

  seLinear_out <- sqrt(diag(COV))[1: lin_len ]
  seRc_out <- sqrt(diag(COV))[-(1:lin_len )]

  ## Waldstatistic
  WaldStatistic<- t(matrix( theta_rc_out )) %*%
    solve(COV[-(1:lin_len), -(1:lin_len) ]) %*%
    matrix( theta_rc_out )

  ## extremum Check
  if( extremumCheck ) {
    hessian <- invisible(hessian(func = gmm_obj, x = res$par,
                                 indices=start_theta2$indices,
                                 blp_results=blp_results,
                                 blp_data=blp_data,
                                 printLevel = 0))
    hessianEig <- eigen(hessian)$values
    isMin_out <- sum(hessianEig > 0) == start_theta2$total_par
    isMin_out <- if(isMin_out) 'positive' else 'negative'
    cat( paste( "Extremum Check:" , isMin_out))
  } else {
    isMin_out <- NA }




  output<- list("theta_rc" = theta_rc_out, # solver results...
                "theta_lin" = theta_lin_out,
                "se_rc" = seRc_out,
                "se_linear" = seLinear_out,
                "local_min" = local_min_out,
                "gradient" = gradient_out,
                "time" = time,
                "outer_it" = outer_it_out,
                "inner_it" = innerItAll_out,
                "delta" = delta_out,
                "xi" = xi_out,

                "#demogrCoef"= blp_data$parameters$total_demogr,
                "#nmkt" = blp_data$parameters$nmkt,
                "#nobs" = nobs,
                "#exoCoef" = length(colnames(blp_data$data$X_exg)),
                "indices" = start_theta2$indices,

                "rand_coef_rownames" = start_theta2$final_row_names_par,
                "rand_coef_colnames" = start_theta2$final_col_names_par,
                "drawsRcMktShape" = blp_data$integration$draws_mktShape,
                "drawsDemMktShape" = blp_data$integration$dD,
                "weights" = blp_data$integration$weights,
                "sij" = sij_out,

                "WaldStatistic" =   WaldStatistic, # Postestimation...
                "IslocalMin" = isMin_out,

                "outerCrit" = solver_reltol,
                "innerCrit" = blp_data$parameters$inner_tol,
                "intMethod" = blp_data$integration$method,
                "intdraws" = length(blp_data$integration$weights),
                "standardErrorMethod" = standardError,
                "call" = call_arguments

  )

  class(output) <- 'blp_est'

  return( output )
}
