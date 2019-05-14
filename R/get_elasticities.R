#' @useDynLib BLPestimatoR
#' @importFrom Rcpp sourceCpp
NULL

#' Calculates elasticities for a given variable and market.
#'
#' @param blp_data data object created by the function \code{BLP_data},
#' @param share_info object with individual and aggregated choice probabilities created by the function \code{getShareInfo},
#' @param theta_lin linear parameter of the variable for which elasticities are calculated for,
#' @param variable character specifying a variable for which elasticities are calculated for,
#' @param products optional: character vector of specific products,
#' @param market character specifying the market in which elasticities are calculated
#' @param printLevel level of output information (default = 1)
#'
#' @return Returns a matrix with elasticities. Value in row j and col i for a variable x,
#' gives the effect of a change in product i's characteristic x on the share of product j.
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
#'
#'
#' @examples
#' K<-2 #number of random coefficients
#' data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
#'                              Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
#'                              Xexo = c("x1", "x2", "x3", "x4", "x5"),
#'                              Xrandom = paste0("x",1:K),instruments = paste0("iv",1:10),
#'                              true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                     Xlin.true.price = -0.2,
#'                                                     Xrandom.true = rep(2,K),
#'                                                     instrument.effects = rep(2,10),
#'                                                     instrument.Xexo.effects = rep(1,5)),
#'                              price.endogeneity = list( mean.xi = -2,
#'                                                        mean.eita = 0,
#'                                                        cov = cbind( c(1,0.7), c(0.7,1))),
#'                              printlevel = 0, seed = 234234 )
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
#' shareObj <- getShareInfo(  blp_data=blp_data,
#'                            par_theta2 = theta_guesses,
#'                            printLevel = 1)
#'
#'
#' get_elasticities(blp_data=blp_data,
#'                  share_info = shareObj ,
#'                  theta_lin = 1,
#'                  variable = "price",
#'                  products = c("4","20"),
#'                  market = 1)
#'
#' @export
get_elasticities <- function( blp_data, share_info, theta_lin, variable, products , market, printLevel = 1){

  if(class(blp_data) != "blp_data")
    stop("blp_data has wrong class. Call BLP_data() first.")

  if(class(share_info) != "shareInfo")
    stop("share_info has wrong class. Call getShareInfo() first.")

  if( missing(variable ))
    stop("Specify variable for which elasticities should be calculated for.")

  if( length(variable) != 1)
    stop("Only one variable can be specified.")

  if( missing(market ))
    stop("Specify mkt in which elasticities should be calculated.")

  if( length(market) != 1)
    stop("Only one market can be specified.")

  X_rand_colNames <- colnames(blp_data$data$X_rand)
  X_lin_colNames <- colnames(blp_data$data$X_lin)

  if ( !( variable %in% c(X_lin_colNames, X_rand_colNames ) ) ){
    stop( paste( "Provided variable", variable , "is not in the set of variables. Constants must be named (Intercept).") )
  }

  ## elasticity calculation----

  original_market_id <- blp_data$parameters$market_id_char_in # market_id_char_in is in order
  market <- as.character(market)

  if( !any(market %in% original_market_id) )
    stop("Market is not available in provided dataset.")

  relevant_Obs_Mkt <- which( market == original_market_id)
  relevant_Mkt <- which( market == unique(original_market_id))

  nobs <- nrow( blp_data$data$X_lin )
  amountDraws <- blp_data$integration$amountDraws
  K <- blp_data$parameters$K
  totalDemographics <- blp_data$parameters$total_demogr
  nprod_Mkt <- length( relevant_Obs_Mkt )
  implShare_Mkt <- share_info$shares[ relevant_Obs_Mkt ]
  product_id_Mkt <- blp_data$parameters$product_id[ relevant_Obs_Mkt ]

  if( missing(products) ){
    product_selector <- 1:nprod_Mkt
  }else{
    product_selector <- products
  }

  if( any(!(product_selector %in% product_id_Mkt ))){
    if(printLevel >0) cat("At least one specified product is not available in the market... selecting all.\n")
    product_selector <- 1:nprod_Mkt
  }

  if( length(theta_lin)!=1 || !is.numeric(theta_lin) || is.na( theta_lin ) )
    stop( "Provide a valid linear parameter theta_lin. If the variable does not enter linerarly, provide zero.")

  # if changingVariable is not modelled as a random coefficient, i.e. use logit:
  if ( !( variable %in% colnames(blp_data$data$X_rand)) ) {
    changingVariable_Mkt <- blp_data$data$X_lin[ relevant_Obs_Mkt , variable, drop = F ]
    # Cols contain the variables that are changed by 1%, and rows contain effects on other products in the choice set.
    EtaMkt <- - matrix( changingVariable_Mkt * theta_lin * implShare_Mkt ,
                        nrow = nprod_Mkt,
                        ncol = nprod_Mkt, byrow = TRUE)
    diag(EtaMkt) <- theta_lin * changingVariable_Mkt * (1 - implShare_Mkt)
  } else {
    # if changingVariable is modelled as a random coefficient:
    changingVariable_Mkt <- blp_data$data$X_rand[ relevant_Obs_Mkt , variable, drop = F ]
    drawsRCMktShape_Mkt <- blp_data$integration$drawsRcMktShape[ relevant_Mkt , , drop = FALSE] # only one line, bec. nodes are used for every product in one market
    drawsRC_tableform_Mkt <- matrix(drawsRCMktShape_Mkt,
                                    nrow = amountDraws,
                                    ncol = K )
    colnames(drawsRC_tableform_Mkt) <- X_rand_colNames #X_lin_colNames, X_rand_colNames

    if ( totalDemographics > 0 ) {
      drawsDemMktShape_Mkt <- blp_data$integration$drawsDemMktShape[ relevant_Mkt, , drop = FALSE]
      demographicReshape_Mkt <- matrix( drawsDemMktShape_Mkt ,
                                        ncol = totalDemographics ,
                                        nrow = amountDraws )
    } else {
      drawsDemMktShape_Mkt <- matrix( NA )
      demographicReshape_Mkt <- matrix( NA )
    }

    sij_Mkt <- share_info$sij[ relevant_Obs_Mkt , ,drop=FALSE]
    weights <- blp_data$integration$weights
    theta2Mat <- share_info$theta2
    sigma <- theta2Mat[ variable , "unobs_sd" ]

    # unobsevered_part contains all individual specific effect parts due to unobs. herterog.
    unobserved_part <- theta_lin +
      sigma * drawsRC_tableform_Mkt[ , variable ]  # get individual price effects

    if ( totalDemographics > 0 ) {
      # extract relevant demographic effects ( ie pick a line of theta2Mat)
      demographic_effects <- matrix( theta2Mat[ variable, -1], # first col is unobs_sd
                                     ncol = totalDemographics,
                                     nrow = amountDraws,
                                     byrow = TRUE )

      # multiply every demographic coefficient with all demogr. draws:
      observed_part <- rowSums( demographic_effects *
                                  demographicReshape_Mkt )

    } else {
      observed_part <- 0
    }

    beta_i <- observed_part + unobserved_part
    Omega <- matrix( NA_real_ ,
                     nrow = nprod_Mkt,
                     ncol = nprod_Mkt)
    scalar <- - matrix( 1/implShare_Mkt ) %*% matrix( changingVariable_Mkt , nrow = 1)
    diag( scalar ) <- - diag( scalar )
    Omega[, ] <- ( t( beta_i )[ rep(1, nprod_Mkt) , ] *
                     t( weights )[ rep( 1 , nprod_Mkt), ] * sij_Mkt ) %*% t( sij_Mkt )
    diag( Omega ) <- c(
      ( t( beta_i )[rep( 1, nprod_Mkt ), ] * sij_Mkt * (1 - sij_Mkt) )
      %*% weights
    )
    EtaMkt <- Omega * scalar

  }

  # after if/esle:
  rownames( EtaMkt ) <- colnames( EtaMkt ) <- product_id_Mkt

  return( EtaMkt[product_selector,product_selector] )

}



# get_elasticities <- function( blp_data, blp_estimation, variable, products , market){
#
#   if(class(blp_data) != "blp_data")
#     stop("blp_data has wrong class. Call BLP_data() first.")
#
#   if(class(blp_estimation) != "blp_est")
#     stop("blp_estimation has wrong class. Call estimate_BLP() first.")
#
#   if( missing(variable ))
#     stop("Specify variable for which elasticities should be calculated for.")
#
#   if( length(variable) != 1)
#     stop("Only one variable can be specified.")
#
#   if( missing(market ))
#     stop("Specify mkt in which elasticities should be calculated.")
#
#   if( length(market) != 1)
#     stop("Only one market can be specified.")
#
#   X_rand_colNames <- colnames(blp_data$data$X_rand)
#   X_lin_colNames <- colnames(blp_data$data$X_lin)
#
#   if ( !( variable %in% c(X_lin_colNames, X_rand_colNames ) ) ){
#     stop( paste( "Provided variable", variable , "is not in the set of variables. Constants must be named (Intercept).") )
#   }
#
#   ## elasticity calculation----
#
#   original_market_id <- blp_data$parameters$market_id_char_in # market_id_char_in is in order
#   market <- as.character(market)
#
#   if( !any(market %in% original_market_id) )
#     stop("Market is not available in provided dataset.")
#
#   relevant_Obs_Mkt <- which( market == original_market_id)
#   relevant_Mkt <- which( market == unique(original_market_id))
#
#   nobs <- nrow( blp_data$data$X_lin )
#   amountDraws <- blp_data$integration$amountDraws
#   K <- blp_data$parameters$K
#   totalDemographics <- blp_data$parameters$total_demogr
#   nprod_Mkt <- length( relevant_Obs_Mkt )
#   obsShare_Mkt <- blp_data$data$shares[ relevant_Obs_Mkt ]
#   product_id_Mkt <- blp_data$parameters$product_id[ relevant_Obs_Mkt ]
#   delta_mkt <- blp_estimation$delta[ relevant_Obs_Mkt ]
#
#   if( missing(products) ){
#     product_selector <- 1:nprod_Mkt
#   }else{
#     product_selector <- products
#   }
#
#   if( any(!(product_selector %in% product_id_Mkt ))){
#     cat("At least one specified product is not available in the market... selecting all.\n")
#     product_selector <- 1:nprod_Mkt
#   }
#
#
#   betabar <- blp_estimation$theta_lin[ variable, ]
#   if( is.na( betabar ) ){
#     # in case that the variable does not enter linerarly
#     betabar <- 0
#   }
#
#   # if changingVariable is not modelled as a random coefficient, ie use logit:
#   if ( !( variable %in% colnames(blp_data$data$X_rand)) ) {
#
#     changingVariable_Mkt <- blp_data$data$X_lin[ relevant_Obs_Mkt , variable, drop = F ]
#
#     # Cols contain the variables that are changed by 1%, and rows contain effects on other products in the choice set.
#     EtaMkt <- - matrix( changingVariable_Mkt * betabar * obsShare_Mkt ,
#                         nrow = nprod_Mkt,
#                         ncol = nprod_Mkt, byrow = TRUE)
#
#     diag(EtaMkt) <- betabar * changingVariable_Mkt * (1 - obsShare_Mkt)
#
#   } else {
#     # if changingVariable is modelled as a random coefficient:
#     changingVariable_Mkt <- blp_data$data$X_rand[ relevant_Obs_Mkt , variable, drop = F ]
#
#     drawsRCMktShape_Mkt <- blp_data$integration$drawsRcMktShape[ relevant_Mkt , , drop = FALSE] # only one line, bec. nodes are used for every product in one market
#     drawsRC_tableform_Mkt <- matrix(drawsRCMktShape_Mkt,
#                                     nrow = amountDraws,
#                                     ncol = K )
#     colnames(drawsRC_tableform_Mkt) <- X_rand_colNames #X_lin_colNames, X_rand_colNames
#
#     if ( totalDemographics > 0 ) {
#       drawsDemMktShape_Mkt <- blp_data$integration$drawsDemMktShape[ relevant_Mkt, , drop = FALSE]
#       demographicReshape_Mkt <- matrix( drawsDemMktShape_Mkt ,
#                                         ncol = totalDemographics ,
#                                         nrow = amountDraws )
#     } else {
#       drawsDemMktShape_Mkt <- matrix( NA )
#       demographicReshape_Mkt <- matrix( NA )
#     }
#
#
#     sij_Mkt <- blp_estimation$sij[ relevant_Obs_Mkt , ]
#     weights <- blp_data$integration$weights
#
#     theta2Mat <- .get.theta2.reshape(theta2.in = blp_estimation$theta_rc,
#                                      totalRC = K,
#                                      total.demogr.in = totalDemographics,
#                                      indices.in = blp_estimation$indices,
#                                      fill = 0)  # NA are replaced by zeros to simplify x * par in getExpMu
#     rownames(theta2Mat) <- blp_estimation$rand_coef_rownames
#     colnames(theta2Mat) <- blp_estimation$rand_coef_colnames
#
#     sigma <- theta2Mat[ variable , "unobs_sd" ]
#
#     # unobsevered_part contains all individual specific effect parts due to unobs. herterog.
#     unobseved_part <- betabar +
#       sigma * drawsRC_tableform_Mkt[ , variable ]  # get individual price effects
#
#     if ( totalDemographics > 0 ) {
#       # extract relevant demographic effects ( ie pick a line of theta2Mat)
#       demographic_effects <- matrix( theta2Mat[ variable, -1], # first col is unobs_sd
#                                      ncol = totalDemographics,
#                                      nrow = amountDraws,
#                                      byrow = TRUE )
#
#       # multiply every demographic coefficient with all demogr. draws:
#       observed_part <- rowSums( demographic_effects *
#                                   demographicReshape_Mkt )
#
#     } else {
#       observed_part <- 0
#     }
#
#     beta_i <- observed_part + unobseved_part
#     Omega <- matrix( NA_real_ ,
#                      nrow = nprod_Mkt,
#                      ncol = nprod_Mkt)
#     scalar <- - matrix( 1/obsShare_Mkt ) %*% matrix( changingVariable_Mkt , nrow = 1)
#     diag( scalar ) <- - diag( scalar )
#     Omega[, ] <- ( t( beta_i )[ rep(1, nprod_Mkt) , ] *
#                      t( weights )[ rep( 1 , nprod_Mkt), ] * sij_Mkt ) %*% t( sij_Mkt )
#     diag( Omega ) <- c(
#       ( t( beta_i )[rep( 1, nprod_Mkt ), ] * sij_Mkt * (1 - sij_Mkt) )
#       %*% weights
#     )
#     EtaMkt <- Omega * scalar
#
#   }
#
#   # after if/esle:
#   rownames( EtaMkt ) <- colnames( EtaMkt ) <- product_id_Mkt
#
#   return( EtaMkt[product_selector,product_selector] )
#
# }
#
#
#
#
