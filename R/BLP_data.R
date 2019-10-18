#' @useDynLib BLPestimatoR
#' @importFrom Rcpp sourceCpp
NULL

#' Prepares data and parameters related to the BLP algorithm for estimation.
#'
#' @param model the model to be estimated in R's formula syntax,
#' @param market_identifier character specifying the market identifier (variable name must be included in \code{productData}),
#' @param product_identifier character specifying the product identifier (variable name must be included in \code{productData}),
#' @param par_delta optional: numeric vector with values for the mean utility (variable name must be included in \code{productData}),
#' @param group_structure optional: character specifying a group structure for clustered standard erros (variable name must be included in \code{productData}),
#' @param additional_variables optional: character vector specifying variables you want to keep for later analysis (variable names must be included in \code{productData})
#' @param productData data.frame with product characteristics,
#' @param demographic_draws optional: list with demographic draws for each market to consider observed heterogeneity (see details),
#' @param integration_accuracy integer specifying integration accuracy,
#' @param integration_method character specifying integration method,
#' @param integration_draws numeric matrix of manually provided integration draws (see details),
#' @param integration_weights numeric vector of manually provided integration weights,
#' @param integration_seed seed for the draws of Monte Carlo based integration,
#' @param blp_inner_tol tolerance for the contraction mapping (default: 1e-9),
#' @param blp_inner_maxit maximum iterations for the contraction mapping (default: 10000)
#'
#' @details For any form of user provided integration draws, i.e. \code{integration_draws} (unobserved heterogeneity)
#' or \code{demographic_draws} (observed heterogeneity), list entries must be named and contain the variable \code{market_identifier} to allow market matching.
#' Each line in these list entries contains the draws for one market.
#' In case of unobserved heterogeneity, list names must match the random coefficients from the model formula.
#' The \code{par_delta} argument provides the variable name for mean utilitys. For example, in the estimation algorithm these values are used as starting guesses in the contraction mapping.
#' Another example is the evaluation of the GMM, which is also based on the provided mean utilitys.
#' If you need to update \code{par_delta} or any other variable in the data object, use \code{update_BLP_data}.
#'
#' @return Returns an object of class \code{blp_data}.
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
#' @export
BLP_data <- function(model,
                     market_identifier, product_identifier, par_delta, group_structure = NULL, additional_variables=NULL,

                     productData , demographic_draws ,

                     integration_accuracy, integration_method, integration_draws, integration_weights, integration_seed,
                     blp_inner_tol = 1e-9, blp_inner_maxit = 10000){

  #### Formula extraction----
  call_arguments <- match.call(expand.dots = TRUE) # capture the call used to create the model
  formula <- Formula::as.Formula(model)

  # length reports the number of parts on the LHS and RHS:
  model_rhs_length <- length(formula)[2]

  stopifnot(length(formula)[1] == 1L,
            model_rhs_length %in% 1:4)

  # NA Check
  tmp <- model.frame(formula, productData, na.action = na.fail)
  tmp <- NULL

  # shares
  f1 <- formula(formula, lhs = 1, rhs = 0)
  shares <- model.response( model.frame(f1, productData),
                            type = "numeric" )
  if( any( shares < 0 ) || any( shares > 1 ) )  stop( "Shares contain values out of [0,1]." )
  nobs <- length(shares)

  # market and product identifyer
  if(  (!is.character(market_identifier)) || (length(market_identifier)!=1) )
    stop( "market_identifier is not valid." )

  if(  !market_identifier %in% names(productData) )
    stop( "market_identifier is not available in the provided data." )

  market_id_char_in <-  vapply(market_identifier,
                               function(x) as.character( get(x, productData) ) ,
                               character(nobs))
  nmkt <- length( unique(market_id_char_in) )
  market_id_numeric <- .indexing.markets(market_id_char_in) # order of numeric values depends on order of input market identifyer
  market_id_numeric_o <- order(market_id_numeric)
  #numbers in market_id_numeric correspond to order of markets in market_id_char_in: all(unique(market_id_char_in[market_id_numeric_o]) == unique(market_id_char_in))


  if(  !product_identifier %in% names(productData) )
    stop( "market_identifier is not available in the provided data." )
  if(  (!is.character(product_identifier)) || (length(product_identifier)!=1) )
    stop( "product_identifier is not valid." )
  product_id_char_in <-  vapply(product_identifier,
                                function(x) as.character( get(x, productData) ) ,
                                character(nobs))


  # uniqueness check
  tmp <- table(paste0(market_id_char_in,"_",product_id_char_in))
  if( any(tmp>1) )
    warning("Combination of market_identifier and product_identifier is not unique.")

  # BLP parameter check

  if( !missing(blp_inner_tol)){
    if(  (!is.finite(blp_inner_tol)) || (length(blp_inner_tol)!=1) ){
      cat("Invalid blp_inner_tol. Set to default (1e-9).\n")
      blp_inner_tol<- 1e-9
    }
  }

  if( !missing(blp_inner_maxit)){
    if(  (!is.finite(blp_inner_maxit)) || (length(blp_inner_maxit)!=1) ){
      cat("Invalid blp_inner_maxit. Set to default (10000).\n")
      blp_inner_maxit<- 10000
    }
  }


  # mean utility
  missing_delta <- missing(par_delta)

  if( !missing_delta){
    if(  (!is.character(par_delta)) || (length(par_delta)!=1) )
      stop( "par_delta is not valid." )
    par_delta_var_name <- par_delta
    par_delta <- get(par_delta, productData)
    delta_error <- !all(is.finite(exp(par_delta))) # checks NA's and NaN's and infinite values
  } else delta_error <- FALSE

  if( missing_delta || delta_error ){
    par_delta_var_name <- "delta"
    par_delta <- rep(0, nobs)
    cat( "Mean utility (variable name: `delta`) is initialized with 0 because of missing or invalid par_delta argument.\n")
  }


  # linear data
  f2 <- formula(formula, lhs = 0, rhs = 1)
  X_lin <- model.matrix( f2 , productData)
  tmp <- apply(X_lin, 2, function(x) round(sum(abs(diff(x))),
                                           3) == 0)
  if (sum(tmp) > 1)
    stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")

  # exogenous data
  if( model_rhs_length >= 2){
    f3 <- formula(formula, lhs = 0, rhs = 2)
    X_exg <- model.matrix( f3 , productData)
    tmp <- apply(X_exg, 2, function(x) round(sum(abs(diff(x))),
                                             3) == 0)
    if (sum(tmp) > 1)
      stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
  } else X_exg <- NULL

  # random coef. data
  if( model_rhs_length >= 3){
    f4 <- formula(formula, lhs = 0, rhs = 3)
    X_rand <- model.matrix( f4 , productData)
    K <- dim(X_rand)[2]
    tmp <- apply(X_rand, 2, function(x) round(sum(abs(diff(x))),
                                              3) == 0)
    if (sum(tmp) > 1)
      stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
  } else X_rand <- NULL

  # IV's
  if( model_rhs_length >= 4 ){
    f5 <- formula(formula, lhs = 0, rhs = 4)
    IV <- model.matrix( f5, productData )
    tmp <- apply(IV, 2, function(x) round(sum(abs(diff(x))),
                                          3) == 0)
    if (sum(tmp) > 1)
      stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
   } else IV <- NULL

  #### Data preparation----

  ### integration I: normaly distributed RC

  has_own_int <- !missing(integration_weights) && !missing(integration_draws)
  has_int_method <-  !missing(integration_accuracy) && !missing(integration_method)

  if( has_own_int == has_int_method )
    stop("Provide either the name and accuracy of a valid integration method or your own weights and draws.")

  final_order_draws <- colnames(X_rand)

  if(has_own_int){
    integration_method <- "provided_by_user"
    weights <- na.fail( as.matrix(as.numeric( integration_weights ), ncol =1) )
    integration_list_names <- names(integration_draws)

    if(!( length(integration_draws) == K ))
      stop("Provided list of integration draws has not enough entries. Number of random coefficients determines length.")

    if( !setequal( integration_list_names , final_order_draws ))
      stop("Names of list entries for draws (unobs. heterogeneity) do not match with names of random coefficients. Remember to name any constant \"(Intercept)\" .\n")

    final_order <- match( final_order_draws, integration_list_names )
    integration_draws <- integration_draws[ final_order ] # list is now in the order of X_rand

    draws_mktShape <- .draws_listToMatrix( drawList = integration_draws,
                                           amountDraws = length(weights),
                                           market_identifier_pd =get( market_identifier , productData),
                                           market_identifier_list_name = market_identifier,
                                           use = "rc")
    # the order of rows in draws_mktShape in determined by the order of markets in productData$market_id

  }

  if(has_int_method){
    ## c) Generate nodes & weights by a specified accuracy and method:
    tmp<- get_integration_input(dim =  K,
                                method = integration_method,
                                accuracy = integration_accuracy,
                                nmkt = nmkt,
                                seed = integration_seed)
    draws_mktShape <- tmp$nodesMktShape
    weights <- tmp$weights
  }


  stopifnot( all( dim(draws_mktShape) == c(nmkt, length(weights) * K))) # final check

  # integration II: demographic data (optional)
  if( !missing(demographic_draws) ){
    stopifnot(is.list(demographic_draws))
    demographic_names <- names(demographic_draws)
    M <- length( demographic_names )
    dD <- .draws_listToMatrix( drawList = demographic_draws,
                               amountDraws = length(weights),
                               market_identifier_pd =get( market_identifier , productData),
                               market_identifier_list_name = market_identifier,
                               use = "demographics")

    stopifnot( all( dim(dD) == c(nmkt, length(weights) * M))) # final check
  } else {
    demographic_names <- NULL
    dD <- NULL
    M <- 0 }

  # reordering all data according to market identifier

  shares <- shares[ market_id_numeric_o ]
  X_rand <- X_rand[ market_id_numeric_o, ,drop=FALSE]
  X_lin <- X_lin[ market_id_numeric_o, ,drop=FALSE]
  X_exg <- X_exg[ market_id_numeric_o, ,drop=FALSE]
  IV <- IV[ market_id_numeric_o, ,drop=FALSE]
  if( is.null(dD) )
    dD <- matrix(NA) # necessary for Rcpp input type

  if( is.null(group_structure)){
    group_structure <- NULL
  } else{
    if(  (!is.character(group_structure)) || (length(group_structure)!=1) )
      stop( "group_structure is not valid." )
    group_structure <- as.character( get(group_structure, productData) )
    group_structure <- group_structure[market_id_numeric_o]
  }

  market_id_numeric <- market_id_numeric[market_id_numeric_o]
  market_id_char_in <- market_id_char_in[market_id_numeric_o]
  product_id_char_in <- product_id_char_in[market_id_numeric_o]
  par_delta <- par_delta[ market_id_numeric_o ]

  #cdindex
  cdindex <- as.numeric( c(0, cumsum( table(market_id_numeric) )) )

  ## exogenous (included and excluded)
  Z <- cbind(X_exg, IV)
  Z <- Z[, unique(colnames(Z))] #duplicate variables are suppressed

  ## additional_variables
  if( !is.null(additional_variables)){
    additional_data <- data.frame("identifier" = paste0(market_id_char_in,
                                                        product_id_char_in))
    for( i in 1:length(additional_variables)){
      vn_i <- additional_variables[i]
      if(  !vn_i %in% names(productData) )
        stop( paste0(vn_i ," is not available in the provided data." ))
      additional_data[[vn_i]] <- productData[[vn_i]][ market_id_numeric_o ]
    }
  }else{
    additional_data <- NULL
  }

  ## Output

  integration<- list('drawsRcMktShape' = draws_mktShape,
                     'drawsDemMktShape' = dD,
                     'weights' = weights,
                     'method' = integration_method ,
                     'amountDraws' = length(weights) )

  parameters <- list( 'inner_tol' = blp_inner_tol,
                      'inner_maxit'= blp_inner_maxit,
                      'nobs' = nobs,
                      'cdindex' = cdindex,
                      'market_id' = market_id_numeric,
                      'nmkt' = nmkt,
                      'K' = K,
                      'total_demogr'= M,
                      'market_id_numeric_o' = market_id_numeric_o,
                      'demographic_names' = demographic_names,
                      'market_id_char_in' =market_id_char_in,
                      'market_id_varname' = market_identifier,
                      'product_id' = product_id_char_in,
                      'product_id_varname' = product_identifier,
                      'par_delta_varname' = par_delta_var_name,
                      'share_varname' = as.character(f1)[2])

  data <- list('X_lin' = X_lin,
               'X_exg' = X_exg,
               'X_rand' = X_rand,
               'shares' = shares,
               'Z' = Z,
               'group_structure'= group_structure,
               'delta' = par_delta,
               'additional_data' = additional_data)

  out <- list( call_arguments=call_arguments,
               integration = integration,
               parameters =  parameters,
               data = data)

  class(out) <- "blp_data"

  return(out)
}




#' Updates the set of linear, exogenous, random coefficient, share or mean utility variable in the data object.
#'
#' @param data_update data.frame with variables to update (must contain the market_identifier and product_identifier variables as in \code{blp_data}),
#' @param blp_data data object created by the function \code{BLP_data}
#'
#' @return Returns an object of class \code{blp_data}.
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
#' new_data <- data.frame(price = seq(1,10,length.out=500),
#'                        x1 =  seq(2,10,length.out=500),
#'                        cdid = sort(rep(1:25,20)),
#'                        prod_id = rep(1:20,25) )
#' blp_data_example_updated <-update_BLP_data(blp_data = blp_data,
#'                                            data_update = new_data)
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
#' @export
update_BLP_data <- function(data_update,
                            blp_data){

  ## BLP_data class
  if(class(blp_data) != "blp_data")
    stop("Input has wrong class. Call BLP_data() first.")

  if(class(data_update) != "data.frame")
    stop("data_update must be a data.frame.")

  product_id_varname_old <- blp_data$parameters$product_id_varname
  market_id_varname_old <- blp_data$parameters$market_id_varname


  if( is.null(product_id_varname_old) )
    stop("Matching of new data not possible, because product_identifier in blp_data is not available.")

  if( !product_id_varname_old %in% names(data_update) )
    stop(paste0(product_id_varname_old, " is not available in data_update."))
  if( !market_id_varname_old %in% names(data_update) )
    stop(paste0(market_id_varname_old, " is not available in data_update."))


  ## reorder new data according to "old" market_identifier and product_identifier in blp_data
  unique_obs_id_old <- paste0(blp_data$parameters$market_id_char_in,"_",
                              blp_data$parameters$product_id)
  tmp <- table( unique_obs_id_old )
  if( any(tmp>1) )
    stop("Matching not possible. Combination of market_identifier and product_identifier in blp_data is not unique.")

  unique_obs_id_new <- paste0(as.character( data_update[[market_id_varname_old]] ),"_",
                              as.character( data_update[[product_id_varname_old]] ))
  tmp <- table( unique_obs_id_new )
  if( any(tmp>1) )
    stop("Matching not possible. Combination of market_identifier and product_identifier in data_update is not unique.")

  neworder <- match( unique_obs_id_old, unique_obs_id_new  )

  if( any( is.na( neworder ) ) )
    stop("Market/product combinations in new and old data are not matching.")

  data_update <- data_update[neworder, ]
  data_update[product_id_varname_old] <- NULL
  data_update[market_id_varname_old] <- NULL

  ## update all related data objects
  new_variables <- names(data_update)

  for(i in new_variables){
    if(i %in% colnames(blp_data$data$X_lin)){
      blp_data$data$X_lin[,i] <- data_update[[i]]
      cat( paste0("Linear variable ", i ," has been updated.\n"))
    } else if(i %in% colnames(blp_data$data$X_exg)){
      blp_data$data$X_exg[,i] <-  blp_data$data$Z[,i] <- data_update[[i]]
      cat( paste0("Exogenous variable ", i ," has been updated.\n"))
    } else if((i %in% colnames(blp_data$data$Z)) && !(i %in% colnames(blp_data$data$X_exg)) ){
      blp_data$data$Z[,i] <- data_update[[i]]
      cat( paste0("Instrument variable ", i ," has been updated.\n"))
    } else if(i %in% colnames(blp_data$data$X_rand)){
      blp_data$data$X_rand[,i] <- data_update[[i]]
      cat( paste0("Random coefficient variable ", i ," has been updated.\n"))
    }else if(i == blp_data$parameters$share_varname){
      blp_data$data$shares <- data_update[[i]]
      cat( paste0("Share variable ", i ," has been updated.\n"))
    }else if(i == blp_data$parameters$par_delta_varname)  {
      blp_data$data$delta<- data_update[[i]]
      cat( paste0("Mean utility variable ", i ," has been updated.\n"))
    }else{
      cat( paste0("No updates performed!\n"))
    }

  }



  return( blp_data )

}
