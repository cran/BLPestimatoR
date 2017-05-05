#' @useDynLib BLPestimatoR
#' @importFrom Rcpp sourceCpp
NULL

#' Performs a BLP demand estimation.
#'
#' @param Xlin character vector specifying the set of linear variables (variable names
#' must be included in \code{productData})
#'
#' @param Xexo character vector specifying the set of exogenous variables (subset of \code{Xlin})
#'
#' @param Xrandom character vector specifying the set of random coefficients (variable names
#' must be included in \code{productData})
#'
#' @param instruments character vector specifying the set of instrumental variables (variable names
#' must be included in \code{productData})
#'
#' @param shares character vector specifying observed market shares (variable name
#' must be included in \code{productData})
#'
#' @param cdid character vector specifying the market identifier (variable name
#' must be included in \code{productData})
#'
#' @param productData dataframe with product characteristics
#'
#' @param demographics optional: character vector specifying the set of demographic variables (must be included as list entries in \code{demographicData})
#'
#' @param demographicData optional: list with demographic data for each market (see details)
#'
#' @param starting.guesses.theta2 matrix with starting values for the optimization routine
#' (NA entries indicate the exclusion from estimation, i.e. the coefficient is assumed to be zero)
#'
#' @param starting.guesses.delta optional: numeric vector with starting guesses for the mean utility
#'
#' @param solver.method specifies the solver method in \code{optim} or \code{ucminf}
#'
#' @param solver.control list of additional arguments for the optimization routines:
#' \describe{
#' \item{\code{solver.reltol}}{tolerance for the optimization routine }
#' \item{\code{solver.maxit}}{maximum iterations for the optimization routine}
#' \item{\code{...}}{further arguments passed to \code{optim} or \code{ucminf}}}
#'
#' @param blp.control list of additional argruments for the BLP algorithm:
#' \describe{
#' \item{\code{inner.tol}}{tolerance for the contraction mapping}
#' \item{\code{inner.maxit}}{maximum iterations for the contraction mapping}}
#'
#' @param integration.control list of parameters for the BLP integration problem:
#' \describe{
#' \item{\code{method}}{integration method}
#' \item{\code{amountNodes}}{integration accuracy for Monte Carlo based integration}
#' \item{\code{accuracyQuad}}{integration accuracy for integration by quadrature rules}
#' \item{\code{seed}}{seed for the draws of Monte Carlo based integration}
#' \item{\code{nodes}}{set of manually provided integration nodes}
#' \item{\code{weights}}{set of manually provided integration weights}
#' \item{\code{output}}{if \code{TRUE}, integration nodes and individual shares (sij) are available as output} }
#'
#' @param postEstimation.control list of parameters specifying post estimation results:
#' \describe{
#' \item{\code{standardError}}{chose \code{robust} (default) or \code{nonRobust} }
#' \item{\code{extremumCheck}}{if \code{TRUE} (default), second derivatives are checked for the existence of minimum at the point estimate}
#' \item{\code{elasticities}}{character vector specifying the set of variables elasticities are calculated for} }
#'
#' @param printLevel level of output information ranges from 1 (no GMM results) to 4 (every norm in the contraction mapping)
#'
#' @return Returns an object of class 'blp'. This object contains, among others,
#' all estimates for preference parameters and standard errors. Necessary information
#' for further post estimation analysis can be included as well.
#'
#' @details The optimization routines are included in the packages \code{optim} and  \code{ucminf}. Only gradient based methods are supported.
#' The \code{ucminf} clones \code{MATLAB}s' standard trust region optimization routine, which turns out to be effective in
#' avoiding overflow problems in the BLP model. Valid arguments are \code{BFGS} , \code{BFGS_matlab}, \code{L-BFGS-B} or \code{CG}.
#'
#' For solver options, the use of \code{solver.maxit} and \code{solver.rel} is recommended.
#' For conflicts of  \code{solver.maxit} and  \code{solver.reltol} and arguments of the respective solvers,
#' priority is given to the latter. Make sure that additionally provided solver control arguments are valid.
#'
#' For demographics variables, list entries of \code{demographicData} must be named according to \code{demographics}.
#' Each list entry contains a dataframe with the draws for different markets and a variable
#' (according to the \code{cdid} argument) that allows to match the draws with the markets.
#'
#' The logit model is used for elasticity calculation, if the variable of interest is not included in \code{Xrandom}.
#' Columns of the elasticity matrix contain the variables that are changed by 1\%, and rows contain effects on other products in the choice set.
#'
#' @importFrom ucminf ucminf
#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats pchisq
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom mvQuad createNIGrid
#' @importFrom mvQuad rescale
#' @importFrom mvQuad getWeights
#' @importFrom mvQuad getNodes
#' @importFrom numDeriv hessian
#' @importFrom randtoolbox halton
#'
#'
#' @export
#'
#' @examples
#'# Parameters
#' i<-1
#' K<-2
#' Xlin_example <-  c("price", "x1", "x2", "x3", "x4", "x5")
#' Xexo_example <- c("x1", "x2", "x3", "x4", "x5")
#' Xrandom_example <- paste0("x",1:K)
#' instruments_example <- paste0("iv",1:10)
#'
#' # Data generation
#' BLP_data <- get.BLP.dataset(nmkt = 25, nbrn = 20,
#'                             Xlin = Xlin_example,
#'                             Xexo = Xexo_example,
#'                             Xrandom = Xrandom_example,
#'                             instruments = instruments_example,
#'                             true.parameters = list(Xlin.true.except.price = rep(0.2,5),
#'                                                    Xlin.true.price = -0.2, Xrandom.true = rep(2,K),
#'                                                    instrument.effects = rep(2,10),
#'                                                    instrument.Xexo.effects = rep(1,5)),
#'                             price.endogeneity = list( mean.xi = -2,
#'                                                       mean.eita = 0,
#'                                                       cov = cbind( c(1,0.7), c(0.7,1))),
#'                             printlevel = 0, seed = 5326 )
#'
#' # Estimation
#' BLP_est<- estimateBLP(Xlin = Xlin_example,
#'                       Xrandom = Xrandom_example,
#'                       Xexo =  Xexo_example,
#'                       instruments = instruments_example,
#'                       shares = "shares",
#'                       cdid = "cdid",
#'                       productData = BLP_data,
#'                       starting.guesses.theta2 = rep(1,K),
#'                       solver.control = list(maxeval = 5000),
#'                       solver.method = "BFGS_matlab",
#'
#'                       starting.guesses.delta =  rep(1, length(BLP_data$cdid)),
#'                       blp.control = list(inner.tol = 1e-6,
#'                                          inner.maxit = 5000),
#'                       integration.control= list(  method="MLHS",
#'                                                   amountNodes= 100,
#'                                                   seed= 3   ),
#'                       postEstimation.control= list(standardError = "robust",
#'                                                    extremumCheck = TRUE,
#'                                                    elasticities = "price"),
#'                       printLevel = 2)
#'
#' # Show results
#' summary(BLP_est)



estimateBLP <- function(Xlin, Xexo, Xrandom, instruments, demographics,
  						           shares, cdid,

  						           productData, demographicData,

  						           starting.guesses.theta2,
  						           starting.guesses.delta,

                         solver.method="BFGS",
  						           solver.control = list(),
                         blp.control = list(), # inner.tol inner.maxit
  						           integration.control=list(), #amountNodes accuracyQuad seed method output
  						           postEstimation.control= list(), #standardError extremumCheck elasticities
                         printLevel = 0) {

  #### Input Check + First Data Preparations----

          # Existence check of necessary arguments
              # collecting all arguments as a list,
              # implicitly requires the arguments to be non-empty:
              toBeTested<- list("Xlin" = Xlin,
                                "Xexo" = Xexo,
                                "Xrandom"= Xrandom,
                               # "instruments" = instruments, # Update 04/06/17 : iv are not necessary to specify anymore
                                "shares" = shares,
                                "cdid" = cdid,
                                "productData" = productData,
                                "starting.guesses.theta2" = starting.guesses.theta2)

              if( missing(solver.control) ) solver.control<- list()
              if( missing(blp.control) ) blp.control<- list()
              if( missing(integration.control) ) integration.control<- list()
              if( missing(postEstimation.control) ) postEstimation.control<- list()

              if( missing(instruments) ) {
                cat("No IV's are in use...")
                names_productData_before<- names( productData )

                instruments<- "NoIVs"
                productData<- cbind(productData, 10000) # assign some numeric value to pass all the input tests

                colnames(productData)<- c(names_productData_before, instruments)
                }

          # data format productData
              productData <- try( as.data.frame( productData ))
              if( class(productData) ==  "try-error" )
                stop( "Provide data in a valid format (e.g. data.frame or matrix)." )


          # product variable checks


              character_indicator <- c ( is.character( Xlin ) ,
                                         is.character( Xexo ) ,
                                         is.character( Xrandom ) ,
                                         is.character( shares ),
                                         is.character( instruments ) )
              #  cannot be checked in one step, because a vector does coerce the object to a vector of strings (making a subsequent check pointless)

              if( ! all ( character_indicator ) )
                stop( "At least one specified variable is not provided as a string.
                      Data must be provided with the *productData* argument." )

              # all variable strings need to be in productData:
              input_data_strings_withoutCDID <- c(Xlin, Xexo, Xrandom, shares, instruments)  #from now on, it is sure that only characters are used

              if( ! all( c( input_data_strings_withoutCDID , cdid) %in%
                         names(productData) ) )
                stop( "At least one specified variable is not available in the data." )

              # all variables need same length:
              nobsTmp <- unique(c(vapply( c( input_data_strings_withoutCDID,
                                             cdid),
                                         function(x) length(get(x,productData)),
                                         numeric(1)),
                                  dim(instruments)[1]))

              if ( length( nobsTmp ) > 1 )
                stop("All variables need to have the same length.")

              # subset checks:
              if (!all(Xexo %in% Xlin))
                stop("Exogenous variable(s) must be a subset of linear variables.")

              if (!all(Xrandom %in% Xlin))
                message("Random coefficient(s) are not a subset of linear variables.")


          # check starting.guesses.delta
              # check if argument is missing:
              if (missing(starting.guesses.delta)) {
                starting.guesses.delta <- rep(0,nrow(productData))
                message('*starting.guesses.delta* is initialized with a value of 0...')
              }else{

                # check length of vector:
                if( length( starting.guesses.delta ) != nobsTmp )
                  stop("Dimension of *starting.guesses.delta* does not match with number of observations.")

                # check for valid values:
                if( any ( !is.finite(starting.guesses.delta) ) )
                  stop("Starting values for *starting.guesses.delta* must be numeric and finite.")
              }



          # handling of unexpected product data
              # check for shares out of [0,1]:
              min_shares <- min ( get(shares, productData) ,na.rm = TRUE)
              max_shares <- max ( get(shares, productData) ,na.rm = TRUE)

              if( (min_shares < 0) |
                  (max_shares > 1) )
                stop("Shares contain values out of [0,1].")


              # Check *relevant* data for NA's or infinite values and delete them:
              data_withoutCDID<- vapply( input_data_strings_withoutCDID,
                                         function(x) get(x, productData), numeric(nrow(productData) ))
              data_CDID <- get( cdid , productData )

              # except for cdid, everything needs to be finite and numeric:
              if( any ( !is.finite( data_withoutCDID ) &
                        !is.na( data_withoutCDID ) )) # NA's are handled separately
                stop("At least one value in *productData* (except cdid) is not finite.")

              # check data_withoutCDID and data_CDID for NA's
              NAindicator_row <- which( is.na( cbind( data_withoutCDID,
                                                      data_CDID ) ), arr.ind = T )[,1]
              NAcounter_prodData<- length( unique(NAindicator_row) ) # multiple NA'S per obs. are unified

              if(NAcounter_prodData>0) {
                warning(paste(NAcounter_prodData, 'missing value(s) were removed.') )
                productData<- productData[ - NAindicator_row, ]
                starting.guesses.delta <- starting.guesses.delta[ - NAindicator_row ]
                # demographics without a cdid info are reordered later
              }


          # Demographic and starting.guesses.theta2 checks
              # existence check:
              if ( ! missing(demographics) & missing(demographicData))
                stop("Specfied demographic variables without data.")

              # initialize demogr. specific variables:
              if (missing(demographics)) {
                total.demogr <- 0
                relevantDemogr.data <- matrix(NA)
                starting.guesses.theta2 <- matrix(starting.guesses.theta2, ncol=1)

                # dimension check of provided starting guesses:
                if( length(Xrandom) != length(starting.guesses.theta2) )
                  stop("Amount of random coefficients does not match with starting values.")

              } else {

                # check if dem. data is a list:
                if( !is.list( demographicData ) )
                  stop("Demographic data must be provided as a list,
                       where each list entry equals a different demographic variable.")

                # check if all variables are in the dataset:
                if( ! all(demographics %in% names(demographicData)) )
                  stop("Specified demographic variable(s) must be included in *demographicData* with the same name.")

                #length of the character vector:
                total.demogr <- length(demographics)

                # overwrite provided demographicData with *relevant* demographicData according demographics argument:
                demographicData <- demographicData[ demographics ]

                # check if at least one dem. has no cdid argument:
                cdidInDemogr.indicator <- unlist(
                  lapply( demographicData,
                          function(i){
                            indicator <- cdid %in% colnames(i)
                            return(indicator)
                          }) )

                if( !all( cdidInDemogr.indicator )){
                  stop("Cannot match demographic data to markets, because at least one demographic variable
                        has no *cdid* information.")
                }

                # dimension check of provided starting guesses:
                if( !is.matrix(starting.guesses.theta2) |
                    !all(dim(starting.guesses.theta2) == c(length(Xrandom) , 1+length(demographics))) ){
                  stop("Starting guesses must be provided as a matrix with
                       dimensions: [#random coef. X (1+ #demographics)]  ")}

                ### further checks are conducted below...
              }



          # solver
              if( ! ( solver.method %in% c( "BFGS" , "BFGS_matlab", "L-BFGS-B" , "CG"  )) )
                stop("Provide a valid opitmizer (*BFGS* , *BFGS_matlab*, *L-BFGS-B* or *CG* ).")


          # solver.control
              if ("solver.maxit" %in% names(solver.control) ) {
                solver.maxit <- solver.control$solver.maxit
              } else {
                solver.maxit <- 10000 }

              if ("solver.reltol" %in% names(solver.control) ) {
                solver.reltol <- solver.control$solver.reltol
              } else {
                solver.reltol <- 1e-6 }

              # multiple arguments unique to different solvers are not allowed:
              if ( ("maxeval" %in% names(solver.control) &
                   "maxit" %in% names(solver.control))     |
                   ("grtol" %in% names(solver.control) &
                    "reltol" %in% names(solver.control))   )
                stop( "Multiple arguments unique to different solvers are not allowed.")

              # prepare final solver.control list:
              if(solver.method == "BFGS_matlab"){ # ucminf optimizer has other names for control arguments than optim
                  # For conflicts of solver.maxit and arguments of the respective solvers, give priority to arguments of the latter
                  if( "maxeval" %in% names(solver.control) ){
                    solver.maxit <- solver.control$maxeval
                    solver.control$maxeval <- NULL # prevents list entries with same names
                  }
                  if( "grtol" %in% names(solver.control) ){
                    solver.reltol <- solver.control$grtol
                    solver.control$grtol <- NULL # prevents list entries with same names
                  }
                  # prevent name conflicts, because other names than known by the solver are not allowed:
                  solver.control$maxit<- NULL
                  solver.control$reltol<- NULL
                  solver.control$solver.maxit<- NULL
                  solver.control$solver.reltol<- NULL

                  solver.control <- c(solver.control,
                                      "maxeval" = solver.maxit,
                                      "grtol" = solver.reltol)

                  } else {

                    if( "maxit" %in% names(solver.control) ){
                      solver.maxit <- solver.control$maxit
                      solver.control$maxit <- NULL # prevents list entries with same names
                    }
                    if( "reltol" %in% names(solver.control) ){
                      solver.reltol <- solver.control$reltol
                      solver.control$reltol <- NULL # prevents list entries with same names
                    }
                    # prevent name conflicts, because other names than known by the solver are not allowed:
                    solver.control$maxeval<- NULL
                    solver.control$grtol<- NULL
                    solver.control$solver.maxit<- NULL
                    solver.control$solver.reltol<- NULL

                    solver.control <- c(solver.control,
                                        "maxit" = solver.maxit,
                                        "reltol" = solver.reltol)
                    }


              gradientIndicator <- TRUE # gradients are used for optimzation; was specified as an argument of estimateBLP in earlier versions

              if(gradientIndicator){
                gradient <- get.gmm.gr  # function 'gradient' is used in optim
              }else{ gradient <- NULL }


          # blp.control
              # unknown list elements are not allowed:
              if( any( ! (names(blp.control) %in%
                          c("inner.tol","inner.maxit" )) ) )
                stop("Unknown *blp.control* list elements are not allowed.")

              # defaults:
              blp.inner.tol <-   if(is.null( blp.control$inner.tol   ) ) 1e-16 else blp.control$inner.tol
              blp.inner.maxit <- if(is.null( blp.control$inner.maxit ) ) 1000  else blp.control$inner.maxit



          # integration.control

              # unknown list elements are not allowed:
              if( any( ! (names(integration.control) %in%
                          c("method", "amountNodes", "accuracyQuad", "seed", "nodes", "weights" , "output") ) ) )
                stop("Unknown *integration.control* list elements are not allowed.")


              # defaults:
              amountNodes  <-       if(is.null( integration.control$amountNodes  ) ) NA else integration.control$amountNodes
              accuracyQuad <-       if(is.null( integration.control$accuracyQuad ) ) NA else integration.control$accuracyQuad
              seedIntegration <-    if(is.null( integration.control$seed         ) ) NA else integration.control$seed
              methodIntegration <-  if(is.null( integration.control$method       ) ) NA else integration.control$method
              outputIntegration <-  if(is.null( integration.control$output       ) ) FALSE else integration.control$output

              # check if only one optional element of nodes and weights is provided
              if ( ('nodes' %in% names(integration.control) & !('weights' %in% names(integration.control)) ) |
                   ( !('nodes' %in% names(integration.control)) & ('weights' %in% names(integration.control) ) ) ){
                stop("If nodes and weights are provided manually, both arguments must be available.")
                ### Further Checks for manually provided nodes and weights are conducted below
              }


              # checks for NOT manually provided nodes and weights:
              if ( ! all( c( 'nodes' , 'weights' ) %in%
                          names(integration.control)    )) {

                if( is.na( methodIntegration ) )
                  stop("Specify the method you want to use for numerical integration.")

                # Check for valid numerical integration method in *get.integration.input*

                if( is.na( amountNodes ) &
                    is.na( accuracyQuad ) ) stop("Specify integration accuracy.")

                if( ! ( is.finite( amountNodes ) |
                        is.finite( accuracyQuad ) )) stop("Provide numerical value for integration accuracy.")

                if( is.finite( amountNodes ) &
                    is.finite( accuracyQuad ) ) stop("Provide only one integration accuracy argument.")

                if( ( !is.na( seedIntegration) & !is.finite(seedIntegration) ) |
                    length(seedIntegration) == 0  ) stop("Provide valid seed.")
              }

              # outputIntegration must be true or false:
                  if( !(outputIntegration==TRUE | outputIntegration==FALSE) )
                    stop('Provide valid parameter for *integration.control$output* ( TRUE or FALSE )' )




          # postEstimation.control
              # unknown list elements are not allowed:
              if( any( ! (names( postEstimation.control ) %in%
                          c("standardError", "extremumCheck" , "elasticities") ) ) )
                stop("Unknown *postEstimation.control* list elements are not allowed.")


              if( is.null( postEstimation.control$standardError )){
                standardError <-  "robust"
                message('Standard errors are calculated as *robust*...')
                } else {
                  standardError <-  postEstimation.control$standardError

                  # Check for correct specified argument:
                  if(  standardError != "robust" &
                       standardError != "nonRobust"  )
                    stop('Provide valid parameter for *standardError* ( robust or nonRobust )' )
                }


              if( is.null( postEstimation.control$extremumCheck )){
                extremumCheck <- TRUE
                } else {
                  extremumCheck <-  postEstimation.control$extremumCheck

                  # extremumCheck must be true or false:
                  if( !(extremumCheck==TRUE | extremumCheck==FALSE) )
                    stop('Provide valid parameter for *extremumCheck* ( TRUE or FALSE )' )
                }


              if( is.null( postEstimation.control$elasticities )){
                elasticities <- NULL
              } else {
                elasticities <-  postEstimation.control$elasticities

                # elasticities must be a valid variable:
                if( ! all(elasticities %in% names(productData) ) )
                  stop('Provide valid parameter for *elasticities* .' )
              }



          # print level
              if( !is.finite(printLevel) )
                stop("Provide a valid printLevel.")

              if( printLevel < 0 )
                printLevel<- 0
              if( printLevel > 4 )
                printLevel<- 4



  ### Input Storage----

        #Sort data by cdid and renew index:
        cdidOld <- get(cdid, productData) # cdidOld is the "input" cdid
        productData <- productData[ order(cdidOld) , ] # sorts data by market
        cdidOrdered <- get(cdid, productData)
        index.objects <- .indexing.markets( cdidOrdered )

        # Prepare data:
        nobs <- length( cdidOrdered )
        cdindex <- index.objects$cdindex
        cdidOrderedNames <- cdidOrdered
        cdidOrdered <- index.objects$cdid # generates a continously increasing numeric id
        nmkt <-   index.objects$nmkt
        K <- length(Xrandom)  # number of random coefficients; length of the string vector

        Xlin.data <- vapply(Xlin, function(x) get(x, productData), numeric(nobs))
        Xexo.data <- vapply(Xexo, function(x) get(x, productData), numeric(nobs))
        Xrandom.data <- vapply(Xrandom, function(x) get(x, productData), numeric(nobs))
        share.data <- get(shares, productData)
        iv.data <- vapply(instruments, function(x) get(x, productData), numeric(nobs))

        if( colnames(iv.data)[1]=="NoIVs" ){
          iv.data<- c()
        }

        Z <- cbind(Xexo.data, iv.data)  #all Instruments (labeled IV in Nevos' Code)

        W <-  try( solve((t(Z) %*% Z)) )
        if (class(W) == "try-error")
          stop("Problems with singular matrizes. This might be caused by (nearly) linear dependent regressors or weak instruments.")
        xzwz <- t(Xlin.data) %*% Z %*% W %*% t(Z)
        xzwzx <- xzwz %*% Xlin.data
        invxzwzx <- try( solve(xzwzx) )
        if (class(invxzwzx) == "try-error")
          stop("Problems with singular matrizes. This might be caused by (nearly) linear dependent regressors or weak instruments.")

        # for data preparation + storage of demographics see the integration part

  ### Integration----
        # provide nodes & weights manually:
        if ( all( c( 'nodes' , 'weights' ) %in% names(integration.control) ) ){

          if( !( ncol(integration.control$nodes) == K |
                 ncol(integration.control$nodes) == ( K*length(integration.control$weights) ))
              ) {

            stop("Provided set of nodes needs dimensions [AmountNodes X Dim]
                 or [nmkt X (Dim*AmountNodes)].  ")

          }else if( ncol(integration.control$nodes) == K  ){

            if( nrow( integration.control$nodes )!=
                length(integration.control$weights) )
              stop("Length of *weights* does not match amount of nodes implied by *nodes*.")

            message("Provided set of nodes is of dimension [AmountNodes X Dim] and is used for all markets.")

            nodesRcMktShape <- matrix( c( integration.control$nodes ), nrow=1)[rep(1,nmkt),]


          } else if( ncol(integration.control$nodes) ==
                     (K*length(integration.control$weights)  ) ){

            if( nmkt != nrow(integration.control$nodes)  )
              stop("Number of markets does not match with dimensions of provided nodes.
                   Provided set of nodes needs dimensions [AmountNodes X Dim]
                   or [nmkt X (Dim*AmountNodes)]")
            message("Provided set of nodes differs between markets.
                Each row of nodes matrix is assigned to markets ordered by an increasing *cdid* index. ")

            nodesRcMktShape <-  integration.control$nodes
          }

          weights <- integration.control$weights
          methodIntegration <- "provided_by_user"

        }else{
        # or generate nodes & weights by a specified accuracy and method:
        tmp<- get.integration.input(dim =  K,
                                    method = methodIntegration,
                                    amount_nodes = amountNodes,
                                    accuracy = accuracyQuad,
                                    nmkt = nmkt,
                                    seed = seedIntegration)
        nodesRcMktShape <- tmp$nodesMktShape
        weights <- tmp$weights
        }

        amountNodes<- length(weights)


    ### Demographics----

      # Data preparation + storage for demographics :
      if (total.demogr > 0) {


        drawsInDemogr<- unlist(
                        lapply( demographicData,
                               function(i){
                                 draws <- ncol(i) -1 # minus 1 because of cdid col.
                                 return(draws)
                                }))


        if( any( drawsInDemogr < amountNodes ) )
          stop("Number of draws for at least one demographic is
               smaller as the provided integration accuracy.
               Include more draws for demographics.")
          #future version might include automatic resampling for this case



        tmp<- lapply( demographics,
                      function(i){
                        tmpDemogr.data <- as.data.frame( demographicData[[ i ]] )
                        tmpDemogr_cdid <- get(cdid,tmpDemogr.data)

                        if( any( table(tmpDemogr_cdid) >1) )
                          stop("Demographic data is not unique for at least one market.")

                        # becomes relevant, if data were reordered by cdid:
                        relevantObs <- match(unique(cdidOrderedNames), # unique keeps the order and uses
                                                                       # names of original CDID argument (and not the renamed one)
                                             tmpDemogr_cdid)
                        if( any( is.na( relevantObs ) ) )
                          stop("No demographic data available for at least one market.")
                        tmpDemogr.data[[cdid]] <- NULL
                        tmpDemogr.data <- as.matrix( tmpDemogr.data[ relevantObs, 1:amountNodes ]) # now, tmpDemogr.data has order as product.data

                        if( ! all( is.finite( tmpDemogr.data ) ) )
                          stop("Demographic data contains non-numeric values.")

                        return(tmpDemogr.data)
                      } )


        # Demographic data is arranged in a matrix, with the draws of
        # different variables next to each other :
        relevantDemogr.data <- do.call(cbind,tmp)

      }else{
        relevantDemogr.data <- matrix(NA)
      }



  ### Initialising and optimisation----


      indices <- which( !is.na( starting.guesses.theta2 ),
                        arr.ind = TRUE ) # NA's are excluded from the estimation algorithm
      starting.guesses.theta2.optim <- na.omit( c( starting.guesses.theta2) )

      if( any (!is.finite( starting.guesses.theta2.optim ) ))
        stop( "Provide finite starting guesses for *starting.guesses.theta2* .")

      total.par <- length(starting.guesses.theta2.optim)


      # Check for overflow problems at start solution in expmu and expdelta
      # (these parts are calculated by getDelta)
          if( any ( is.infinite (exp ( starting.guesses.delta ) ))){
            warning("Too large starting values for *starting.guesses.delta* ... values are resetted to 0.")
            starting.guesses.delta <- rep(0,nrow(productData))
          }

          theta2Mat_precheck<- .get.theta2.reshape(theta2.in       = starting.guesses.theta2.optim,
                                                   totalRC         = K,
                                                   total.demogr.in = total.demogr,
                                                   indices.in      = indices,
                                                   fill            = 0 )

          InfCheckMatrix <- getExpMu( theta2Matrix = theta2Mat_precheck,
                                      qv = nodesRcMktShape,
                                      Xrandom = Xrandom.data,
                                      cdid = cdidOrdered,
                                      demographics = relevantDemogr.data )
          if( any( is.infinite( InfCheckMatrix )))
            stop("Overflow for starting values of *starting.guesses.theta2* . Provide other starting values or rescale your data.")





      if (printLevel > 0) {
        cat("Starting a BLP demand estimation with ", nobs, " observations in ", nmkt, " markets...\n")
        cat("[integration::method", methodIntegration, " integration::nodes", amountNodes, "]\n")
        cat("[blp::inner.tol", blp.inner.tol, " blp::inner.maxit", blp.inner.maxit, "]\n")
        cat("[solver::method", solver.method, " solver::maxit", solver.maxit, "]\n")
      }


      # making use of global variables (environments),
      # because optim allows just a scalar as output of get.gmm.obj:
      blp.results <- new.env( parent = emptyenv())
        blp.results$deltaOld <- starting.guesses.delta
        blp.results$innerItAll <- c()
        blp.results$negShares<- FALSE
        blp.results$gradient <- rep(NA_real_,
                                    length( starting.guesses.theta2.optim ) )

      blp.integration<- list('nodesRcMktShape' = nodesRcMktShape,
                             'weights' = weights,
                             'method' = methodIntegration ,
                             'amountNodes' = amountNodes)


      blp.parameters <- list( 'inner.tol' = blp.inner.tol,
                              'inner.maxit'= blp.inner.maxit,
                              'gradientIndicator' = gradientIndicator,
                              'nobs' = nobs,
                              'cdindex' = cdindex,
                              'cdid' = cdidOrdered,
                              'nmkt' = nmkt,
                              'K' = K,
                              'total.demogr'= total.demogr,
                              'indices' = indices,
                              'total.par' = total.par
                              )


      blp.data <- list('Xlin' = Xlin.data,
                       'Xexo' = Xexo.data,
                       'Xrandom' = Xrandom.data,
                       'obsshare' = share.data,
                       'Z' = Z,
                       'W' = W,
                       'xzwz' = xzwz,
                       'invxzwzx' = invxzwzx,
                       'demographics' = relevantDemogr.data)


      start.time <- Sys.time()

      if( solver.method== "BFGS_matlab"){ # In case the Matlab Solver should be used
        res <- ucminf(par = starting.guesses.theta2.optim,
                      fn = get.gmm.obj, gr = gradient,

                      control = solver.control,

                      blp.integration = blp.integration,
                      blp.parameters =  blp.parameters,
                      blp.data = blp.data,
                      blp.results = blp.results,

                      printLevel = printLevel)
        solverMessage<- if( res$convergence>=0 ) "Successful convergence" else paste("See error code (ucminf package)", res$convergence )
        outer.it.out <-  res$info["neval"]

         }else{ # Standard R optimizers

      res <- optim(par = starting.guesses.theta2.optim,
                   fn = get.gmm.obj, gr = gradient,

                   method = solver.method,
                   control = solver.control,

                   blp.integration = blp.integration,
                   blp.parameters =  blp.parameters,
                   blp.data = blp.data,
                   blp.results = blp.results,

                   printLevel = printLevel)
      solverMessage<- if( res$convergence==0 ) "Successful convergence" else paste("See error code (optim package)", res$convergence )
      outer.it.out <-  res$counts[1]

      }
      end.time <- Sys.time()
      time <- end.time - start.time


  ### Post estimation----

      innerItAll.out <- blp.results$innerItAll

      cat("------------------------------------------ \n")
      cat(paste("Solver message:", solverMessage ,"\n") )
          if( !( solverMessage=="Successful convergence" ) ) stop( "Cannot compute post estimation results due to failed minimization routine." )
      cat("------------------------------------------ \n")
      cat("Final GMM evaluation at optimal parameters: \n")

      # the next call ensures that values that are written to environments
      # and are used for post estimation analysis are based on the optimal
      # set of parameters and not just the last step of the solver:
      finalTmp <- get.gmm.obj(theta2 = res$par,
                              blp.integration,
                              blp.parameters,
                              blp.data,
                              blp.results,
                              printLevel=3)


      delta.out<- blp.results$deltaOld
      sij.out <- blp.results$sij

      theta.lin.out <- blp.results$bet
      names(theta.lin.out) <- Xlin

      theta.rc.out <- res$par
      names(theta.rc.out)[1:K] <- Xrandom # full naming is done in summary

      gradient.out <- blp.results$gradient
      jacob.out <- blp.results$jacobian
      xi.out <- blp.results$xi


      ### se (robust and non robust)

      a <- t(cbind(Xlin.data, jacob.out)) %*% Z
      IVres <- Z * xi.out[ , rep(1, dim(Z)[2])]
      b <- t(IVres) %*% IVres
      tmpSE <- try( solve(a %*% W %*% t(a)) )
      if (class(tmpSE) == "try-error"){
        stop("Standard errors cannot be computed due to singular matrizes.")
      } else{
        if( standardError == "robust") { # default
          cat("Using robust standard errors... \n")
          COV <- tmpSE %*% a %*% W %*% b %*% W %*% t(a) %*% tmpSE
        } else {
          cat("Using non robust standard errors... \n")
          COV <- c( (t(xi.out) %*% xi.out)/nobs ) * tmpSE
          }
        seLinear.out <- sqrt(diag(COV))[1:length(Xlin)]
        seRc.out <- sqrt(diag(COV))[-(1:length(Xlin))]

      }


      ### Waldstatistic
      WaldStatistic<- t(matrix( theta.rc.out )) %*%
        solve(COV[-(1:length(Xlin)), -(1:length(Xlin)) ]) %*%
        matrix( theta.rc.out )

      ### extremum Check
      if( extremumCheck ) {
        #from now on, dont use blp.results environment anymore:
        blp.parameters$gradientIndicator <- FALSE # Hessian does not need gradient
        hessian <- invisible(hessian(func = get.gmm.obj, x = res$par,
                           blp.integration = blp.integration,
                           blp.parameters =  blp.parameters,
                           blp.data = blp.data,
                           blp.results = blp.results,

                           printLevel = 0))
        hessianEig <- eigen(hessian)$values
        isMin.out <- sum(hessianEig > 0) == (K)
        isMin.out <- if(isMin.out){'positive'}else{'negative'}
        cat( paste( "Extremum Check:" , isMin.out))
      } else {
        isMin.out <- NA }


      ### priceElasticities

      if( ! is.null(elasticities ) ) {
        elasticitiy_results <- lapply( elasticities,
                                       function(i){
                                         res_tmp <- list()

                                         # use orginal names for the outputlist
                                         for( j in unique( cdidOrderedNames )){
                                           res_tmp[[j]] <- get.Elasticities( changingVariable_name = i, #parameters...
                                                                             market = j,

                                                                             productData = productData ,           # data...
                                                                             Xrandom_name = Xrandom,
                                                                             delta = delta.out,

                                                                             cdid_name = cdid,
                                                                             share_name = shares ,

                                                                             all_randomCoef = theta.rc.out,
                                                                             all_linearCoef = theta.lin.out,
                                                                             nodesRCMktShape = nodesRcMktShape,
                                                                             nodesDemMktShape = relevantDemogr.data ,
                                                                             Integration_Weights = weights,
                                                                             totalDemographics = total.demogr,
                                                                             indices_InParameters = indices,
                                                                             sij = sij.out )
                                         }
                                         names(res_tmp) <- unique( cdidOrderedNames )
                                         return(res_tmp)
                                       }        )

        names( elasticitiy_results ) <- elasticities

      } else {
        elasticitiy_results <- NA
      }

      # Delete integration information (can be very large objects), if specified
      if( !outputIntegration ) {
        nodesRcMktShape <- "suppressed_by_default"
        relevantDemogr.data <- "suppressed_by_default"
        weights <- "suppressed_by_default"
        sij.out <- "suppressed_by_default"
      }

### Output list ----
      output<- list("theta.rc" = theta.rc.out, # solver results...
                    "theta.lin" = theta.lin.out,
                    "se.rc" = seRc.out,
                    "se.linear" = seLinear.out,
                    "local.minimum" = res$value,
                    "gradient" = gradient.out,
                    "time" = time,
                    "outer.it" = outer.it.out,
                    "inner.it" = innerItAll.out,
                    "delta" = delta.out,
                    "xi" = xi.out,

                    "demographicNames" = if( total.demogr > 0 ) demographics , # Names and amounts...
                    "RcNames" = Xrandom,
                    "#demogrCoef"= total.demogr ,
                    "#nmkt" = nmkt,
                    "#nobs" = nobs,
                    "#exoCoef" = length(Xexo),
                    "indices" = indices,

                    "nodesRcMktShape" = nodesRcMktShape, # optional integration information (can be large!)...
                    "nodesDemMktShape" = relevantDemogr.data,
                    "weights" = weights,
                    "sij" = sij.out,

                    "WaldStatistic" =   WaldStatistic, # Postestimation...
                    "IslocalMin" = isMin.out,
                    "elasticities" = elasticitiy_results,

                    "solver" = solver.method, # Parameters important for other functions
                    "outerCrit" = solver.reltol,
                    "innerCrit" =   blp.inner.tol,
                    "intMethod" = methodIntegration ,
                    "intNodes" = amountNodes,
                    "startingGuesses" =  starting.guesses.theta2,
                    "standardErrorMethod" = standardError)




  class(output) <- 'blp'

  return( output )

}








get.gmm.obj <- function(theta2 ,
                        blp.integration,
                        blp.parameters,
                        blp.data,
                        blp.results,
                        printLevel){

      theta2Mat<- .get.theta2.reshape(theta2.in       = theta2,
                                      totalRC         = blp.parameters$K,
                                      total.demogr.in = blp.parameters$total.demogr,
                                      indices.in      = blp.parameters$indices,
                                      fill            = 0 ) # NA are replaced by zeros to simplify x * par in getExpMu

      deltaOld <- blp.results$deltaOld # delta vector from previous run (without the blp.results environment)

      #Call C++ function:
      tmp <- getDelta(  theta2    = theta2Mat,
                        cdid      = blp.parameters$cdid,
                        cdindex   = blp.parameters$cdindex,
                        innerCrit = blp.parameters$inner.tol,
                        indices   = blp.parameters$indices,
                        innerMaxit= blp.parameters$inner.maxit,
                        Xrandom          = blp.data$Xrandom,
                        obsshare         = blp.data$obsshare,
                        deltaOld         = deltaOld,
                        nodesDemMktShape = blp.data$demographics ,
                        nodesRcMktShape  = blp.integration$nodesRcMktShape,
                        weights          = blp.integration$weights,
                        printLevel = printLevel)


      delta <- tmp$delta
      counter <- tmp$counter

      bet <- NA
      xi <- NA
      sij<- NA
      jacobian <- NA
      gradient <- NA

      if (any(is.nan(delta))) {
        gradient <- rep(Inf, blp.parameters$total.par)
        f <- Inf
      } else {
        #save delta as start solution for the next run,
        #only if contraction mapping converged in the previous step:
        blp.results$deltaOld <- delta
        ## Objective
        bet <- blp.data$invxzwzx %*% (blp.data$xzwz %*% delta)
        xi <- delta - blp.data$Xlin %*% bet
        tmp2 <- t(xi) %*% blp.data$Z
        f <- c(tmp2 %*% blp.data$W %*% t(tmp2))

        ## Gradient
        if ( blp.parameters$gradientIndicator ) {
          sij <- matrix(NA_real_, nrow = blp.parameters$nobs , ncol = blp.integration$amountNodes )
          # init tmp3 ?
          sij[ , ] <- tmp$sijMod * exp(delta)
          jacobian <- get.jacob(sij = sij,
                            blp.data = blp.data,
                            blp.parameters = blp.parameters,
                            blp.integration = blp.integration,
                            printLevel)

          if (any(is.na(jacobian))) {
          if (printLevel > 0) { cat("\t gradient contains Na's --> objective value and gradient replaced by Inf \n")}
            gradient <- rep(Inf, blp.parameters$total.par)
            f <- Inf
            }else{

           gradient <- 2 * t(jacobian) %*% blp.data$Z %*% blp.data$W %*%
            t(blp.data$Z) %*% xi

          } #is.na(jacobian)
          } #end gradient

      } #end is.nan(delta)



      if (printLevel >= 1) {
        cat("gmm objective:", round(f, 4))
        if ( any(is.nan(delta)))  cat(" [delta contains NaN's] ")
        cat("\n")
      }
      if (printLevel >= 2) {
        cat("\t theta (RC): ")
        cat( round(theta2Mat[ ,1] , 2), "\n")
        if( blp.parameters$total.par>0 ){
          cat("\t theta (demogr.): ")
          cat( round(c(theta2Mat[ ,-1]) , 2), "\n")
        }
        cat("\t inner iterations: ")
        cat( counter, "\n")
        cat("\t gradient: " )
        cat( round(gradient,4), "\n")
      }

      # save to environment (not provided as return object, because optim only accepts a single number as output)
      blp.results$bet <- bet
      blp.results$gradient <- gradient
      blp.results$jacobian <- jacobian
      blp.results$xi <- xi
      blp.results$innerItAll <- c(blp.results$innerItAll, tmp$counter)
      blp.results$sij <- sij
      if(tmp$negShares==TRUE) blp.results$negShares <- TRUE

      return(f)
  }




get.integration.input <- function(dim, method, amount_nodes,
                                  accuracy, nmkt, seed) {

  if (!is.na(seed)) set.seed(seed)

  if (!(method %in% c("MLHS","Halton","randHalton","MC","sgGH","sgNH"))){
    stop("Integration method not available! Choose *MLHS*, *Halton*, *randHalton*, *MC*, *sgGH*, *sgNH*.")
  }

  ## Calc. Nodes and Weights ----
  output <- list()
  attributes(output) <- list(method = method, amountNodes = NA)


  # MLHS
  if (method == "MLHS") {
    if (is.na(amount_nodes)){
      stop("For method 'MLHS' the argument 'amountNodes' in the 'integration.control' needs to be specified")
    }

    nodes <- replicate(nmkt, .MLHS(dim, amount_nodes))
    nodes_allMkt <- qnorm(t(sapply(1:nmkt,function(x){ c(nodes[,,x]) } )))
    weights <- rep(1/amount_nodes, amount_nodes)
  }


  # Halton
  if (method == "Halton") {
    if (is.na(amount_nodes)){
      stop("For method 'Halton' the argument 'amountNodes' in the 'integration.control' needs to be specified")
    }

    nodes <- replicate(nmkt, .Halton(dim, amount_nodes, randomized = FALSE))
    nodes_allMkt <- qnorm(t(sapply(1:nmkt,function(x){ c(nodes[,,x]) } )))
    weights <- rep(1/amount_nodes, amount_nodes)
  }

  # randomized Halton
  if (method == "randHalton") {
    if (is.na(amount_nodes)){
      stop("For method 'Halton' the argument 'amountNodes' in the 'integration.control' needs to be specified")
    }

    nodes <- replicate(nmkt, .Halton(dim, amount_nodes, randomized = TRUE))
    nodes_allMkt <- qnorm(t(sapply(1:nmkt,function(x){ c(nodes[,,x]) } )))
    weights <- rep(1/amount_nodes, amount_nodes)
  }

  # MC
  if (method == "MC") {
    if (is.na(amount_nodes)){
      stop("For method 'MC' the argument 'amountNodes' in the 'integration.control' needs to be specified")
    }

    nodes_allMkt <- matrix(rnorm(nmkt * amount_nodes * dim), nrow=nmkt)
    weights <- rep(1/amount_nodes, amount_nodes)
  }

  # Gauss Hermite
  if (method == "sgGH") {
    if (is.na(accuracy)){
      stop("For method 'sgGH' the argument 'accuracyQuad' in the 'integration.control' needs to be specified")
    }

    grid <- createNIGrid(dim, type = "GHN", level = accuracy, ndConstruction = "sparse")
    nodes <- c(getNodes(grid))
    nodes_allMkt <- t(replicate(nmkt, nodes))
    weights<-c(getWeights(grid))
  }

  # Gauss Hermite nested
  if (method == "sgNH") {
    if (is.na(accuracy)){
      stop("For method 'sgNH' the argument 'accuracyQuad' in the 'integration.control' needs to be specified")
    }

    grid <- createNIGrid(dim, "nHN", accuracy, ndConstruction = "sparse")
    nodes <- c(getNodes(grid))
    nodes_allMkt <- t(replicate(nmkt, nodes))
    weights<-c(getWeights(grid))
  }

  # Sparse Grids Trapezoidal
  if (method == "sgTr1") {
    if (is.na(accuracy)){
      stop("For method 'sgTr1' the argument 'accuracyQuad' in the 'integration.control' needs to be specified")
    }

    grid <- createNIGrid(dim, 'cNC1', level = accuracy, ndConstruction = "sparse", level.trans = TRUE)
    rescale(grid,cbind(rep(-3.5,dim),rep(+3.5,dim)))
    nodes <- getNodes(grid)
    nodes_allMkt <- t(replicate(nmkt, c(nodes)))
    weights <- getWeights(grid)
    tmp <- apply(nodes, MARGIN = 1,
                 FUN = function(x){
                   prod(dnorm(x))
                 })
    weights <- weights * tmp
  }


  if (method == "sgTr") {
    if (is.na(accuracy)){
      stop("For method 'sgTr' the argument 'accuracyQuad' in the 'integration.control' needs to be specified")
    }

    grid <- createNIGrid(dim, 'Trapez', level = accuracy, ndConstruction = "sparse", level.trans = TRUE)
    rescale(grid,cbind(rep(-3.5,dim),rep(+3.5,dim)))
    nodes <- getNodes(grid)
    nodes_allMkt <- t(replicate(nmkt, c(nodes)))
    weights <- getWeights(grid)
    tmp <- apply(nodes, MARGIN = 1,
                     FUN = function(x){
                       prod(dnorm(x))
                     })
    weights <- weights * tmp
  }


  output$nodesMktShape <- nodes_allMkt
  output$weights <- weights
  attributes(output)$amountNodes <- length(weights)

  return(output)
}




get.jacob <- function(sij,
                      blp.data,
                      blp.parameters,
                      blp.integration,
                      printLevel ) {

  nodesRcMktShape <- blp.integration$nodesRcMktShape

  qvt <- matrix(NA_real_, ncol = blp.integration$amountNodes * blp.parameters$K, nrow = 1)
  dt <- matrix(NA_real_, ncol = blp.integration$amountNodes * blp.parameters$total.dem, nrow = 1)

  # jacobian

  jacob.out <- matrix(NA_real_, nrow = blp.parameters$nobs, ncol = blp.parameters$total.par)

  for (i in 1:blp.parameters$nmkt) {

    market.identifier <- (blp.parameters$cdindex[i] + 1): blp.parameters$cdindex[i + 1]
    nprodt <- length(market.identifier)
    x2t <- matrix(NA_real_, nrow = nprodt, ncol = blp.parameters$K)
    sijt <- matrix(NA_real_, nrow = nprodt, ncol = blp.integration$amountNodes)
    dsdxi <- matrix(NA_real_, nrow = nprodt, ncol = nprodt)
    dsdtheta <- matrix(NA_real_, nrow = nprodt, ncol = blp.parameters$total.par)
    #tmp <- matrix(NA_real_, nrow = nprodt, ncol = total.par)
    qvt[ , ] <- nodesRcMktShape[i, ]
    if (blp.parameters$total.demogr > 0) {
       dt[ , ] <- blp.data$demographics[i, ]}

    x2t[ , ] <- blp.data$Xrandom[market.identifier, , drop = FALSE]  # w?hle die werte x2
    sijt[ , ] <- sij[market.identifier, ]
    dsdxi[ , ] <- .dstdxit(sijt,
                          blp.integration$weights)
    dsdtheta[ , ] <- .dstdtheta(sijt = sijt,
                               xt = x2t,
                               qvt = qvt,
                               dt = dt,
                               weights = blp.integration$weights,
                               blp.parameters = blp.parameters)


    tmp <- try(-solve(dsdxi) %*% dsdtheta, silent = T) # ;% tmp darf nicht initialisiert werden, da der Klassenbefehl mit try error nicht mehr geht

    if (class(tmp) == "try-error") {
      if (printLevel >= 2) {
        warning(paste("Error in jacobian (market ", i,
          ") : Singular matrix occured")) }
      return(NA)
    } else {
      jacob.out[market.identifier, ] <- tmp
    }

  }  # end market loop
  return(jacob.out)
}



get.gmm.gr <- function(theta2 ,
                       blp.integration,
                       blp.parameters,
                       blp.data,
                       blp.results,
                       printLevel) {
  if( blp.parameters$gradientIndicator == TRUE ) {
  return(blp.results$gradient) }else{
    return( NULL ) }

}




.dstdxit <- function(sijt, weights) {
  numProdt <- dim(sijt)[1]
  weightsMat <- matrix( weights, ncol=1)
  dsdxi <- matrix( NA_real_, nrow = numProdt, ncol = numProdt)

  dsdxi[ , ] <- -(t(weightsMat)[rep(1, numProdt), ]* sijt) %*%  t(sijt)
  diag(dsdxi) <- c( ( sijt * (1 - sijt) ) %*% weightsMat)
  return(dsdxi)
}


.dstdtheta <- function(sijt, xt, qvt, dt, weights, blp.parameters ) {
  numProdt <- dim(xt)[1]
  K <- dim(xt)[2]
  total.demogr <- blp.parameters$total.demogr
  total.parameters <- blp.parameters$total.par
  amountNodes <- dim(qvt)[2]/K
  weightsMat <- matrix( weights, ncol=1)

  # allocating matrices
  dsdtheta.out <- matrix(NA_real_, nrow = numProdt, ncol = total.parameters)
  sumterm <- vector(mode = "numeric", length = amountNodes)
  bracket <- matrix(NA_real_, nrow = numProdt, ncol = amountNodes)
  scalar <- matrix(NA_real_, nrow = numProdt, ncol = amountNodes)

  for (i in 1:K) {  # i iterates over Random Coefficients describing unobs. het.

    coef_for_i <- blp.parameters$indices[ ,1] == i # coefficients that belong to i (RC and demographics) are in line k of the parameter indices matrix

    with_unobs_het_i <- 1 %in% blp.parameters$indices[ coef_for_i, 2 ] # check if RC for unobs. het. is not NA (if so, there is no 1 as a col index)

    # the following 2 lines are needed anyway (even if RC for unobs. het = NA)
    sumterm[ ] <- t(xt[, i, drop = F]) %*% sijt  # sum over x_mt^k * s_mti (m=product)  for every person i (in cols)
    bracket[ , ] <- xt[, i] - matrix(sumterm, nrow = 1)[rep(1, numProdt), ]  # for a given k: substract that sum from a product characteristic (for every person i)

    if( with_unobs_het_i ){
        scalar[ , ] <- qvt[rep(1, numProdt), ((i - 1) * amountNodes +
                                                1):(i * amountNodes)] * sijt # calc draw_i*s_jti  for every person i (in cols)
        RCRow <- which( coef_for_i )[1] # first element decribes unobserved heterogeneity
        dsdtheta.out[, RCRow ] <- (scalar * bracket) %*% weightsMat
    }

    # Demographics (all interactions with Random Coefficient k)
    if( total.demogr > 0 ){

        if( with_unobs_het_i ){ # check if RC for unobs. het. is not NA (if so, there is no 1 as a col index)
            demographicsRow <- which( coef_for_i )[-1] # first one is for unobs. het.
        }else{
          demographicsRow <- which( coef_for_i )
        }

        relevantDemographicsForI <- blp.parameters$indices[ demographicsRow, 2] - 1 # Dem starts in sec. col (always)
        if( length( relevantDemographicsForI ) >= 1 ){
          for ( j in 1: length( relevantDemographicsForI ) ) {
          # sumterm and bracket can be used from former loop
            relevantDemographic <- relevantDemographicsForI[j]
            scalar[, ] <- dt[rep(1, numProdt), ((relevantDemographic - 1) * amountNodes + 1):
                                  (relevantDemographic * amountNodes)] * sijt  # calc draw_i*s_jti  for every person i (in cols)
            dsdtheta.out[, demographicsRow[j] ] <- (scalar * bracket) %*% weightsMat
          }}
        }
  }
  return(dsdtheta.out)
}



get.Elasticities <- function( changingVariable_name, #parameters...
                              market,

                              productData,           # data...
                              Xrandom_name,
                              delta,

                              cdid_name,
                              share_name,

                              all_randomCoef,
                              all_linearCoef,
                              nodesRCMktShape ,
                              nodesDemMktShape ,
                              Integration_Weights,
                              totalDemographics ,
                              indices_InParameters,
                              sij) {

  ## input checks
    cdid_all <- productData[ , cdid_name ]

    if ( !( changingVariable_name %in% names( productData ) ) )
      stop( paste( "Provided variable", changingVariable_name , "is not in the set of variables.") )

    if ( !( market %in% cdid_all ) )
      stop( paste( "Provided market", market , "is not included in cdid.") )


  ## cross elasticities calculation

  # VariableExtraction
  relevant_Obs_Mkt <- which( market ==  cdid_all )  # cdidOld contains the input from BLPestimator
  relevant_Mkt <- which( market ==  unique(cdid_all) )  # determines a relevant numeric market id

  nobs <- nrow( productData )
  nprod_Mkt <- length( relevant_Obs_Mkt )
  obsShare_Mkt <- productData[ relevant_Obs_Mkt , share_name ]

  betabar <- all_linearCoef[ changingVariable_name ]
  if( is.na( betabar ) ){
    # in case that the variable does not enter linerarly
    betabar <- 0
  }

  changingVariable_Mkt <- unlist( productData[ relevant_Obs_Mkt , changingVariable_name,
                                               drop = F ] )

  # if changingVariable is not modelled as a random coefficient, ie use logit:
  if ( !( changingVariable_name %in%
          Xrandom_name) ) {
    # Use Logit elasticities, because '" , changingVariable_name , "' is not modelled as random coefficient...

    # Cols contain the variables that are changed by 1%, and rows contain effects on other products in the choice set.

    EtaMkt <- - matrix( changingVariable_Mkt * betabar * obsShare_Mkt ,
                        nrow = nprod_Mkt,
                        ncol = nprod_Mkt, byrow = TRUE)

    diag(EtaMkt) <- betabar * changingVariable_Mkt * (1 - obsShare_Mkt)


  } else {
  # if changingVariable is modelled as a random coefficient:

    amountNodes <- ncol( sij )
    K <- ncol( nodesRCMktShape ) / amountNodes
    XrandomData_Mkt <- as.matrix( productData[ relevant_Obs_Mkt ,
                                               Xrandom_name ,
                                               drop = F ] )
    nodesRCMktShape_Mkt <- nodesRCMktShape[ relevant_Mkt , ,
                                            drop = FALSE] # only one line, bec. nodes are used for every product in one market
    nodesRC_tableform_Mkt <- matrix(nodesRCMktShape_Mkt,
                                    nrow = amountNodes,
                                    ncol = K )
    colnames(nodesRC_tableform_Mkt) <- Xrandom_name


    if ( totalDemographics > 0 ) {
      nodesDemMktShape_Mkt <- nodesDemMktShape[ relevant_Mkt, ,
                                                drop = FALSE]
      demographicReshape_Mkt <- matrix( nodesDemMktShape_Mkt ,
                                        ncol = totalDemographics,
                                        nrow = amountNodes )
    } else {
      nodesDemMktShape_Mkt <- matrix( NA )
      demographicReshape_Mkt <- matrix( NA )
    }

    cdid_Mkt <- rep(1, nprod_Mkt)

    sij_Mkt <- sij[ relevant_Obs_Mkt , ]
    weights <- matrix( Integration_Weights )
    expdelta_Mkt<- exp( delta[ relevant_Obs_Mkt ] )



    theta2Mat <- .get.theta2.reshape(theta2.in = c( all_randomCoef ),
                                     totalRC = K,
                                     total.demogr.in = totalDemographics,
                                     indices.in = indices_InParameters,
                                     fill = 0)  # NA are replaced by zeros to simplify x * par in getExpMu
    rownames(theta2Mat) <- Xrandom_name

    expmu_Mkt <- getExpMu(theta2Mat,
                          nodesRCMktShape_Mkt,
                          XrandomData_Mkt,
                          cdid_Mkt,
                          nodesDemMktShape_Mkt)



    sigma <- all_randomCoef[ changingVariable_name ]

    # unobsevered_part contains all individual specific effect parts due to unobs. herterog.
    unobseved_part <- betabar +
      sigma * nodesRC_tableform_Mkt[ , changingVariable_name ]  # get individual price effects

    if ( totalDemographics > 0 ) {
      # extract relevant demographic effects ( ie pick a line of theta2Mat)
      demographic_effects <- matrix( theta2Mat[ changingVariable_name, -1],
                                     ncol = totalDemographics,
                                     nrow = amountNodes,
                                     byrow = TRUE )

      # multiply every demographic coefficient with all demogr. draws:
      observed_part <- rowSums( demographic_effects *
                                    demographicReshape_Mkt)

    } else {
      observed_part <- 0
    }

    beta_i <- observed_part + unobseved_part

    Omega <- matrix( NA_real_ ,
                     nrow = nprod_Mkt,
                     ncol = nprod_Mkt)

    scalar <- - matrix( 1/obsShare_Mkt ) %*% matrix( changingVariable_Mkt , nrow = 1)

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
  rownames( EtaMkt ) <- colnames( EtaMkt ) <- paste0("product", 1:nprod_Mkt)

  return( EtaMkt )


}

