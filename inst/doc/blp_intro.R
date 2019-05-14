## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(BLPestimatoR)

## ------------------------------------------------------------------------
nevos_model <- as.formula("share ~  price + productdummy |
    0+ productdummy |
    price + sugar + mushy |
    0+ IV1 + IV2 + IV3 + IV4 + IV5 + IV6 + IV7 + IV8 + IV9 + IV10 + 
    IV11 + IV12 + IV13 + IV14 + IV15 + IV16 + IV17 + IV18 + IV19 + IV20")

## ------------------------------------------------------------------------
head(productData)

## ------------------------------------------------------------------------
demographicData$income[1:4, 1:5]

demographicData$incomesq[1:4, 1:5]

demographicData$age[1:4, 1:5]

demographicData$child[1:4, 1:5]

## ------------------------------------------------------------------------
originalDraws$constant[1:4,1:5]

# renaming constants:
names(originalDraws)[1] <- "(Intercept)"

originalDraws$price[1:4,1:5]

originalDraws$sugar[1:4,1:5]

originalDraws$mushy[1:4,1:5]

## ------------------------------------------------------------------------
productData$startingGuessesDelta  <- c( log(w_guesses)) # include orig. draws in the product data

nevo_data <- BLP_data(model = nevos_model,
                      market_identifier="cdid",
                      par_delta = "startingGuessesDelta",
                      product_identifier = "product_id",
                      productData = productData,
                      demographic_draws = demographicData,
                      blp_inner_tol = 1e-6, blp_inner_maxit = 5000,
                      integration_draws = originalDraws, 
                      integration_weights= rep(1/20,20))

## ------------------------------------------------------------------------
#before:
theta_guesses
theta_guesses[theta_guesses==0] <- NA 
colnames(theta_guesses) <- c("unobs_sd", "income", "incomesq", "age", "child" )
rownames(theta_guesses) <- c("(Intercept)", "price" , "sugar", "mushy")

# correctly named:
theta_guesses

## ------------------------------------------------------------------------

blp_est <- estimateBLP( blp_data=nevo_data,
                        par_theta2 = theta_guesses,
                        solver_method = "BFGS", solver_maxit = 1000, solver_reltol = 1e-6,
                        standardError = "heteroskedastic",
                        extremumCheck = FALSE ,
                        printLevel = 1 )

summary(blp_est)

## ------------------------------------------------------------------------
nevo_data2 <- BLP_data(model = nevos_model,
                      market_identifier="cdid",
                     
                      product_identifier = "product_id",
                      productData = productData,
                      integration_method = "MLHS", 
                      integration_accuracy = 20, integration_seed = 213)

blp_est2 <- estimateBLP( blp_data=nevo_data2, printLevel = 1 )

summary(blp_est2) 

## ------------------------------------------------------------------------
# extract parameters from output
theta1_price <- blp_est$theta_lin["price",]
theta2 <- matrix(NA , nrow = 4, ncol = 5)
colnames(theta2) <- c("unobs_sd", "income", "incomesq", "age", "child" )
rownames(theta2) <- c("(Intercept)", "price" , "sugar", "mushy")
for( i in 1:13){
theta2[blp_est$indices[i,1], blp_est$indices[i,2]] <- blp_est$theta_rc[i]
}

delta_data <- data.frame( "product_id" = nevo_data$parameters$product_id,
                     "cdid" = nevo_data$parameters$market_id_char_in,
                     "startingGuessesDelta" = blp_est$delta )
# always use update_BLP_data() to update data object to maintain consistent data 
nevo_data <- update_BLP_data(data_update = delta_data, 
                             blp_data = nevo_data)

shareObj <- getShareInfo(  blp_data=nevo_data, 
                            par_theta2 = theta2,
                            printLevel = 1)


get_elasticities(blp_data=nevo_data, 
                 share_info = shareObj, 
                 theta_lin = theta1_price,
                 variable = "price",
                 products = c("cereal_1","cereal_4"),
                 market = "market_2")


## ------------------------------------------------------------------------

delta_eval <- getDelta_wrap(
  blp_data = nevo_data,
  par_theta2 = theta_guesses,
  printLevel = 4
)

productData$startingGuessesDelta[1:6]
delta_eval$delta[1:6]
delta_eval$counter

gmm <- gmm_obj_wrap(
  blp_data = nevo_data,
  par_theta2 = theta_guesses,
  printLevel = 2
)

gmm$local_min

## ------------------------------------------------------------------------
shareObj <- getShareInfo(  blp_data=nevo_data,
                           par_theta2 = theta_guesses,
                           printLevel = 4)
shareObj$shares[1:6]

## ------------------------------------------------------------------------
# market 2:
derivatives1 <- dstdtheta_wrap(  blp_data=nevo_data,
                                par_theta2 = theta_guesses,
                                 market = "market_2")
derivatives2 <- dstddelta_wrap(  blp_data=nevo_data,
                                par_theta2 = theta_guesses,
                                 market = "market_2")

jac_mkt2 <- -solve(derivatives2)%*%derivatives1

jac_mkt2[1:5,1:4]

# all markets
jacobian_nevo <- getJacobian_wrap(blp_data=nevo_data,
                                 par_theta2 = theta_guesses,
                                 printLevel = 2)

jacobian_nevo[25:29,1:4] # compare to jac_mkt2


