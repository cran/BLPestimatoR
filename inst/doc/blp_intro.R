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
head(productData_cereal)

## ------------------------------------------------------------------------
demographicData_cereal$income[1:4, 1:5]

demographicData_cereal$incomesq[1:4, 1:5]

demographicData_cereal$age[1:4, 1:5]

demographicData_cereal$child[1:4, 1:5]

## ------------------------------------------------------------------------
originalDraws_cereal$constant[1:4, 1:5]

# renaming constants:
names(originalDraws_cereal)[1] <- "(Intercept)"

originalDraws_cereal$price[1:4, 1:5]

originalDraws_cereal$sugar[1:4, 1:5]

originalDraws_cereal$mushy[1:4, 1:5]

## ------------------------------------------------------------------------
productData_cereal$startingGuessesDelta <- c(log(w_guesses_cereal)) # include orig. draws in the product data

cereal_data <- BLP_data(
  model = nevos_model,
  market_identifier = "cdid",
  par_delta = "startingGuessesDelta",
  product_identifier = "product_id",
  productData = productData_cereal,
  demographic_draws = demographicData_cereal,
  blp_inner_tol = 1e-6, blp_inner_maxit = 5000,
  integration_draws = originalDraws_cereal,
  integration_weights = rep(1 / 20, 20)
)

## ------------------------------------------------------------------------
# before:
theta_guesses_cereal
theta_guesses_cereal[theta_guesses_cereal == 0] <- NA
colnames(theta_guesses_cereal) <- c("unobs_sd", "income", "incomesq", "age", "child")
rownames(theta_guesses_cereal) <- c("(Intercept)", "price", "sugar", "mushy")

# correctly named:
theta_guesses_cereal

## ------------------------------------------------------------------------
cereal_est <- estimateBLP(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  solver_method = "BFGS", solver_maxit = 1000, solver_reltol = 1e-6,
  standardError = "heteroskedastic",
  extremumCheck = FALSE,
  printLevel = 1
)

summary(cereal_est)

## ------------------------------------------------------------------------
cereal_data2 <- BLP_data(
  model = nevos_model,
  market_identifier = "cdid",

  product_identifier = "product_id",
  productData = productData_cereal,
  integration_method = "MLHS",
  integration_accuracy = 20, integration_seed = 213
)

cereal_est2 <- estimateBLP(blp_data = cereal_data2, printLevel = 1)

summary(cereal_est2)

## ------------------------------------------------------------------------
# extract parameters from output
theta1_price <- cereal_est$theta_lin["price", ]
theta2 <- matrix(NA, nrow = 4, ncol = 5)
colnames(theta2) <- c("unobs_sd", "income", "incomesq", "age", "child")
rownames(theta2) <- c("(Intercept)", "price", "sugar", "mushy")
for (i in 1:13) {
  theta2[cereal_est$indices[i, 1], cereal_est$indices[i, 2]] <- cereal_est$theta_rc[i]
}

delta_data <- data.frame(
  "product_id" = cereal_data$parameters$product_id,
  "cdid" = cereal_data$parameters$market_id_char_in,
  "startingGuessesDelta" = cereal_est$delta
)
# always use update_BLP_data() to update data object to maintain consistent data
cereal_data <- update_BLP_data(
  data_update = delta_data,
  blp_data = cereal_data
)

shareObj <- getShareInfo(
  blp_data = cereal_data,
  par_theta2 = theta2,
  printLevel = 1
)

get_elasticities(
  blp_data = cereal_data,
  share_info = shareObj,
  theta_lin = theta1_price,
  variable = "price",
  products = c("cereal_1", "cereal_4"),
  market = "market_2"
)

## ------------------------------------------------------------------------
delta_eval <- getDelta_wrap(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  printLevel = 4
)

productData_cereal$startingGuessesDelta[1:6]
delta_eval$delta[1:6]
delta_eval$counter

gmm <- gmm_obj_wrap(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  printLevel = 2
)

gmm$local_min

## ------------------------------------------------------------------------
shareObj <- getShareInfo(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  printLevel = 4
)

shareObj$shares[1:6]

## ------------------------------------------------------------------------
# market 2:
derivatives1 <- dstdtheta_wrap(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  market = "market_2"
)
derivatives2 <- dstddelta_wrap(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  market = "market_2"
)

jac_mkt2 <- -solve(derivatives2) %*% derivatives1

jac_mkt2[1:5, 1:4]

# all markets
jacobian_nevo <- getJacobian_wrap(
  blp_data = cereal_data,
  par_theta2 = theta_guesses_cereal,
  printLevel = 2
)

jacobian_nevo[25:29, 1:4] # compare to jac_mkt2

## ------------------------------------------------------------------------
# add owner matix to productData
own_pre <- dummies_cars
colnames(own_pre) <- paste0("company", 1:26)
productData_cars <- cbind(productData_cars, own_pre)

# construct instruments
nobs <- nrow(productData_cars)
X <- data.frame(
  productData_cars$const, productData_cars$hpwt,
  productData_cars$air, productData_cars$mpg, productData_cars$space
)

sum_other <- matrix(NA, nobs, ncol(X))
sum_rival <- matrix(NA, nobs, ncol(X))
sum_total <- matrix(NA, nobs, ncol(X))

for (i in 1:nobs) {
  other_ind <- productData_cars$firmid == productData_cars$firmid[i] &
    productData_cars$cdid == productData_cars$cdid[i] &
    productData_cars$id != productData_cars$id[i]
  rival_ind <- productData_cars$firmid != productData_cars$firmid[i] &
    productData_cars$cdid == productData_cars$cdid[i]
  total_ind <- productData_cars$cdid == productData_cars$cdid[i]

  sum_other[i, ] <- colSums(X[other_ind == 1, ])
  sum_rival[i, ] <- colSums(X[rival_ind == 1, ])
  sum_total[i, ] <- colSums(X[total_ind == 1, ])
}

colnames(sum_other) <- paste0("IV", 1:5)
colnames(sum_rival) <- paste0("IV", 6:10)
productData_cars <- cbind(productData_cars, sum_other, sum_rival)
head(productData_cars)

# To show similarities between implementations of other authors,
# the variable "const" is used, although constants are considered by default.
blps_model <- as.formula("share ~  0 + const + price + hpwt + air + mpg + space |
                        0 + const + hpwt + air + mpg + space |
                        0 + price + const + hpwt + air + mpg |
                        0 + IV1 + IV2 + IV3 + IV4 + IV5 + IV6 + IV7 + IV8 + IV9 + IV10")

car_data <- BLP_data(
  model = blps_model,
  market_identifier = "cdid",
  product_identifier = "id",
  additional_variables = paste0("company", 1:26), # check reordering works
  productData = productData_cars,
  blp_inner_tol = 1e-9,
  blp_inner_maxit = 5000,
  integration_method = "MLHS",
  integration_accuracy = 50, integration_seed = 48
)

## ------------------------------------------------------------------------
set.seed(121)
theta_guesses <- matrix(rnorm(5))
rownames(theta_guesses) <- c("price", "const", "hpwt", "air", "mpg")
colnames(theta_guesses) <- "unobs_sd"


car_est <- estimateBLP(
  blp_data = car_data,
  par_theta2 = theta_guesses,
  solver_method = "BFGS", solver_maxit = 1000, solver_reltol = 1e-6,
  extremumCheck = FALSE, printLevel = 0
)

summary(car_est)

## ------------------------------------------------------------------------
## Pre-Merger data
own_pre <- as.matrix(car_data$data$additional_data[, paste0("company", 1:26)])
delta_pre <- car_est$delta
theta1_price <- car_est$theta_lin["price", ]
theta2_price <- car_est$theta_rc["unobs_sd*price"]
theta2_all <- matrix(car_est$theta_rc)
rownames(theta2_all) <- c("price", "const", "hpwt", "air", "mpg")
colnames(theta2_all) <- "unobs_sd"

## update mean utility in data ( always use update_BLP_data() to update data object to maintain consistent data )
delta_data <- data.frame(
  "id" = car_data$parameters$product_id,
  "cdid" = car_data$parameters$market_id,
  "delta" = delta_pre
)
car_data_updated <- update_BLP_data(
  data_update = delta_data,
  blp_data = car_data
)

## ------------------------------------------------------------------------
## calculate sij
shareObj <- getShareInfo(
  blp_data = car_data_updated,
  par_theta2 = theta2_all,
  printLevel = 0
)

## computation of marginal costs
market_id <- car_data$parameters$market_id
nmkt <- length(unique(market_id))
markups <- numeric(length(market_id))

sh <- shareObj$shares
prices_pre <- car_data$data$X_rand[, "price"]


for (i in 1:nmkt) {
  mkt_ind <- market_id == i
  share_i <- sh[ mkt_ind ]
  price_pre_i <- prices_pre[ mkt_ind ]
  scalar_i <- matrix(1 / share_i) %*% matrix(price_pre_i, nrow = 1)
  elasticities_i <- get_elasticities(
    blp_data = car_data_updated,
    share_info = shareObj,
    theta_lin = theta1_price,
    variable = "price",
    market = i,
    printLevel = 0
  )

  derivatives_i <- elasticities_i / scalar_i # partial derivatives of shares wrt price
  own_pre_i <- own_pre[ mkt_ind, ]
  own_prod_pre_i <- own_pre_i %*% t(own_pre_i) # if element (i,j) equals 1, that means that prod i and j are produced by same firm
  markups[mkt_ind] <- c(-solve(t(derivatives_i) * own_prod_pre_i) %*% share_i)
}
marg_cost <- prices_pre - markups

## ------------------------------------------------------------------------
# Merger between company 16 and 19 (i.e. GM and Chrysler)
prices_post <- numeric(2217)
own_post <- cbind(
  own_pre[, 1:15],
  own_pre[, 16] + own_pre[, 19],
  own_pre[, 17:18],
  own_pre[, 20:26]
)

## ---- eval = FALSE-------------------------------------------------------
#  foc_bertrand_mkt <- function(par, own_prod, blp_data, mkt, marg_cost, theta_lin, theta_rc) {
#    # argument par: candidate for post merger prices
#    # arguments own_prod, blp_data, mkt, marg_cost, theta_lin, theta_rc: see previous code blocks
#  
#    # post merger updates: update the BLP_data object for market i
#    tmp <- data.frame(
#      "id" = blp_data$parameters$product_id,
#      "cdid" = blp_data$parameters$market_id,
#      "delta" = blp_data$data$delta,
#      "price" = blp_data$data$X_rand[, "price"]
#    )
#  
#    market_ind <- blp_data$parameters$market_id == mkt
#    delta_old <- blp_data$data$delta
#    prices_pre <- blp_data$data$X_rand[, "price"]
#    tmp$price[ market_ind ] <- par
#    tmp$delta[ market_ind ] <- delta_old[market_ind] - prices_pre[market_ind] * theta_lin + par * theta_lin
#  
#  
#    new_blp_data <- update_BLP_data(
#      blp_data = blp_data,
#      data_update = tmp
#    )
#  
#    ShareObj <- getShareInfo(
#      blp_data = new_blp_data,
#      par_theta2 = theta_rc,
#      printLevel = 0
#    )
#  
#    implied_shares <- as.matrix(ShareObj$shares[market_ind])
#  
#    elasticities_post_mkt <- get_elasticities(
#      blp_data = new_blp_data,
#      share_info = ShareObj,
#      theta_lin = theta_lin,
#      variable = "price",
#      market = mkt,
#      printLevel = 0
#    )
#  
#    scalar_mkt <- matrix(1 / implied_shares) %*% matrix(par, nrow = 1)
#    derivatives_mkt <- elasticities_post_mkt / scalar_mkt
#  
#    markups_post <- c(-solve(t(derivatives_mkt) * own_prod) %*% implied_shares)
#    differences <- par - marg_cost[market_ind] - markups_post
#  
#    return(differences)
#  }

## ----eval=FALSE----------------------------------------------------------
#  library(nleqslv) # to solve non linear first order conditions
#  for (i in 1:nmkt) {
#    mkt_ind <- market_id == i
#    own_post_i <- own_post[ mkt_ind, ]
#    own_prod_post_i <- own_post_i %*% t(own_post_i)
#    price_pre_i <- prices_pre[ mkt_ind ]
#  
#    solution <- nleqslv(
#      x = price_pre_i, foc_bertrand_mkt, # startingguesses: price_pre_i
#      own_prod = own_prod_post_i,
#      blp_data = car_data_updated,
#      mkt = i,
#      marg_cost = marg_cost,
#      theta_lin = theta1_price,
#      theta_rc = theta2_all
#    )
#  
#    prices_post[ market_id == i ] <- solution$x
#  }

