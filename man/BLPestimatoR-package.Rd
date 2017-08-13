\name{BLPestimatoR-package}
\alias{BLPestimatoR-package}
\alias{BLPestimatoR}
\docType{package}
\title{
  BLP demand estimation for differentiated products
}
\description{
  Provides the estimation algorithm to perform the demand estimation described in Berry, Levinsohn and Pakes (1995) <DOI:10.2307/2171802>. The routine uses analytic gradients and offers a large number of implemented integration methods and optimization routines.
}

\details{
  \tabular{ll}{
    Package: \tab BLPestimatoR\cr
    Type: \tab Package\cr
    Version: \tab 0.1.5\cr
    Date: \tab 2017-08-12\cr
    License: \tab GPL-3\cr
  }
}
\author{
  Daniel Brunner (HHU of Duesseldorf / Germany)

  Constantin Weiser (HHU of Duesseldorf / Germany)

  Andre Romahn (HHU of Duesseldorf / Germany)

  Maintainer: Daniel Brunner <daniel.brunner@hhu.de>
}
\references{
  Steven Berry, James Levinsohn, Ariel Pakes (1995): Automobile Prices in Market Equilibrium <DOI:10.2307/2171802>


  Christopher R. Knittel,  Konstantinos Metaxoglou (2014): Estimation of Random-Coefficient Demand Models: Two Empiricists Perspective <DOI:10.1162/REST_a_00394>
}
\keyword{ package }
\keyword{ math }
\examples{
# Parameters
 i<-1
 K<-2
 Xlin_example <-  c("price", "x1", "x2", "x3", "x4", "x5")
 Xexo_example <- c("x1", "x2", "x3", "x4", "x5")
 Xrandom_example <- paste0("x",1:K)
 instruments_example <- paste0("iv",1:10)

 # Data generation
 BLP_data <- get.BLP.dataset(nmkt = 25, nbrn = 20,
                             Xlin = Xlin_example,
                             Xexo = Xexo_example,
                             Xrandom = Xrandom_example,
                             instruments = instruments_example,
                             true.parameters = list(Xlin.true.except.price = rep(0.2,5),
                                                    Xlin.true.price = -0.2, Xrandom.true = rep(2,K),
                                                    instrument.effects = rep(2,10),
                                                    instrument.Xexo.effects = rep(1,5)),
                             price.endogeneity = list( mean.xi = -2,
                                                       mean.eita = 0,
                                                       cov = cbind( c(1,0.7), c(0.7,1))),
                             printlevel = 0, seed = 5326 )

 # Estimation
 BLP_est<- estimateBLP(Xlin = Xlin_example,
                       Xrandom = Xrandom_example,
                       Xexo =  Xexo_example,
                       instruments = instruments_example,
                       shares = "shares",
                       cdid = "cdid",
                       productData = BLP_data,
                       starting.guesses.theta2 = rep(1,K),
                       solver.control = list(maxeval = 5000),
                       solver.method = "BFGS_matlab",

                       starting.guesses.delta =  rep(1, length(BLP_data$cdid)),
                       blp.control = list(inner.tol = 1e-6,
                                          inner.maxit = 5000),
                       integration.control= list(  method="MLHS",
                                                   amountNodes= 100,
                                                   seed= 3   ),
                       postEstimation.control= list(standardError = "robust",
                                                    extremumCheck = TRUE,
                                                    elasticities = "price"),
                       printLevel = 2)

 # Show results
 summary(BLP_est)
}
