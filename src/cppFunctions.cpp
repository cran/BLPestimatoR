#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// Includes/namespaces
#include <iostream>
#include <string>


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix  getSij(const NumericMatrix &expmu,
                      const NumericVector &expdelta,
                      const IntegerVector &cdindex){

  int nmkt = cdindex.size()-1;
  int amountNodes = expmu.ncol();
  int nobs = expmu.nrow();
  double sumMktNodes = 0;
  int nprodt;
  int startpos;
  NumericMatrix sij(nobs,amountNodes);

  for(int i=0;i<nmkt;i++){  // markets

    nprodt = cdindex[i+1]-cdindex[i] ;
    startpos = cdindex[i]+1;

    for( int z=0; z<amountNodes;z++) { // draws
      sumMktNodes = 0;
      for( int j=0;j<nprodt;j++) { // products
        sumMktNodes += expmu(startpos+(j-1),z)*expdelta[startpos+(j-1)];
      }

      for( int j=0;j<nprodt;j++) { // products
        sij(startpos+(j-1),z)= expmu(startpos+(j-1),z)*expdelta[startpos+(j-1)]/(1+sumMktNodes);
      }
    }
  }


  return sij;
}



// [[Rcpp::export]]
NumericVector  getSjtMod( const NumericMatrix &expmu,
                         const NumericVector &expdelta,
                         const int &nprodt,
                         const int &startpos,
                         const NumericVector &weights){


  int amountNodes = expmu.ncol();

  double sumMktNodes;

  double si0_mod;
  NumericVector sj_mod(nprodt);

  for( int z=0; z<amountNodes;z++) { //iterate over draws
      si0_mod = 0 ;
      sumMktNodes = 1; // one comes from logit denominator
      for( int j=0;j<nprodt;j++) { // iterate over products
        // for a given market and individual: sum over utility components
        sumMktNodes += expmu(startpos+(j-1),z)*expdelta[startpos+(j-1)];
      }
      // basically the individual prob. for the outside option times a weight:
      si0_mod =  weights[z]/sumMktNodes;

      for( int j=0;j<nprodt;j++) { // iterate over products
        sj_mod[j] += expmu(startpos+(j-1),z)*si0_mod;
      }

  }

  return sj_mod;
}



// [[Rcpp::export]]
NumericMatrix  getExpMu(const NumericMatrix &theta2Matrix,
                        const NumericMatrix &qv,
                        const NumericMatrix &Xrandom,
                        const IntegerVector &cdid ,
                        const NumericMatrix &demographics){
  int K = Xrandom.ncol();
  int totalDem = theta2Matrix.ncol()-1; //1 is RC column
  int nobs = Xrandom.nrow();
  int amountNodes = qv.ncol()/K;
  int startpointRC ;
  int startpointDem ;
  int mktInd;
  double demPart = 0;
  NumericMatrix expmu(nobs,amountNodes);


  for( int i=0; i<K;i++){ // iteriere ueber RC
    startpointRC= i*amountNodes;

    for( int z=0; z<nobs; z++){
      mktInd=  cdid[z] - 1;
      for( int r=0; r<amountNodes; r++){

        if(totalDem > 0){
          demPart=0;
        for ( int d=0; d<totalDem; d++){
        startpointDem = d*amountNodes;
        demPart += demographics(mktInd, startpointDem + r) * theta2Matrix(i,d+1) ;
        }}

        expmu(z,r) += (qv(mktInd,startpointRC+r) * theta2Matrix(i,0) + demPart) *Xrandom(z,i);
      }

   }
  }


  for( int z=0; z<nobs; z++){
    for( int r=0; r<amountNodes; r++){
      expmu(z,r) = exp( expmu(z,r) );
   }
  }


  return expmu;
}


// [[Rcpp::export]]
List getDelta( const NumericMatrix &theta2,
               const NumericVector &deltaOld,

               const IntegerVector &cdid,
               const IntegerVector &cdindex,
               const NumericMatrix &Xrandom,
               const NumericVector &obsshare,
               const double &innerCrit,
               const int &innerMaxit,
               const int &printLevel,
               const NumericMatrix &indices,
               const NumericMatrix &nodesRcMktShape,
               const NumericMatrix &nodesDemMktShape,
               const NumericVector &weights ){

  // initialising objects
  int counter = 0;
  int nobs = cdid.size();
  int nmkt = cdindex.size()-1;
  double dist = 100.0;
  NumericVector dist_m(nmkt,100.0);
  NumericVector market_convergence(nmkt);

  int amountRC = theta2.nrow();
  int amountNodes = nodesRcMktShape.ncol() / amountRC; // also possible with demopgraphics


  NumericMatrix expMu( nobs , amountNodes );
  NumericMatrix sij( nobs , amountNodes );
  NumericVector expDeltaStart = exp( deltaOld );
  NumericVector expDeltaNew( nobs );
  bool negShares = FALSE;
  bool negWeights = min(weights) < 0 ;

  int nprodt;
  int startpos;

  // input check delta:
  if( any( is_na( deltaOld ) ) | any( is_nan( deltaOld ) ) ) {
    for(int i = 0; i < nobs; ++i) {
      expDeltaStart[i] = 1;
    }
  }

  // calc. first part of utility
  expMu = getExpMu( theta2,
                    nodesRcMktShape,
                    Xrandom,
                    cdid,
                    nodesDemMktShape) ;

  // start of the contraction mapping
  if( printLevel == 4 )  Rcpp::Rcout << "----------------------" << std::endl ;

  while( (dist > innerCrit) & (counter < innerMaxit) ){
    counter++;
    for( int i=0; i<nmkt; i++){

      if( market_convergence[i] == 1 )  continue ;
      nprodt = cdindex[i+1]-cdindex[i] ;
      startpos = cdindex[i]+1;

      // calc. sjt_mod: (_mod because of skipped delta multiplication)
      NumericVector sjt_mod(nprodt);
      sjt_mod = getSjtMod(  expMu,
                            expDeltaStart,
                            nprodt,
                            startpos,
                            weights );

      // update delta(I):
      for( int j=0; j<nprodt; j++){
        expDeltaNew[startpos+(j-1)] = obsshare[startpos+(j-1)] / sjt_mod[j];
      }

      //Handling of negative shares for quadrature
      if( negWeights ) {
        if( Rcpp::min(sjt_mod) < 0 ){ //if( any( shMod < 0 ).is_true() ) {
          if( !negShares ){
            Rcpp::Rcout << "Integration rule produced negative shares --> reset delta to 0" ;
            negShares = TRUE ;
          }
          Rcpp::Rcout << "." ;

          for(int negInd = 0; negInd < nprodt; negInd++) {
            if( sjt_mod[negInd] < 0 ){
              expDeltaNew[startpos + (negInd-1) ] = 1;
            }
          }
        }
      } //end negWeights

      // update delta(II):
      dist_m[i] = Rcpp::max(  Rcpp::abs( expDeltaNew  - expDeltaStart ));

      // 1 indicates that market converged
      if( dist_m[i] < innerCrit) market_convergence[i] = 1;

      // avoids pointers and creates a deep copy:
      for( int j=0; j<nprodt; j++){
        expDeltaStart[startpos+(j-1)] = expDeltaNew[startpos+(j-1)];
      }
    } // end market loop

    dist = Rcpp::max(dist_m);

    if( printLevel == 4 ) Rcpp::Rcout << "\t dist: "  << dist << std::endl ;
  }// end while


  if( negShares ) Rcpp::Rcout << std::endl ;

  // check if maxit was reached:
  if( counter == innerMaxit){
    List ret;
    ret["delta"] = NAN;
    ret["sij"] = NAN;
    ret["counter"] = counter;
    ret["negShares"] = negShares;
    return ret ;
  }

  // calc sij for gradient (one single time)
  sij = getSij(expMu,
               expDeltaStart,
               cdindex );

  List ret;
  ret["delta"] = log(expDeltaNew);
  ret["sij"] = sij;
  ret["counter"] = counter;
  ret["negShares"] = negShares;
  return ret ;
}

// [[Rcpp::export]]
arma::mat dstddelta_c( arma::mat &sijt,
                     arma::mat &weights ) { // weight are provided as a [amountNodes X 1] matrix

  int numProdt = sijt.n_rows;
  arma::mat Dsdxi(numProdt, numProdt);
  arma::mat weightsMat = repmat( weights,
                                 1, //num_copies_per_row
                                 numProdt); //num_copies_per_col
  Dsdxi = -(( weightsMat.t() % sijt ) * sijt.t() );
  Dsdxi.diag() = (sijt % (1 - sijt)) * weights ;

  return ( Dsdxi );
}




// [[Rcpp::export]]
arma::mat dstdtheta_c(arma::mat &sijt_arma,
                      const NumericMatrix &indices,
                      arma::mat &xt_arma,
                      arma::mat &qvt_arma,
                      arma::mat &dt_arma,
                      arma::mat &weights_arma ){ // weight are provided as a [amountNodes X 1] matrix
  int numProdt = xt_arma.n_rows;
  int K = xt_arma.n_cols;
  int total_parameter = indices.nrow();
  //int total_demographics = max( indices(_,1) ) -1;
  int amountNodes = qvt_arma.n_cols / K;


  arma::mat dsdtheta_out(numProdt, total_parameter);
  arma::vec sumterm(amountNodes);
  arma::mat bracket_arma(numProdt, amountNodes);
  arma::mat scalar_arma(numProdt, amountNodes);

  arma::mat xti_arma(numProdt,
                     1);
  arma::mat xtiLarge_arma(numProdt,
                          amountNodes);
  arma::mat dsdtheta_arma( numProdt,
                           total_parameter );


  for (int i = 0; i < K; i++) {

    // At the C-level, all R objects are stored in a common datatype,
    // the SEXP, or S-expression. All R objects are S-expressions so every C function
    // that you create must return a SEXP as output and take SEXPs as inputs. (Technically,
    // this is a pointer to a structure with typedef SEXPREC.) A SEXP is a variant type,
    // with subtypes for all R’s data structures.

    // coef_for_i <- blp.parameters$indices[ ,1] == i
    // // NumericMatrix indices =  blp_parameters["indices"] ;
    //NumericVector row_indices = indices(_, 0);
    arma::vec row_indices = indices(_, 0);
    arma::vec col_indices = indices(_, 1);

    // get rows of "indices" that belong to random coefficient i
    arma::uvec row_indices_i = arma::find( row_indices == (i+1) );

    // check column indizes for random coef. i for 1 (=first column index in R)
    bool with_unobs_het_i = min( col_indices( row_indices_i ) ) == 1;
    // check column indizes for random coef. i for values > 1: interaction with demogr.
    bool with_obs_het_i = max( col_indices( row_indices_i ) ) > 1;

    xti_arma =  xt_arma.col(i);
    arma::mat sumterm_arma = trans( xti_arma ) * sijt_arma;    // (m=product)  for every person i (in cols)

    xtiLarge_arma = repmat( xti_arma,
                            1, //num_copies_per_row
                            amountNodes); //num_copies_per_col
    bracket_arma = xtiLarge_arma.each_row() - sumterm_arma;

    if (with_unobs_het_i) {
      scalar_arma = sijt_arma.each_row() % qvt_arma.submat(  0,                          // first_row
                                       i  * amountNodes,           // first_col
                                       0,                          // last_row
                                       (i+1) * amountNodes - 1  ); // last_col

      dsdtheta_arma.col( row_indices_i(0) ) = (scalar_arma % bracket_arma) * weights_arma;
    }

    // Demographics:
    if( with_obs_het_i ){
      int relevantDemographic_ij;
      int startDem;

      if( with_unobs_het_i ){
        startDem=1;
      }else{
        startDem=0;
      }

      arma::uvec demographicsRow = row_indices_i.subvec(startDem,
                                                        row_indices_i.size() -1 );

      arma::vec relevantDemographics_i = col_indices.elem(demographicsRow) - 1 ; //Dem starts in sec. col (always)

      int amountDemogr_i = relevantDemographics_i.size();

      if( amountDemogr_i >= 1 ){

        for (int j = 0; j < amountDemogr_i ; j++) {

          relevantDemographic_ij = relevantDemographics_i( j ) - 1 ;

          scalar_arma = sijt_arma.each_row() % dt_arma.submat(  0,     // first_row
                                           (relevantDemographic_ij)  * amountNodes,           // first_col
                                           0,                          // last_row
                                           (relevantDemographic_ij+1) * amountNodes -1   );

          dsdtheta_arma.col(demographicsRow(j)) = (scalar_arma % bracket_arma) * weights_arma;

        }
      } //end relevantDemographicsForI loop
    }
  } // end i loop

  return ( dsdtheta_arma );
}




// [[Rcpp::export]]
NumericMatrix jacob_c( NumericMatrix &sij,
                       const NumericMatrix &indices,
                       const List &blp_data,
                       const List &blp_parameters,
                       const List &blp_integration,
                       const int &printLevel) {


  int nobs = (int)blp_parameters["nobs"];
  int K = (int)blp_parameters["K"];
  int total_par = indices.nrow();
  int total_dem =  max( indices(_,1) ) -1;
  int nmkt = (int)blp_parameters["nmkt"];
  Rcpp::NumericVector cdindex((SEXP)blp_parameters["cdindex"]);

  NumericMatrix Xrandom((SEXP)blp_data["X_rand"]);

  int amountDraws = (int)blp_integration["amountDraws"];
  NumericMatrix nodesRcMktShape((SEXP)blp_integration["drawsRcMktShape"]);
  NumericMatrix nodesDemMktShape((SEXP)blp_integration["drawsDemMktShape"]);
  Rcpp::NumericVector weights((SEXP)blp_integration["weights"]);

  arma::mat Xrandom_arma( Xrandom.begin(),
                          Xrandom.nrow(),
                          Xrandom.ncol(),
                          false);
  arma::mat nodesRcMktShape_arma(nodesRcMktShape.begin(),
                                 nodesRcMktShape.nrow(),
                                 nodesRcMktShape.ncol(),
                                 false);
  arma::mat nodesDemMktShape_arma(nodesDemMktShape.begin(),
                                  nodesDemMktShape.nrow(),
                                  nodesDemMktShape.ncol(),
                                  false);
  arma::mat sij_arma(sij.begin(), // Fehler wenn const vor sij Argument??
                     sij.nrow(),
                     sij.ncol(),
                     false);

  arma::mat weights_arma(weights.begin(),
                         weights.size(), //row
                         1, //col
                         false);

  arma::mat qvt_arma(1, K * amountDraws);
  arma::mat dt_arma(1, amountDraws * total_dem);
  arma::mat jacob_arma( nobs ,
                        total_par);

  // iterate over markets:
  for (int t = 0; t < nmkt; t++) {

    arma::mat x2t_arma = Xrandom_arma.submat(  cdindex(t),               // first_row
                                               0,                        // first_col
                                               cdindex(t+1)-1,             // last_row
                                               K - 1  );                 // last_col

    arma::mat sijt_arma = sij_arma.submat(  cdindex(t),               // first_row
                                            0,                        // first_col
                                            cdindex(t+1)-1,             // last_row
                                            amountDraws - 1  );                 // last_col


    arma::mat qvt_arma = nodesRcMktShape_arma.row( t );


    if (total_dem >0) {
      dt_arma = nodesDemMktShape_arma.row( t ) ;
    }


    arma::mat dsdxi_arma = dstddelta_c( sijt_arma,
                                         weights_arma);

    arma::mat dsdtheta_arma = dstdtheta_c( sijt_arma,
                                           indices,
                                          x2t_arma,
                                          qvt_arma,
                                          dt_arma,
                                          weights_arma);

    arma::mat invertible(1,1);
    //  inverse check: inv(B,A) resets B and returns a bool set to false (exception is not thrown)

    if( arma::inv(invertible, dsdxi_arma ) ) {
      jacob_arma.submat(  cdindex(t),               // first_row
                          0,                        // first_col
                          cdindex(t+1)-1,             // last_row
                          total_par - 1  ) = - dsdxi_arma.i() * dsdtheta_arma;

    } else {
      if (printLevel >= 2) {
        Rcpp::Rcout << "Error in jacobian (market " << t
                    << ") : Singular matrix occured" << std::endl;
      }
      return 0;
    }

  } // end of t loop

  return ( wrap( jacob_arma ) );
} // end of function

