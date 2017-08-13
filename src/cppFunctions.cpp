//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;

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
  int K = theta2Matrix.nrow();
  int totalDem = theta2Matrix.ncol()-1; //1 is RC column
  int nobs = Xrandom.nrow();
  int amountNodes = qv.ncol()/K;
  int startpointRC ;
  int startpointDem ;
  int mktInd;
  double demPart = 0;
  NumericMatrix xvPart(nobs,amountNodes);
  NumericMatrix expmu(nobs,amountNodes);

  for( int i=0; i<K;i++){ // iteriere ueber RC
    startpointRC= i*amountNodes;

    for( int z=0; z<nobs; z++){
      mktInd= cdid[z]-1;
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


//[[Rcpp::export]]
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
  double current_max;
  double current_dist;

  // input check delta:
  if( any( is_na( deltaOld ) ) | any( is_nan( deltaOld ) ) ) {
    for(int i = 0; i < nobs; ++i) {
      expDeltaStart[i] = 1;
    }
  }


  // calc. the individual part of utility
  expMu = getExpMu( theta2,
                    nodesRcMktShape,
                    Xrandom,
                    cdid,
                    nodesDemMktShape) ;


  // start of the contraction mapping
  if( printLevel == 4 ){
    Rcpp::Rcout << "----------------------" << std::endl ;
  }

  while( (dist > innerCrit) & (counter < innerMaxit) ){
    counter++;

   for( int i=0; i<nmkt; i++){

    if( market_convergence[i] == 1 ){
      continue ;
      }

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
    //end negWeights

    // update delta(II):
      // init maximum value with first difference in a market:
      current_max = std::abs( expDeltaNew[startpos-1] - expDeltaStart[startpos-1] );

      // search max distance between old and new delta
      for( int j=0; j<nprodt; j++){
        current_dist = std::abs(expDeltaNew[startpos+(j-1)]  - expDeltaStart[startpos+(j-1)]) ;
        if( current_dist > current_max ){
          current_max = current_dist;
        }
      }

      if( current_max < innerCrit){
        market_convergence[i] = 1; // 1 indicates that market converged
      }

      dist_m[i] = current_max ;


     // avoids pointers and creates a deep copy:
    for( int j=0; j<nprodt; j++){
      expDeltaStart[startpos+(j-1)] = expDeltaNew[startpos+(j-1)];
    }

   } // end market loop
   dist = Rcpp::max(dist_m);

    if( printLevel == 4 ){
      Rcpp::Rcout << "\t dist: "  << dist << std::endl ;
    }


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



