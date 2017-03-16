//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix  getSijMod(const NumericMatrix &expmu,  const NumericVector  &expdelta,  const IntegerVector  &cdindex){

                int nmkt = cdindex.size()-1;
                int amountNodes = expmu.ncol();
                int nobs = expmu.nrow();
                double sumMktNodes = 0;
                int nprodt;
                int startpos;

                NumericMatrix numer(nobs,amountNodes);
                 for(int i=0;i<nobs;i++){
                      for(int j=0;j<amountNodes;j++){
                       numer(i,j)=expmu(i,j)*expdelta[i];
                       }
                 }


                NumericMatrix sij(nobs,amountNodes);
                for(int i=0;i<nmkt;i++){  // iteriere über Märkte

                  nprodt = cdindex[i+1]-cdindex[i] ;
                  startpos = cdindex[i]+1;

                  for( int z=0; z<amountNodes;z++) { //iteriere über draws
                    sumMktNodes = 0;
                    for( int j=0;j<nprodt;j++) { // iteriere über Produkte
                      sumMktNodes += numer(startpos+(j-1),z);
                    }

                    for( int j=0;j<nprodt;j++) { // iteriere über Produkte
                      sij(startpos+(j-1),z)= expmu(startpos+(j-1),z)/(1+sumMktNodes);
                    }
                  }
                }


                return sij;
}




// [[Rcpp::export]]
NumericMatrix  getExpMu(const NumericMatrix &theta2Matrix, const NumericMatrix &qv,const NumericMatrix &Xrandom,const IntegerVector &cdid , const NumericMatrix &demographics){
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
List getDelta( const NumericMatrix & theta2,
               const NumericVector & deltaOld,

               const IntegerVector & cdid,
               const IntegerVector & cdindex,
               const NumericMatrix & Xrandom,
               const NumericVector & obsshare,
               const double & innerCrit,
			   const int & innerMaxit,
			   const int & printLevel,
               const NumericMatrix indices,
               const NumericMatrix nodesRcMktShape,
               const NumericMatrix nodesDemMktShape,
               const NumericVector weights ){

  // initialising objects
  double crit = 100.0;
  int counter = 0;
  int nobs = cdid.size();

  int amountRC = theta2.nrow();
  int amountNodes = nodesRcMktShape.ncol() / amountRC; // also possible with demopgraphics

  NumericMatrix sijMod( nobs , amountNodes );
  NumericMatrix expMu( nobs , amountNodes );
  NumericVector shMod( nobs );
  NumericVector expDeltaStart = exp( deltaOld );
  double shareTmp;
  bool negShares = FALSE;
  bool negWeights = min(weights) < 0 ;

  if( any( is_na( deltaOld ) ) | any( is_nan( deltaOld ) ) ) {
    for(int i = 0; i < nobs; ++i) {
      expDeltaStart[i] = 1;
    }
  }

  NumericVector expDeltaNew( nobs );


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

  while( crit > innerCrit ){
    counter++;

  	if( counter == innerMaxit){
  		  List ret;
  		  ret["delta"] = NAN;
  		  ret["sijMod"] = NAN;
  		  ret["counter"] = counter;
  		  ret["negShares"] = negShares;
  		  return( ret );
  	}

    sijMod = getSijMod( expMu,
                        expDeltaStart,
                        cdindex );

    for( int i=0; i<nobs;i++){
      shareTmp=0;
      for( int j=0; j<amountNodes;j++){
        shareTmp += sijMod(i,j)*weights[j];
      }
      shMod[i] =  shareTmp;
    }

    expDeltaNew = obsshare / shMod;

    //Handling of negative shares for quadrature
    if( negWeights ) {
    if( min(shMod) < 0 ){ //if( any( shMod < 0 ).is_true() ) {
      if( !negShares ){
        Rcpp::Rcout << "Integration rule produced negative shares --> reset delta to 0" ;
        negShares = TRUE ;
      }
      Rcpp::Rcout << "." ;

      for(int negInd = 0; negInd < nobs; negInd++) {
        if( shMod[negInd] < 0 ){
          expDeltaNew[negInd] = 1;
        }
      }
    }
    } //end negWeights

    crit = max(abs(expDeltaNew - expDeltaStart));
    expDeltaStart = clone(expDeltaNew);

	if( printLevel == 4 ){
      Rcpp::Rcout << "\t crit: "  << crit << std::endl ;
    }


  }

  if( negShares ) Rcpp::Rcout << std::endl ;




  List ret;
  ret["delta"] = log(expDeltaNew);
  ret["sijMod"] = sijMod;
  ret["counter"] = counter;
  ret["negShares"] = negShares;
  return( ret );
}






