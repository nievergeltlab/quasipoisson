#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
mat matprod (mat amat , mat bmat) {
 mat pro = amat * bmat;
 return(pro);
} //matrix product function

// [[Rcpp::export]]
mat mattridotprod (mat amat , mat bmat, mat cmat) {
 mat pro = amat % bmat % cmat;
 return(pro);
} //matrix hadamard product function

// [[Rcpp::export]]
mat mattridotprod2 (vec amat , mat bmat, mat cmat) {
 mat pro = ones(bmat.n_rows,bmat.n_cols);
 
  for (int i =0; i < bmat.n_cols; i++){
  pro.col(i) = bmat.col(i) % amat % cmat.col(i) ;
  // cout << "The value of i : " << i << "\n";
  }
  
 return(pro);
} //matrix hadamard product function


// [[Rcpp::export]]
mat lambdakfactorial_alt3 (vec L, vec eL, int maxy) {
 mat Lk = ones(L.n_elem,maxy+1); 
 for (int i =1; i <= maxy; i++){
  vec ii = vec(L.n_elem);
  ii.fill(i);
  Lk.col(i) = (L%Lk.col(i-1))%(1/ii);
  //cout << "The value of i : " << i << "\n";
 }
  for (int i =0; i <= maxy; i++){
  Lk.col(i) = eL % Lk.col(i);
  //cout << "The value of i : " << i << "\n";
 }
   //cout << "The rows of lK : " << lK.n_rows << "\n";
   //cout << "The cols of lK : " << lK.n_cols << "\n";
   
 return (Lk);
}

// [[Rcpp::export]]
mat obs_min_expected (vec L, int maxy) {
 mat Lk = ones(L.n_elem,maxy+1); 
 for (int i =0; i <= maxy; i++){
  vec ii = vec(L.n_elem);
  ii.fill(i);
  Lk.col(i) = L - ii;
  //cout << "The value of i : " << i << "\n";
 }

 return (Lk);
}




