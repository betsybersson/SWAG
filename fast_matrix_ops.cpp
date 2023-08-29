// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>


// [[Rcpp::export]]
SEXP csolve(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd C = A.inverse();
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP multAB(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP multABC(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd D = A * B * C;
  
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP multABAprime(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd D = A * B * A.transpose();
  
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP multAprimeBA(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd D = A.transpose() * B * A;
  
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP multAprimeA(const Eigen::Map<Eigen::MatrixXd> A ){
  const int n(A.cols());
  Eigen::MatrixXd AtA(Eigen::MatrixXd(n, n).setZero().
                 selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
  return Rcpp::wrap(AtA);
}

// [[Rcpp::export]]
SEXP multAAprime(const Eigen::Map<Eigen::MatrixXd> A){
  const int m(A.rows());
  Eigen::MatrixXd AAt(Eigen::MatrixXd(m,m).setZero().
                        selfadjointView<Eigen::Lower>().rankUpdate(A));
  return Rcpp::wrap(AAt);
}

// [[Rcpp::export]]
SEXP ccholU (const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  Eigen::MatrixXd U(Eigen::MatrixXd(llt.matrixU()));
  return Rcpp::wrap(U);
}

// [[Rcpp::export]]
SEXP ccholL (const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  Eigen::MatrixXd U(Eigen::MatrixXd(llt.matrixL()));
  return Rcpp::wrap(U);
}

// [[Rcpp::export]]
SEXP crwish (const Eigen::Map<Eigen::MatrixXd> A,Eigen::Map<Eigen::MatrixXd> Z){
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  const int m(A.rows());
  Eigen::MatrixXd U(Eigen::MatrixXd(llt.matrixU()));
  Eigen::MatrixXd C = Z * U;
  Eigen::MatrixXd wish(Eigen::MatrixXd(m, m).setZero().
        selfadjointView<Eigen::Lower>().rankUpdate(C.adjoint()));
  return Rcpp::wrap(wish);
}    

// [[Rcpp::export]]
SEXP crinvwish (const Eigen::Map<Eigen::MatrixXd> A,Eigen::Map<Eigen::MatrixXd> Z){
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  const int m(A.rows());
  Eigen::MatrixXd U(Eigen::MatrixXd(llt.matrixU()));
  Eigen::MatrixXd C = Z * U;
  Eigen::MatrixXd wish(Eigen::MatrixXd(m, m).setZero().
                         selfadjointView<Eigen::Lower>().rankUpdate(C.adjoint()));
  Eigen::MatrixXd iwish = wish.inverse();
  
  return Rcpp::wrap(iwish);
} 

// [[Rcpp::export]]
SEXP crmatnorm (const Eigen::Map<Eigen::MatrixXd> M,Eigen::Map<Eigen::MatrixXd> U,Eigen::Map<Eigen::MatrixXd> V,Eigen::Map<Eigen::MatrixXd> Z){
  const int N(Z.rows()),P(Z.cols());
  const Eigen::LLT<Eigen::MatrixXd> Ullt(U);
  const Eigen::LLT<Eigen::MatrixXd> Vllt(V);
  Eigen::MatrixXd VU(Eigen::MatrixXd(Vllt.matrixU()));
  Eigen::MatrixXd UL(Eigen::MatrixXd(Ullt.matrixL()));
  Eigen::MatrixXd nor = M + UL * Z * VU;
  return Rcpp::wrap(nor);
}    

// [[Rcpp::export]]
SEXP crmvnorm (const Eigen::Map<Eigen::MatrixXd> M,Eigen::Map<Eigen::MatrixXd> V,Eigen::Map<Eigen::MatrixXd> Z){
  const int N(Z.rows()),P(Z.cols());
  const Eigen::LLT<Eigen::MatrixXd> Vllt(V);
  Eigen::MatrixXd VU(Eigen::MatrixXd(Vllt.matrixU()));
  Eigen::MatrixXd nor = M + Z * VU;
  return Rcpp::wrap(nor);
}    