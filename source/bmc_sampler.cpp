#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::List ProdRes(const List& Y, const List& X, const List& beta_ij,
                   IntegerMatrix Start, IntegerMatrix End, IntegerVector m_j) {
  int J = Y.size();
  int p;
  List Res(J);
  NumericVector Resj;
  double res;
  arma::vec Yvec;
  arma::mat Xmat;
  arma::mat Bmat;
  for (int j = 0; j < J; j++) {
    Resj = NumericVector(max(End(_,j)));
    Yvec = as<arma::vec>(Y[j]);
    Xmat = as<arma::mat>(X[j]);
    p = Xmat.n_cols;
    Bmat = as<arma::mat>(beta_ij[j]);
   for (int i = 0; i < m_j[j]; i++) {
     for (int k = Start(i,j)-1; k < End(i,j); k++) {
       res = as_scalar(Yvec[k] - Xmat.submat(k,0,k,p-1)*Bmat.submat(0,i,p-1,i)); // (startrow: k, startcolumn: 0, endrow: k, endcolumn: p-1)
       Resj[k] = res;
     }
   }
   Res[j] = Resj;
  }
  return Res;
}

// [[Rcpp::export]]
Rcpp::List lapplycpp(const List& C, NumericVector thetah) {
  int J = C.size();
  List Thetah(J);
  IntegerVector idx;
  IntegerVector idx2;
  for (int j = 0; j < J; j++) {
    idx = C[j];
    IntegerVector ones (idx.length(), 1);
    idx2 = idx - ones;
    Thetah[j] = thetah[idx2];
  }
  return Thetah;
}

// [[Rcpp::export]]
Rcpp::List DvidD(const List& Y, const List& orgX, const arma::mat& d_ij, const List& idx_j,
                   IntegerMatrix Start, IntegerMatrix End, IntegerVector m_j) {
  int J = Y.size();
  List Res(J);
  NumericVector Resj;
  double res;
  arma::vec Yvec;
  arma::vec Xvec;
  arma::vec idxj;
  for (int j = 0; j < J; j++) {
    Resj = NumericVector(max(End(_,j)));
    Yvec = as<arma::vec>(Y[j]);
    Xvec = as<arma::vec>(orgX[j]);
    idxj = as<arma::vec>(idx_j[j]);
    for (int i = 0; i < m_j[j]; i++) {
      for (int k = Start(i,j)-1; k < End(i,j); k++) {
        res = Yvec[k]/exp(Xvec[k]*d_ij(idxj[i]-1,j)/2);
        Resj[k] = res;
      }
    }
    Res[j] = Resj;
  }
  return Res;
}

// [[Rcpp::export]]
Rcpp::List DvidD_mat(const List& X, const List& orgX, const arma::mat& d_ij, const List& idx_j,
                 IntegerMatrix Start, IntegerMatrix End, IntegerVector m_j) {
  int J = X.size();
  int p; 
  List Res(J);
  arma::mat Resj;
  arma::rowvec res;
  arma::mat Ymat;
  arma::vec Xvec;
  arma::vec idxj;
  arma::rowvec dum;
  for (int j = 0; j < J; j++) {
    Ymat = as<arma::mat>(X[j]);
    p = Ymat.n_cols;
    Resj = arma::mat(max(End(_,j)), p);
    dum = arma::rowvec(p);
    Xvec = as<arma::vec>(orgX[j]);
    idxj = as<arma::vec>(idx_j[j]);
    for (int i = 0; i < m_j[j]; i++) {
      for (int k = Start(i,j)-1; k < End(i,j); k++) {
        res = Ymat.submat(k,0,k,p-1) % dum.fill(1/exp(Xvec[k]*d_ij(idxj[i]-1,j)/2));
        Resj.submat(k,0,k,p-1) = res;
      }
    }
    Res[j] = Resj;
  }
  return Res;
}

// [[Rcpp::export]]
arma::mat loglikecpp(const List& e, const arma::mat& d_ij, const List& orgX, const arma::vec& sigj_sq, 
                     const List& idx_j, IntegerMatrix Start, IntegerMatrix End, IntegerVector m_j) {
  int J = d_ij.n_cols;
  int m = d_ij.n_rows; 
  arma::mat loglike (m, J);
  double D;
  double res;
  arma::vec Xvec;
  arma::vec evec;
  arma::vec idxj;
  
  for (int j = 0; j < J; j++) {
    Xvec = as<arma::vec>(orgX[j]);
    evec = as<arma::vec>(e[j]);
    idxj = as<arma::vec>(idx_j[j]);
    for (int i = 0; i < m_j[j]; i++) {
      res = 0;
      for (int k = Start(i,j)-1; k < End(i,j); k++) {
        D = exp(Xvec[k]*d_ij(idxj[i]-1,j)/2);
        res = res + R::dnorm(evec[k], 0, D*sqrt(sigj_sq[j]), TRUE);
      }
      loglike(idxj[i]-1, j) = res;
    }
  }
  return loglike;
}

// [[Rcpp::export]]
Rcpp::List gbcpp_post(const List& Y, const List& X,
                      const arma::mat& pi_ij, const arma::mat& invSigj, const arma::mat& Sigj, const arma::vec& sigj_sq,
                      const List& idx_j, IntegerVector m_j, IntegerMatrix Start, IntegerMatrix End) {
  int J = pi_ij.n_cols;
  int m = pi_ij.n_rows;
  int p;
  arma::vec Yvec;
  arma::mat Xmat;
  arma::vec idxj;
  arma::vec iSigvec;
  arma::vec Sigvec;
  arma::mat gamma_ij;
  List beta_ij(J);
  arma::mat betajmat;
  arma::colvec dum;
  
  gamma_ij = arma::mat(m,J).fill(NA_REAL);
  for (int j = 0; j < J; j++) {
    Yvec = as<arma::vec>(Y[j]);
    Xmat = as<arma::mat>(X[j]);
    p = Xmat.n_cols;
    idxj = as<arma::vec>(idx_j[j]);
    iSigvec = invSigj.submat(0,j,p*p-1,j);
    Sigvec = Sigj.submat(0,j,p*p-1,j);
    betajmat = arma::mat(p,m_j[j]);
    dum = arma::colvec(p);
    
    for (int i = 0; i < m_j[j]; i++) {
      arma::mat invV = pinv(reshape(iSigvec,p,p) + Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1).t()*
        (arma::eye(End(i,j)-Start(i,j)+1,End(i,j)-Start(i,j)+1)/sigj_sq[j])*
        Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1));
      arma::vec mu0 = (Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1).t() *
        Yvec.subvec(Start(i,j)-1,End(i,j)-1)) % dum.fill(1/sigj_sq[j]);
      
      double loggam1 = -0.5*log(det(reshape(Sigvec,p,p)*
                                Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1).t()*
                                (arma::eye(End(i,j)-Start(i,j)+1,End(i,j)-Start(i,j)+1)/sigj_sq[j])*
                                Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1) + arma::eye(p,p))) +
                                0.5*as_scalar(mu0.t()*invV*mu0) + log(pi_ij(idxj[i]-1,j));
      double loggam0 = log(1-pi_ij(idxj[i]-1,j));
      
      NumericVector gamma = rbinom(1,1,1/(1+exp(-loggam1+loggam0)));
      gamma_ij(idxj[i]-1,j) = gamma[0];
      
      if (gamma[0] == 1) {
        arma::vec mu = invV*mu0;
        arma::mat Z = arma::randn(p,1);
        betajmat.submat(0,i,p-1,i) = mu + arma::chol(arma::symmatu(invV),"lower")*Z;
      } else {
        betajmat.submat(0,i,p-1,i) = arma::zeros(p);
      }
    }
    beta_ij[j] = betajmat;
  }
  
  List out;
  out["beta_ij"] = beta_ij;
  out["gamma_ij"] = gamma_ij;
  return out;
}

// temporary
// [[Rcpp::export]]
Rcpp::List bcpp_post(const List& Y, const List& X,
                     const arma::mat& gamma_ij, const arma::mat& invSigj, const arma::mat& Sigj, const arma::vec& sigj_sq,
                     const List& idx_j, IntegerVector m_j, IntegerMatrix Start, IntegerMatrix End) {
  int J = gamma_ij.n_cols;
  int p;
  arma::vec Yvec;
  arma::mat Xmat;
  arma::vec idxj;
  arma::vec iSigvec;
  arma::vec Sigvec;
  List beta_ij(J);
  arma::mat betajmat;
  arma::colvec dum;
  
  for (int j = 0; j < J; j++) {
    Yvec = as<arma::vec>(Y[j]);
    Xmat = as<arma::mat>(X[j]);
    p = Xmat.n_cols;
    idxj = as<arma::vec>(idx_j[j]);
    iSigvec = invSigj.submat(0,j,p*p-1,j);
    Sigvec = Sigj.submat(0,j,p*p-1,j);
    betajmat = arma::mat(p,m_j[j]);
    dum = arma::colvec(p);
    
    for (int i = 0; i < m_j[j]; i++) {
      arma::mat invV = pinv(reshape(iSigvec,p,p) + Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1).t()*
        (arma::eye(End(i,j)-Start(i,j)+1,End(i,j)-Start(i,j)+1)/sigj_sq[j])*
        Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1));
      arma::vec mu0 = (Xmat.submat(Start(i,j)-1,0,End(i,j)-1,p-1).t() *
        Yvec.subvec(Start(i,j)-1,End(i,j)-1)) % dum.fill(1/sigj_sq[j]);
      
      if (gamma_ij(idxj[i]-1,j) == 1) {
        arma::vec mu = invV*mu0;
        arma::mat Z = arma::randn(p,1);
        betajmat.submat(0,i,p-1,i) = mu + arma::chol(arma::symmatu(invV),"lower")*Z;
      } else {
        betajmat.submat(0,i,p-1,i) = arma::zeros(p);
      }
    }
    beta_ij[j] = betajmat;
  }
  
  return beta_ij;
}
