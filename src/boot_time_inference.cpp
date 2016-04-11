// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include "RcppArmadillo.h"
#include <math.h>
using namespace Rcpp;

NumericMatrix rbindCpp(NumericMatrix a, NumericMatrix b){
  if( a.ncol() != b.ncol() )
    stop("rbind failed due to mismatch in number of columns");

  int arow = a.nrow();
  int brow = b.nrow();
  NumericMatrix out(arow + brow, a.ncol());
  for(int j = 0; j < arow+brow; j++){
    if(j < arow){
      out(j,_) = a(j,_);
    } else {
      out(j,_) = b(j-arow,_);
    }
  }
  return(out);
}

NumericMatrix rbindCpp(NumericVector a, NumericVector b){
  if( a.length() != b.length() )
    stop("rbind failed due to mismatch in length");

  NumericMatrix out(2, a.length());
  out(0,_) = a;
  out(1,_) = b;
  return(out);
}


NumericMatrix cbindCpp(NumericVector a, NumericVector b) {
  int acoln = 1;
  int bcoln = 1;

  NumericMatrix out = no_init_matrix(a.length(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a;
    } else {
      out(_, j) = b;
    }
  }
  return out;
}


NumericMatrix cbindCpp(NumericMatrix b, NumericVector a) {
  int acoln = 1;
  int bcoln = b.ncol();

  NumericMatrix out(a.length(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a;
    } else {
      out(_, j) = b(_,j-acoln);
    }
  }
  return out;
}


NumericMatrix cbindCpp(NumericVector a, NumericMatrix b) {
  int acoln = 1;
  int bcoln = b.ncol();

  NumericMatrix out(a.length(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a;
    } else {
      out(_, j) = b(_,j-acoln);
    }
  }
  return out;
}

NumericMatrix cbindCpp(NumericMatrix a, NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();

  NumericMatrix out(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_,j);
    } else {
      out(_, j) = b(_,j-acoln);
    }
  }
  return out;
}

IntegerVector C(IntegerVector a, IntegerVector b){
  int lena = a.length();
  int lenb = b.length();
  IntegerVector out(lena + lenb);
  for(int i = 0; i < lena+lenb; i++){
    if(i < lena){
      out(i) = a(i);
    } else {
      out(i) = b(i-lena);
    }
  }
  return(out);
}

List fastLm(NumericVector yr, NumericMatrix Xr) {

  int k = Xr.ncol();
  int n = Xr.nrow();

  arma::mat X(Xr.begin(), n, k, false); // reuses memory and avoids extra copy
  arma::colvec y(yr.begin(), yr.size(), false);

  arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
  arma::colvec resid = y - X*coef;            // residuals

  double sig2 = arma::as_scalar( arma::trans(resid)*resid/(n-k) );
  // std.error of estimate
  arma::colvec stderrest = arma::sqrt(
    sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );

  return List::create(Named("coef")         = coef,
                      Named("stderr")       = stderrest,
                      Named("resid")        = resid
  );
}


IntegerVector Csample( IntegerVector x,
                       int size = 1,
                       bool replace = false,
                       NumericVector prob = NumericVector::create()) {
  IntegerVector sam = RcppArmadillo::sample(x, size, replace, prob);
  return sam;
}

IntegerVector sbSequenceCpp(int T, int b_av, int length){
  IntegerVector index_sequence(2*T);

  for(int j = 0; j < 2*T; j++){
    int k;
    if(j >= T){
      k = j-T;
    } else {
      k = j;
    }
    index_sequence[j] = k+1;
  }

  IntegerVector sequence = rep(IntegerVector::get_na(), length + T);
  IntegerVector temp = seq(1,T); // to sample from

  int current = -1;
  while(current < (length-1)){
    int start = as<int>(Csample(temp));
    // int start = as<int>(sample_(temp, 1, false));
    NumericVector b_num = rgeom(1, 1.0/b_av) + 1.0;
    int b_int = as<int>(b_num);

    IntegerVector _to_set = seq(current+1,current+b_int);
    IntegerVector _idx = seq(start-1, start-1 + b_int - 1);
    // replace with NA (when index out of bound error)
    for(int i = 0; i < _idx.length(); i++){
      if(_idx[i] >= index_sequence.length() || _to_set[i] >= sequence.length()){
        // do nothing
      } else {
        sequence[_to_set[i]] = index_sequence[_idx[i]];
      }
    }

    current = current + b_int;
  }

  return sequence[Range(0,length-1)];
}

double sharpeRatioDiff(NumericMatrix ret) {

  NumericVector ret1 = ret(_,0);
  NumericVector ret2 = ret(_,1);
  double mu1_hat = mean(ret1);
  double mu2_hat = mean(ret2);
  double sig1_hat = sd(ret1);
  double sig2_hat = sd(ret2);
  double SR1_hat = mu1_hat/sig1_hat;
  double SR2_hat = mu2_hat/sig2_hat;
  double diff = SR1_hat - SR2_hat;

  return diff;
}

NumericVector pow2(NumericVector v){
  int n = v.length();
  NumericVector out(n);
  for(int j = 0; j < n; j++)
    out(j) = v(j) * v(j);

  return( out );
}

NumericVector pow15(NumericVector v){
  int n = v.length();
  NumericVector out(n);
  for(int j = 0; j < n; j++)
    out(j) = pow(v(j), 1.5);

  return( out );
}

NumericMatrix computeVhat(NumericMatrix ret){
  NumericVector ret1 = ret(_,0);
  NumericVector ret2 = ret(_,1);

  NumericMatrix V_hat1 = cbindCpp(ret1 - mean(ret1), ret2 - mean(ret2));
  NumericMatrix V_hat2 = cbindCpp(pow2(ret1) - mean(pow2(ret1)), pow2(ret2) - mean(pow2(ret2)));

  return(cbindCpp(V_hat1,V_hat2));
}

List baseSVD(NumericMatrix Xnm) {
  arma::mat X = as<arma::mat>(wrap(Xnm));
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");
  return List::create(_["d"] = S,
                      _["u"] = U,
                      _["v"] = V);
}

inline static double sqrt_double( double x ){ return ::sqrt( x ); }
NumericVector sqrtCpp(const NumericVector & orig) {
  NumericVector vec(orig.size());     // create a target vector of the same size
  std::transform(orig.begin(), orig.end(), vec.begin(), sqrt_double);

  return(vec);
}

double computeAlphaHat(NumericMatrix V_hat){
  int p = V_hat.ncol();
  double numerator = 0.0;
  double denominator = 0.0;

  Environment env1("package:stats");
  Function ar_ols = env1["ar.ols"];

  for(int i = 0; i < p; i++){
    List fit = ar_ols(V_hat(_,i), Named("aic",false), Named("order.max",1));
    double rho_hat = as<double>(fit["ar"]);
    double sig_hat = as<double>(sqrtCpp(fit["var.pred"]));
    numerator = numerator + 4 * pow(rho_hat,2) * pow(sig_hat,4) / pow(1-rho_hat,8);
    denominator = denominator + pow(sig_hat,4) / pow(1-rho_hat,4);
  }
  return( numerator/denominator );
}

NumericMatrix computeGammaHat(NumericMatrix V_hat, int j){
  int T = V_hat.nrow();
  int p = V_hat.ncol();
  NumericMatrix temp(p,p);
  arma::mat Gamma_hat = as<arma::mat>(wrap(temp));
  arma::mat arma_V_hat = as<arma::mat>(wrap(V_hat));

  if(j >= T)
    stop("j must be smaller than the row dimension!");

  for(int i = j; i < T; i++){
    Gamma_hat = Gamma_hat + arma_V_hat.row(i).t() * arma_V_hat.row(i-j);
  }
  arma::mat result = Gamma_hat/T;
  return(wrap(result));
}

double kernelParzen(double x){
  double result;
  double abs_x = fabs(x);

  if(abs_x <= 0.5){
    result = 1.0 - 6.0*pow(x,2.0) + 6.0*pow(abs_x,3.0);
  } else if(abs_x <= 1){
    result = 2.0 * pow(1-abs_x,3);
  } else {
    result = 0.0;
  }
  return(result);
}

NumericMatrix computePsiHat(NumericMatrix V_hat){

  int T = V_hat.nrow();
  double alpha_hat = computeAlphaHat(V_hat);
  double S_star = 2.6614 * pow(alpha_hat * T, 0.2);
  NumericMatrix Psi_hat = computeGammaHat(V_hat, 0);
  arma::mat arma_Psi_hat = as<arma::mat>(Psi_hat);

  int j = 1;
  while( j < S_star ){
    arma::mat Gamma_hat = as<arma::mat>(wrap(computeGammaHat(V_hat,j)));
    arma_Psi_hat = arma_Psi_hat + kernelParzen(j/S_star) * (Gamma_hat + Gamma_hat.t());
    j = j+1;
  }

  arma_Psi_hat =  arma_Psi_hat * ((double)T/(T-4));
  return(wrap(arma_Psi_hat));
}

double computeSeParzenPw(NumericMatrix ret){
  NumericVector ret1 = ret(_,0);
  NumericVector ret2 = ret(_,1);
  double mu1_hat = mean(ret1);
  double mu2_hat = mean(ret2);
  double gamma1_hat = mean(pow2(ret1));
  double gamma2_hat = mean(pow2(ret2));

  NumericVector gradient(4);
  gradient(0) = gamma1_hat / ( pow(gamma1_hat - pow(mu1_hat,2), 1.5) );
  gradient(1) = (-1) * gamma2_hat / (pow(gamma2_hat - pow(mu2_hat,2), 1.5));
  gradient(2) = (-0.5) * mu1_hat / (pow(gamma1_hat - pow(mu1_hat,2), 1.5));
  gradient(3) = 0.5 * mu2_hat / (pow(gamma2_hat - pow(mu2_hat,2),1.5));

  int T = ret1.length();
  NumericMatrix V_hat = computeVhat(ret);
  NumericMatrix A_ls(4,4);
  NumericMatrix V_star(T-1,4);

  NumericVector reg1 = as<NumericVector>(wrap(V_hat(Range(0,T-2),Range(0,0))));
  NumericVector reg2 = as<NumericVector>(wrap(V_hat(Range(0,T-2),Range(1,1))));
  NumericVector reg3 = as<NumericVector>(wrap(V_hat(Range(0,T-2),Range(2,2))));
  NumericVector reg4 = as<NumericVector>(wrap(V_hat(Range(0,T-2),Range(3,3))));

  for(int j = 0; j < 4; j++){
    NumericVector y = as<NumericVector>(wrap(V_hat(Range(1,T-1),Range(j,j))));
    NumericMatrix X = cbindCpp(reg1, cbindCpp(reg2, cbindCpp(reg3, reg4)));
    List fit = fastLm(y, X);
    A_ls(j,_) = as<NumericVector>(wrap(fit["coef"]));
    V_star(_,j) = as<NumericVector>(wrap(fit["resid"]));
  }

  List svd_A = baseSVD(A_ls);
  NumericVector d = as<NumericVector>(wrap(svd_A["d"]));
  NumericVector d_adj = Rcpp::clone(d);
  for(int i = 0; i < 4; i++){
    if(d(i) > 0.97){
      d_adj(i) = 0.97;
    } else if(d(i) < -0.97) {
      d_adj(i) = -0.97;
    }
  }
  arma::mat v = svd_A["v"];
  arma::mat u = svd_A["u"];
  arma::mat A_hat = u * arma::diagmat(as<arma::vec>(d_adj)) * v.t();
  arma::mat D = arma::inv(arma::eye(4,4) - A_hat);

  NumericMatrix reg_mat = rbindCpp(rbindCpp(reg1,reg2),rbindCpp(reg3,reg4));

  for(int j = 0; j < 4; j++){
    arma::mat dotProduct = A_hat.row(j) * as<arma::mat>(reg_mat);
    V_star(_,j) = as<NumericVector>(wrap(V_hat(Range(1,T-1),Range(j,j)))) -
      as<NumericVector>(wrap(dotProduct));
  }

  NumericMatrix Psi_hat = computePsiHat(V_star);
  arma::mat arma_Psi_hat = D * as<arma::mat>(wrap(Psi_hat)) * D.t();
  arma::vec arma_gradient = as<arma::vec>(wrap(gradient));
  double se = arma::as_scalar(arma::sqrt(arma_gradient.t() * arma_Psi_hat * (arma_gradient/T)));
  return(se);
}

IntegerVector cbbSequence(int T, int b){
  int l = T/b;
  IntegerVector index_sequence = C(seq(1,T),seq(1,b));

  IntegerVector sequence(T);
  IntegerVector start_points = Csample(seq(1,T),l,true);

  for(int j = 1; j <= l; j++){
    int start = as<int>(wrap(start_points(j-1)));
    IntegerVector to_write = Range((j-1)*b + 1,j*b)-1;
    IntegerVector to_read = Range(start-1,start+b-2);
    sequence[to_write] = index_sequence[to_read];
  }

  return(sequence[sequence > 0]);
}

template <class T>
inline double do_mean( T& x ) {
  return mean(x);
}

NumericVector col_means( NumericMatrix& X ) {

  int nCols = X.ncol();
  NumericVector out = no_init(nCols);

  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = do_mean( tmp );
  }

  return out;
}

//' calculate bootstrap inference
//'
//' @param ret refers to a returns matrix (columns = 2)
//' @param b blocksize to use. Choose the optimal blocksize by
//' \code{\link{blockSizeCalibrate}}
//' @param Delta_null refers to null hypothesis, where Delta means the difference
//' in sharpe ratio
//' @param M number of bootstrap resamples (see Eq. 9 in vignette)
//'
//' @return a list with estimated difference in sharpe ratio and p-value of the
//' test with \code{Delat_null} as null hypothesis
//'
//' @examples
//' \dontrun{
//' DATA <- bootTimeInference:::ret.hedge
//' opt <- blockSizeCalibrate(DATA, b_vec = c(1,2,4,6,8,10), K = 5000)
//' bootTimeInference(DATA, b = opt$b.optimal, M = 4999)
//' }
//' DATA <- bootTimeInference:::ret.hedge
//' # small example, please choose an appropiate K (see vignette)
//' opt <- blockSizeCalibrate(DATA, b_vec = c(1,2,4,6), K = 10)
//' bootTimeInference(DATA, b = opt$b.optimal, M = 4999)
//' @export
// [[Rcpp::export]]
List bootTimeInference(NumericMatrix ret,
                       int b,
                       int M,
                       double Delta_null = 0.0){

  int T = ret.nrow();
  int l = floor(T/b);
  double Delta_hat = sharpeRatioDiff(ret);
  double d = fabs(Delta_hat - Delta_null) / computeSeParzenPw(ret);
  double p_value = 1.0;

  // NumericMatrix ret_star;
  NumericMatrix ret_star(T - (T % b),ret.ncol());
  NumericVector gradient(4);
  NumericMatrix temp(4,4);
  for(int m = 0; m < M; m++){
    IntegerVector to_extract = cbbSequence(T,b);
    for(int i = 0; i < to_extract.length(); i++){
      ret_star(i,_) = ret(to_extract[i]-1,_);
    }
    double Delta_hat_star = sharpeRatioDiff(ret_star);
    NumericVector ret1_star = ret_star(_,0);
    NumericVector ret2_star = ret_star(_,1);
    double mu1_hat_star = mean(ret1_star);
    double mu2_hat_star = mean(ret2_star);
    double gamma1_hat_star = mean(pow2(ret1_star));
    double gamma2_hat_star = mean(pow2(ret2_star));

    gradient(0) = gamma1_hat_star/pow(gamma1_hat_star - pow(mu1_hat_star,2),1.5);
    gradient(1) = (-1)*gamma2_hat_star/pow(gamma2_hat_star - pow(mu2_hat_star,2),1.5);
    gradient(2) = (-0.5) * mu1_hat_star / pow(gamma1_hat_star - pow(mu1_hat_star,2),1.5);
    gradient(3) = 0.5 * mu2_hat_star / pow(gamma2_hat_star - pow(mu2_hat_star,2),1.5);

    NumericMatrix y_star = cbindCpp(
      cbindCpp(ret1_star - mu1_hat_star, ret2_star - mu2_hat_star),
      cbindCpp(pow2(ret1_star) - gamma1_hat_star, pow2(ret2_star) - gamma2_hat_star)
    );

    arma::mat Psi_hat_star = as<arma::mat>(wrap(temp));
    for(int j = 1; j <= l; j++){
      NumericMatrix sub = y_star(Range((j-1)*b,j*b-1),Range(0,y_star.ncol()-1));
      NumericVector zeta_star = pow(b,0.5) * col_means(sub);
      arma::vec arma_zeta_star = as<arma::vec>(wrap(zeta_star));
      Psi_hat_star = Psi_hat_star + arma_zeta_star * arma_zeta_star.t();
    }

    Psi_hat_star = Psi_hat_star/(double)l;
    arma::vec arma_gradient = as<arma::vec>(wrap(gradient));
    double se_star = arma::as_scalar(arma::sqrt(arma_gradient.t() * Psi_hat_star * arma_gradient/T));
    double d_star = fabs(Delta_hat_star - Delta_hat) / se_star;
    if(d_star >= d)
      p_value = p_value + 1.0;
  }

  p_value = p_value/(M + 1);
  return(List::create(
      _["Difference"] = Delta_hat,
      _["p.Value"] = p_value
  ));
}

//' calculate the optimal block size
//'
//' by Algorithm 3.1 in vignette on page 5
//'
//' @param ret refers to a returns matrix (columns = 2)
//' @param b_vec see vignette
//' @param alpha see vignette
//' @param M see vignette
//' @param K see vignette
//' @param b_av see vignette
//' @param T_start see vignette
//'
//' @return a list with probabilities and optimal block size
//'
//' @export
// [[Rcpp::export]]
List blockSizeCalibrate(NumericMatrix ret,
                        IntegerVector b_vec = IntegerVector::create(1, 3, 6, 10),
                        double alpha = 0.05,
                        int M = 199,
                        int K = 1000,
                        int b_av = 5,
                        int T_start = 50) {

  int b_len = b_vec.length();
  NumericVector emp_reject_probs = rep(0.0, b_len);
  double Delta_hat = sharpeRatioDiff(ret);
  NumericVector ret1 = ret(_,0);
  NumericVector ret2 = ret(_,1);
  int T = ret1.length();
  NumericMatrix Var_data(T_start + T, 2);
  Var_data(0,_) = ret(0,_);

  IntegerVector range1 = seq(1, T-1);
  IntegerVector range2 = seq(0, T-2);
  NumericVector intercept = rep(1.0, ret1.length());
  List fit1 = fastLm(ret1[range1], cbindCpp(intercept[range2], cbindCpp(ret1[range2],ret2[range2])));
  List fit2 = fastLm(ret2[range1], cbindCpp(intercept[range2], cbindCpp(ret1[range2],ret2[range2])));

  NumericVector coef1 = as<NumericVector>(wrap(fit1["coef"]));
  NumericVector coef2 = as<NumericVector>(wrap(fit2["coef"]));
  NumericMatrix resid_mat = cbindCpp(as<NumericVector>(fit1["resid"]),
                                     as<NumericVector>(fit2["resid"]));

  for(int k = 0; k < K; k++){
    // create resid_mat_star
    NumericMatrix resid_mat_star(T_start + T, resid_mat.cols());
    // fill with NA by default
    int xsize = resid_mat_star.nrow() * resid_mat_star.ncol();
    for (int i = 0; i < xsize; i++) {
      resid_mat_star[i] = NumericMatrix::get_na();
    }
    // fill first row with 0
    for(int c = 0; c < resid_mat_star.ncol(); c++)
      resid_mat_star(0,c) = 0.0;
    IntegerVector index = sbSequenceCpp(T-1, b_av, T_start + T - 1) - 1;
    for(int j = 0; j < index.length(); j++){
      // handling NA values
      if( IntegerVector::is_na(index(j)) ){
        // do nothing
      } else if( (j+1) > resid_mat_star.nrow() ){
        // do nothing
      } else if( index[j] > resid_mat.nrow() ){
        // do nothing
      } else {
        resid_mat_star(j+1,_) = resid_mat(index[j],_);
      }
    }

    for(int t = 1; t < (T_start + T); t++){
      Var_data(t,0) = coef1[0] + coef1[1]*Var_data(t-1,0) +
        coef1[2]*Var_data(t-1,1) + resid_mat_star(t,0);
      Var_data(t,1) = coef2[0] + coef2[1]*Var_data(t-1,0) +
        coef2[2]*Var_data(t-1,1) + resid_mat_star(t,1);
    }

    NumericMatrix Var_data_trunc = Var_data(Range(T_start,T_start+T-1), Range(0,Var_data.ncol()-1));
    for(int j = 0; j < b_len; j++){
      List bTI = bootTimeInference(Var_data_trunc, b_vec[j], M, Delta_hat);
      double p_value = as<double>(wrap(bTI["p.Value"]));
      if(p_value <= alpha)
        emp_reject_probs[j] = emp_reject_probs[j] + 1;
    }
  }
  emp_reject_probs = emp_reject_probs/(double)K;
  Environment env = Environment::base_namespace();
  Function order = env["order"];
  IntegerVector b_order = order(abs(emp_reject_probs - alpha));
  int b_opt = as<int>(wrap(b_vec(b_order[0]-1)));
  NumericMatrix b_vec_with_probs = rbindCpp(as<NumericVector>(wrap(b_vec)), emp_reject_probs);

  return(List::create(
      _["Empirical.Rejection.Probs"] = b_vec_with_probs,
      _["b.optimal"] = b_opt
  ));
}
