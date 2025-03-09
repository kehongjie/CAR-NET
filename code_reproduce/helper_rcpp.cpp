#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// Function to calculate the standard deviation of an Eigen vector
double calculate_sd(const VectorXd &v) {
  double mean = v.mean();
  double sq_sum = (v.array() - mean).square().sum();
  return std::sqrt(sq_sum / v.size());
}

// Helper function to calculate residuals
VectorXd cal_res(const VectorXd &y, const MatrixXd &X) {
  MatrixXd XtX = X.transpose() * X;
  MatrixXd XtX_inv = XtX.completeOrthogonalDecomposition().pseudoInverse();
  VectorXd beta = XtX_inv * X.transpose() * y;
  VectorXd res = y - X * beta;
  
  if (calculate_sd(res) == 0) {
    res = res.array() + ArrayXd::Random(res.size()) * 1e-50;
  }
  
  return res;
}

// Optimized correlation function
double cor_optimized(const VectorXd &x, const VectorXd &y) {
  int n = x.size();
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_x_sq = 0.0;
  double sum_y_sq = 0.0;
  double sum_xy = 0.0;

  for (int i = 0; i < n; ++i) {
    sum_x += x[i];
    sum_y += y[i];
    sum_x_sq += x[i] * x[i];
    sum_y_sq += y[i] * y[i];
    sum_xy += x[i] * y[i];
  }

  double mean_x = sum_x / n;
  double mean_y = sum_y / n;

  double numerator = sum_xy - n * mean_x * mean_y;
  double denominator = std::sqrt((sum_x_sq - n * mean_x * mean_x) * (sum_y_sq - n * mean_y * mean_y));
  
  if (denominator == 0) {
    stop("Denominator is zero, indicating zero variance in one of the vectors.");
  }
  
  return numerator / denominator;
}

// [[Rcpp::export]]
double partial_cor(NumericVector x, NumericVector y, NumericMatrix z) {
  Map<VectorXd> x_map(as<Map<VectorXd>>(x));
  Map<VectorXd> y_map(as<Map<VectorXd>>(y));
  Map<MatrixXd> z_map(as<Map<MatrixXd>>(z));
  
  VectorXd res1 = cal_res(x_map, z_map);
  VectorXd res2 = cal_res(y_map, z_map);
  
  // Calculate the correlation between the residuals using the optimized function
  double cor = cor_optimized(res1, res2);
  
  return cor;
}
