
#include <Rcpp.h>
#include <math.h>   
using namespace Rcpp;

//compute a wilcoxon sum rank test, over the ranks given by X. summing all ranks where Y=1
double Compute_Wilcoxon_Stat(NumericVector X, IntegerVector Y){
  double stat=0.0;
  for(int i=0;i<X.length();i++){
    if(Y(i) == 1)
      stat += X(i);
  }
  return stat;
}

//function for computing the wilcoxon rank sum test, over a matrix of given permutations. Each column in Y is a different permutation of the data
// [[Rcpp::export]]
List rcpp_Wilcoxon_PermTest_Given_Permutations(NumericVector X, IntegerMatrix Y) {
  
  NumericVector stats(Y.ncol());
  for(int i=0;i<Y.ncol();i++){
    stats(i) = Compute_Wilcoxon_Stat(X,Y(_,i));
  }
  List z  = List::create( stats);
  return z ;
}


//compute a wilcoxon sum rank test, over the ranks given by X. summing all ranks where Y=1
double Compute_Spearman_Stat(NumericVector X, NumericVector Y){
  double stat=0.0;
  for(int i=0;i<X.length();i++){
      stat += X(i) * Y(i);
  }
  return stat;
}

//function for computing the wilcoxon rank sum test, over a matrix of given permutations. Each column in Y is a different permutation of the data
// [[Rcpp::export]]
List rcpp_Spearman_PermTest_Given_Permutations(NumericVector X, NumericMatrix Y) {
  NumericVector stats(Y.ncol());
  for(int i=0;i<Y.ncol();i++){
    stats(i) = Compute_Spearman_Stat(X,Y(_,i));
  }
  List z  = List::create( stats);
  return z ;
}


// two part test for counts data
double Compute_TwoPartTest(NumericVector X, IntegerVector Y, NumericVector additional_return_values){
  double total_counts_0 = 0.0;
  double total_counts_1 = 0.0;
  double counts_zeros_0 = 0.0;
  double counts_zeros_1 = 0.0;
  double nr_samples_0 = 0.0;
  double nr_samples_1 = 0.0;
  double nr_samples_0_non_zero = 0.0;
  double nr_samples_1_non_zero = 0.0;
  
  //double subsample_size = SubsampleSize(0);
  double eps = 0.01;
  
  for(int i=0 ; i<X.length(); i++){
    if(Y(i) == 1){
      if(X(i)<= eps){
        counts_zeros_1 += 1.0;
      }else{
        total_counts_1 += X(i);
      }
      nr_samples_1 = nr_samples_1 + 1.0;
    }else{
      if(X(i)<= eps){
        counts_zeros_0 += 1.0;
      }else{
        total_counts_0 += X(i);  
      }
    }
    
  }
  nr_samples_0 = ((double)(Y.length())) - nr_samples_1;
  nr_samples_1_non_zero = nr_samples_1 - counts_zeros_1;
  nr_samples_0_non_zero = nr_samples_0 - counts_zeros_0;
  
  double N_dbl_Wilcoxon = nr_samples_1_non_zero + nr_samples_0_non_zero;
    
  double Mean_Wilcoxon = nr_samples_1_non_zero * (N_dbl_Wilcoxon + 1) / 2.0;
  
  double Var_Wilcoxon = nr_samples_0_non_zero * nr_samples_1_non_zero * (N_dbl_Wilcoxon + 1) / 12.0;
  
  double Z_wilcoxon = (total_counts_1 - Mean_Wilcoxon)/std::sqrt(Var_Wilcoxon);
  
  double phat_pooled = (counts_zeros_0 + counts_zeros_1) / (nr_samples_0 + nr_samples_1);
  
  double phat_0 = counts_zeros_0 / nr_samples_0;
  
  double phat_1 = counts_zeros_1 / nr_samples_1;
  
  double stat = Z_wilcoxon * Z_wilcoxon;
  
  double Z_proportion = 0;
  
  if(counts_zeros_0 + counts_zeros_1 >0){
    Z_proportion = (phat_0 - phat_1) / std::sqrt(phat_pooled * (1.0 - phat_pooled)*(1.0/nr_samples_1 + 1.0/nr_samples_0));
    stat +=Z_proportion * Z_proportion;
  }
    
  additional_return_values[0] = Z_proportion;
  additional_return_values[1] = Z_wilcoxon;
  additional_return_values[2] = total_counts_1;
  additional_return_values[3] = counts_zeros_0;
  additional_return_values[4] = counts_zeros_1;
  additional_return_values[5] = nr_samples_0_non_zero;
  additional_return_values[6] = nr_samples_1_non_zero;
  
  return stat;
}

// function for computing the two part test, given a matrix of permutations
// [[Rcpp::export]]
List rcpp_TwoPartTest_Given_Permutations(NumericVector X, IntegerMatrix Y) {
  NumericVector stats(Y.ncol());
  NumericVector AdditionalReturnValues(7);
  double stat = Compute_TwoPartTest(X,Y,AdditionalReturnValues);
  for(int i=0;i<Y.ncol();i++){
    stats(i) = Compute_TwoPartTest(X,Y(_,i),AdditionalReturnValues);
  }
  List z  = List::create( stats);
  return z ;
}

//function for comuting the wilcoxon rank sum test
// [[Rcpp::export]]
NumericVector rcpp_Compute_Wilcoxon_Signed_Rank_Stat(NumericVector RankedDifferences, IntegerVector SecondGroupBigger){
  double stat=0.0;
  for(int i=0;i<RankedDifferences.length();i++){
    if(SecondGroupBigger(i) == 1)
      stat += RankedDifferences(i);
  }
  NumericVector ret(1);
  ret(0) = stat;
  return ret;
}


//compute a welch t-test,
double Compute_Welch_Stat(NumericVector X, IntegerVector Y){
  double _n0 = 0.0;
  double _n1 = 0.0;
  double _sum_x_0 = 0.0;
  double _sum_x_1 = 0.0;
  double _sum_x2_0 = 0.0;
  double _sum_x2_1 = 0.0;
  
  for(int i=0;i<X.length();i++){
    if(Y(i) == 1){
      _n1 = _n1+1.0;
      _sum_x_1 += X(i);
      _sum_x2_1 += X(i)*X(i);
    }else{
      _sum_x_0 += X(i);
      _sum_x2_0 += X(i)*X(i);
      _n0 = _n0+1.0;
    }
  }
  
  double _mean_0 = _sum_x_0/_n0;
  double _mean_1 = _sum_x_1/_n1;
  double _S2_0 = (_sum_x2_0 - _n0*_mean_0*_mean_0)/(_n0 - 1.0);
  double _S2_1 = (_sum_x2_1 - _n1*_mean_1*_mean_1)/(_n1 - 1.0);
  
  double stat = (_mean_0 - _mean_1)/(std::sqrt( _S2_0/_n0 + _S2_1/_n1 ));
  return stat;
}

//function for computing the wilcoxon rank sum test, over a matrix of given permutations. Each column in Y is a different permutation of the data
// [[Rcpp::export]]
List rcpp_Welch_PermTest_Given_Permutations(NumericVector X, IntegerMatrix Y) {
  
  NumericVector stats(Y.ncol());
  for(int i=0;i<Y.ncol();i++){
    stats(i) = Compute_Welch_Stat(X,Y(_,i));
  }
  List z  = List::create( stats);
  return z ;
}

