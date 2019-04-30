
#include <Rcpp.h>
#include <math.h>   
using namespace Rcpp;

//RNGScope scope; // ensure RNG gets set/reset


double Compute_Wilcoxon_Stat(NumericVector X, IntegerVector Y){
  double stat=0.0;
  for(int i=0;i<X.length();i++){
    if(Y(i) == 1)
      stat += X(i);
  }
  return stat;
}

// [[Rcpp::export]]
List rcpp_Wilcoxon_PermTest(NumericVector X, IntegerVector Y,IntegerVector B,IntegerVector DoWald,IntegerVector ReportPerm) {
  
  int _perm_vec_length = 1;
  int reportperm = (ReportPerm(0) == 1);
  if(reportperm)
    _perm_vec_length = B(0);
  NumericVector perms(_perm_vec_length);
  
  double stat = Compute_Wilcoxon_Stat(X,Y);
  double pval = 0;
  double dist_bigger_or_equal_to = 0;
  int b = 0;
  int dowald = (DoWald(0) == 1);
  bool flag = true;
  IntegerVector sampled_Y;
  double current_stat;
  double Perm_Sum = 0;
  double Perm_Sum_Squared = 0;
  int Perm_Nr_performed = 0;
  //we compute the sum of the vector, and the mean of the distribution
  double _sum = 0.0;
  int _nr_1 = 0;
  double _dist_H0_mean = 0.0;
  for(int i=0;i<X.length();i++){
    _sum += X(i);
    if(Y(i)==1)
      _nr_1++;
  }
  _dist_H0_mean = _sum*((double)_nr_1)/((double)X.length());
  double stat_dist_from_mean = std::abs(stat - _dist_H0_mean);
  double _current_distance = stat_dist_from_mean;
  //NumericVector Perm_Stats(B(0)+2);
  
  if(B(0) == 0){
    flag = false;
    pval = 1.0;
  }
    
  
  while(flag){

    sampled_Y = sample(Y,Y.length(),false);
    current_stat = Compute_Wilcoxon_Stat(X,sampled_Y);
    //Rprintf("perm_stat: %lf \n\r",current_stat);
    _current_distance = std::abs(current_stat - _dist_H0_mean);

    if(_current_distance >= stat_dist_from_mean)
      dist_bigger_or_equal_to += 1.0;
    pval = (dist_bigger_or_equal_to + 1.0)/ (((double)b) + 2.0); // note that b starts at zero
    //perms(b) = current_stat;
    b++;
    if(b>=B(0))
      flag = false;
    if(dowald && b >= 1000 && pval>0.1)
      flag = false;
    if(dowald && b >= 10000 && pval>0.001)
      flag = false;
    if(dowald && b >= 20000 && pval>0.0005)
      flag = false;
    if(b>=1){
      Perm_Sum += current_stat;
      Perm_Sum_Squared += current_stat*current_stat;
      Perm_Nr_performed ++;
      if(reportperm)
        perms(b-1) = current_stat;
    }
  }
  
  //IntegerVector sampled_Y = sample(Y,Y.length(),false);
  List z  = List::create( stat, pval,b,Perm_Sum,Perm_Sum_Squared,Perm_Nr_performed,perms);
  return z ;
}


// [[Rcpp::export]]
List rcpp_Wilcoxon_PermTest_Given_Permutations(NumericVector X, IntegerMatrix Y) {
  
  NumericVector stats(Y.ncol());
  //IntegerVector _current_Y;
  for(int i=0;i<Y.ncol();i++){
    stats(i) = Compute_Wilcoxon_Stat(X,Y(_,i));
  }
  List z  = List::create( stats);
  return z ;
}

// two part test
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

// [[Rcpp::export]]
List rcpp_TwoPartTest(NumericVector X, IntegerVector Y,IntegerVector B,IntegerVector DoWald,IntegerVector ReportPerm) {
  int _perm_vec_length = 1;
  int reportperm = (ReportPerm(0) == 1);
  if(reportperm)
    _perm_vec_length = B(0);
  NumericVector perms(_perm_vec_length);
  
  NumericVector AdditionalReturnValues(7);
  NumericVector AdditionalReturnValues_Perms(7); // this is just for needing to pass a reference
  double stat = Compute_TwoPartTest(X,Y,AdditionalReturnValues);
  double pval = 0;
  double dist_bigger_or_equal_to = 0;
  int b = 0;
  int dowald = (DoWald(0) == 1);
  bool flag = true;
  IntegerVector sampled_Y;
  double current_stat;
  
  if(B(0) == 0){
    flag = false;
    pval = 1.0;
  }

  while(flag){
    
    sampled_Y = sample(Y,Y.length(),false);
    current_stat = Compute_TwoPartTest(X,sampled_Y,AdditionalReturnValues_Perms);
    //Rprintf("perm_stat: %lf \n\r",current_stat);
    if(current_stat >= stat)
      dist_bigger_or_equal_to += 1.0;
    pval = (dist_bigger_or_equal_to + 1.0)/ (((double)b) + 2.0); // note that b starts at zero
    //perms(b) = current_stat;
    b++;
    if(b>=B(0))
      flag = false;
    if(dowald && b >= 1000 && pval>0.1)
      flag = false;
    if(dowald && b >= 10000 && pval>0.001)
      flag = false;
    if(dowald && b >= 20000 && pval>0.0005)
      flag = false;
    if(reportperm)
      perms(b-1) = current_stat;
  }
  
  //IntegerVector sampled_Y = sample(Y,Y.length(),false);
  List z  = List::create( stat, pval,b,perms,AdditionalReturnValues);
  return z ;
}


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
