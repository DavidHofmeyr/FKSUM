#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]

NumericVector cbin_alloc(NumericVector x, int nbin, double min, double max){
  int n = x.size();
  NumericVector output(n);
  double skip = (nbin-1.0)/(max-min);
  for(int i=0; i<n; i++) output[i] = floor((x[i]-min)*skip) + 1;
  return output;
}


// [[Rcpp::export]]

NumericVector sm_bin_wts(NumericVector x, NumericVector y, int nbin, double min, double max){
  int n = x.size();
  NumericVector output(nbin);
  double skip = (nbin-1.0)/(max-min);
  int fl = 0;
  double xe;
  for(int i=0; i<n; i++){
    xe = (x[i]-min)*skip;
    fl = floor(xe);
    if(fl<(nbin-1) && fl>=0){
      output[fl+1] += (xe-fl)*y[i];
      output[fl] += (fl+1-xe)*y[i];
    }
    else if(fl>=(nbin-1)) output[nbin-1] += y[i];
    else output[0] += y[i];
  }
  return output;
}



// [[Rcpp::export]]

NumericVector bin_wts(NumericVector x, NumericVector y, int nbin, double min, double max){
  int n = x.size();
  NumericVector output(nbin);
  double skip = (nbin-1.0)/(max-min);
  int fl = 0;
  double xe;
  for(int i=0; i<n; i++){
    xe = (x[i]-min)*skip;
    fl = floor(xe);
    if(fl<(nbin-1) && fl>=0){
      output[fl] += y[i];
    }
    else if(fl>=(nbin-1)) output[nbin-1] += y[i];
    else output[0] += y[i];
  }
  return output;
}



/*
pair (x, y) should be sorted according to values in x (non-decreasing). x_eval short be sorted (non-decreasing).
*/

// [[Rcpp::export]]

NumericVector ksum(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, NumericVector Counts = NumericVector(1)){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericVector output(n_eval);
  double denom;
  double exp_mult;
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  for(int i=0; i<=ord; i++) Ly(i,0) = pow(-x[0], i)*y[0];
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      Ly(j,i) = pow(-x[i],j)*y[i] + exp((x[i-1]-x[i])/h)*Ly(j,i-1);
      Ry(j,n-i-1) = exp((x[n-i-1]-x[n-i])/h)*(pow(x[n-i],j)*y[n-i]+Ry(j,n-i));
    }
  }
  if(Counts.size()==1){
    NumericVector counts(n_eval);
    int count = 0;
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]>=x[n-1]){
        for(int j=i; j<n_eval; j++) counts[j] = n;
        break;
      }
      else{
        while(x[count]<=x_eval[i]){
          count += 1;
          if(count==n) break;
        }
        counts[i] = count;
      }
    }
    for(int orddo=0; orddo<=ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j=2; j<=orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num/orddo;
        for(int i=2; i<=orddo; i++){
          coefs[i-1] = num/denom1/denom2;
          denom1 *= i;
          denom2 /= (orddo-i+1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i=0; i<n_eval; i++){
        ix = round(counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i]-x[0])/h);
          output[i] += betas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult+pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
  }
  else{
    for(int orddo=0; orddo<=ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j=2; j<=orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num/orddo;
        for(int i=2; i<=orddo; i++){
          coefs[i-1] = num/denom1/denom2;
          denom1 *= i;
          denom2 /= (orddo-i+1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i=0; i<n_eval; i++){
        ix = round(Counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i]-x[0])/h);
          output[i] += betas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) output[i] += betas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult+pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
  }
  return output;
}



/*
pair (x, y) should be sorted according to values in x (non-decreasing). x_eval short be sorted (non-decreasing).
*/

// [[Rcpp::export]]

NumericVector dksum(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, NumericVector Counts = NumericVector(0)){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericVector output(n_eval);
  double denom;
  double exp_mult;
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k+1)*betas[k+1]-betas[k];
  tbetas[ord] = -betas[ord];
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  for(int i=0; i<=ord; i++) Ly(i,0) = pow(-x[0], i)*y[0];
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      Ly(j,i) = pow(-x[i],j)*y[i] + exp((x[i-1]-x[i])/h)*Ly(j,i-1);
      Ry(j,n-i-1) = exp((x[n-i-1]-x[n-i])/h)*(pow(x[n-i],j)*y[n-i]+Ry(j,n-i));
    }
  }
  if(Counts.size()==0){
    int count = 0;
    NumericVector counts(n_eval);
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]>=x[n-1]){
        for(int j=i; j<n_eval; j++) counts[j] = n;
        break;
      }
      else{
        while(x[count]<=x_eval[i]){
          count += 1;
          if(count>=(n-.0000001)) break;
        }
        counts[i] = count;
      }
    }
    for(int orddo=0; orddo<=ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j=2; j<=orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num/orddo;
        for(int i=2; i<=orddo; i++){
          coefs[i-1] = num/denom1/denom2;
          denom1 *= i;
          denom2 /= (orddo-i+1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i=0; i<n_eval; i++){
        ix = round(counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i]-x[0])/h);
          output[i] += tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          for(int j=0; j<=orddo; j++) output[i] += tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) output[i] -= tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
  }
  else{
    for(int orddo=0; orddo<=ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j=2; j<=orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num/orddo;
        for(int i=2; i<=orddo; i++){
          coefs[i-1] = num/denom1/denom2;
          denom1 *= i;
          denom2 /= (orddo-i+1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i=0; i<n_eval; i++){
        ix = round(Counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i]-x[0])/h);
          output[i] += tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          for(int j=0; j<=orddo; j++) output[i] += tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) output[i] -= tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
  }
  return output;
}



/*
pair (x, y) should be sorted according to values in x (non-decreasing). x_eval short be sorted (non-decreasing).
*/

// [[Rcpp::export]]

NumericMatrix kndksum(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, NumericVector Counts = NumericVector(0)){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericMatrix output(n_eval, 2);
  double denom;
  double exp_mult;
  NumericMatrix Ly(ord + 1, n);
  NumericMatrix Ry(ord + 1, n);
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k+1)*betas[k+1]-betas[k];
  tbetas[ord] = -betas[ord];
  for(int i=0; i<=ord; i++) Ly(i,0) = pow(-x[0], i)*y[0];
  for(int i=1; i<n; i++){
    for(int j=0; j<=ord; j++){
      Ly(j,i) = pow(-x[i],j)*y[i] + exp((x[i-1]-x[i])/h)*Ly(j,i-1);
      Ry(j,n-i-1) = exp((x[n-i-1]-x[n-i])/h)*(pow(x[n-i],j)*y[n-i]+Ry(j,n-i));
    }
  }
  if(Counts.size()==0){
    int count;
    NumericVector counts(n_eval);
    count = 0;
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]>=x[n-1]){
        for(int j=i; j<n_eval; j++) counts[j] = n;
        break;
      }
      else{
        while(x[count]<=x_eval[i]){
          count += 1;
          if(count==n) break;
        }
        counts[i] = count;
      }
    }
    for(int orddo=0; orddo<=ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j=2; j<=orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num/orddo;
        for(int i=2; i<=orddo; i++){
          coefs[i-1] = num/denom1/denom2;
          denom1 *= i;
          denom2 /= (orddo-i+1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i=0; i<n_eval; i++){
        ix = round(counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i]-x[0])/h);
          output(i,0) += betas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          output(i,1) += tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          for(int j=0; j<=orddo; j++){
            output(i,0) += betas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
            output(i,1) += tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
          }
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++){
            output(i,0) += betas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult+pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
            output(i,1) -= tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
          }
        }
      }
    }
  }
  else{
    for(int orddo=0; orddo<=ord; orddo++){
      NumericVector coefs(orddo + 1);
      coefs[0] = coefs[orddo] = 1;
      if(orddo>1){
        double num = 1;
        for(int j=2; j<=orddo; j++) num *= j;
        double denom1 = 1;
        double denom2 = num/orddo;
        for(int i=2; i<=orddo; i++){
          coefs[i-1] = num/denom1/denom2;
          denom1 *= i;
          denom2 /= (orddo-i+1);
        }
      }
      denom = pow(h, orddo);
      int ix;
      for(int i=0; i<n_eval; i++){
        ix = round(Counts[i]);
        if(ix==0){
          exp_mult = exp((x_eval[i]-x[0])/h);
          output(i,0) += betas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          output(i,1) += tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult*y[0];
          for(int j=0; j<=orddo; j++){
            output(i,0) += betas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
            output(i,1) += tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
          }
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++){
            output(i,0) += betas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult+pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
            output(i,1) -= tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
          }
        }
      }
    }
  }
  return output;
}
