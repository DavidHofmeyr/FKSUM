#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


double d_abs(double x){
  if(x>0) return x;
  else return -x;
}



// [[Rcpp::export]]

double fk_md(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k+1)*betas[k+1]-betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0;
  double stdv;
  NumericVector df(n_eval);
  if(al>.0000000001){
    double var = 0;
    for(int i=0; i<n; i++) var += x[i]*x[i];
    var /= (n-1);
    stdv = pow(var,.5);
    std::sort(x.begin(),x.end());
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
          df[i] -= tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult;
          for(int j=0; j<=orddo; j++) df[i] -= tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) df[i] += tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]<(-al*stdv)) df[i] -= 2.0*C*(-al*stdv-x_eval[i]);
      if(x_eval[i]>(al*stdv)) df[i] += 2.0*C*(-al*stdv+x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 10000000000000;
    int pos = 0;
    double eps = pow(0.1, 10);
    while(pos<(n_eval-1)){
      if(df[pos]<0 && df[pos+1]>0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos+1];
        int ix = round(counts[pos]);
        while((hi-lo)>eps && d_abs(df_at_min)>eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5*lo+0.5*hi;
          exp_mult = exp((x[ix-1]-mid)/h);
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
            for(int j=0; j<=orddo; j++){
              df_at_min += tbetas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult-pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
              f_at_min += betas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult+pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
            }
          }
          if(mid<(-al*stdv)){
            df_at_min -= 2.0*C*(-al*stdv-mid);
            f_at_min += C*pow(-al*stdv-mid,2);
          }
          if(mid>(al*stdv)){
            df_at_min += 2.0*C*(-al*stdv+mid);
            f_at_min += C*pow(-al*stdv+mid,2);
          }
          if(df_at_min<(-eps)) lo = mid;
          if(df_at_min>eps) hi = mid;
        }
        if(f_at_min<miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else{
    for(int i=1; i<n; i++){
      double add = 0;
      for(int j=0; j<=ord; j++) add += betas[j]*pow(d_abs(x[i]/h),j);
      miny += add*exp(-d_abs(x[i])/h);
    }
  }
  return miny;
}



// [[Rcpp::export]]

NumericVector fk_md_dp(NumericVector xo, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  int n = xo.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k+1)*betas[k+1]-betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0;
  double stdv;
  NumericVector df(n_eval);
  if(al>.0000000001){
    double var = 0;
    NumericVector x(n);
    for(int i=0; i<n; i++){
      x[i] = xo[i];
      var += x[i]*x[i];
    }
    var /= (n-1);
    stdv = pow(var,.5);
    std::sort(x.begin(),x.end());
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
          df[i] -= tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult;
          for(int j=0; j<=orddo; j++) df[i] -= tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) df[i] += tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]<(-al*stdv)) df[i] -= 2.0*C*(-al*stdv-x_eval[i]);
      if(x_eval[i]>(al*stdv)) df[i] += 2.0*C*(-al*stdv+x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 10000000000000;
    int pos = 0;
    double eps = pow(0.1, 10);
    while(pos<(n_eval-1)){
      if(df[pos]<0 && df[pos+1]>0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos+1];
        int ix = round(counts[pos]);
        while((hi-lo)>eps && d_abs(df_at_min)>eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5*lo+0.5*hi;
          exp_mult = exp((x[ix-1]-mid)/h);
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
            for(int j=0; j<=orddo; j++){
              df_at_min += tbetas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult-pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
              f_at_min += betas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult+pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
            }
          }
          if(mid<(-al*stdv)){
            df_at_min -= 2.0*C*(-al*stdv-mid);
            f_at_min += C*pow(-al*stdv-mid,2);
          }
          if(mid>(al*stdv)){
            df_at_min += 2.0*C*(-al*stdv+mid);
            f_at_min += C*pow(-al*stdv+mid,2);
          }
          if(df_at_min<(-eps)) lo = mid;
          if(df_at_min>eps) hi = mid;
        }
        if(f_at_min<miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else minx = 0;
  NumericVector dp(n);
  double add;
  for(int i=0; i<n; i++){
    add = 0;
    if(xo[i]>minx){
      for(int j=0; j<=ord; j++) add += tbetas[j]*pow((xo[i]-minx)/h,j);
      dp[i] = add*exp((minx-xo[i])/h);
    }
    else{
      for(int j=0; j<=ord; j++) add -= tbetas[j]*pow((minx-xo[i])/h,j);
      dp[i] = add*exp((xo[i]-minx)/h);
    }
  }
  if(minx<(-al*stdv)){
    double cnst = 2.0*al*C/stdv/(n-1.0)*(minx+al*stdv);
    for(int i=0; i<n; i++) dp[i] += cnst*xo[i];
  }
  if(minx>(al*stdv)){
    double cnst = 2.0*al*C/stdv/(n-1.0)*(al*stdv-minx);
    for(int i=0; i<n; i++) dp[i] += cnst*xo[i];
  }
  return(dp);
}



// [[Rcpp::export]]

double fk_is_minim_md(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k+1)*betas[k+1]-betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0;
  double stdv;
  double modef;
  double mode1 = 10000000000;
  NumericVector df(n_eval);
  NumericVector dfpen(n_eval);
  if(al>.0000000001){
    double var = 0;
    for(int i=0; i<n; i++) var += x[i]*x[i];
    var /= (n-1);
    stdv = pow(var,.5);
    std::sort(x.begin(),x.end());
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
          df[i] -= tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult;
          for(int j=0; j<=orddo; j++) df[i] -= tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) df[i] += tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]<(-al*stdv)) dfpen[i] -= 2.0*C*(-al*stdv-x_eval[i]);
      if(x_eval[i]>(al*stdv)) dfpen[i] += 2.0*C*(-al*stdv+x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 10000000000000;
    int pos = 0;
    double eps = pow(0.1, 10);
    while(pos<(n_eval-1)){
      if(df[pos]>0 && df[pos+1]<0){
        if(x_eval[pos]<mode1){
          mode1 = x_eval[pos];
        }
        modef = x_eval[pos];
      }
      if((df[pos]+dfpen[pos])<0 && (df[pos+1]+dfpen[pos+1])>0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos+1];
        int ix = round(counts[pos]);
        while((hi-lo)>eps && d_abs(df_at_min)>eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5*lo+0.5*hi;
          exp_mult = exp((x[ix-1]-mid)/h);
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
            for(int j=0; j<=orddo; j++){
              df_at_min += tbetas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult-pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
              f_at_min += betas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult+pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
            }
          }
          if(mid<(-al*stdv)){
            df_at_min -= 2.0*C*(-al*stdv-mid);
            f_at_min += C*pow(-al*stdv-mid,2);
          }
          if(mid>(al*stdv)){
            df_at_min += 2.0*C*(-al*stdv+mid);
            f_at_min += C*pow(-al*stdv+mid,2);
          }
          if(df_at_min<(-eps)) lo = mid;
          if(df_at_min>eps) hi = mid;
        }
        if(f_at_min<miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else minx = 0;
  double ret = 0;
  if(minx>mode1 && minx<modef) ret = 1;
  return(ret);
}


// [[Rcpp::export]]

double fk_md_b(NumericVector x, NumericVector y, NumericVector x_eval, double h, NumericVector betas, double al, double C){
  int n = x.size();
  int n_eval = x_eval.size();
  int ord = betas.size()-1;
  NumericVector tbetas(ord+1);
  for(int k=0; k<ord; k++) tbetas[k] = (k+1)*betas[k+1]-betas[k];
  tbetas[ord] = -betas[ord];
  double miny;
  double minx = 0;
  double stdv;
  NumericVector df(n_eval);
  if(al>.0000000001){
    double var = 0;
    for(int i=0; i<n; i++) var += x[i]*x[i];
    var /= (n-1);
    stdv = pow(var,.5);
    std::sort(x.begin(),x.end());
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
          df[i] -= tbetas[orddo]*pow(x[0]-x_eval[i], orddo)/denom*exp_mult;
          for(int j=0; j<=orddo; j++) df[i] -= tbetas[orddo]*coefs[j]*pow(-x_eval[i],orddo-j)*Ry(j,0)/denom*exp_mult;
        }
        else{
          exp_mult = exp((x[ix-1]-x_eval[i])/h);
          for(int j=0; j<=orddo; j++) df[i] += tbetas[orddo]*coefs[j]*(pow(x_eval[i], orddo-j)*Ly(j,ix-1)*exp_mult-pow(-x_eval[i],orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
        }
      }
    }
    for(int i=0; i<n_eval; i++){
      if(x_eval[i]<(-al*stdv)) df[i] -= 2.0*C*(-al*stdv-x_eval[i]);
      if(x_eval[i]>(al*stdv)) df[i] += 2.0*C*(-al*stdv+x_eval[i]);
    }
    double f_at_min, df_at_min, lo, hi, mid;
    miny = 10000000000000;
    int pos = 0;
    double eps = pow(0.1, 10);
    while(pos<(n_eval-1)){
      if(df[pos]<0 && df[pos+1]>0){
        df_at_min = 1.0;
        lo = x_eval[pos];
        hi = x_eval[pos+1];
        int ix = round(counts[pos]);
        while((hi-lo)>eps && d_abs(df_at_min)>eps){
          df_at_min = 0.0;
          f_at_min = 0;
          mid = 0.5*lo+0.5*hi;
          exp_mult = exp((x[ix-1]-mid)/h);
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
            for(int j=0; j<=orddo; j++){
              df_at_min += tbetas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult-pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
              f_at_min += betas[orddo]*coefs[j]*(pow(mid, orddo-j)*Ly(j,ix-1)*exp_mult+pow(-mid,orddo-j)*Ry(j,ix-1)/exp_mult)/denom;
            }
          }
          if(mid<(-al*stdv)){
            df_at_min -= 2.0*C*(-al*stdv-mid);
            f_at_min += C*pow(-al*stdv-mid,2);
          }
          if(mid>(al*stdv)){
            df_at_min += 2.0*C*(-al*stdv+mid);
            f_at_min += C*pow(-al*stdv+mid,2);
          }
          if(df_at_min<(-eps)) lo = mid;
          if(df_at_min>eps) hi = mid;
        }
        if(f_at_min<miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else{
    for(int i=1; i<n; i++){
      double add = 0;
      for(int j=0; j<=ord; j++) add += betas[j]*pow(d_abs(x[i]/h),j);
      miny += add*exp(-d_abs(x[i])/h);
    }
  }
  return minx;
}
