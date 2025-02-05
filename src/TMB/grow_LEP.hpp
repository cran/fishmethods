#include<iostream>
#include <Eigen/Dense>
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type> Type calc_dlnorm(Type x, Type mean, Type sd) {
  Type resid=log(x)-log(mean);
  Type denom=2*sd*sd;
  Type constant=Type(sqrt(2*M_PI));
  Type ans=log((1/(x*constant*sd))*exp(-(resid*resid)/denom));
  return ans;
}
template<class Type> vector<Type> dlnorm(vector<Type> x, Type mean, Type sd) {
  vector<Type> result(x.size());
  for (int i = 0; i < x.size(); i++) {
    result(i) = calc_dlnorm(x(i),mean,sd);
  }
  return result;
}

template<class Type>Type grow_LEP(objective_function<Type>* obj)
{
  //Version 4
  DATA_VECTOR(l1);
  DATA_VECTOR(l2);
  DATA_VECTOR(dt);
  DATA_VECTOR(measure);
  DATA_INTEGER(model1);
  DATA_INTEGER(fix_beta);
  DATA_INTEGER(model_measure);
  int n = l1.size();
  PARAMETER(log_K1);
  PARAMETER(log_K2);
  PARAMETER(log_mu_Linf);
  PARAMETER(log_mean_age);
  PARAMETER(log_sigma_mu_Linf);
  PARAMETER(log_sigma_age);
  PARAMETER(log_sigma_resid);
  PARAMETER_VECTOR(A); // random effect for age
  PARAMETER(log_alpha);
  PARAMETER(log_beta);
  PARAMETER(log_measure);

  vector<Type>mu1(n);
  vector<Type>mu2(n);
  vector<Type>sigma2_l1(n);
  vector<Type>sigma2_l2(n);
  vector<Type>cov(n);
  vector<Type>residuals(n);
  vector<Type>lresids1(n);
  vector<Type>lresids2(n);
  vector<Type>pred_Linf(n);
  Type f1=0.0;
  Type f2=0.0;
  vector<Type>hl1l2a(n);
  Type K1=exp(log_K1);
  Type K2=exp(log_K2);
  Type mu_Linf=exp(log_mu_Linf);
  Type sigma_mu_Linf=exp(log_sigma_mu_Linf);
  Type mean_age=exp(log_mean_age);
  Type sigma_age=exp(log_sigma_age);
  Type sigma_resid=exp(log_sigma_resid);
  Type beta=exp(log_beta);
  Type alpha=exp(log_alpha);
  Type sigma_measure=exp(log_measure);
  Type p1;
  Type p2;
  Type eq11;
  Type eq12;
  Type eq21;
  Type rho;
  Type q1;
  Type q2;
  Type q3;
  Type q12;
  Type nll = 0.0;
  using namespace std;
  using namespace density;
  //contribution of the A random effects
  nll -= sum(dlnorm(A,mean_age,sigma_age)); //mean_age is log-transformed in function
 
  for (int i=0;i<n;i++){
    if(model1==1){
      f1=1-exp(-K1*A(i));
      f2=1-exp(-K1*(A(i)+dt(i)));
    }
    if(model1==2){
      eq11=(1+exp(-beta*(A(i)-alpha)))/(1+exp(beta*alpha));
      eq12=(-(K2-K1))/beta;
      eq21=((1+exp(-beta*(A(i)+dt(i)-alpha)))/(1+exp(beta*alpha)));
      f1=1-exp(-K2*(A(i)))*pow(eq11,eq12);
      f2=1-exp(-K2*(A(i)+dt(i)))*pow(eq21,eq12);
    }
    mu1(i) = mu_Linf*f1;
    mu2(i) = mu_Linf*f2;
    sigma2_l1(i)=square(sigma_mu_Linf)*square(f1)+square(sigma_resid);
    sigma2_l2(i)=square(sigma_mu_Linf)*square(f2)+square(sigma_resid)+square(sigma_measure*measure(i));
    cov(i)=square(sigma_mu_Linf)*f1*f2;
    rho=cov(i)/(sqrt(sigma2_l1(i))*sqrt(sigma2_l2(i)));
    q1=square(l1(i)-mu1(i))/sigma2_l1(i);
    q2=2*rho*(((l1(i)-mu1(i))*(l2(i)-mu2(i)))/(sqrt(sigma2_l1(i))*sqrt(sigma2_l2(i))));
    q3=square(l2(i)-mu2(i))/sigma2_l2(i);
    q12=q1-q2+q3;
    lresids1(i)=l1(i)-mu1(i);
    lresids2(i)=l2(i)-mu2(i);
    hl1l2a(i)=log((1/(Type(2*M_PI)*sqrt(sigma2_l1(i))*sqrt(sigma2_l2(i))*sqrt(1-square(rho))))*exp(-q12/(2*(1-square(rho)))));
    p1=square(sigma_mu_Linf)/(square(sigma_resid)+square(sigma_mu_Linf)*(square(f1)+square(f2)));
    p2=f1*(l1(i)-mu_Linf*f1)+f2*(l2(i)-mu_Linf*f2);
    pred_Linf(i)=mu_Linf+p1*p2;
  }
  nll-=sum(hl1l2a);
  ADREPORT(mu_Linf);
  ADREPORT(sigma_mu_Linf);
  ADREPORT(K1);
  if(model1==2){
    ADREPORT(K2);
    ADREPORT(alpha);
    if(fix_beta==0) ADREPORT(beta);
  }
  ADREPORT(mean_age);
  ADREPORT(sigma_age);
  ADREPORT(sigma_resid);
  if(model_measure==1) ADREPORT(sigma_measure);
  REPORT(pred_Linf);
  REPORT(lresids1);
  REPORT(lresids2);
  REPORT(mu1);
  REPORT(mu2);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
