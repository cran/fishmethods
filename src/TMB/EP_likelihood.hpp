#include <iostream>
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
/** \brief  Approximate inverse normal cumulative distribution function, similar to R's qnorm (one-argument case only).
* \details
To be replaced by more accurate version based on Rmath library.
*/
#
template<class Type>
Type EP_likelihood(objective_function<Type>* obj)
{  
  DATA_VECTOR(age);  
  DATA_VECTOR(vlength);
  DATA_VECTOR(vn);
  DATA_VECTOR(rlength);  
  DATA_VECTOR(vp);
  DATA_VECTOR(age_ind);   
  DATA_VECTOR(len_ind);   
  DATA_SCALAR(sample_size);   
  DATA_INTEGER(bin_size);   

  PARAMETER(log_Linf);  
  PARAMETER(log_k);    
  PARAMETER(log_CV); 
  PARAMETER(ao); 
  PARAMETER_VECTOR(lambda);

  Type Linf = exp(log_Linf);  
  Type k = exp(log_k);      
  Type CV = exp(log_CV);  

  int n = age.size(); 
  int na = age_ind.size(); 
  int nr = len_ind.size(); 
  int nq = rlength.size();

  vector<Type> pl(nr);
  vector<Type> norma(na);
  vector<Type> Qh(nq);

  pl = 0;
  Qh = 0;
  norma = 0;

  for(int i = 0;i < na;++i){
    Type age_now = age_ind(i);
    Type mu_now = Linf*(Type(1)-exp(-k*(age_now-ao)));
    Type sig_now = CV*mu_now;
    norma(i) = dnorm(len_ind,mu_now,sig_now).sum();    
  }
  
  vector <Type> pa(na); 
  vector <Type> log_p(na);
  Type lterm = 0.0;
  
  for(int i = 0;i <(na-1);++i){
    lterm += log(1 + exp(lambda(i)));
    log_p(i) = lambda(i) - lterm;
  }
  log_p(na-1) = -lterm;
  pa = exp(log_p);
  
  for(int i = 0;i < nr;++i){
    Type l_now = len_ind(i);
    vector<Type> pla(na);
    for(int j = 0;j < na;++j){
      Type age_now = age_ind(j);
      Type mu_now = Linf*(Type(1)-exp(-k*(age_now-ao)));
      Type sig_now = CV*mu_now;
      pla(j) = dnorm(l_now,mu_now,sig_now);
    }
    pl(i) = ((pla/norma)*pa).sum();

  }
  
  for(int i = 0;i < nq;++i){
    for(int j = 0;j < bin_size;++j){
      Qh(i) += pl(i*bin_size+j);
    }
  }
  
  Type nll = sample_size*log((Qh*vp).sum()); 

  for(int i = 0;i < n;++i){
    Type age_now = age(i);
    Type mu_now = Linf*(Type(1)-exp(-k*(age_now-ao)));
    Type sig_now = CV*mu_now;
    Type l_now = vlength(i);
    int i_now=0;
    for(int j = 0;j < na;++j){
      if(age_ind(j)==age_now){i_now=j;break;}
    }
    Type pla = dnorm(l_now,mu_now,sig_now)/norma(i_now);
    nll -= vn(i)*log(pla*pa(i_now)+Type(0.0000000000000001));
  }
  
  REPORT(Linf);  
  REPORT(k);    
  REPORT(CV);   
  REPORT(ao);   
  REPORT(pa); 


  ADREPORT(Linf);  
  ADREPORT(k);    
  ADREPORT(CV);   
  ADREPORT(ao);   
  ADREPORT(pa); 

    
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
