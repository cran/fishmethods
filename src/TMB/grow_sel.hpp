 #include <iostream>
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type> Type square(Type x){return x*x;}
template<class Type>
Type grow_sel(objective_function<Type>* obj)
{
 DATA_VECTOR(age);
 DATA_VECTOR(size);
 DATA_VECTOR(weights);
 DATA_VECTOR(minlimit);
 DATA_VECTOR(maxlimit);
 DATA_VECTOR(minmax);
 DATA_INTEGER(switch_varpar);
 DATA_VECTOR(linf_prior);
 DATA_INTEGER(linf_pdf);
 DATA_VECTOR(k_prior);
 DATA_INTEGER(k_pdf);
 DATA_VECTOR(t0_prior);
 DATA_INTEGER(t0_pdf);
 DATA_VECTOR(varpar_prior);
 DATA_INTEGER(varpar_pdf);
 PARAMETER(linf);
 PARAMETER(k);
 PARAMETER(t0);
 PARAMETER(varpar);
   
 int n=age.size();//number of elements
 int i;
 Type sqrt_2pi=sqrt(2.0*3.14159265358979);
 Type dzero=0.00000001;
 vector<Type>size_pred(n);
 vector<Type>sigma(n);
 vector<Type>numer(n); //numerator of likelihood fcn
 vector<Type>denom1(n);  //denominator of likelihood fcn
 vector<Type>denom2(n);//denominator of likelihood fcn
 vector<Type>zscore1(n); //interim calc in likelihood fcn 
 vector<Type>zscore2(n); //interim calc in likelihood fcn 
 vector<Type>zscore3(n); //interim calc in likelihood fcn 
 Type switch_filter=0.0;//if turned on, alert user that data were filtered;switch turned off to start
 Type f_fit=0;
 Type f_prior=0;
 Type fval=0;

 //Program
 size_pred=linf*(1.-exp(-k*(age-t0)));
 switch(switch_varpar){
       case 1: //option to estimate sigma
          sigma=varpar;
          break;
       case 2: //option to estimate sigma proportional to mean (CV)
          sigma=varpar*size_pred;
          break;
       case 3: //option to estimate variance proportional to mean
          sigma=sqrt(varpar*size_pred);
          break;                
   }

   for (i=0;i<=n-1;i++){
     zscore1(i)=(size(i)-size_pred(i))/sigma(i);
     zscore2(i)=(minlimit(i)-size_pred(i))/sigma(i);
     zscore3(i)=(maxlimit(i)-size_pred(i))/sigma(i);  
   }
    
   for (i=0; i<=n-1; i++){
          numer(i)= (1.0/(sqrt_2pi*sigma(i)))*exp(-0.5*square(zscore1(i)));
	 denom1(i)=1.0-pnorm(zscore2(i));
	 denom2(i)=pnorm(zscore3(i));
	 if(minmax(i)<2){	
	   if (denom1(i)<=dzero) {denom1(i)=dzero;} //floor to avoid possibility of taking log(zero) during search process 
           if (size(i)>=minlimit(i)) {f_fit+=-weights(i)* (log(numer(i))-log(denom1(i)));} //removes fish smaller than the minimum size limit
              else {switch_filter=1.0;}
         }
        if(minmax(i)>1){    
          if (denom2(i)<=dzero) {denom2(i)=dzero;} //floor to avoid possibility of taking log(zero) during search process 
          if (size(i)<=maxlimit(i)) {f_fit+=-weights(i)* (log(numer(i))-log(denom2(i)));} //removes fish larger than the maximum size limit
             else {switch_filter=1.0;}
       }
    }
  // REPORT(numer);
  // REPORT(denom1);
  // REPORT(denom2);
   
    fval+=f_fit;
       
   //Linf prior
   Type pred=linf;
   Type prior=linf_prior(0);
   Type var=linf_prior(1);
   int pdf=linf_pdf;
   Type LkvalTmp;
   Type alpha, beta, ab_iq;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
           if(pred<=0) LkvalTmp=1000000000;
           else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
          break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-lgamma(alpha+beta)+lgamma(alpha)+lgamma(beta);
          else LkvalTmp=1000000000;
          break;
       }
   
     f_prior+=LkvalTmp; 
     //K prio
   pred=k;
   prior=k_prior(0);
   var=k_prior(1);
   pdf=k_pdf;
   
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
           if(pred<=0) LkvalTmp=1000000000;
           else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
          break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-lgamma(alpha+beta)+lgamma(alpha)+lgamma(beta);
          else LkvalTmp=1000000000;
          break;
       }
  

     f_prior+=LkvalTmp; 

   //t0 prior
   pred=t0;
   prior=t0_prior(0);
   var=t0_prior(1);
   pdf=t0_pdf;
   
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
           if(pred<=0) LkvalTmp=1000000000;
           else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
          break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-lgamma(alpha+beta)+lgamma(alpha)+lgamma(beta);
          else LkvalTmp=1000000000;
          break;
       }
      f_prior+=LkvalTmp; 

  //varpar prior
   pred=varpar;
   prior=varpar_prior(0);
   var=varpar_prior(1);
   pdf=varpar_pdf;
   
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
           if(pred<=0) LkvalTmp=1000000000;
           else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
          break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-lgamma(alpha+beta)+lgamma(alpha)+lgamma(beta);
          else LkvalTmp=1000000000;
          break;
       }

     f_prior+=LkvalTmp;
  
     fval+=f_prior;
    //  REPORT(fval);

   REPORT(switch_filter);
   REPORT(size_pred);
   REPORT(f_prior);
   REPORT(f_fit);
   REPORT(fval);
    return fval;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
