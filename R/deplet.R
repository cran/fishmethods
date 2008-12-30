deplet<-function(catch=NULL,effort=NULL,method=c("l","d","ml","hosc","hesc","hemqle","sch"),nboot=500
         ){ 
 if(any(method %in% c("l","b","ml","hosc","hesc","hemqle"))){            
     if(is.null(catch)) 
         stop ("catch vector does not exist")
     if(is.null(effort)) 
         stop ("effort vector does not exist")
      if(length(catch)!=length(effort)) 
        stop("unequal vector lengths")
   }

  if(any(method=="sch")){
     if(is.null(catch)) 
         stop ("catch vector does not exist") 
     }
    x<-as.data.frame(cbind(catch,effort))
    names(x)<-c("catch","effort")
    x$samp<-seq(1,length(x$catch),1)
    x$cpue<-x$catch/x$effort
    x$K[1]<-0;x$F[1]<-0
    nsam<-length(x$samp)
    x$K[c(seq(2,nsam,1))]<-cumsum(x$catch[c(seq(1,nsam-1,1))])
    x$F[c(seq(2,nsam,1))]<-cumsum(x$effort[c(seq(1,nsam-1,1))])

  if(length(x$catch)<3) stop("only two observations")
         
#estimate N and q using Leslie and Davis method
  if(length(x$catch)>=3){
     if(any(method=="l")){
       ld<-lm(cpue~K,data=x)
       names(ld$coefficients)[2]<-"q"
       Nld<--coef(summary(ld))[1]/coef(summary(ld))[2]
       qld<--coef(summary(ld))[2]
       s2ld<-summary(ld)$sigma^2
       SEnld<-sqrt((s2ld/qld^2)*((1/nsam)+(Nld-mean(x$K))^2/sum((x$K-mean(x$K))^2)))
       CInld<-c(Nld-qt(0.975,nsam-2)*SEnld,Nld+qt(0.975,nsam-2)*SEnld)
       SEqld<-coef(summary(ld))[4]
       CIqld<-c(qld-qt(0.975,nsam-2)*SEqld,qld+qt(0.975,nsam-2)*SEqld)
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Nld,0),round(SEnld,1),round(CInld[1],1),round(CInld[2],1)),
                  cbind(as.numeric(qld),as.numeric(SEqld),CIqld[1],CIqld[2]))
         dimnames(ans$results)<-list(cbind("N","q"),c("Estimate","SE","95% LCI","95% UCI"))
         ans$gof<-matrix(NA,2L,1L)
         ans$gof<-rbind(round(sum((x$catch-predict(ld)*x$effort)^2/(predict(ld)*x$effort)),1), 
                   as.numeric(summary(ld)$fstatistic[3]))
         dimnames(ans$gof)<-list(c("Chi-square","df"),"Value")
         ans$anova<-anova(ld)   
         ans$summary<-summary(ld) 
         ans$residuals<-as.vector(residuals(ld))
         l.out<<-ans
     }
#estimate N and q using effort-corrected Delury method
    if(any(method=="d")){
  	x$E<-x$F+x$effort*0.5
   	cd<-lm(log(cpue)~E,data=x)
      names(cd$coefficients)[2]<-"q"
   	Ncd<-exp(coef(summary(cd))[1])/-coef(summary(cd))[2]
   	qcd<--coef(summary(cd))[2]
   	s2cd<-summary(cd)$sigma^2
   	SEncd<-sqrt((s2cd*Ncd^2)*(1/nsam+((qcd*mean(x$E)-1)/qcd)^2*(1/sum((x$E-mean(x$E))^2))))
   	CIncd<-c(Ncd-qt(0.975,nsam-2)*SEncd,Ncd+qt(0.975,nsam-2)*SEncd)
   	SEqcd<-coef(summary(cd))[4]
   	CIqcd<-c(qcd-qt(0.975,nsam-2)*SEqcd,qcd+qt(0.975,nsam-2)*SEqcd)
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Ncd,0),round(SEncd,1),round(CIncd[1],1),round(CIncd[2],1)),
                  cbind(as.numeric(qcd),as.numeric(SEqcd),CIqcd[1],CIqcd[2]))
         dimnames(ans$results)<-list(cbind("N","q"),c("Estimate","SE","95% LCI","95% UCI"))
  
         ans$gof<-matrix(NA,2L,1L)
         ans$gof<-rbind(round(sum((x$catch-exp(predict(cd))*x$effort)^2/(exp(predict(cd))*x$effort)),1), 
                   as.numeric(summary(cd)$fstatistic[3]))
         dimnames(ans$gof)<-list(c("Chi-square","df"),"Value")
         ans$anova<-anova(cd)   
         ans$summary<-summary(cd) 
         ans$residuals<-as.vector(residuals(cd))
         d.out<<-ans
     }

 if(any(method=="sch")){
     	nsam<-length(catch)
	pcatch1<-c(rep(0,nsam))
	T<-c(rep(0,nsam))
	T[c(seq(1,nsam,1))]<-cumsum(catch[c(seq(1,nsam,1))])
	TC<-sum(catch)
	N<-TC*1.2
	q1<-T[nsam]/(nsam*N-sum(T[1:nsam-1]))
	parms1<-c(N,q1)                                   
      model1<-function(x){					
   		N<-x[1]
  		q1<-x[2]
  		pcatch1[seq(1,nsam,1)]<-N*q1*(1-q1)^(seq(1,nsam,1)-1)
   		pT<-sum(pcatch1)
 		L1<-sum(log((seq(1,T[nsam],1)+N-T[nsam])/seq(1,T[nsam],1)))
  		H<-sum(catch[seq(1,nsam,1)]*log(catch[seq(1,nsam,1)]/pcatch1[seq(1,nsam,1)]))  
  		G<-N*log(N)-T[nsam]*log(T[nsam])-(N-T[nsam])*log(N-pT)-L1
   		G+H
	   }   
	lower<-c(sum(catch),0)                  
	upper<-c(sum(catch)*20,2)                 #upper limits of N and q
	results1<-optim(parms1, model1, gr = NULL,lower=lower,upper=upper,method=c("L-BFGS-B"), 
      control=list(maxit=100000),hessian=TRUE)
	cov1<-solve(results1$hessian)
      pcatch1[seq(1,nsam,1)]<-results1$par[1]*results1$par[2]*(1-results1$par[2])^(seq(1,nsam,1)-1)
      resid1<-catch-pcatch1
         t<-1.96
         ans<-NULL  
         ans$Model_1<-matrix(NA,2L,4L)
         ans$Model_1<-rbind(cbind(round(results1$par[1],1),
                              round(sqrt(cov1[1,1]),2),
                              round(results1$par[1],1)-t*round(sqrt(cov1[1,1]),2),
                              round(results1$par[1],1)+t*round(sqrt(cov1[1,1]),2)),
                            cbind(results1$par[2],sqrt(cov1[2,2]),
                               results1$par[2]-t*sqrt(cov1[2,2]),
                               results1$par[2]+t*sqrt(cov1[2,2])))
         dimnames(ans$Model_1)<-list(cbind("N","q1"),c("Estimate","SE","95% LCI","95% UCI"))
         ans$Model1_Parameter_Correlation<-cov2cor(cov1)
      #Model 2  
	pcatch2<-c(rep(0,nsam))
      T<-c(rep(0,nsam))
	T[c(seq(1,nsam,1))]<-cumsum(catch[c(seq(1,nsam,1))])
	TC<-sum(catch)
	N<-TC*1.2
	q12<-catch[1]/N
	q2<-(T[nsam]-catch[1])/
      	((nsam-1)*(N-catch[1])-sum(T[1:nsam-1]-catch[1]))
  	parms2<-c(N,q12,q2)
  	model2<-function(x){	
   		N<-x[1]
  		q12<-x[2]
   		q2<-x[3]
   		pcatch2[1]<-q12*N
   		pcatch2[seq(2,nsam,1)]<-N*q2*(1-q12)*(1-q2)^(seq(2,nsam,1)-2)
   		pT<-sum(pcatch2)
  		L1<-sum(log((seq(1,T[nsam],1)+N-T[nsam])/seq(1,T[nsam],1)))
  		H<-sum(catch[seq(1,nsam,1)]*log(catch[seq(1,nsam,1)]/pcatch2[seq(1,nsam,1)]))  
  		G<-N*log(N)-T[nsam]*log(T[nsam])-(N-T[nsam])*log(N-pT)-L1
   		G+H
		}
	lower<-c(sum(catch),0,0)                   #lower limits of N and q
	upper<-c(sum(catch)*20,2,2)                 #upper limits of N and q
	results2<-optim(parms2, model2, gr = NULL,lower=lower,upper=upper,method=c("L-BFGS-B"), 
      control=list(maxit=100000),hessian=TRUE)
	cov2<-solve(results2$hessian)
      	pcatch2[1]<-results2$par[2]*results2$par[1]
   		pcatch2[seq(2,nsam,1)]<-results2$par[1]*results2$par[3]*(1-results2$par[2])*(1-results2$par[3])^(seq(2,nsam,1)-2)
   		resid2<-catch-pcatch2

         t<-1.96
         ans$Model_2<-matrix(NA,3L,4L)
         ans$Model_2<-rbind(cbind(round(results2$par[1],1),round(sqrt(cov2[1,1]),2),
                              round(results2$par[1],1)-t*round(sqrt(cov2[1,1]),2),
                              round(results2$par[1],1)+t*round(sqrt(cov2[1,1]),2)),
                            cbind(results2$par[2],sqrt(cov2[2,2]),
                              results2$par[2]-t*sqrt(cov2[2,2]),
                              results2$par[2]+t*sqrt(cov2[2,2])),
                            cbind(results2$par[3],sqrt(cov2[3,3]),
                              results2$par[3]-t*sqrt(cov2[3,3]),
                              results2$par[3]+t*sqrt(cov2[3,3])))
         dimnames(ans$Model_2)<-list(cbind("N","q1","q2"),c("Estimate","SE","95% LCI","95% UCI"))
          ans$Model2_Parameter_Correlation<-cov2cor(cov2)

         ans$gof<-matrix(NA,2L,3L)
         ans$gof<-rbind(cbind(sum((catch-pcatch1)^2/pcatch1),nsam-2),
                        cbind(sum((catch-pcatch2)^2/pcatch2),nsam-3)) 
                   
         dimnames(ans$gof)<-list(c("Model 1","Model 2"),c("Chi-square","df"))

	   ans$comparison<-matrix(NA,3L,1L)
         ans$comparison<-rbind(round(results1$value,2),round(results2$value,2),
                          2*(round(results1$value,2)-round(results2$value,2)))
          dimnames(ans$comparison)<-list(c("Model 1 F","Model 2 F","Chi-square"),c("Value"))
          ans$residuals<-data.frame(Model_1=resid1,Model_2=resid2)
          sch.out<<-ans
  }
    
     if(any(method=="ml")){
        dd<-x
        parms<-c(k=0.0005)
        ml<-function(y){
          k<-y[1] 
          dd$p<-1-exp(-k*dd$effort)
          dd$q<-1-dd$p
            for(i in 1:as.numeric(length(dd$catch))){
              dd$nf[i]<-sum(log(seq(1,dd$catch[i],1)))
              if(i==1) dd$qp[i]<-dd$p[i]
              if(i>1)  dd$qp[i]<-dd$p[i]*prod(dd$q[1:i-1])
            }
            Q<<-1-sum(dd$qp)
            dd$pn<-dd$catch*log(dd$qp/(1-Q))
            SP<-sum(dd$pn)
            XS<-sum(log(seq(1,sum(dd$catch),1)))
            LL2<-((XS-sum(dd$nf))+SP)
            LL2*-1 
       }
       upper<--log(0.0005)/max(dd$effort)
       Kout<-optimize(ml,lower=0,upper=upper,tol=0.0000000001)
       Korig<-Kout$minimum
       L<-Kout$objective
       Norig<-sum(dd$catch)/(1-Q)
       
        Bout<-data.frame(k=NULL,N=NULL)
        mlb<-function(y){
          k<-y[1]  
          Bdata$p<-1-exp(-k*Bdata$effort)
          Bdata$q<-1-Bdata$p
            for(i in 1:as.numeric(length(Bdata$catch))){
              dd$nf[i]<-sum(log(seq(1,Bdata$catch[i],1)))
              if(i==1) Bdata$qp[i]<-Bdata$p[i]
              if(i>1)  Bdata$qp[i]<-Bdata$p[i]*prod(Bdata$q[1:i-1])
            }
            Qb<<-1-sum(Bdata$qp)
            Bdata$pn<-Bdata$catch*log(Bdata$qp/(1-Qb))
            SP<-sum(Bdata$pn)
            XS<-sum(log(seq(1,sum(Bdata$catch),1)))
            LL2<-((XS-sum(Bdata$nf))+SP)
            LL2*-1 
         }
       for(i in 1:nboot){
       newc<-rmultinom(1, Norig, c(dd$catch/Norig,1-(sum(dd$catch)/Norig)))
       newc<-ifelse(newc==0,1,newc)
       Bdata<-data.frame(catch=newc[-length(newc)],effort=dd$effort)
       upper<--log(0.0005)/max(Bdata$effort)
       Kb<-optimize(mlb,lower=0,upper=upper,tol=0.0000000001)$minimum
       Nb<-sum(Bdata$catch)/(1-Qb)
       Bout[i,1]<-Kb;Bout[i,2]<-Nb
       }
         sek<-sd(Bout[,1]);seN<-sd(Bout[,2])
         t<-qt(0.975,length(dd$catch)-2)
         LCn<-Norig-t*seN
         UCn<-Norig+t*seN
         LCq<-Korig-t*sek
         UCq<-Korig+t*sek 
         ans<-NULL  
         ans$results<-matrix(NA,2L,6L)
         ans$results<-rbind(cbind(round(Norig,0),round(seN,1),round(LCn,1),round(UCn,1)),
                  cbind(Korig,sek,LCq,UCq))
         dimnames(ans$results)<-list(cbind("N","q"),c("Estimate","SE","95% LCI","95% UCI"))
          dd$p<-1-exp(-Korig*dd$effort)
          dd$q<-1-dd$p
            for(i in 1:as.numeric(length(dd$catch))){
              if(i==1) dd$qp[i]<-dd$p[i]
              if(i>1)  dd$qp[i]<-dd$p[i]*prod(dd$q[1:i-1])
            }
          dd$pred<-Norig*dd$qp
         
         ans$gof<-matrix(NA,2L,1L)
         ans$gof<-rbind(round(sum((dd$catch-dd$pred)^2/dd$pred),1), 
                   length(dd$catch)-2)
         dimnames(ans$gof)<-list(c("Chi-square","df"),"Value")
         ans$residuals<-as.vector(dd$catch-dd$pred)
         ml.out<<-ans
    }  
 
#Homogeneous Sample Coverage
 if(any(method=="hosc")){
     dd<-x
      nsam<-length(x[,1])
      dd$Dk<-cumsum(dd$catch)
      dd$ue<-dd$catch/dd$effort
      for(i in 1:as.numeric(nsam-1)){dd$Cp[i]<-1-(dd$ue[i+1]/(dd$ue[1]))} 
        dd$check<-ifelse(dd$Cp<=0,1,0)
        dd$check<-ifelse(dd$check[1]==1,0,ifelse(dd$check[nsam]==1,0,
               dd$check))
        tau<-ifelse(any(dd$check==TRUE),which(dd$check==1)-1,nsam-1)
        Nsc0<-dd$Dk[tau]/dd$Cp[tau]
        prob<-c(dd$catch/Nsc0,1-(sum(dd$catch)/Nsc0))
        BootNsc0<-data.frame(rep=NULL,Nscb=NULL)


     #Bootstrap
        kk<-x
     for(j in 1:nboot){
       kk$catch<-rmultinom(1,Nsc0,prob)[1:nsam]
       kk$Dk<-cumsum(kk$catch)
       kk$ue<-kk$catch/kk$effort
       for(i in 1:as.numeric(nsam-1)){kk$Cp[i]<-1-(kk$ue[i+1]/(kk$ue[1]))} 
        kk$check<-ifelse(kk$Cp<=0,1,0)
        kk$check<-ifelse(kk$check[1]==1,0,ifelse(kk$check[nsam]==1,0,
               kk$check))
       tauk<-ifelse(any(kk$check==TRUE),which(kk$check==1)-1,nsam-1)
       Nscb<-kk$Dk[tauk]/kk$Cp[tauk]
       BootNsc0[j,1]<-1
       BootNsc0[j,2]<-Nscb
      }
         seN<-sd(BootNsc0[,2])
         C<-exp(1.96*sqrt(log(1+((seN^2)/((Nsc0-sum(dd$catch))^2)))))
         LCn<-sum(dd$catch)+(Nsc0-sum(dd$catch))/C
         UCn<-sum(dd$catch)+(Nsc0-sum(dd$catch))*C
       
         ans<-NULL  
         ans$results<-matrix(NA,1L,6L)
         ans$results<-rbind(cbind(round(Nsc0,0),round(seN,1),round(LCn,1),round(UCn,1)))   
         dimnames(ans$results)<-list(cbind("N"),c("Estimate","SE","95% LCI","95% UCI"))  
         hosc.out<<-ans
   }  

 #Heterogenous Sample Coverage
   if(any(method=="hesc")){
      dd<-x
      nsam<-length(dd[,1])
      dd$Dk<-cumsum(dd$catch)
      dd$ue<-dd$catch/dd$effort
      for(i in 1:as.numeric(nsam-1)){dd$Cp[i]<-1-(dd$ue[i+1]/(dd$ue[1]))} 
      dd$check<-ifelse(dd$Cp<=0,1,0)
      dd$check<-ifelse(dd$check[1]==1,0,ifelse(dd$check[nsam]==1,0,
               dd$check))
     tau<-ifelse(any(dd$check==TRUE),which(dd$check==1)-1,nsam-1)
     Nsc0<-dd$Dk[tau]/dd$Cp[tau]
     epsilon<-max(Nsc0*((dd$catch[1]-dd$catch[2]*(dd$effort[1]/dd$effort[2]))/(dd$catch[1]^2))-1,0)
     dd$Ak<-cumsum(dd$effort)/dd$effort
     Nsc1<-(dd$Dk[tau]/dd$Cp[tau])+(dd$Ak[tau]*dd$catch[tau]*epsilon)/dd$Cp[tau]
     probb<-c(dd$catch/Nsc1,1-(sum(dd$catch)/Nsc1))
     BootNsc1<-data.frame(rep=NULL,Nscb=NULL)
     CVNsc1<-sqrt(epsilon)

     #Bootstrap
        kk<-x
     for(j in 1:nboot){
       kk$catch<-rmultinom(1,Nsc1,probb)[1:nsam]
       kk$Dk<-cumsum(kk$catch)
       kk$ue<-kk$catch/kk$effort
       for(i in 1:as.numeric(nsam-1)){kk$Cp[i]<-1-(kk$ue[i+1]/(kk$ue[1]))} 
        kk$check<-ifelse(kk$Cp<=0,1,0)
        kk$check<-ifelse(kk$check[1]==1,0,ifelse(kk$check[nsam]==1,0,
               kk$check))
       tauk<-ifelse(any(kk$check==TRUE),which(kk$check==1)-1,nsam-1)
       Nscb0<-kk$Dk[tauk]/kk$Cp[tauk]
       epsilonb<-max(Nscb0*((kk$catch[1]-kk$catch[2]*(kk$effort[1]/kk$effort[2]))/(kk$catch[1]^2))-1,0)
       kk$Akb<-cumsum(kk$effort)/kk$effort
       Nsc1b<-(kk$Dk[tauk]/kk$Cp[tauk])+(kk$Akb[tauk]*kk$catch[tauk]*epsilonb)/kk$Cp[tauk]
       BootNsc1[j,1]<-1
       BootNsc1[j,2]<-Nsc1b
      }
         seN<-sd(BootNsc1[,2])
         C<-exp(1.96*sqrt(log(1+((seN^2)/((Nsc1-sum(dd$catch))^2)))))
         LCn<-sum(dd$catch)+(Nsc1-sum(dd$catch))/C
         UCn<-sum(dd$catch)+(Nsc1-sum(dd$catch))*C
       
         ans<-NULL  
         ans$results<-matrix(NA,1L,7L)
         ans$results<-rbind(cbind(round(Nsc1,0),round(seN,1),round(LCn,1),round(UCn,1),round(CVNsc1,3)))   
         dimnames(ans$results)<-list(cbind("N"),c("Estimate","SE","95% LCI","95% UCI"," lambda CV"))  
         hesc.out<<-ans
   }
  #Heterogenous MQLE
   if(any(method=="hemqle")){
     dd<-x
      nsam<-length(dd[,1])
      dd$Dk<-cumsum(dd$catch)
      dd$ue<-dd$catch/dd$effort
       for(i in 1:as.numeric(nsam-1)){dd$Cp[i]<-1-(dd$ue[i+1]/(dd$ue[1]))} 
       dd$check<-ifelse(dd$Cp<=0,1,0)
       dd$check<-ifelse(dd$check[1]==1,0,ifelse(dd$check[nsam]==1,0,
               dd$check))
       tau<-ifelse(any(dd$check==TRUE),which(dd$check==1)-1,nsam-1)
       Nsc0<-dd$Dk[tau]/dd$Cp[tau]
       epsilon<-max(Nsc0*((dd$catch[1]-dd$catch[2]*(dd$effort[1]/dd$effort[2]))/(dd$catch[1]^2))-1,0)
       dd$Ak<-cumsum(dd$effort)/dd$effort
       dd$I<-ifelse(dd$Cp>0,1,0)
     for(i in 1:nsam){
        if(i==1){ dd$sql[i]<-0;dd$sql2<-0}
        if(i>1){
           dd$sql[i]<-(dd$effort[i]^2/dd$catch[i])*(dd$Dk[i-1]+(dd$Ak[i-1]*dd$catch[i-1]*epsilon))*dd$I[i-1]
           dd$sql2[i]<-(dd$effort[i]^2/dd$catch[i])*dd$Cp[i-1]*dd$I[i-1]
          }
     }
     Nmqle<-sum(dd$sql)/sum(dd$sql2)
     probb<-c(dd$catch/Nmqle,1-(sum(dd$catch)/Nmqle))
     BootNmq<-data.frame(rep=NULL,Nscb=NULL)
     CVNmqle<-sqrt(epsilon)

     #Bootstrap
        kk<-x
     for(j in 1:nboot){     
       kk$catch<-rmultinom(1,Nmqle,probb)[1:nsam]
       kk$catch<-ifelse(kk$catch==0,1,kk$catch)
       kk$Dk<-cumsum(kk$catch)
       kk$ue<-kk$catch/kk$effort
       for(i in 1:as.numeric(nsam-1)){kk$Cp[i]<-1-(kk$ue[i+1]/(kk$ue[1]))} 
        kk$check<-ifelse(kk$Cp<=0,1,0)
        kk$check<-ifelse(kk$check[1]==1,0,ifelse(kk$check[nsam]==1,0,
               kk$check))
       tauk<-ifelse(any(kk$check==TRUE),which(kk$check==1)-1,nsam-1)
       Nscb0<-kk$Dk[tauk]/kk$Cp[tauk]
       epsilonb<-max(Nscb0*((kk$catch[1]-kk$catch[2]*(kk$effort[1]/kk$effort[2]))/(kk$catch[1]^2))-1,0)
       kk$Akb<-cumsum(kk$effort)/kk$effort
        kk$I<-ifelse(kk$Cp>0,1,0)
       for(i in 1:nsam){
        if(i==1){ kk$sql[i]<-0;kk$sql2[i]<-0}
        if(i>1){
           kk$sql[i]<-(kk$effort[i]^2/kk$catch[i])*(kk$Dk[i-1]+(kk$Ak[i-1]*kk$catch[i-1]*epsilonb))*kk$I[i-1]
           kk$sql2[i]<-(kk$effort[i]^2/kk$catch[i])*kk$Cp[i-1]*kk$I[i-1]
         }
       }
       Nmqb<-sum(kk$sql)/sum(kk$sql2)
       BootNmq[j,1]<-i
       BootNmq[j,2]<-Nmqb
    }
         seN<-sd(BootNmq[,2])
         C<-exp(1.96*sqrt(log(1+((seN^2)/((Nmqle-sum(dd$catch))^2)))))  
         LCn<-sum(dd$catch)+(Nmqle-sum(dd$catch))/C
         UCn<-sum(dd$catch)+(Nmqle-sum(dd$catch))*C
         ans<-NULL  
         ans$results<-matrix(NA,1L,7L)
         ans$results<-rbind(cbind(round(Nmqle,0),round(seN,1),round(LCn,1),round(UCn,1),round(CVNmqle,3)))   
         dimnames(ans$results)<-list(cbind("N"),c("Estimate","SE","95% LCI","95% UCI"," lambda CV"))  
         hemqle.out<<-ans
  }
 }
}

