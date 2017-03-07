###############################################################################
#
#                  Survival Rate Estimators
#
###############################################################################

agesurv<-function(age=NULL,full=NULL,last=NULL,estimate=c("s","z"),method=c("lr","he","cr","crcb","ripois","wlr","pois")){
     if(is.null(age)) 
         stop ("vector does not exist")
      if(!is.numeric(age)) 
         stop ("vector is not numeric")
      if(is.null(full)) 
         stop ("fully-recruited age not specified") 
      age<-age[!is.na(age)]
      st_obj<-as.data.frame(table(age))
      st_obj[,1]<-as.numeric(levels(st_obj[,1])[as.integer(st_obj[,1])])
      if(is.null(last)) last<-max(age) else last<-last 
      d<<-subset(st_obj,st_obj[,1]>=full & st_obj[,1]<=last)
        names(d)<-c("age","number")
      if(d[1,1]!=full) 
         stop ("Age specified as fully-recruited does not exist.")
      cnt<-1
     if(length(d[,1])<=2){
        print(paste("warning: only", length(d[,1]),"ages!!!"))
        rown<-length(method)*length(estimate)-1
      }
     if(length(d[,1])>2) 
        rown<-length(method)*length(estimate)  
     results<-data.frame(Method=rep("NA",rown),
                    Parameter=rep("NA",rown),
                    Estimate=rep(as.numeric(NA),rown),
                    SE=rep(as.numeric(NA),rown),stringsAsFactors=F)


    if(any(method=="lr")){ 
      if(length(d[,1])>2){ 
        extd1<-d[d[,2]>0,]
	       cc<-lm(log(extd1[,2])~extd1[,1])
         Zcc<-coef(summary(cc))[2,1]*-1
         SEZcc<-round(coef(summary(cc))[2,2],3)
         Scc<-exp(coef(summary(cc))[2,1])
         SEScc<-Scc*SEZcc

      if(any(estimate=="s")){
         results[cnt,1]<-"Linear Regression"
         results[cnt,2]<-"S"
         results[cnt,3]<-round(Scc,2)
         results[cnt,4]<-round(SEScc,3)
         cnt<-cnt+1
       }
      if(any(estimate=="z")){
         results[cnt,1]<-"Linear Regression"
         results[cnt,2]<-"Z"
         results[cnt,3]<-round(Zcc,2)
         results[cnt,4]<-round(SEZcc,3)
         cnt<-cnt+1
       }
     }
   } #lr
     
     if(any(method=="wlr")){ 
       if(length(d[,1])>2){   
         extd1<-d[d[,2]>0,] 
         LRglm_out1<-try(lm(log(extd1[,2])~extd1[,1]),silent=TRUE)
         if(!any(class(LRglm_out1)=="try-error")){
           preds<-predict(LRglm_out1)
           preds[preds<0]<-0
           preds<-preds/sum(preds)
           LRglm_out2<-try(lm(log(extd1[,2])~extd1[,1],weights=preds),silent=TRUE)
           if(!any(class(LRglm_out2)=="try-error")){
             Zmb<-round(abs(summary(LRglm_out2)$coefficients[2,1]),4)
             SEZmb<-round(summary(LRglm_out2)$coefficients[2,2],4)
             Smb<-exp(-Zmb)
             SESmb<-Smb*SEZmb
           }
         }  
         if(any(class(LRglm_out1)=="try-error")) {Zmb<-NA;SEZmb<-NA;Smb<-NA;SESmb<-NA}
         if(any(class(LRglm_out2)=="try-error")) {Zmb<-NA;SEZmb<-NA;Smb<-NA;SESmb<-NA}
         if(any(estimate=="s")){
           results[cnt,1]<-"Weighted Linear Regression"
           results[cnt,2]<-"S"
           results[cnt,3]<-round(Smb,2)
           results[cnt,4]<-round(SESmb,3)
           cnt<-cnt+1
         }
         if(any(estimate=="z")){
           results[cnt,1]<-"Weighted Linear Regression"
           results[cnt,2]<-"Z"
           results[cnt,3]<-round(Zmb,2)
           results[cnt,4]<-round(SEZmb,3)
           cnt<-cnt+1
         }
       }
     } 
     
   if(any(method=="he")){
     	Sh<-1-d[,2][d[,1]==full]/sum(d[,2])
	SESh<-sqrt(Sh*(1-Sh)/sum(d[,2]))
	Zh<--log(1-d[,2][d[,1]==full]/sum(d[,2]))
	SEZh<-sqrt(((1-Sh)^2)/(sum(d[,2])*Sh))

      if(any(estimate=="s")){
         results[cnt,1]<-"Heincke"
         results[cnt,2]<-"S"
         results[cnt,3]<-round(Sh,2)
         results[cnt,4]<-round(SESh,3)
         cnt<-cnt+1
       }
      if(any(estimate=="z")){
         results[cnt,1]<-"Heincke"
         results[cnt,2]<-"Z"
         results[cnt,3]<-round(Zh,2)
         results[cnt,4]<-round(SEZh,3)
         cnt<-cnt+1
       }
    }
   if(any(method=="cr")){
       j<-seq(0,length(d[,1])-1,1)
       Scr<-sum(j*d[,2])/(sum(d[,2])+sum(j*d[,2])-1)
       SEScr<-sqrt(Scr*(Scr-(sum(j*d[,2])-1)/(sum(d[,2])+
                 sum(j*d[,2])-2)))
 	      Zcr<--log(Scr)
	      SEZcr<-sqrt(((1-Scr)^2)/(sum(d[,2])*Scr))          
       if(any(estimate=="s")){
          results[cnt,1]<-"Chapman-Robson"
          results[cnt,2]<-"S"
          results[cnt,3]<-round(Scr,2)
          results[cnt,4]<-round(SEScr,3)
          cnt<-cnt+1
        }
       if(any(estimate=="z")){
          results[cnt,1]<-"Chapman-Robson"
          results[cnt,2]<-"Z"
          results[cnt,3]<-round(Zcr,2)
          results[cnt,4]<-round(SEZcr,3)
          cnt<-cnt+1
        }
    }
     if(any(method=="crcb")){
       j<-seq(0,length(d[,1])-1,1)
       Scr<-sum(j*d[,2])/(sum(d[,2])+sum(j*d[,2])-1)
       SEScr<-sqrt(Scr*(Scr-(sum(j*d[,2])-1)/(sum(d[,2])+
            sum(j*d[,2])-2)))
       X<-sum(j*d[,2])
       n<-sum(d[,2])
       bc<-((n-1)*(n-2))/(n*(n+X-1)*(X+1))
       Zcr<--log(Scr)
       pN<-sum(d[,2])*(1-Scr)*Scr^j
       chi<-sum((d[,2]-pN)^2/pN)
       chat<-sqrt(chi/(length(d[,2])-1))
       SEZcr<-sqrt(((1-Scr)^2)/(sum(d[,2])*Scr))          
       if(any(estimate=="s")){
         results[cnt,1]<-"Chapman-Robson CB"
         results[cnt,2]<-"S"
         results[cnt,3]<-round(Scr,2)
         results[cnt,4]<-round(SEScr*chat,3)
         cnt<-cnt+1
       }
       if(any(estimate=="z")){
         results[cnt,1]<-"Chapman-Robson CB"
         results[cnt,2]<-"Z"
         results[cnt,3]<-round(-log(Scr)-bc,2)
         results[cnt,4]<-round(SEZcr*chat,3)
         cnt<-cnt+1
       }
     }
     
     
     if(any(method=="pois")){
       max.age<-max(d[,1])
       extd<-rbind(d,cbind(age=(max.age+1):(3*max.age),number=rep(0,max.age)))
       poiglm_out1<-try(glm(number~age,family=poisson(link=log),data=extd,
                            control=list(maxit=1000)),silent=TRUE)
       if(!any(class(poiglm_out1)=="try-error")){
         if(any(estimate=="s")){
           results[cnt,1]<-"Poisson Model"
           results[cnt,2]<-"S"
           dodo<-summary(poiglm_out1)$coefficients[2,1] #Z
           pr<-as.numeric(predict(poiglm_out1,type="response"))
           X2<-(extd[,2]-pr)^2/pr
           chat<-sum(X2[which(pr>=1)])/(poiglm_out1$df.res-length(which(pr<1)))
           chat<-ifelse(chat<1,1,chat)
           dodo1<-summary(poiglm_out1)$coefficients[2,2]*sqrt(chat)
           tf<-bt.log(dodo,dodo1,sum(d[,2]))
           results[cnt,3]<-round(as.numeric(tf[1]),2)
           results[cnt,4]<-round(as.numeric(tf[4]),3)
           cnt<-cnt+1
         }
         if(any(estimate=="z")){
           results[cnt,1]<-"Poisson Model"
           results[cnt,2]<-"Z"
           results[cnt,3]<-round(abs(summary(poiglm_out1)$coefficients[2,1]),2)
           pr<-as.numeric(predict(poiglm_out1,type="response"))
           X2<-(extd[,2]-pr)^2/pr
           chat<-sum(X2[which(pr>=1)])/(poiglm_out1$df.res-length(which(pr<1)))
           chat<-ifelse(chat<1,1,chat)
           results[cnt,4]<-round(summary(poiglm_out1)$coefficients[2,2]*sqrt(chat),3)
           cnt<-cnt+1
         }
       }
     } 
     
  if(any(method=="ripois")){
      max.age<-max(d[,1])
      my.control<-glmerControl(optCtrl=list(maxfun=10000),optimizer="bobyqa")
    extd<-rbind(d,cbind(age=(max.age+1):(3*max.age),number=rep(0,max.age)))
    wer<-glmer(number~age+(1|age),family=poisson,data=extd,control=my.control)
    if(any(estimate=="s")){
      results[cnt,1]<-"Random-Intercept Poisson Model"
      results[cnt,2]<-"S"
      tf<-bt.log(summary(wer)$coefficients[2,1],summary(wer)$coefficients[2,2],sum(d[,2]))
      results[cnt,3]<-round(as.numeric(tf[1]),2)
      results[cnt,4]<-round(as.numeric(tf[4]),3)
      cnt<-cnt+1
    }
    if(any(estimate=="z")){
     results[cnt,1]<-"Random-Intercept Poisson Model"
     results[cnt,2]<-"Z"
     results[cnt,3]<-round(abs(summary(wer)$coefficients[2,1]),2)
     results[cnt,4]<-round(summary(wer)$coefficients[2,2],3) 
     cnt<-cnt+1
    }
  }
     out<-list(results,d);names(out)<-c("results","data")
   return(out)
}


