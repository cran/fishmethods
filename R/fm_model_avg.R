
fm_model_avg<-function(...,global=NULL,chat=1.0){
             foo.call <- as.character(match.call())[-1] 
             all.x <- list(...)
             n<-length(all.x)
             names(all.x) <- unique(foo.call)
            
#Check model type
       check<-NULL
       chyrs<-NULL
	 for(i in 1:n){ 
             check[i]<-all.x[[i]]$type
             chyrs[i]<-length(all.x[[i]]$parameter_estimates$Occasion)
          }
        if(length(unique(check))>1) stop ("Model outputs are not of the same type ")
	  if(length(unique(chyrs))>1) stop("Model Occasions are not the same length")
        
   #Create Statistics Table
        outpt<-array(NA,dim=c(n,11))
        rownames(outpt)<-c(unique(foo.call))
        colnames(outpt)<-c("Log-Likelihood","No. Parms","AIC","AICc","N","QAIC"
                 ,"QAICc","dQAICc","e(-0.5*dQAICc)","QAICc Wgts","Global chat")
         cnt<-0
        for(x in list(...)){  
          cnt<-cnt+1
          outpt[cnt,1]<-x$statistics[1]
          outpt[cnt,2]<-x$statistics[2]
          outpt[cnt,3]<-x$statistics[3]
          outpt[cnt,4]<-x$statistics[4]
          outpt[cnt,5]<-x$statistics[5]
          outpt[cnt,6]<-((-2*outpt[cnt,1])/chat)+2*(outpt[cnt,2]+1)
          outpt[cnt,7]<-outpt[cnt,6]+(2*(outpt[cnt,2]+1)*(outpt[cnt,2]+2))/(outpt[cnt,5]-outpt[cnt,2])
         }
         outpt[,8]<-outpt[,7]-min(outpt[,7])
         outpt[,9]<-exp(-0.5*outpt[,8])
         outpt[,10]<-outpt[,9]/sum(outpt[,9])
         outpt[,11]<-chat

  #################Generate Adjusted SE          
       for(i in 1:n){
         all.x[[i]]$parameter_estimates$adjFSE<-sqrt(all.x[[i]]$parameter_estimates$FSE^2*chat)
         all.x[[i]]$parameter_estimates$adjMSE<-sqrt(all.x[[i]]$parameter_estimates$MSE^2*chat)
         all.x[[i]]$parameter_estimates$adjpSE<-sqrt(all.x[[i]]$parameter_estimates$pSE^2*chat)
        }
 #Wgt SE
	  yrlen<-unique(chyrs)
	       FF<-data.frame(Occasion=all.x[[1]]$parameter_estimates$Occasion,avgF=NA,Wgt_SE=NA,Uncond_SE=NA )
         M<-data.frame(Occasion=all.x[[1]]$parameter_estimates$Occasion,avgM=NA,Wgt_SE=NA,Uncond_SE=NA )
	       P<-data.frame(Occasion=all.x[[1]]$parameter_estimates$Occasion,avgP=NA,Wgt_SE=NA,Uncond_SE=NA )
	      tempF<-NULL;tempFWSE<-NULL;tempFUSE<-NULL
        tempM<-NULL;tempMWSE<-NULL;tempMUSE<-NULL
        tempP<-NULL;tempPWSE<-NULL;tempPUSE<-NULL
	      for(i in 1:n){
		        tempF<-cbind(tempF,all.x[[i]]$parameter_estimate$F*outpt[i,10])
		        tempFWSE<-cbind(tempFWSE,all.x[[i]]$parameter_estimate$adjFSE*outpt[i,10])
            tempM<-cbind(tempM,all.x[[i]]$parameter_estimates$M*outpt[i,10])
		        tempMWSE<-cbind(tempMWSE,all.x[[i]]$parameter_estimates$adjMSE*outpt[i,10])
            tempP<-cbind(tempP,all.x[[i]]$parameter_estimates$p*outpt[i,10])
		        tempPWSE<-cbind(tempPWSE,all.x[[i]]$parameter_estimates$adjpSE*outpt[i,10])
		      }
      
       FF$avgF<-rowSums(tempF)
       FF$Wgt_SE<-rowSums(tempFWSE)
       M$avgM<-rowSums(tempM)
       M$Wgt_SE<-rowSums(tempMWSE)
     
       P$avgP<-rowSums(tempP)
       P$Wgt_SE<-rowSums(tempPWSE)
#calculate unconditional SE
  	 for(i in 1:n){
	    tempFUSE<-cbind(tempFUSE,outpt[i,10]*(all.x[[i]]$parameter_estimates$adjFSE^2+
                 (all.x[[i]]$parameter_estimates$F-FF$avgF)^2))
   	    tempMUSE<-cbind(tempMUSE,outpt[i,10]*(all.x[[i]]$parameter_estimates$adjMSE^2+
                 (all.x[[i]]$parameter_estimates$M-M$avgM)^2))
   	   
	    tempPUSE<-cbind(tempPUSE,outpt[i,10]*(all.x[[i]]$parameter_estimates$adjpSE^2+
                 (all.x[[i]]$parameter_estimates$p-P$avgP)^2))
       
       }
 	FF$Uncond_SE<-sqrt(rowSums(tempFUSE))
  	M$Uncond_SE<-sqrt(rowSums(tempMUSE))
  	P$Uncond_SE<-sqrt(rowSums(tempPUSE))
      ans<-NULL
      ans$statistics<-outpt
	ans$model_averaged_F<-FF
	ans$model_averaged_M<-M
	ans$model_averaged_P<-P
    return(ans)  
} 

