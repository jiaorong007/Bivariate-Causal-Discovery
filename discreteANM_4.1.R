# conduct permutation test to determine direction
# same version as 2.7.2
# remove constant function
# consider Fisher Exact Test
# permutation
# by Rong Jiao
# dir=1, X causes Y; dir=-1, Y causes X; dir=0, no causation; dir=2, both directions, very rare;
# treat genotype values 0,1,2 as non-cyclic, i.e. cycX=0
# If X/Y is non-cyclic, set cycX/cycY=0
library(descr)
permdANM<-function(X,Y,number_of_permutations=5000,level=0.05,cycX=1,cycY=1){
  output=fit_both_dir_discrete(X,cycX=cycX,Y,cycY=cycY,level=level)
  P_X2Y=output$p_val_fw
  P_Y2X=output$p_val_bw
  diff_estimated=P_Y2X-P_X2Y
  abs_diff=abs(diff_estimated)
  
  # P value of causal test by permutation test
  perm_X2Y<-perm_Y2X<-c()
  sizeY=length(Y)
  permutation_output=lapply(1:number_of_permutations, function(permutation){
    Y_perm=sample(Y,sizeY)
    fit_both_dir_discrete(X,cycX=cycX,Y=Y_perm,cycY=cycY,level=level)})
  matrix_perm=matrix(unlist(permutation_output),nrow=number_of_permutations,byrow=TRUE)
  ncol_mat=ncol(matrix_perm)
  perm_X2Y=matrix_perm[,ncol_mat-1]
  perm_Y2X=matrix_perm[,ncol_mat]
  
  perm_diff<-perm_Y2X-perm_X2Y
  abs_perm_diff=abs(perm_diff)
  
  Pc=sum(abs_perm_diff>=abs_diff)/number_of_permutations
  dir=ifelse((Pc>=level)|(P_X2Y<level & P_Y2X<level),0,ifelse(P_X2Y>P_Y2X,1,ifelse(P_X2Y<P_Y2X,-1,2)))
  
  p_val_ind=ifelse((length(unique(Y))==1|length(unique(X))==1),1,association_test_pvalue(Y,X))
  
  list(dir=dir,P_no_causation=Pc,P_X2Y=P_X2Y,P_Y2X=P_Y2X,p_val_ind=p_val_ind,fct_fw=output$fct_fw,fct_bw=output$fct_bw)
}

########################### source script ##################################
fit_both_dir_discrete<-function(X,cycX,Y,cycY,level){
  if (cycY==0){
    fit.output=fit_discrete(X,Y,level)
  }else if(cycY==1){
    fit.output=fit_discrete_cyclic(X,Y,level)
  }
  fct_fw=fit.output$fct
  p_val_fw=fit.output$p_val
  
  if (cycX==0){
    fit.output=fit_discrete(Y,X,level)
  }else if(cycX==1){
    fit.output=fit_discrete_cyclic(Y,X,level)
  }
  fct_bw=fit.output$fct
  p_val_bw=fit.output$p_val
  
  # options(warn=-1)
  # p_val_ind=ifelse((length(unique(Y))==1|length(unique(X))==1),1,association_test_pvalue(Y,X))
  # 
  # options(warn=0)
  
  list(fct_fw=fct_fw,fct_bw=fct_bw,p_val_fw=p_val_fw,p_val_bw=p_val_bw)
}

#####################
# cyclic
fit_discrete_cyclic<-function(X,Y,level){
  options(warn=-1)
  # parameter
  num_iter=10
  num_pos_fct=min(max(Y)-min(Y),10)
  
  # rescaling
  # X_new takes values from 1 ... X_new_max
  # Y_values are everything between Y_min and Y_max
  X_values=unique(X)
  Y_values=seq(min(Y),max(Y),by=1)
  
  if (length(X_values)==1|length(Y_values)==1){
    fct=rep(1,length(X_values))*Y_values[1]
    p_val=1
  }else{
    p<-CrossTable(c(X,rep(NA,length(Y_values))),c(Y,Y_values),prop.chisq = FALSE)$t
    
    fct=c()
    cand=list()
    for (i in 1:length(X_values)){
      b=order(p[i,])
      for (k in 1:ncol(p)){
        p[i,k]=ifelse(k==b[length(b)],p[i,k]+1,p[i,k]+1/(2*abs(k-b[length(b)])))
      }
      b=order(p[i,])
      cand[[i]]=b
      fct=c(fct,Y_values[b[length(b)]])
    }
    
    X_new=X
    for (i in 1:nrow(p)){
      X_new[X==rownames(p)[i]]=i
    }
    
    yhat=fct[X_new]
    eps=(Y-yhat)%%(max(Y)-min(Y)+1)
    p_val=ifelse((length(unique(eps))==1),1,association_test_pvalue(eps,X))
    
    # correct=TRUE as default; if correct=FALSE, completely consistant to original MATLAB scripts
    
    # remove constant fitting
    p_val=ifelse(length(unique(fct))==1,0,p_val)
    
    i=0

    while(p_val<level & i<num_iter){
      for (j_new in sample.int(length(X_values))){  # for each X, take the Y value with largest Prob
        pos_fct=list()
        p_val_comp<-p_val_comp2<-c()
        for (j in 1:(num_pos_fct+1)){
          pos_fct[[j]]=fct
          pos_fct[[j]][j_new]=Y_values[cand[[j_new]][length(cand[[j_new]])-(j-1)]]
          yhat=pos_fct[[j]][X_new]
          eps=(Y-yhat)%%(max(Y)-min(Y)+1)
          if (length(unique(eps))==1){
            p_val_comp=c(p_val_comp,1)
            p_val_comp2=c(p_val_comp2,0)
          }else{
            options(warn=-1)
            chi_sq=chisq.test(eps,X)
            options(warn=0)
            #p_val_comp=c(p_val_comp,chi_sq$p.value)
            p_val_comp=c(p_val_comp,association_test_pvalue(eps,X))
            p_val_comp2=c(p_val_comp2,chi_sq$statistic)
          }
        }
        aa=max(p_val_comp)
        j_max=which(p_val_comp==aa)
        if (aa<level){
          j_max=which(p_val_comp2==min(p_val_comp2))
        }
        fct=pos_fct[[min(j_max)]]
        yhat=fct[X_new]
        eps=(Y-yhat)%%(max(Y)-min(Y)+1)
        p_val=ifelse((length(unique(eps))==1),1,association_test_pvalue(eps,X))
      }
      i=i+1
      
      # remove constant fitting, and take the Y value with the second largest Prob or second
      # minimal statistic
      if (length(unique(fct))==1){
        if (length(j_max)>1){
          j_max=min(j_max[j_max!=min(j_max)])
        }else if(aa<level){
          second_min_p_val_comp2=min(p_val_comp2[p_val_comp2!=min(p_val_comp2)])
          j_max=which(p_val_comp2==second_min_p_val_comp2)
        }else{
          second_max_p_val_comp=max(p_val_comp[p_val_comp!=max(p_val_comp)])
          if (second_max_p_val_comp<level){
            old_j_max=j_max
            j_max=which(p_val_comp2==min(p_val_comp2))
            if (j_max==old_j_max){
              second_min_p_val_comp2=min(p_val_comp2[p_val_comp2!=min(p_val_comp2)])
              j_max=which(p_val_comp2==second_min_p_val_comp2)
            }
          }else{
            j_max=which(p_val_comp==second_max_p_val_comp)
          }
        }
        
        fct=pos_fct[[min(j_max)]]
        yhat=fct[X_new]
        eps=(Y-yhat)%%(max(Y)-min(Y)+1)
        p_val=ifelse((length(unique(eps))==1),1,association_test_pvalue(eps,X))
      }
    }
  }
  options(warn=0)
  list(fct=fct,p_val=p_val)
}

###################
# non_cyclic
fit_discrete<-function(X,Y,level){
  options(warn=-1)
  # parameter
  num_iter=10
  num_pos_fct=min(max(Y)-min(Y),20)
  
  # rescaling
  # X_new takes values from 1 ... X_new_max
  # Y_values are everything between Y_min and Y_max
  X_values=unique(X)
  Y_values=seq(min(Y),max(Y),by=1)
  
  if (length(X_values)==1|length(Y_values)==1){
    fct=rep(1,length(X_values))*Y_values[1]
    p_val=1
  }else{
    p<-CrossTable(c(X,rep(NA,length(Y_values))),c(Y,Y_values),prop.chisq = FALSE)$t
    
    fct=c()
    cand=list()
    for (i in 1:length(X_values)){
      b=order(p[i,])
      for (k in 1:ncol(p)){
        p[i,k]=ifelse(k==b[length(b)],p[i,k]+1,p[i,k]+1/(2*abs(k-b[length(b)])))
      }
      b=order(p[i,])
      cand[[i]]=b
      fct=c(fct,Y_values[b[length(b)]])
    }
    # the following script more convenient compared to MATLAB
    X_new=X
    for (i in 1:nrow(p)){
      X_new[X==rownames(p)[i]]=i
    }
    
    yhat=fct[X_new]
    eps=Y-yhat
    p_val=ifelse((length(unique(eps))==1),1,association_test_pvalue(eps,X))

    # remove constant fitting
    p_val=ifelse(length(unique(fct))==1,0,p_val)
    
    i=0
    
    while(p_val<level & i<num_iter){
      for (j_new in sample.int(length(X_values))){
        pos_fct=list()
        p_val_comp<-p_val_comp2<-c()
        for (j in 1:(num_pos_fct+1)){
          pos_fct[[j]]=fct
          pos_fct[[j]][j_new]=Y_values[cand[[j_new]][length(cand[[j_new]])-(j-1)]]
          yhat=pos_fct[[j]][X_new]
          eps=Y-yhat
          
          if (length(unique(eps))==1){
            p_val_comp=c(p_val_comp,1)
            p_val_comp2=c(p_val_comp2,0)
          }else{
            options(warn=-1)
            chi_sq=chisq.test(eps,X)
            options(warn=0)
            p_val_comp=c(p_val_comp,association_test_pvalue(eps,X))
            p_val_comp2=c(p_val_comp2,chi_sq$statistic)
          }
        }
        aa=max(p_val_comp)
        j_max=which(p_val_comp==aa)
        if (aa<level){
          j_max=which(p_val_comp2==min(p_val_comp2))
        }
        fct=pos_fct[[min(j_max)]]
        yhat=fct[X_new]
        eps=Y-yhat
        p_val=ifelse((length(unique(eps))==1),1,association_test_pvalue(eps,X))
      }
      i=i+1
      
      # remove constant fitting, and take the Y value with the second largest Prob or second
      # minimal statistic
      if (length(unique(fct))==1){
        if (length(j_max)>1){
          j_max=min(j_max[j_max!=min(j_max)])
        }else if(aa<level){
          second_min_p_val_comp2=min(p_val_comp2[p_val_comp2!=min(p_val_comp2)])
          j_max=which(p_val_comp2==second_min_p_val_comp2)
        }else{
          second_max_p_val_comp=max(p_val_comp[p_val_comp!=max(p_val_comp)])
          if (second_max_p_val_comp<level){
            old_j_max=j_max
            j_max=which(p_val_comp2==min(p_val_comp2))
            if (j_max==old_j_max){
              second_min_p_val_comp2=min(p_val_comp2[p_val_comp2!=min(p_val_comp2)])
              j_max=which(p_val_comp2==second_min_p_val_comp2)
            }
          }else{
            j_max=which(p_val_comp==second_max_p_val_comp)
          }
        }
        
        fct=pos_fct[[min(j_max)]]
        yhat=fct[X_new]
        eps=Y-yhat
        p_val=ifelse((length(unique(eps))==1),1,association_test_pvalue(eps,X))
      }
    }
    fct=fct+round(mean(eps))
  }
  options(warn=0)
  list(fct=fct,p_val=p_val)
}

association_test_pvalue<-function(x,y){
  options(warn=-1)
  chisq.association=chisq.test(x,y)
  p.value=try(fisher.test(x,y,hybrid=TRUE)$p.value,silent=TRUE)
  options(warn=0)
  p.value=ifelse(class(p.value)=="try-error",chisq.association$p.value,p.value)
  return(p.value)
}
