# use gamma approximation based test rather than permutation test for HSIC association.
library(assist)
library(psych)
library(dHSIC)
library(fANCOVA)
library(entropy)

# conduct permutation test to determine direction
permANM<-function(X,Y,number_of_permutations=5000,level=0.05,
                  fit.method=c("ssr","B.spline","loess"),
                  measurement=c("difference","ratio"),score=c("HSIC","Entropy-empirical")){

  if (length(fit.method)==3){
    fit.method="ssr"
  }
  if (length(measurement)==2){
    measurement="difference"
  }
  if (length(score)==2){
    score="HSIC"
  }
  fit.output=fit_both_dir_continuous(X,Y,fit.method=fit.method,score=score,level=level)
  P_X2Y=fit.output$p_val_fw
  P_Y2X=fit.output$p_val_bw
  constant_fit_X2Y=fit.output$constant_fit_fw
  constant_fit_Y2X=fit.output$constant_fit_bw
  
  if (measurement=="difference"){
    diff_obs=P_Y2X-P_X2Y
    obs_measure=abs(diff_obs)
  }else{
    ratio_obs=P_Y2X/P_X2Y
    obs_measure=ifelse(ratio_obs>1,ratio_obs,1/ratio_obs)
  }
  
  # P value of causation test by permutation test based on difference
  perm_X2Y<-perm_Y2X<-c()
  sizeY=length(Y)
  
  permutation_output=lapply(1:number_of_permutations, function(permutation){
    Y_perm=sample(Y,sizeY)
    fit_both_dir_continuous(X,Y_perm,fit.method=fit.method,score=score,level=level)})
  matrix_perm=matrix(unlist(permutation_output),nrow=number_of_permutations,byrow=TRUE)
  ncol_mat=ncol(matrix_perm)
  perm_X2Y=matrix_perm[,ncol_mat-4]
  perm_Y2X=matrix_perm[,ncol_mat-3]
  
  if (measurement=="difference"){
    diff_perm<-perm_Y2X-perm_X2Y
    perm_measures=abs(diff_perm)
  }else{
    perm_measures<-perm_Y2X/perm_X2Y
    perm_measures[perm_measures<1]=(1/perm_measures)[perm_measures<1]
  }
  
  Pc=sum(perm_measures>=obs_measure)/number_of_permutations
  dir=ifelse((Pc>=level)|(P_X2Y<level & P_Y2X<level),0,ifelse(P_X2Y>P_Y2X,1,ifelse(P_X2Y<P_Y2X,-1,2)))
  dir=ifelse(dir==2,dir-3*constant_fit_X2Y-1*constant_fit_Y2X+2*constant_fit_X2Y*constant_fit_Y2X,
             ifelse(dir==1,dir-constant_fit_X2Y,ifelse(dir==-1,dir+constant_fit_Y2X,dir)))
  
  list(dir=dir,P_X2Y=P_X2Y,P_Y2X=P_Y2X,P_no_causation=Pc,P_ind=fit.output$p_val_ind,
       constant_fit_X2Y=constant_fit_X2Y,constant_fit_Y2X=constant_fit_Y2X)
}

# fit both direction for contiANM
fit_both_dir_continuous<-function(x,y,m=c(2,3,4),N=40,limnla=c(-10,3),
                                  fit.method,score,level){
  if (length(m)==3){
    m=2
  }
  fit.output=contiANM(x,y,m=m,N=N,limnla=limnla,fit.method=fit.method,score=score,level=level)
  statistic_fw=fit.output$statistic
  p_val_fw=fit.output$pvalue
  constant_fit_fw=fit.output$constant_fit
  fit.output=contiANM(y,x,m=m,N=N,limnla=limnla,fit.method=fit.method,score=score,level=level)
  statistic_bw=fit.output$statistic
  p_val_bw=fit.output$pvalue
  constant_fit_bw=fit.output$constant_fit
  
  p_val_ind=dhsic.test(x,y,method="gamma")$p.value
  
  list(statistic_fw=statistic_fw,statistic_bw=statistic_bw,p_val_fw=p_val_fw,p_val_bw=p_val_bw,
       constant_fit_fw=constant_fit_fw,constant_fit_bw=constant_fit_bw,p_val_ind=p_val_ind)
}

# continuous ANM
contiANM<-function(x,y,m,N,limnla,fit.method,score,level){
  # m: m=2, cubic spline
  #    m=3, quintic spline
  #    m=4, septic spline
  
  # limnla: a vector of length one or two, specifying a 
  #         search range for log10(n*lambda), where lambda is the 
  #         smoothing parameter and n is the sample size. If it is 
  #         a single value, the smoothing parameter will be fixed at 
  #         this value. This option is only applicable to spline 
  #         smoothing with a single smoothing parameter.
  # 
  # require(assist)
  # require(psych)
  # require(dHSIC)
  # require(fANCOVA)
  # require(entropy)
  
  n.sub=length(x)
  if (n.sub!=length(y)){
    stop("lengths of x and y do not match")
  }else{
    test.index=seq(1, n.sub, 5)
    if (fit.method=="ssr"){
      x=x-min(x)
      y=y-min(y)
    }
    x.train=x[-test.index]
    y.train=y[-test.index]
    x.test=x[test.index]
    y.test=y[test.index]
    m.test=length(x.test)
    
    ##### B spline #####
    if (fit.method=="B.spline"){
      BS=smooth.spline(x.train,y.train,nknots=N)
      fitted <- predict(BS,x.test)$y
      eps=y.test-fitted
      
      new.x.train=x.train[order(x.train)]
      fitted.new.train=predict(BS,new.x.train)$y
    }
    
    ##### smoothing splines #####
    if (fit.method=="ssr"){
      if (m==2){
        B<-ssr(y.train~x.train,rk=cubic2(x.train),limnla=limnla) # based on classic polynomials
      }else if(m==3){
        B<-ssr(y.train~x.train+I(x.train^2),rk=quintic2(x.train),limnla=limnla)
      }else if(m==4){
        B<-ssr(y.train~x.train+I(x.train^2)+I(x.train^3),rk=septic2(x.train),limnla=limnla)
      }
      fitted=predict(B,data.frame(x.train=x.test),pstd=FALSE)$fit
      eps=y.test-fitted
      
      new.x.train=x.train[order(x.train)]
      fitted.new.train=B$fit[order(x.train)]
    }
    
    ##### LOESS #####
    if (fit.method=="loess"){
      new.x.train=x.train[order(x.train)]
      new.y.train=y.train[order(x.train)]
      # training span parameter based on AICC
      span=(loess.as(new.y.train,new.x.train)$pars)$span
      cars.lo <- loess(new.y.train ~ new.x.train,span=span)
      fitted <- predict(cars.lo,x.test)
      eps=y.test-fitted
      
      fitted.new.train=cars.lo$fitted
    }
    
    ########################
    options(warn=-1)
    p_constant=summary(lm(fitted.new.train~new.x.train))$coefficients[2,4]
    options(warn=0)
    constant_fit=ifelse(p_constant>level,1,0)
    
    x.test <- x.test[!is.na(eps)]
    eps <- eps[!is.na(eps)]
    
    # two scores
    if (score=="HSIC"){
      pvalue=dhsic.test(x.test,eps,method="gamma")$p.value
      #pvalue
      statistic=dhsic.test(x.test,eps,method="gamma")$statistic
      #statistic
    }else if(score=="Entropy-empirical"){
      statistic=EdS_entropy(x.test)+EdS_entropy(eps)
      pvalue=1/statistic
    }
    list(statistic=statistic,pvalue=pvalue,constant_fit=constant_fit)
  }
}

EdS_entropy<-function(x){
  x=sort(x)
  N=length(x)
  return(mean(sapply(1:(N-1), function(i){
    dx=x[i+1]-x[i]
    dx=ifelse(dx!=0,log(abs(dx)),0)}))-digamma(1)+digamma(N))
}
