# TACC
args = commandArgs(trailingOnly=TRUE)
sim_N1 = as.numeric(args[1])
sim_N2 = as.numeric(args[2])
###########################################
library(parallel)
setwd("/work/04956/rjiao/stampede2/Frontiers/first_author/simulation/simulation2_power")
source("contiANM_2.1.R")
Total_N=100*1000
N_sim=500
#Nsample_vec=c(200,500,1000,2000,5000)
Nsample_vec=c(8000)
level_vec=c(.05,.01)
len_Nsample=length(Nsample_vec)

N=3

n_core=sim_N2-sim_N1+1
cl = makeCluster(getOption("cl.cores", n_core))
clusterExport(cl,varlist=ls(),envir=environment()) # make data available to different cores
temp_rlt = mclapply(sim_N1:sim_N2,mc.cores=n_core,function(j){
  set.seed(j)
  X1=rnorm(Total_N)
  #X2=rnorm(Total_N)
  #X=cbind(X1,X2)
  X1n=rnorm(N)
  #X2n=rnorm(N)
  r=rnorm(1)
  w <- runif(N)
  w <- w / sum(w)
  output_list=lapply(1:N,function(i){
    #w[i]*exp(-r*((X1-X1n[i])^2+(X2-X2n[i])^2))
    w[i]*exp(-r*((X1-X1n[i])^2))
  })
  output_matrix=matrix(unlist(output_list),nrow=Total_N,byrow=FALSE)
  # change effect of noise
  Y=apply(output_matrix,1,sum)+rnorm(Total_N,sd=.1)
  #Y=apply(output_matrix,1,sum)+rnorm(Total_N)

  #sink(paste("output/power_sim_",j,".txt",sep=""))
  sink(paste("output8000/power_sim_",j,".txt",sep=""))
  cat(c("Nsamples","level","dir","P_X2Y","P_Y2X","P_no_causation","P_ind"),"\n")
  for (Nsamples in Nsample_vec){
      #index=sample.int(Nsamples)
      index=sample(1:Total_N,Nsamples)
      subX=X1[index]
      subY=Y[index]
      for (level in level_vec){
        N_sim=ifelse(level==0.05,100,500)
        output=try(permANM(X=subX,Y=subY,number_of_permutations=N_sim,level=level))
        if (class(output)!="try-error"){
          text=c(output$dir,output$P_X2Y,output$P_Y2X,output$P_no_causation,output$P_ind)
          cat(c(Nsamples,level,text),"\n")
        }
      }
    }
  sink()
})
stopCluster(cl)
