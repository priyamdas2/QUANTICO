rm(list=ls())
#install.packages("pracma")
#install.packages("rqPen")
#install.packages("lme4")
library(pracma)
library(rqPen)
library(parallel)
library(lme4)
setwd("Y:/QUANTICO/QUANTICO REPRODUCIBLE CODES/Sample Size 100 Codes/LASSO QR Method")

no_quantiles <- 9
num_rep <- 25
n <- 100         # Sample Size
num_of_p <- 10   # Number of Level 1 covariates : Can be set = 5,10...

f <- function(rep) {
  library(pracma)
  library(rqPen)
  library(parallel)
  library(lme4)
  print(rep)
  no_quantiles <- 9
  n <- 100 
  num_of_p <- 10



MSE <- array(0,no_quantiles)
AUC <- array(0,no_quantiles)
TP <- TPR <- FP <- FPR <- TN <- FN <- MCC <- array(0,no_quantiles)
presence <- matrix(0,no_quantiles,num_of_p)  


AUC_iters <- 501
AUC_center <- (AUC_iters+1)/2


for(i in 1:no_quantiles)
{if(i<=5)
{true_coords <- c(1,2)
null_coords <- c(3:num_of_p)
}else
{true_coords <- c(1,2,3)
null_coords <- c(4:num_of_p)}
  
  tau_here <- i*0.1
  filename <- paste("./Data/QR_P_sim_n_",n,"_seed_",rep,".csv", sep = "")
  QR_P_all <- read.csv(filename, header=FALSE, sep=",")
  
  filename <- paste("./Data/QR_Y_sim_n_",n,"_seed_",rep,".csv", sep = "")
  QR_Y <- read.csv(filename, header=FALSE, sep=",")
  
  filename <- paste("./Data/QR_Y_sim_TRUEs_n_",n,"_seed_",rep,".csv", sep = "")
  QR_y_TRUE  <- read.csv(filename, header=FALSE, sep=",")
  
  set.seed(i)
  QR_P <- QR_P_all[,2:(num_of_p+1)]
  x <- as.matrix(QR_P)
  y <- as.matrix(QR_Y)
  
  
  cv_model <- cv.rq.pen(x,y, tau = tau_here, penalty = "LASSO")
  
  lambda_opt <- cv_model$cv[which.min(cv_model$cv[,2]),1]
  
  cv_model_final <- cv.rq.pen(x,y, tau = tau_here, penalty = "LASSO",lambda = lambda_opt)
  
  aaaa <- coefficients(cv_model_final)
  
  for(ii in 1:num_of_p)
  {if(abs(aaaa[ii+1])>10^-7)
  {presence[i,ii] <- 1
  }else
  {presence[i,ii] <- 0}}
  
  TP[i] <- sum(presence[i,true_coords])
  TPR[i] <- TP[i]/length(true_coords)
  
  FP[i] <- sum(presence[i,null_coords])
  FPR[i] <- FP[i]/length(null_coords)
  
  TN[i] <- length(true_coords) - TP[i]
  FN[i] <- length(null_coords) - FP[i]
  
  MCC[i] <- (TP[i]*TN[i] - FP[i]*FN[i])/sqrt((TP[i]+FP[i])*(TP[i]+FN[i])*(TN[i]+FP[i])*(TN[i]+FN[i]))
  
  est <- as.matrix(cbind(matrix(1,n,1),QR_P))%*%as.matrix(aaaa)
  
  MSE[i] <- sum((est - QR_y_TRUE[,i])*(est - QR_y_TRUE[,i]))/n
  
  ### AUC calculation
  
  log_lambda_center <- log(lambda_opt)
  
  log_lambda_array <- array(AUC_iters,1)
  
  for(jj in 1:(AUC_center-1))
  {log_lambda_array[AUC_center-jj] <- log_lambda_center - jj*0.02
  log_lambda_array[AUC_center+jj] <- log_lambda_center + jj*0.02}
  log_lambda_array[AUC_center] <- log_lambda_center
  
  
  TP_here <- TPR_here <- array(0,AUC_iters)
  FP_here <- FPR_here <- array(0,AUC_iters)
  presence_here <- matrix(0,AUC_iters,num_of_p)
  
  for(jj in 1:AUC_iters)
  {print(c(rep,i,jj))
    cv_model_final_case <- cv.rq.pen(x,y, tau = tau_here, penalty = "LASSO",lambda = exp(log_lambda_array[jj]))
    bbbb <- coefficients(cv_model_final_case)
    
    
    
    for(ii in 1:num_of_p)
    {if(abs(bbbb[ii+1])>10^-7)
    {presence_here[jj,ii] <- 1
    }else
    {presence_here[jj,ii] <- 0}}
    
    TP_here[jj] <- sum(presence_here[jj,true_coords])
    TPR_here[jj] <- TP_here[jj]/length(true_coords)
    
    FP_here[jj] <- sum(presence_here[jj,null_coords])
    FPR_here[jj] <- FP_here[jj]/length(null_coords)}
  
  areas <- unique(cbind(as.matrix(TPR_here), as.matrix(FPR_here)))
  areas_reverse <- areas[dim(areas)[1]:1,]
  if(dim(areas)[1]>1)
  {AUC[i] <- trapz(areas_reverse[,2], areas_reverse[,1])}
  else{AUC[i] <- 1}
}

Summary <- rbind(TPR,FPR,MCC,AUC,MSE)


return(Summary)
}

no_cores <- detectCores() - 1
#cl <- makeCluster(min(no_cores,num_rep))
cl <- makeCluster(25)
save2 <- parLapply(cl, 1:num_rep,f)

Summary_all <- array(0,c(5,no_quantiles,num_rep))

for(kkkk in 1:num_rep)
{Summary_all[,,kkkk] <- save2[[kkkk]]}



mean_collection <- matrix(0,5,no_quantiles)
sd_collection <- matrix(0,5,no_quantiles)



for(iiii in 1:5)
{for(jjjj in 1:no_quantiles)
{mean_collection[iiii,jjjj] <- mean(Summary_all[iiii,jjjj,1:num_rep])
  sd_collection[iiii,jjjj] <- sd(Summary_all[iiii,jjjj,1:num_rep])}}


write.table(mean_collection, file=paste("ZZZZ_MEAN_n_",n,"_p_",num_of_p,"_rep_",num_rep,".csv",
                                    sep=""),sep=",",row.names=FALSE, col.names=FALSE)

write.table(sd_collection, file=paste("ZZZZ_SD_n_",n,"_p_",num_of_p,"_rep_",num_rep,".csv",
                                      sep=""),sep=",",row.names=FALSE, col.names=FALSE)
