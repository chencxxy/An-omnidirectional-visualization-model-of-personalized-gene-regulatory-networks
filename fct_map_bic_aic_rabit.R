setwd("/Users/chixiangchen/Desktop/phd in psu/research/Prof. Wu/rabits")

#install.packages("mvtnorm")
library(mvtnorm)
library(MASS)
set.seed(1324)

###############################################
#analysis 2: do log transformation for the data
###############################################

#######################################
#useful functions in functional mapping
#######################################
#calculate P_ij in the paper
P_ij<-function (y_i,x_i,beta,sigma,omega)
{
  beta_sigma<-rbind(beta,sigma)
  whole<-omega*apply(beta_sigma,2,function (x) exp(-t(y_i-x_i%*%x[1:2])%*%(y_i-x_i%*%x[1:2])/2/x[3])/(sqrt(2*pi*x[3]))^(length(y_i))  )     
  #as.matrix(whole/sum(whole),ncol=1)
  whole
}

#calcuate P matrix
P_matrix<-function (x) {
  Y<-x
  x_i<-cbind(1,X[!is.na(Y)])
  y_i<-Y[!is.na(Y)]
  p_ij_1<-P_ij(y_i,x_i,beta,sigma,omega)
  p_ij_1[p_ij_1==0]<-1e-100
  p_ij_1
}

#calculate one part in mle
fct_A<-function (x)
{
  Y<-x[-1]
  P_j_i<-x[1]
  x_i<-cbind(1,as.numeric(X[!is.na(Y)]))
  y_i<-Y[!is.na(Y)]
  t(P_j_i*t(y_i)%*%x_i)
}

#calculate another part in mle
fct_B<-function (x)
{
  Y<-x[-1]
  P_j_i<-x[1]
  x_i<-cbind(1,as.numeric(X[!is.na(Y)]))
  y_i<-Y[!is.na(Y)]
  P_j_i*t(x_i)%*%x_i
}

#calculate sigma_square
fct_sigma<-function (x)
{
  Y<-x[-1]
  P_j_i<-x[1]
  x_i<-cbind(1,as.numeric(X[!is.na(Y)]))
  y_i<-Y[!is.na(Y)]
  numerator<-P_j_i*t(y_i-x_i%*%beta_new[,j])%*%(y_i-x_i%*%beta_new[,j])
  numerator
}

#recalculate the row mean
mean_nonzero<-function (x)
{
  xx<-x[!is.na(x)]
  mean(xx)
}

#read all data
data_rabit<-read.csv("rabits_delete.csv",header = F,row.names = 1)
head(data_rabit)
dim(data_rabit)
str(data_rabit)

#read gene expression 
data_gene_exp<-read.csv("genes.csv",header = F,row.names = 1)
head(data_gene_exp)
str(data_gene_exp)
#data_gene_exp<-as.data.frame(sapply(data_gene_exp, unlist))
#data_gene_exp<-as.data.frame(sapply(data_gene_exp, as.numeric))
#head(data_gene_exp)
#str(data_gene_exp)
dim(data_gene_exp)

#read outcome
data_outcome<-unlist(data_rabit[2,])
data_outcome<-ifelse(data_outcome=="HIGH",1,0)

#create patient and time
times<-unlist(data_rabit[1,])
times<-ifelse(times=="0.08",0.08,ifelse(times=="1",1,ifelse(times=="3",3,ifelse(
  times=="7",7,ifelse(times=="14",14,ifelse(times=="30",30,ifelse(times=="90",90,180)))
))))
data_phenotype<-rbind(times,data_outcome)
head(data_phenotype)
head(data_gene_exp)
colnames(data_gene_exp)<-colnames(data_phenotype)<-rep(1:73)
data_all<-rbind(data_phenotype,data_gene_exp)
head(data_all)
dim(data_all)


#read the data
rownames(data_all)<-rep(1:nrow(data_all))
dim(data_all)

#make every value positive
data_gene_exp<-data_gene_exp-min(data_gene_exp)+0.1

#summary statistics for row means:
row_mean<-apply(data_gene_exp,1, mean)
summary(row_mean)
#par(mfrow=c(1,1))
#sum(row_mean>=7)/length(row_mean)
qqnorm(y=row_mean)

# order the expression based on the mean value
bind<-cbind(row_mean,rep(1:length(row_mean)))  #add index and sort it
order_bind<-bind[order(bind[,1],decreasing = T),] 
dim(order_bind)

#order the data by expression index:
total<-apply(data_gene_exp,2,sum)
exp_index<-total[order(total)]
X<-log(exp_index)
data_gene_exp_order<-data_gene_exp[,order(total)]
dim(data_gene_exp_order)
Y<-data_gene_exp_order[order_bind[,2],]
Y[Y!=0]<-log(Y[Y!=0])
Y[Y==0]<-NA
dim(Y)
dev.off()
plot(y=as.numeric(Y[140,]),x=X)
length(X)
length(unique(X))
summary(Y)
row_mean<-matrix(apply(Y,1, mean_nonzero),ncol=1)
dim(row_mean)
summary(row_mean)
k_cluster<-c(30,45,50,55,60,65,80)
#k_cluster<-c(80,100,120,135,145,150)
AIC<-rep()
BIC<-rep()
kk<-0

timestart<-Sys.time()
for (iii in k_cluster)
{
  #use k means to find initial:beta, sigma, omega
  cluster<-iii
  kmeans_cluster<-kmeans(row_mean,cluster)
  cluster_index<-matrix(kmeans_cluster$cluster,ncol=1)
  omega<-kmeans_cluster$size
  omega<-omega/sum(omega)
  
  
  Y_matrix<-as.matrix(Y,ncol=ncol(Y))
  Y_matrix_data<-data.frame(t(Y_matrix))
  Y_resp<-as.matrix(stack(Y_matrix_data)[,1],ncol=1)
  X_cov<-as.matrix(rep(X,time=nrow(Y)))
  Cluster_index<-as.matrix(rep(cluster_index,each=ncol(Y)),ncol=1)
  data_initial<-cbind(Y_resp,X_cov,Cluster_index)
  
  coefficients<-rep()
  red_var<-rep()
  for (i in unique(cluster_index))
  {
    fit<-lm(Y_resp[Cluster_index %in% i,]~X_cov[Cluster_index %in% i,])
    red_var<-c(red_var,var(fit$residuals))
    coefficients<-cbind(coefficients,fit$coefficients)
  }
  plot(red_var)
  summary(red_var)
  plot(coefficients[1,],coefficients[2,])
  summary(coefficients[1,])
  summary(coefficients[2,])
  total_number<-nrow(Y)
  beta<-coefficients
  sigma<-red_var
  
  #initial
  total<-0
  sum_t<-apply(Y[,],1,function (x) length(x[!is.na(x)]))
  unique(sum_t)
  
  ################
  #EM algorithm
  ################
  #update omega_new
  
  repeat
  {
    P_1<-apply(Y,1,P_matrix)
    de_sum<-apply(P_1,2,sum)
    P_1_comb<-rbind(de_sum,P_1)
    P<-apply(P_1_comb,2,function (x) {
      deno<-x[1]  
      del<-x[-1]
      del/deno
    }
    ) 
    omega_new<-as.vector(apply(P,1,mean))
    omega_old<-omega
    omega<-omega_new
    
    #update beta_new
    beta_new<-rep()
    for (j in 1:nrow(P))
    {
      new_data<-cbind(P[j,],Y)
      A<-apply(new_data,1,fct_A)
      B<-apply(new_data,1,fct_B)
      A_sum<-apply(A,1,sum)
      B_sum<-matrix(apply(B,1,sum),ncol=2)
      beta_new<-cbind(beta_new,ginv(B_sum)%*%A_sum)
    }
    
    #update sigma
    #numerator_sum<-rep(0,length=nrow(P))
    sigma_new<-rep()
    for (j in 1:nrow(P))
    {
      new_data<-cbind(P[j,],Y)
      numerator<-apply(new_data,1,fct_sigma)
      sigma_update<-sum(numerator)/(sum(sum_t*P[j,]))
      sigma_new<-c(sigma_new,sigma_update)
    }
    
    beta_old<-beta
    beta<-beta_new
    sigma_old<-sigma
    sigma<-sigma_new
    
    check_beta<-mean(abs(beta-beta_old))/ncol(beta)
    check_sigma<-mean(abs(sigma-sigma_old))
    check_omega<-mean(abs(omega-omega_old))
    total<-total+1
    print(total)
    if (check_beta<0.01 & check_sigma<0.01 & check_omega<0.001)
    {
      index_z<-matrix(apply(P,2,function (x) ifelse(x==max(x),1,0)),nrow(P))
      #P_1[P_1==0]<-9.681414e-320
      likelihood<-sum(log((index_z*P_1)[index_z*P_1>0]))
      break
    }
  }
  
  select_index<-which(apply(index_z,1,sum)!=0)
  length(select_index)
  kk<-kk+1
  p<-(2+1+1)*cluster
  AIC[kk]<--2*likelihood+2*p
  BIC[kk]<--2*likelihood+p*log(total_number*ncol(Y))
  
  print(kk)
}

timeend<-Sys.time()
timeend-timestart
plot(AIC)
plot(BIC)
#plot each cluster curve in one plot
par(mfrow=c(1,1))
fitted<-exp(cbind(1,X) %*% beta[,1])
plot(x=exp_index,y=fitted,type="l",ylim = c(0,10))
cl<-rainbow(cluster)
for (i in select_index)
{
  fitted<-exp(cbind(1,X) %*% beta[,i])
  lines(x=exp_index,y=fitted,col=cl[i])
}

plot(x=exp_index,y=data_gene_exp_order[order_bind[,2],][1250,],type="l",ylim = c(0,10))

