Diff_mat <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  X_expanded <- array(rep(X, each = n), dim = c(n, n, p))
  X_transposed <- aperm(X_expanded, c(2, 1, 3))
  difference_matrix <- X_transposed - X_expanded 
  
  return(difference_matrix)
}

Gramma_Mat<-function(X,Kernel){
  Dis_Mat<- abs(Diff_mat(X))
  if (Kernel=='Euclid'){
    X_gram=-Dis_Mat
  }
  if (Kernel=='Gauss'){
    #sigmas <- apply(Dis_Mat^2, 3, median)
    #Dis_Mat <- sweep(Dis_Mat, MARGIN=3, STATS=sqrt(0.5*sigmas), FUN="/")
    sigma <- sqrt(0.5*median(dist(X)^2))
    #X_gram<-array(0,c(n,n,p))
    #for(i in 1:p){
    #  X_gram[,,i]<-exp(-(Dis_Mat[,,i])^2/sigmas[i] )
    #}
    X_gram=dnorm(Dis_Mat, mean=0, sd=sigma)
  }
  if (Kernel=='Laplace'){
    sigma = sqrt(median(Dis_Mat^2))
    X_gram = exp(-Dis_Mat/sigma)
  }
  return(X_gram)
}



T_n<-function(X,Y,D,kernel){
  n <- nrow(X)
  p <- ncol(X)
  D_Mat<-array(rep(D, each = n), dim = c(n, n))
  Y_Mat<-array(rep(Y, each = n), dim = c(n, n))
  Risk_Mat<-Y_Mat<=t(Y_Mat)
  Risk_Mat_T<-t(Risk_Mat)
  Psi<-t(D_Mat)*Risk_Mat_T-D_Mat*(Risk_Mat)
  Bar_Psi1<-rowMeans(Psi)
  Bar_Psi1_M<-t(array(rep(Bar_Psi1, each = n), dim = c(n, n)))
  Bar_Psi12_M<-Bar_Psi1_M*t(Bar_Psi1_M)
  Bar_Psi1_Psi_M<-Bar_Psi1_M*Psi
  #Bar_PsiPsi_M<-rep(0,n^3)
  #dim(Bar_PsiPsi_M)<-c(n,n,n)
  #for(i in 1:n){
  #  Bar_PsiPsi_M[i,,]=t((D_Mat[,i]*(Y>=Y[i])-D_Mat[i,]*(Y[i]<=Y))*
  #                        (D_Mat*(t(Y_Mat)>=Y_Mat)-t(D_Mat)*(Y_Mat>=t(Y_Mat))))
  #}
  # Cal Bar_PsiPsi_M
  Bar_PsiPsi_M <- array(NA, dim = c(n, n, n))
  for (k in 1:n) {
    Bar_PsiPsi_M[,,k] <- Psi[,k] %o% Psi[,k]
  }
  
  ###Estimation of T_n###
  Bar_PsiPsi_M1<-apply(Bar_PsiPsi_M, c(1, 2), mean)
  K_X_Gram<-Gramma_Mat(X,kernel)
  Temp1<-sapply(1:p, function(l) 
    sum(K_X_Gram[,,l]*(n^2*Bar_Psi12_M+2*n*Bar_Psi1_Psi_M-n*Bar_PsiPsi_M1-Psi^2)))
  K_ii<-array(apply(K_X_Gram, 3, diag), dim = c(n, p))
  
  #K_ii<-array(rep(dnorm(0, mean=0, sd=sqrt(0.5*sigmas)), each = n), dim = c(n, p))
  Temp2<-as.vector((n*rowMeans(Psi^2)-n^2*Bar_Psi1^2)%*%K_ii)
  Temp<-Temp1+Temp2
  n_4<-(n*(n-1)*(n-2)*(n-3))
  T_np<-sum(Temp)/n_4
  #T_np
  ###Estimation of S^2###
  mean_dim1 <- apply(K_X_Gram, c(2, 3), mean)
  mean_dim3 <- apply(K_X_Gram, 3, mean)
  mean_dim1_expanded <- array(rep(mean_dim1, each = n), dim = c(n, n, p))
  mean_dim2_expanded<-aperm(mean_dim1_expanded, c(2,1,3))
  mean_dim3_expanded <- array(rep(mean_dim3, each = n * n), dim = c(n, n, p))
  L1_s<-K_X_Gram - mean_dim1_expanded - mean_dim2_expanded + mean_dim3_expanded
  L1_pn<-(apply(L1_s, c(1,2), sum))^2
  L2_n<- Bar_Psi12_M^2
  c_n<-1#((n-1)^2)/(n-3)^2
  S2_np<-c_n*sum(L1_pn*L2_n-diag(diag(L1_pn*L2_n)))/(n*(n-1))
  
  
  ###Scale###
  n_2<-(n*(n-1)/2)
  Z_alpha<-sqrt(n_2)*T_np/sqrt(S2_np)
  list(T_np=T_np,S2_np=S2_np,Z_alpha=Z_alpha)
}


#wild_boot_logrank<-function(X,Y,D,kernel,num_bootstrap=999){
#  n <- nrow(X)
#  p <- ncol(X)
#  D_Mat<-array(rep(D, each = n), dim = c(n, n))
#  Y_Mat<-array(rep(Y, each = n), dim = c(n, n))
#  Risk_Mat<-Y_Mat<=t(Y_Mat)
#  Risk_Mat_T<-t(Risk_Mat)
#  Psi<-t(D_Mat)*Risk_Mat_T-D_Mat*(Risk_Mat)
#  Bar_Psi1<-rowMeans(Psi)
#  Bar_Psi1_M<-t(array(rep(Bar_Psi1, each = n), dim = c(n, n)))
#
#  Bar_Psi1_Psi_M<-Bar_Psi1_M*Psi
#  Bar_PsiPsi_M <- array(NA, dim = c(n, n, n))
#  for (k in 1:n) {
#    Bar_PsiPsi_M[,,k] <- Psi[,k] %o% Psi[,k]
#  }
#  
#  ###Estimation of T_n###
#  Bar_PsiPsi_M1<-apply(Bar_PsiPsi_M, c(1, 2), mean)
#  K_X_Gram<-Gram_K(X,kernel)
#  K_ii<-diag(K_X_Gram)
#  
#  Temp1<-(n^2*Bar_Psi12_M+2*n*Bar_Psi1_Psi_M-n*Bar_PsiPsi_M1-Psi^2)
#  Temp2<-(Psi^2-n*Bar_Psi1_M*Psi)

#  Temp<-K_X_Gram*Temp1+Temp2*K_ii
#  n_4<-(n*(n-1)*(n-2)*(n-3))
#  T_np<-sum(Temp)/n_4
  
#  #samp <- lapply(1:num_bootstrap, function(x) sample(1:n))
#  #statistic_list <- sapply(1:num_bootstrap,
#  #                         function(x) {
#  #                           distX_samp <- K_X_Gram[samp[[x]],samp[[x]]]
#  #                           K_ii_samp<-diag(distX_samp)
#  #                           Temp3<-distX_samp*Temp1+Temp2*K_ii_samp
#  #                           res<- sum(Temp_3)/n_4
#  #                           return(res)
#  #               })
#  
#  statistic_list <- rep(0,num_bootstrap)
#  for (b in 1:num_bootstrap){
#    W <- rbinom(n, 1, 1/2)*2-1
#    WM <- outer(W, W)
#    diag(WM)<-0
#    bootstrap_statistic <- sum(WM * Temp)/n_4
#    statistic_list[b] <- bootstrap_statistic
#  }
#  pval <- (1 + length(which(statistic_list>T_np))) / (1 + num_bootstrap)
  
  
#  list(KLR=T_np,pvalue=pval)
#}

T_n_bootstrap<-function(X,Y,D,kernel,num_bootstrap){
  n <- nrow(X)
  p <- ncol(X)
  D_Mat<-array(rep(D, each = n), dim = c(n, n))
  Y_Mat<-array(rep(Y, each = n), dim = c(n, n))
  Risk_Mat<-Y_Mat<=t(Y_Mat)
  Risk_Mat_T<-t(Risk_Mat)
  Psi<-t(D_Mat)*Risk_Mat_T-D_Mat*(Risk_Mat)
  Bar_Psi1<-rowMeans(Psi)
  Bar_Psi1_M<-t(array(rep(Bar_Psi1, each = n), dim = c(n, n)))
  Bar_Psi12_M<-Bar_Psi1_M*t(Bar_Psi1_M)
  Bar_Psi1_Psi_M<-Bar_Psi1_M*Psi
  #Bar_PsiPsi_M<-rep(0,n^3)
  #dim(Bar_PsiPsi_M)<-c(n,n,n)
  #for(i in 1:n){
  #  Bar_PsiPsi_M[i,,]=t((D_Mat[,i]*(Y>=Y[i])-D_Mat[i,]*(Y[i]<=Y))*
  #                        (D_Mat*(t(Y_Mat)>=Y_Mat)-t(D_Mat)*(Y_Mat>=t(Y_Mat))))
  #}
  # Cal Bar_PsiPsi_M
  Bar_PsiPsi_M <- array(NA, dim = c(n, n, n))
  for (k in 1:n) {
    Bar_PsiPsi_M[,,k] <- Psi[,k] %o% Psi[,k]
  }
  
  ###Estimation of T_n###
  Bar_PsiPsi_M1<-apply(Bar_PsiPsi_M, c(1, 2), mean)
  K_X_Gram<-Gram_K(X,kernel)
  Temp1<-K_X_Gram*(n^2*Bar_Psi12_M+2*n*Bar_Psi1_Psi_M-
                          n*Bar_PsiPsi_M1-Psi^2)
  dim(Temp1)<-c(n,n,p)
  sum_Temp1_dim3 <- apply(Temp1, c(1, 2), sum)
  sum(sum_Temp1_dim3)
  
  K_ii<-array(apply(K_X_Gram, 3, diag), dim = c(n, p))
  rowS_K_ii<-rowSums(K_ii)
  
  
  sum_Temp2_dim3<-(Psi^2-n*Bar_Psi1_M*Psi)*rowS_K_ii
  sum(sum_Temp2_dim3)
  
  Temp<-sum_Temp1_dim3+sum_Temp2_dim3
  n_4<-(n*(n-1)*(n-2)*(n-3))
  T_np<-sum(Temp)/n_4

  statistic_list <- rep(0,num_bootstrap)
  for (b in 1:num_bootstrap){
    W <- rnorm(n)
    WM <- outer(W, W)
    diag(WM)<-0
    bootstrap_statistic <- sum(WM * Temp)/n_4
    statistic_list[b] <- bootstrap_statistic
  }
  pval <- (1 + length(which(statistic_list>T_np))) / (1 + num_bootstrap)
  
  S2_np<-var(c(statistic_list,T_np))
  Z_alpha<-T_np/sqrt(S2_np)
  
  list(T_np=T_np,S2_np=S2_np,Z_alpha=pval)
}




T_new<-function(X,Y,D){
  n <- nrow(X)
  p <- ncol(X)
  D_Mat<-array(rep(D, each = n), dim = c(n, n))
  Y_Mat<-array(rep(Y, each = n), dim = c(n, n))
  Risk_Mat<-Y_Mat<=t(Y_Mat)
  Risk_Mat_T<-t(Risk_Mat)
  Psi<-t(D_Mat)*Risk_Mat_T-D_Mat*(Risk_Mat)
 
  K_X_Gram<-Gramma_Mat(X,'Gauss')
  Temp_new <- rep(0,p)
  for(s in 1:p){
    Temp_3<-0
    count<-0
    for(i in 1:n){
      for(j in 1:n){
        if (j != i) {
          for(k in 1:n){
            if (k != i && k != j) {
              for(l in 1:n){
                if (l != i && l != j && l != k) {
                  Temp_3 <- Temp_3 + K_X_Gram[i,j,s]*Psi[i,k]*Psi[j,l]
                  count<-count+1
                  }
                }
              }
            }
          }
        }
    }
    Temp_new[s]<-Temp_3
  }
  
  n_4<-(n*(n-1)*(n-2)*(n-3))
  T_n<-sum(Temp_new)/n_4
  T_n
  
  
  Temp5<-0
  Temp4<-0
  for(s in 1:p){
    for(i in 1:n){
      for(j in 1:n){
        Temp5<-Temp5+K_X_Gram[i,i,s]*((Psi[i,j])^2-n*Bar_Psi1[i]*Psi[i,j])
        Temp4<-Temp4+K_X_Gram[i,j,s]*(n^2*Bar_Psi1[i]*Bar_Psi1[j]+2*n*Bar_Psi1[i]*Psi[i,j]-n*Bar_PsiPsi_M1[i,j]-Psi[i,j]^2)
      }
    }
  }
  
  
  
}  
T_new<-function(X,Y,D){
  n <- nrow(X)
  p <- ncol(X)
  D_Mat<-array(rep(D, each = n), dim = c(n, n))
  Y_Mat<-array(rep(Y, each = n), dim = c(n, n))
  Risk_Mat<-Y_Mat<=t(Y_Mat)
  Risk_Mat_T<-t(Risk_Mat)
  Psi<-t(D_Mat)*Risk_Mat_T-D_Mat*(Risk_Mat)
  
  K_X_Gram<-Gram_K(X,kernel)
    Temp_3<-0
    count<-0
    for(i in 1:n){
      for(j in 1:n){
        if (j != i) {
          for(k in 1:n){
            if (k != i && k != j) {
              for(l in 1:n){
                if (l != i && l != j && l != k) {
                  Temp_3 <- Temp_3 + K_X_Gram[i,j]*Psi[i,k]*Psi[j,l]
                  count<-count+1
                }
              }
            }
          }
        }
      }
    }
  n_4<-(n*(n-1)*(n-2)*(n-3))
  T_n<-sum(Temp_3)/n_4
  T_n
  
  
  Temp5<-0
  Temp4<-0
    for(i in 1:n){
      for(j in 1:n){
        Temp5<-Temp5+K_X_Gram[i,i]*((Psi[i,j])^2-n*Bar_Psi1[i]*Psi[i,j])
        Temp4<-Temp4+K_X_Gram[i,j]*(n^2*Bar_Psi1[i]*Bar_Psi1[j]+2*n*Bar_Psi1[i]*Psi[i,j]-n*Bar_PsiPsi_M1[i,j]-Psi[i,j]^2)
      }
    }
}  






