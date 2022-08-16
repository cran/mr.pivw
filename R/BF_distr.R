
BF_distr=function(data,beta_hat,tau2,lambda,n_boot=1000,seed_boot=1){
  i=0
  set.seed(seed_boot)
  z_b=numeric()

  while(i < n_boot){
    gamma= data[,2]
    Gamma_se  = data[,3]
    gamma_se  = data[,4]
    alpha = rnorm(length(gamma),0, sqrt(tau2))
    Gamma_hat = rnorm(length(gamma),gamma*beta_hat+alpha, Gamma_se)
    gamma_hat = rnorm(length(gamma),gamma, gamma_se)

    v1 = sum(Gamma_se^(-4)*(Gamma_hat^2*gamma_hat^2-(Gamma_hat^2-Gamma_se^2-tau2)*(gamma_hat^2-gamma_se^2)))
    v2 = sum(Gamma_se^(-4)*(4*(gamma_hat^2-gamma_se^2)*gamma_se^2 + 2*gamma_se^4))
    v12= sum(2*Gamma_se^(-4)*Gamma_hat*gamma_hat*gamma_se^2)
    mu1= sum(Gamma_se^(-2)*Gamma_hat*gamma_hat)
    mu2= sum(Gamma_se^(-2)*(gamma_hat^2-gamma_se^2))
    mu2p = mu2/2+sign(mu2)*sqrt(mu2^2/4+lambda*v2)
    mu1p = mu1+(v12/v2)*(mu2p-mu2)
    w = mu2p/(2*mu2p-mu2)
    v1p = v1
    v2p = w^2*v2
    v12p = w*v12

    temp=(mu1p-beta_hat*mu2p)^2/(v1p-2*beta_hat*v12p+beta_hat^2*v2p)
    if(temp<0){
      next
    }else{
      z_b=c(temp,z_b)
      i=i+1
    }
  }
  z_b=sort(z_b)
  return(z_b)
}





































