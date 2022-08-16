mr.pivw <- function(data,lambda=1,plei=TRUE,sel.pval=NULL,delta=0,Boot.Fieller=TRUE,sig=0.05) {

  if(delta!=0 & plei!=TRUE){
    sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
    data = data[sel,]
  }

  Gamma_hat = data[,1]
  gamma_hat = data[,2]
  Gamma_se  = data[,3]
  gamma_se  = data[,4]

  v2 = sum(Gamma_se^(-4)*(4*(gamma_hat^2-gamma_se^2)*gamma_se^2 + 2*gamma_se^4))
  v12= sum(2*Gamma_se^(-4)*Gamma_hat*gamma_hat*gamma_se^2)
  mu1= sum(Gamma_se^(-2)*Gamma_hat*gamma_hat)
  mu2= sum(Gamma_se^(-2)*(gamma_hat^2-gamma_se^2))
  mu2p = mu2/2 + sign(mu2)*sqrt(mu2^2/4+lambda*v2)
  mu1p = mu1 + (v12/v2)*(mu2p-mu2)
  beta_hat = mu1p/mu2p

  if(plei==TRUE){
    tau2 = sum(((Gamma_hat-beta_hat*gamma_hat)^2 - Gamma_se^2 -beta_hat^2*gamma_se^2) * Gamma_se^(-2))/sum(Gamma_se^(-2))
    if (tau2 < 0){tau2 = 0}
  }else{
    tau2 = 0
  }

  if(delta!=0 & plei==TRUE){
    sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
    data = data[sel,]
    Gamma_hat = data[,1]
    gamma_hat = data[,2]
    Gamma_se  = data[,3]
    gamma_se  = data[,4]

    v2 = sum(Gamma_se^(-4)*(4*(gamma_hat^2-gamma_se^2)*gamma_se^2 + 2*gamma_se^4))
    v12= sum(2*Gamma_se^(-4)*Gamma_hat*gamma_hat*gamma_se^2)
    mu1= sum(Gamma_se^(-2)*Gamma_hat*gamma_hat)
    mu2= sum(Gamma_se^(-2)*(gamma_hat^2-gamma_se^2))
    mu2p = mu2/2 + sign(mu2)*sqrt(mu2^2/4+lambda*v2)
    mu1p = mu1 + (v12/v2)*(mu2p-mu2)
    beta_hat = mu1p/mu2p
  }
  beta_se = 1/abs(mu2p)*sqrt(sum((gamma_hat^2/Gamma_se^2)*(1+tau2*Gamma_se^(-2))
                                 + beta_hat^2*(gamma_se^2/Gamma_se^4)*(gamma_hat^2+gamma_se^2)))
  pval = 2*pnorm(abs(beta_hat),0,beta_se,lower.tail = FALSE)
  CI_ll = beta_hat+qnorm(sig/2)*beta_se
  CI_ul = beta_hat-qnorm(sig/2)*beta_se
  CI = cbind(CI_ll,CI_ul)
  colnames(CI)=c("lower bound","upper bound")

  if(Boot.Fieller==TRUE){

    z_b=BF_distr(data,beta_hat,tau2,lambda)
    v1p = sum(Gamma_se^(-4)*(Gamma_hat^2*gamma_hat^2-(Gamma_hat^2-Gamma_se^2-tau2)*(gamma_hat^2-gamma_se^2)))
    w = mu2p/(2*mu2p-mu2)
    v2p = w^2*v2
    v12p = w*v12
    z0 = mu1p^2/v1p
    pval_BF=sum(z0<z_b)/length(z_b)

    qt=z_b[round(length(z_b)*(1-sig))]
    A = mu2p^2-qt*v2p
    B = 2*(qt*v12p - mu1p*mu2p)
    C = mu1p^2 -qt*v1p
    D = B^2-4*A*C

    if(A>0){
      r1=(-B-sqrt(D))/2/A
      r2=(-B+sqrt(D))/2/A
      CI_BF = cbind(r1,r2)
    }else if (D>0){
      r1=(-B-sqrt(D))/2/A
      r2=(-B+sqrt(D))/2/A
      CI1=c(-Inf,r1)
      CI2=c(r2,Inf)
      CI_BF = rbind(CI1,CI2)
    }else{
      CI_BF=cbind(-Inf,Inf)
    }

    colnames(CI_BF)=c("lower bound","upper bound")
  }

  out=list(beta_hat,beta_se,pval,CI)
  names(out)=c("beta","se","pval (Normal)","CI (Normal)")
  if(Boot.Fieller==TRUE){
    out[[5]]=pval_BF
    out[[6]]=CI_BF
    names(out)[5]="pval (Bootstrap Fieller)"
    if(nrow(CI_BF)==2){names(out)[6]="CI (Bootstrap Fieller): the union of CI1 and CI2"
    }else{names(out)[6]="CI (Bootstrap Fieller)"}
  }
  if(plei==TRUE){out[[7]]=tau2;names(out)[7]="tau2"}

  return(out)
}





