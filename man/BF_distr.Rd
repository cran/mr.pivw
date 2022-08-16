\name{BF_distr}
\alias{BF_distr}
\title{Generate bootstrap samples for the bootstrapping Fieller's confidence interval}
\usage{
  BF_distr(data,beta_hat,tau2,lambda,n_boot=1000,seed_boot=1)
}
\description{
  This function is used to generate bootstrap samples for the bootstrapping Fieller's confidence interval.
}

\arguments{
  \item{data}{A matrix or data frame consists of four columns: the 1st (2nd) column contains the estimated genetic effects on the outcome (exposure); the 3rd (4th) column contains the estimated standard errors of the estimated genetic effects on the outcome (exposure).}
  \item{beta_hat}{The estimated causal effect.}
  \item{tau2}{The estimated variance of the horizontal pleiotropy.}
  \item{lambda}{The penalty parameter in the pIVW estimator. The penalty parameter plays a role in the bias-variance trade-off of the estimator.}
  \item{n_boot}{The sample size of the bootstrap samples. By default, n_boot=1000.}
  \item{seed_boot}{The seed for random sampling in the bootstrap method. By default, seed_boot=1.}
}

\value{
  \item{z_b}{A vector containing the bootstrap samples for the bootstrapping Fieller's confidence interval.}
}



