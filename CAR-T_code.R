library(dplyr)
library(rjags); library(coda)
library(dfcrm)
library(MASS)

# Simulation function
CAR_sim <- function(dose, time, beta1_true, beta2_true, sigma_true, rho_true, PI, TTL, alpha, n.cohort, 
                    priorMean, priorCov, a, b, delta, c_overdose, c1=0.1, c2=0.1, c3=0.1, E_scen, omega=NULL, z_obs=NULL,
                    T_threshold=NULL, Z_threshold=NULL, C1=FALSE, C2=FALSE, AUC=FALSE, x0=1, n.chains=1, n.adapt=1e+4) {
  # Inputs:
  # ----------------------------------------
  # n.cohort: each entry gives the number of patients in each cohort
  # E_scen  : 1: efficacy data simulated from the Bi-Exp model
  #           2: Bi-Exp with random effect
  #           3: Two-population 
  #           4: observed data + noise
  # z_obs   : KxJ matrix, each column contains the observed expansion for a dose
  
  # rjags to simulate from the posterior
  jags.script <- "
  model {
    # Likelihood
    for(i in 1:length(x)) {
      for(k in 1:length(time)) {
        f[i,k] <- x[i]*(exp(beta1*time[k]*x[i]+log(2))- exp(beta2*time[k]*x[i]))
      }
      Z[i,1:length(time)] ~ dmnorm.vcov(f[i,1:length(time)], Sigma)
    }
    
    # Prior
    beta1 ~ dnorm(mu[1], (sigma2E*priorCov[1])^(-1)) T(,0)
    beta2 ~ dnorm(mu[2], (sigma2E*priorCov[2])^(-1)) T(,beta1)
    sigma_inv ~ dgamma(a, b)
    sigma2E <- sigma_inv^(-1)
    rho ~ dunif(0,1)
    for (k in 1:length(time)) {
      Sigma[k,k] <- sigma2E
    }
    for (AR in 1:(length(time)-1)) {
      for (k in (AR+1):length(time)) {
        Sigma[k-AR,k] <- sigma2E*rho^AR
        Sigma[k,k-AR] <- sigma2E*rho^AR
      }
    }
    
    # escalation crit
    for(j in 1:length(d)) {
      C1[j] <- d[j]*(exp(beta1*d[j]*T_threshold+log(2))-exp(beta2*d[j]*T_threshold))
      for(t in 1:time[length(time)]) {
        Z2[j,t] <- d[j]*(exp(beta1*d[j]*t+log(2))-exp(beta2*d[j]*t))
        C[j,t] <- (Z2[j,t]>=Z_threshold)
      }
      C2[j] <- sum(C[j,1:time[length(time)]])
      C3[j] <- (exp(beta1*d[j]*T_max+log(2))-2)/beta1 - (exp(beta2*d[j]*T_max)-1)/beta2
    }
  }
  "
  Calc_D <- Vectorize(function(b1, b2, d, Z_threshold, T_max=91) {
    mod.Z <- function(b1, b2, d, t) d * (exp(b1*t*d+log(2)) - exp(b2*t*d))
    f <- function(t, d, b1, b2) {
      Z <- d * (exp(b1*t*d+log(2)) - exp(b2*t*d))
      (Z - Z_threshold)^2
    }
    t_max <- (log(b2/b1)-log(2))/((b1-b2)*d)
    D <- L <- U <- numeric()
    for (j in 1:length(d)) {
      if(mod.Z(b1, b2, d[j], t_max[j])<Z_threshold) {
        D[j] <- 0
      } else {
        if(d[j]>=Z_threshold) {
          L[j] <- 0
        } else {
          L[j] <- optimise(f, c(0,t_max[j]), d=d[j], b1=b1, b2=b2)$minimum
        }
        U[j] <- optimise(f, c(max(0,t_max[j]), T_max), d=d[j], b1=b1, b2=b2)$minimum
        D[j] <- U[j]-L[j]
      }
    }
    D
  }, c("b1", "b2"))
  Calc_AUC <- Vectorize(function(b1, b2, d, T_max=91) {
    AUC_fit <- (exp(b1*d*T_max+log(2))-2)/b1 - (exp(b2*d*T_max)-1)/b2
    sel <- AUC_fit >= 0.9*AUC_fit[length(d)]
    AUC_fit
  }, c("b1", "b2"))
  # set up
  K <- length(time)
  x <- pop <- u <- mu_E <- numeric(sum(n.cohort))  # dose assigned to patient i
  x[1:n.cohort[1]] <- dose[x0]
  Z <- matrix(nrow=sum(n.cohort), ncol=K)
  y <- numeric(sum(n.cohort))
  Sigma <- matrix(nrow=K, ncol=K)
  for (AR in 0:(K-1)) {
    for (k in (AR+1):K) {
      Sigma[k-AR,k] <- rho_true^AR
      Sigma[k,k-AR] <- rho_true^AR
    }
  }
  Sigma <- sigma_true^2*Sigma
  i <- 1
  Continue <- TRUE
  # simulate the trial
  while (i<=length(n.cohort) && Continue) {
    # simulate pseudo data
    if(i==1) {
      for(j in 1:sum(n.cohort[1:i])) {
        # Safety data
        y[j] <- rbinom(1, 1, PI[which(dose==x[j])])
        # Efficacy data
        if (E_scen==1) {
          pop[j] <- 1
          Z[j,1:K] <- x[j]*( exp(beta1_true*x[j]*time+pop[j]*log(2)) - exp(beta2_true*x[j]*time)*pop[j] ) +
            mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==2) {
          u[j] <- rnorm(1, 0, omega)
          Z[j,1:K] <- x[j]*( exp(beta1_true*x[j]*time+log(2)) - exp(beta2_true*x[j]*time) ) + u[j] +
            mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==3) {
          pop[j] <- rbinom(1, 1, prob=0.2*which(dose==x[j]))
          Z[j,1:K] <- x[j]*( exp(beta1_true*x[j]*time+pop[j]*log(2)) - exp(beta2_true*x[j]*time)*pop[j] ) +
            mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==4) {
          d_ast <- which(d==x[j])
          Z[j,1:K] <- z_obs[1:K,d_ast] + mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==5) {
          u[j] <- rnorm(1, 0, omega)
          mu_E <- x[j]*( exp(beta1_true*x[j]*time+log(2)) - exp(beta2_true*x[j]*time) ) + u[j]
          Z[j,1:K] <- rnbinom(K, mu=mu_E, size=sigma_true^2)
        }
      }
    } else {
      for(j in (sum(n.cohort[1:(i-1)])+1):sum(n.cohort[1:i])) {
        y[j] <- rbinom(1, 1, PI[which(dose==x[j])])
        if (E_scen==1) {
          pop[j] <- 1
          Z[j,1:K] <- x[j]*( exp(beta1_true*x[j]*time+pop[j]*log(2)) - exp(beta2_true*x[j]*time)*pop[j] ) +
            mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==2) {
          u[j] <- rnorm(1, 0, omega)
          Z[j,1:K] <- x[j]*( exp(beta1_true*x[j]*time+log(2)) - exp(beta2_true*x[j]*time) ) + u[j] +
            mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==3) {
          pop[j] <- rbinom(1, 1, prob=0.2*which(dose==x[j]))
          Z[j,1:K] <- x[j]*( exp(beta1_true*x[j]*time+pop[j]*log(2)) - exp(beta2_true*x[j]*time)*pop[j] ) +
            mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==4) {
          d_ast <- which(d==x[j])
          Z[j,1:K] <- z_obs[1:K,d_ast] + mvrnorm(1, mu=rep(0, K), Sigma)
        }
        if (E_scen==5) {
          u[j] <- rnorm(1, 0, omega)
          mu_E <- x[j]*( exp(beta1_true*x[j]*time+log(2)) - exp(beta2_true*x[j]*time) ) + u[j]
          Z[j,1:K] <- rnbinom(K, mu=mu_E, size=sigma_true^2)
        }
      }
    }
    # fit model
    # Safety model
    level <- sapply(x[1:sum(n.cohort[1:i])], function(i) which(dose==i))
    tox.fit <- crm(alpha, TTL, tox=y[1:sum(n.cohort[1:i])], level=level, model="empiric",
                   conf.level = 1-2*c_overdose)
    d.accept <- dose[tox.fit$ptoxU<=TTL+2*delta]  # safe doses
    if(length(d.accept)==0){
      Continue <- FALSE
    }
    if(Continue){
      # Efficacy model
      model.fit <- jags.model(textConnection(jags.script), quiet = TRUE,
                              data=list(Z=Z[1:sum(n.cohort[1:i]),], mu=priorMean, 
                                        priorCov=priorCov, x=x[1:sum(n.cohort[1:i])], 
                                        time=time, a=a, b=b, d=d.accept, 
                                        T_threshold=T_threshold, T_max=t[K], Z_threshold=Z_threshold),
                              n.chains=chains, n.adapt=n.adapt)
      update(model.fit, n.adapt, progress.bar="none")
      tt<-jags.samples(model.fit, 
                       c('beta1','beta2', 'sigma2E', 'rho', 'C1', 'C2', 'C3'), 
                       n.adapt, progress.bar="none")
      beta1 <- mean(tt$beta1[1,,])
      beta2 <- mean(tt$beta2[1,,])
      sigma2E <- mean(tt$sigma2E[1,,])
      sigma <- sqrt(sigma2E)
      rho <- mean(tt$rho[1,,])
      
      # recommend next dose level
      if(i<length(n.cohort)) {
        time_fit <- Z_fit <- AUC_fit <- numeric(length(d.accept))
        if(C1) {
          for (j in 1:length(d.accept)) {
            time_fit[j] <- mean(tt$C1[j,,])
          }
          j_ast <- which.max(time_fit)
          indiff <- (1-c1)*time_fit[j_ast]
          x[(sum(n.cohort[1:i])+1):sum(n.cohort[1:(i+1)])] <- d.accept[which(time_fit>=indiff)[1]]
        }
        if(C2) {
          for (j in 1:length(d.accept)) {
            Z_fit[j] <- mean(tt$C2[j,,])
          }
          j_ast <- which.max(Z_fit)
          indiff <- (1-c2)*Z_fit[j_ast]
          x[(sum(n.cohort[1:i])+1):sum(n.cohort[1:(i+1)])] <- d.accept[which(Z_fit>=indiff)[1]]
        }
        if(AUC==TRUE) {
          for (j in 1:length(d.accept)) {
            AUC_fit[j] <- mean(tt$C3[j,,])
          }
          j_ast <- which.max(AUC_fit)
          indiff <- (1-c3)*AUC_fit[j_ast]
          x[(sum(n.cohort[1:i])+1):sum(n.cohort[1:(i+1)])] <- d.accept[which(AUC_fit>=indiff)[1]]
        }
      }
    }
    
    i <- i+1
  }
  d.select <- rep(0, length(dose))
  if(Continue) {
    d.select[which(dose==x[length(x)])] <- 1
  }
  return(list(d.sel=d.select, x=x, pop=pop, y=y, ptox=tox.fit$ptox, a=tox.fit$estimate, Z=Z,
              b1=beta1, b2=beta2, sigma2E=sigma, rho=rho, C1=time_fit, C2=Z_fit, C3=AUC_fit))
}


# non-parametric benchmark
Benchmark <- function(n, S, K, rho, d, PI, nu, delta, c_overdose, b1, b2, sigmaE, 
                      criterion, T_threshold=NULL, Z_threshold=NULL, 
                      E_scen, omega=1, z_obs=NULL) {
  # Inputs:
  # -------------------------------------------------
  # n        : number of patients.
  # S        : number of simulated trials.
  # K        : number of measurements.
  # d        : dose levels.
  # PI       : true toxicity probabilities 
  # nu       : TTL
  # delta    : width of acceptable toxicity interval
  # c_overdose: threshold for overdose control
  # b1, b2, sigmaE: true parameter values.
  # criterion: escalation/de-escalation criteria.
  #            1 for choosing the highest at T_threshold;
  #            2 for choosing the longest duration above Z_threshold;
  #            3 for choosing the lowest dose within 10% highest AUC.
  select <- dplyr::select
  t <- seq(0, by=7, len=K)
  n <- nrow(z_obs)
  Ben <- function(s, criterion) {
    d.sel <- numeric(length(d))
    Continue <- TRUE
    # Safety model
    u_T <- runif(n)
    Y <- sapply(1:n, function(i) as.integer(u_T[i]<PI)) %>% t()
    p_hat <- colMeans(Y)
    p_overdose <- 1 - pnorm(sqrt(n)*(nu+2*delta-p_hat)/sqrt(p_hat*(1-p_hat)))
    d.safe <- d[p_overdose<c_overdose]
    if(length(d.safe)==0){
      Continue <- FALSE
    }
    if(Continue) {
      # Efficacy model
      mu_x <- rep(0, K)
      Sigma_x <- matrix(0, nrow=K, ncol=K)
      for (AR in 0:(K-1)) {
        for (k in AR:K) {
          Sigma_x[k-AR,k] <- rho^AR
          Sigma_x[k,k-AR] <- rho^AR
        }
      }
      u_E <- sapply(1:n, function(i) {
        x_i <- MASS::mvrnorm(1, mu_x, Sigma_x)
        pnorm(x_i)
      }) %>% t()
      if(E_scen==1) {
        Z <- array(dim=c(n, K, length(d.safe)))
        for (j in 1:length(d.safe)) {
          Z[,,j] <- sapply(1:n, function(i) {
            mu_j <- d.safe[j]*(exp(b1*t*d.safe[j]+log(2))-exp(b2*t*d.safe[j]))
            qnorm(u_E[i,], mu_j, sigmaE)
          }) %>% t()
        }
      }
      if(E_scen==2) {
        Z <- array(dim=c(n, K, length(d.safe)))
        for (j in 1:length(d.safe)) {
          Z[,,j] <- sapply(1:n, function(i) {
            # find sample quantile
            z_quan <- numeric(K)
            mu_j <- d.safe[j]*(exp(b1*t*d.safe[j]+log(2))-exp(b2*t*d.safe[j]))
            rand_i <- rnorm(1000, 0, omega)
            for (k in 1:K) {
              z_quan[k] <- quantile(mu_j[k]+rnorm(1000,0,sigmaE)+rand_i, probs=u_E[i,k])
            }
            z_quan
          }) %>% t()
        }
      }
      if(E_scen==3) {
        Z <- array(dim=c(n, K, length(d.safe)))
        for (j in 1:length(d.safe)) {
          Z[,,j] <- sapply(1:n, function(i) {
            z_quan <- numeric(K)
            #pop <- rbinom(length(d.safe), 1, prob=0.2*(1:length(d.safe)))
            #mu_j <- d.safe[j]*(exp(b1*t*d.safe[j]+pop[j]*log(2))-pop[j]*exp(b2*t*d.safe[j]))
            for (k in 1:K) {
              pop <- rbinom(1000, 1, prob=0.2*j)
              mu_jk <- d.safe[j]*(exp(b1*t[k]*d.safe[j]+pop*log(2))-pop*exp(b2*t[k]*d.safe[j]))
              z_quan[k] <- quantile(mu_jk+rnorm(1000,0,sigmaE), probs=u_E[i,k])
            }
            z_quan
          }) %>% t()
        }
      }
      if(E_scen==4) {
        Z <- array(dim=c(n, K, length(d.safe)))
        for (j in 1:length(d.safe)) {
          Z[,,j] <- sapply(1:n, function(i) {
            mu_j <- z_obs[,j]
            qnorm(u_E[i,], mu_j, sigmaE)
          }) %>% t()
        }
      }
      if(E_scen==5){
        Z <- array(dim=c(n, K, length(d.safe)))
        for (j in 1:length(d.safe)) {
          Z[,,j] <- sapply(1:n, function(i) {
            # find sample quantile
            z_quan <- numeric(K)
            mu_j <- d.safe[j]*(exp(b1*t*d.safe[j]+log(2))-exp(b2*t*d.safe[j]))
            rand_i <- rnorm(1000, 0, omega)
            for (k in 1:K) {
              z_quan[k] <- quantile(rnbinom(length(mu_j), mu=mu_j+rand_i, size=sigmaE^2), probs=u_E[i,k])
            }
            z_quan
          }) %>% t()
        }
      }
      mu_hat <- apply(Z, c(2,3), mean)
      if(criterion==1) {
        mu_new <- numeric(ncol(mu_hat))
        for (j in 1:ncol(mu_hat)) {
          mu_new[j] <- approx(t, mu_hat[,j], xout=T_threshold, rule=2)$y
        }
        j_ast <- which.max(mu_new)
        indiff <- 0.9*mu_new[j_ast]
        j_1 <- which(mu_new>=indiff)[1]
        d.sel[which(d==d.safe[j_1])] <- 1
      }
      if(criterion==2) {
        mu_new <- matrix(nrow=t[K], ncol=ncol(mu_hat))
        for (j in 1:ncol(mu_hat)) {
          mu_new[,j] <- approx(t, mu_hat[,j], xout=1:t[K], rule=2)$y
        }
        D_hat <- apply(mu_new, 2, function(v) sum(v>=Z_threshold))
        j_ast <- which.max(D_hat)[1]
        indiff <- 0.9*D_hat[j_ast]
        j_2 <- which(D_hat>=indiff)[1]
        d.sel[which(d==d.safe[j_2])] <- 1
      }
      if(criterion==3) {
        mu_new <- matrix(nrow=(1+t[K]), ncol=ncol(mu_hat))
        for (j in 1:ncol(mu_hat)) {
          mu_new[,j] <- approx(t, mu_hat[,j], xout=0:t[K], rule=2)$y
        }
        AUC_hat <- sapply(1:length(d.safe), function(j) sum(mu_new[(2:(t[K]+1)),j]))
        j_ast <- which.max(AUC_hat)[1]
        indiff <- 0.9*AUC_hat[j_ast]
        j_3 <- which(AUC_hat>=indiff)[1]
        d.sel[which(d==d.safe[j_3])] <- 1
      }
    }
    d.sel
  }
  d.sel <- sapply(1:S, Ben, criterion) %>% rowMeans()
  d.sel
}