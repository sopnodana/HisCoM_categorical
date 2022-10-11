mu_prob <- function(cumprob, nobs, ncat1) {
  mat <- matrix(cumprob, nobs, ncat1, TRUE)
  mat <- rbind(mat[, 1], diff(t(mat)), 1 - mat[, ncat1])
  mat <- c(mat)
  return(mat)
}

norm    <- function(u) sqrt(sum(u^2)) 
PHARAOH_multinomial_ISLR_nosd <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(nobs/(nobs-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  w_Kro_idx        = which(t(W2) == 99)
  F_mat <- X_mat %*% W2
  
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  
  family <- make.link("logit") 
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  exclude <- seq(ncat, nobs * ncat, ncat)
  ncat1 <- ncat-1
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  
  ans <- diag(ncat1)
  ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
  
  while(iter < maxiter){
    W_old            <- W_new
    W2[w_star_99]    <- W_old
    beta_old         <- beta_new
    est_old          <- c(W_old, beta_old)
    F_mat            <- X_mat %*% W2
    eta              <- drop(F_mat %*% beta_old)   
    
    fitproball       <- mu_prob(linkinv(eta), nobs, ncat1)  ###linkinv(eta) produce cumulative probability
    fitprob          <- fitproball[-exclude]
    dummy            <- mu.eta(eta)
    pi               <- fitprob
    resids           <- y_mat - pi
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector         <- Nsubs == i
        mueta_i1         <- dummy[selector]
        mueta_i          <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
        dpi_eta          <- ans*mueta_i
        B_mat_i          <- B_mat[selector,]
        Pi_vct           <- pi[selector]
        V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        V_inv            <- ginv(V_mat)  #inversematrix
        res_i            <- resids[selector]
        G_i              <- t(dpi_eta) %*% V_inv %*% dpi_eta
        Z_i              <- eta[selector] + (ginv(dpi_eta)%*%res_i)
        Z_1i             <- Z_i - z2[selector]
        Sw_300           <- (t(B_mat_i) %*% G_i %*% Z_1i)   
        Sw_mat           <- Sw_mat + Sw_300
        Hw_300           <- (t(B_mat_i) %*% G_i %*% B_mat_i)
        Hw_mat           <- Hw_mat + Hw_300
      }
      ##
      w_j                <-   ginv(Hw_mat + lambda1*diag(Nj)) %*% Sw_mat 
      w_j                <-   sqrt(nobs)*w_j/ norm(X_j%*%w_j) 
      W2[k:kk,(ncat1+j)] <- w_j
      F_mat[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    
    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      mueta_i1 <- dummy[selector]
      mueta_i <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
      dpi_eta <- ans*mueta_i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      G_i <- t(dpi_eta) %*% V_inv %*% dpi_eta
      Z_i <- eta[selector] + (ginv(dpi_eta)%*%res_i)   
      Sb_300 <- (t(F_mat_i) %*% G_i %*% Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% G_i %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    p_beta <- rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    p_beta[c(1:ncat1, indx)] <- 0
    beta_new <- ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    est_new <- c(W_new, beta_new)
    crit <- sum(abs(est_new - est_old))
    iter <- iter + 1
    if(iter%%5==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  
  ## deviance calculate
  #eta_F <- drop(X_mat %*% W2 %*% beta_coef)
  #fitproball <- mu_prob(linkinv(eta_F), nobs, ncat1)  
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001)) 
  ##
  fit <- list()
  fit$xnames <- xnames
  fit$nvar <- nvar
  fit$beta_coef <- as.numeric(beta_coef)
  fit$weight_coef <- weight_coef
  #fit$sd <- sdP
  fit$W <- W2
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  #fit$beta <- beta_rslt
  fit
}

norm    <- function(u) sqrt(sum(u^2)) 
PHARAOH_multinomial_ISLR <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(nobs/(nobs-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  w_Kro_idx        = which(t(W2) == 99)
  #set.seed(2)
  W2[w_star_99] <- runif(num_windex)
  
  F_mat <- X_mat %*% W2
  
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  
  family <- make.link("logit") 
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  exclude <- seq(ncat, nobs * ncat, ncat)
  ncat1 <- ncat-1
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  
  ans <- diag(ncat1)
  ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
  
  while(iter < maxiter){
    W_old            <- W_new
    W2[w_star_99]    <- W_old
    beta_old         <- beta_new
    est_old          <- c(W_old, beta_old)
    F_mat            <- X_mat %*% W2
    eta              <- drop(F_mat %*% beta_old)   
    
    fitproball       <- mu_prob(linkinv(eta), nobs, ncat1)  ###linkinv(eta) produce cumulative probability
    fitprob          <- fitproball[-exclude]
    dummy            <- mu.eta(eta)
    pi               <- fitprob
    resids           <- y_mat - pi
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector         <- Nsubs == i
        mueta_i1         <- dummy[selector]
        mueta_i          <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
        dpi_eta          <- ans*mueta_i
        B_mat_i          <- B_mat[selector,]
        Pi_vct           <- pi[selector]
        V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        V_inv            <- ginv(V_mat)  #inversematrix
        res_i            <- resids[selector]
        G_i              <- t(dpi_eta) %*% V_inv %*% dpi_eta
        Z_i              <- eta[selector] + (ginv(dpi_eta)%*%res_i)
        Z_1i             <- Z_i - z2[selector]
        Sw_300           <- (t(B_mat_i) %*% G_i %*% Z_1i)   
        Sw_mat           <- Sw_mat + Sw_300
        Hw_300           <- (t(B_mat_i) %*% G_i %*% B_mat_i)
        Hw_mat           <- Hw_mat + Hw_300
      }
      ##
      w_j                <-   ginv(Hw_mat + lambda1*diag(Nj)) %*% Sw_mat 
      w_j                <-   sqrt(nobs)*w_j/ norm(X_j%*%w_j) 
      W2[k:kk,(ncat1+j)] <- w_j
      F_mat[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    
    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      mueta_i1 <- dummy[selector]
      mueta_i <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
      dpi_eta <- ans*mueta_i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      G_i <- t(dpi_eta) %*% V_inv %*% dpi_eta
      Z_i <- eta[selector] + (ginv(dpi_eta)%*%res_i)   
      Sb_300 <- (t(F_mat_i) %*% G_i %*% Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% G_i %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    p_beta <- rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    p_beta[c(1:ncat1, indx)] <- 0
    beta_new <- ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    est_new <- c(W_new, beta_new)
    crit <- sum(abs(est_new - est_old))
    iter <- iter + 1
    if(iter%%1==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  
  ####Information Matrix######
  
  Phi_mat    <-  kronecker(X_mat, t(beta_coef))[,w_Kro_idx]
  
  nw <- ncol(Phi_mat)
  nb <- ncol(F_mat)
  
  I11_mat  <- matrix(0, nw, nw, FALSE)  
  I12_mat  <- matrix(0, nw, nb, FALSE) 
  I22_mat  <- matrix(0, nb, nb, FALSE) 
  for(i in 1:nobs){
    selector         <- Nsubs == i
    mueta_i1         <- dummy[selector]
    mueta_i          <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
    dpi_eta          <- ans*mueta_i
    G_i              <- t(dpi_eta) %*% V_inv %*% dpi_eta
    
    
    Pi_vct           <- pi[selector]
    V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
    V_inv            <- ginv(V_mat)  #inversematrix
    res_i            <- resids[selector]
    
    
    Phi_mat_i          <- Phi_mat[selector,]
    F_mat_i          <- F_mat[selector,]
    
    I11_300           <- (t(Phi_mat_i ) %*% G_i %*% Phi_mat_i )
    I11_mat           <- I11_mat + I11_300
    
    I22_300           <- (t(F_mat_i ) %*% G_i %*% F_mat_i )
    I22_mat           <- I22_mat + I22_300
    
    I12_300           <- (t(Phi_mat_i ) %*% G_i %*% F_mat_i )
    I12_mat           <- I12_mat + I12_300
  }
  
  
  I_F <- rbind(cbind(I11_mat , I12_mat), cbind(t(I12_mat) , I22_mat))
  I <- I_F + diag(c(rep(lambda1,nw), rep(0, ncat1), rep(lambda2,(nb-ncat1-length(indx))), rep(0, length(indx))))
  
  Vcov <- ginv(I) %*%I_F %*%ginv(I)
  beta.std <- sqrt(diag(Vcov))[-(1:nw)]
  ##
  rr <- rankMatrix(Vcov)
  
  z.beta <- beta_coef/beta.std
  p.beta <- 2*pnorm(abs(z.beta),lower.tail=FALSE)#, scientific=TRUE, digits = 5)
  #p.beta <- format(pchisq(z.beta**2, df = df1,lower.tail=FALSE), scientific=TRUE, digits = 5)
  beta_rslt <- data.frame(round(beta_coef,5), beta.std, 
                          round(z.beta,5), p.beta)
  colnames(beta_rslt) <- c( "Estimate", "Std.Error", "z value", "Pr(>|z|)")
  
  ## deviance calculate
  #eta_F <- drop(X_mat %*% W2 %*% beta_coef)
  #fitproball <- mu_prob(linkinv(eta_F), nobs, ncat1)  
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001)) 
  ##
  fit <- list()
  fit$xnames <- xnames
  fit$nvar <- nvar
  fit$beta_coef <- as.numeric(beta_coef)
  fit$weight_coef <- weight_coef
  #fit$sd <- sdP
  fit$W <- W2
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  fit$beta <- beta_rslt
  fit
}


########
Hiscom_GLM_update_CV <- function(y, path, path_var, data = dat, lambda1 ,
                                 lambda2, fold, miter, tol1){
  n <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  grid1 <- lambda1
  grid2 <- lambda2
  d1 <- expand.grid(grid1, grid2)
  nlambda <- nrow(d1)
  #nlambda <- length(lambda1)
  #dat_cv <- data.frame(y, dat)
  set.seed(2)
  data.idx <- sample(1:nrow(dat),nrow(dat))
  CV_dat <- NULL
  for (i in 1:nlambda){
    tryCatch({
      cat("NCV", i, "\n")
      lam1 <- as.numeric(d1[i,])[1]
      lam2 <- as.numeric(d1[i,])[2]
      #lam1 <- lam2 <- lambda1[i]
      like_Li <- 0
      dev <- 0
      for( k in 1:fold){
        cat("fold=", k, "\n")
        testindex	<- (floor((k-1) * n/fold)+1) : floor(k * n/fold)
        trainindex	= setdiff(1:n, testindex)
        data_test	= dat[data.idx[testindex],]
        data_train	= dat[data.idx[trainindex],]
        
        ###
        # Specify Y 
        ytest		= y[data.idx[testindex]]
        ntest <- length(ytest)
        ytrain		= y[data.idx[trainindex]]
        
        ###
        Ridge_model <- PHARAOH_multinomial_ISLR(ytrain, path, path_var, data =  data_train,
                                                maxiter = 500,  lambda1 = lam1, lambda2 = lam2,
                                                tol = 0.001)
        
        w_train <- Ridge_model$W
        beta_train <- Ridge_model$beta_coef
        ###
        ncase <- length(ytest)
        Y1 <- rep(ytest, each = ncat)
        Intercept <- rep.int(seq(ncat), length(ytest))
        
        ###
        X.all <- data_test[,match(path_var, colnames(data_test))]
        xnames <- colnames(X.all)
        X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
        X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
        X_mat_1 <- scale(X_mat_1)*sqrt(ncase/(ncase-1))
        #diag(t(X_mat_1)%*%X_mat_1)
        Intercept1 <- rep.int(seq(ncat-1), length(ytest))
        X_mat_0 <- model.matrix(~factor(Intercept1)-1 )
        X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
        
        X_mat <- cbind(X_mat_0, X_mat_1)
        
        #####
        eta_F <- drop(X_mat %*% w_train %*% beta_train)
        family <- make.link("logit") 
        linkinv <- family$linkinv
        mu.eta <- family$mu.eta
        nobs <- length(ytest)
        
        
        fitprob <- mu_prob(linkinv(eta_F), nobs, ncat1)  
        Y1 <- rep(ytest, each = ncat)
        Intercept <- rep.int(seq(ncat), length(ytest))
        Y1_mat <- as.numeric(Y1 == Intercept)
        like_Li[k] <- sum(Y1_mat * log(fitprob))
        dev[k] <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 1e-08)) 
      } ###end K
      cv.like_Li <- sum(like_Li)/fold
      cv.dev <- sum(dev)/fold
      ##
      CV_rslt <- c(lam1, lam2, cv.like_Li, cv.dev)
      cat("CV_Value",  CV_rslt, "\n" )
      CV_dat <- rbind(CV_dat, CV_rslt)
      colnames(GCV_dat) <- c("lambda1", "lambda2", "cv.like_Li", "cv.dev")
    }, error = function(e){})
  } ##end i
  return(CV_dat)
}

##############################################
##### Same lambda#####
Hiscom_Ord_CV <- function(y, path, path_var, data, lambda.vec, fold){
  n <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  nlambda <- length(lambda.vec)
  set.seed(2)
  data.idx <- sample(1:nrow(data),nrow(data))
  CV_dat <- NULL
  for (i in 1:nlambda){
    tryCatch({
      cat("NCV", i, "\n")
      lam1 <- lam2 <- lambda.vec[i]
      like_Li <- c()
      dev <- c()
      for( k in 1:fold){
        tryCatch({
        cat("fold=", k, "\n")
        testindex	<- (floor((k-1) * n/fold)+1) : floor(k * n/fold)
        trainindex	= setdiff(1:n, testindex)
        data_test	= data[data.idx[testindex],] #data[testindex,] #
        data_train	= data[data.idx[trainindex],]  #data[-testindex,] # 
        
        ###
        # Specify Y 
        ytest		= y[data.idx[testindex]]  #y[testindex]  #
        ntest <- length(ytest)
        ytrain		= y[data.idx[trainindex]] #y[-testindex]   #
        #
        ###
        Ridge_model <- PHARAOH_multinomial_ISLR_CV(ytrain, path, path_var, 
                                                   indx= c(68,69,70), data =  data_train,
                                                   maxiter = 500,  lambda1 = lam1, lambda2 = lam2,
                                                   tol = 0.0001)
        w_train <- Ridge_model$W
        beta_train <- Ridge_model$beta_coef
        crit <- Ridge_model$crit
        cat("crit=", crit, "\n")
        ###
        ncase <- length(ytest)
        
        ###
        X.all <- data_test[,match(path_var, colnames(data_test))]
        xnames <- colnames(X.all)
        X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
        X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
        X_mat_1 <- scale(X_mat_1)*sqrt(ncase/(ncase-1))
        #diag(t(X_mat_1)%*%X_mat_1)
        Intercept1 <- rep.int(seq(ncat-1), length(ytest))
        X_mat_0 <- model.matrix(~factor(Intercept1)-1 )
        X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
        
        X_mat <- cbind(X_mat_0, X_mat_1)
        
        #####
        eta_F <- drop(X_mat %*% w_train %*% beta_train)
        family <- make.link("logit") 
        linkinv <- family$linkinv
        mu.eta <- family$mu.eta
        nobs <- length(ytest)
        
        
        fitprob <- mu_prob(linkinv(eta_F), nobs, ncat1)  
        Y1 <- rep(ytest, each = ncat)
        Intercept <- rep.int(seq(ncat), length(ytest))
        Y1_mat <- as.numeric(Y1 == Intercept)
        like_Li2 <- sum(Y1_mat * log(fitprob))
        dev2 <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 1e-08)) 
        
        tol <- 0.0001
        if(crit < tol){
          dev1 <- dev2
        }else{
          dev1 <- "NA"
        }
        
        if(crit < tol){
          like_Li1 <- like_Li2
        }else{
          like_Li1 <- NA
        }
        cat("dev1=", dev1, "like_Li1 =", like_Li1, "\n")
        like_Li[k]    <- like_Li1 
        dev[k]        <- dev1
       }, error = function(e){})   
      } ###end K
      cv.like_Li <- sum(like_Li)/fold
      cv.dev <- sum(dev)/fold
      ##
      CV_rslt <- c(lam1, lam2, cv.like_Li, cv.dev)
      cat("CV_Value",  CV_rslt, "\n" )
      CV_dat <- rbind(CV_dat, CV_rslt)
      colnames(CV_dat) <- c("lambda1", "lambda2", "cv.like_Li", "cv.dev")
    }, error = function(e){})
  } ##end i
  return(CV_dat)
}

PHARAOH_multinomial_ISLR_CV <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(nobs/(nobs-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  
  w_Kro_idx        = which(t(W2) == 99)
  #set.seed(2)
  W2[w_star_99] <- runif(num_windex)
  
  
  F_mat <- X_mat %*% W2
  
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  
  family <- make.link("logit") 
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  exclude <- seq(ncat, nobs * ncat, ncat)
  ncat1 <- ncat-1
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  
  ans <- diag(ncat1)
  ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
  
  
  while(iter < maxiter){
    W_old <- W_new
    W2[w_star_99] <- W_old
    beta_old <- beta_new
    est_old <- c(W_old, beta_old)
    F_mat <- X_mat %*% W2
    eta <- drop(F_mat %*% beta_old)   #drop(X_mat %*% coef_new)
    
    fitproball <- mu_prob(linkinv(eta), nobs, ncat1)  ###linkinv(eta) produce cumulative probability
    fitprob <- fitproball[-exclude]
    dummy <- mu.eta(eta)
    pi <- fitprob
    resids <- y_mat - pi
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector <- Nsubs == i
        #eta_i <- eta[selector]
        mueta_i1 <- dummy[selector]
        mueta_i<- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
        dpi_eta <- ans*mueta_i
        B_mat_i <- B_mat[selector,]
        Pi_vct <- pi[selector]
        V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        V_inv <- ginv(V_mat)  #inversematrix
        res_i <- resids[selector]
        G_i <- t(dpi_eta) %*% V_inv %*% dpi_eta
        Z_i <- eta[selector] + (ginv(dpi_eta)%*%res_i)
        Z_1i <- Z_i - z2[selector]
        Sw_300 <- (t(B_mat_i) %*% G_i %*% Z_1i)   
        Sw_mat <- Sw_mat + Sw_300
        Hw_300 <-  (t(B_mat_i) %*% G_i %*% B_mat_i)
        Hw_mat <- Hw_mat + Hw_300
      }
      ##
      w_j             =  ginv(Hw_mat + lambda1*diag(Nj)) %*% Sw_mat 
      w_j             = sqrt(nobs)*w_j/ norm(X_j%*%w_j) 
      W2[k:kk,(ncat1+j)]       = w_j
      F_mat[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    
    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      mueta_i1 <- dummy[selector]
      mueta_i <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
      dpi_eta <- ans*mueta_i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      G_i <- t(dpi_eta) %*% V_inv %*% dpi_eta
      Z_i <- eta[selector] + (ginv(dpi_eta)%*%res_i)   
      Sb_300 <- (t(F_mat_i) %*% G_i %*% Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% G_i %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    p_beta <- rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    p_beta[c(1:ncat1, indx)] <- 0
    
    beta_new <- ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    est_new <- c(W_new, beta_new)
    crit <- sum(abs(est_new - est_old))
    iter <- iter + 1
    if(iter%%5==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  ###
  fit <- list()
  fit$beta_coef <- as.numeric(beta_coef)
  fit$W <- W2
  fit$crit <- crit
  fit
}

###################Penalty Function#####################
########################################################
QPenalty <- function(theta, lambda, penalty){
  ab.theta <- abs(theta)
  sn.theta <- sign(theta)
  if(penalty == "SCAD"){
    a=3.7
    tem0.vec = (lambda) * (ab.theta < lambda)
    tem1.vec = ((a*lambda-ab.theta)/(a-1)) * (ab.theta>=lambda)*(ab.theta<a*lambda)
    penalty_val <- (tem0.vec + tem1.vec)*sn.theta
  }
  if(penalty == "MCP"){
    a <- 2.7
    b1 <- (lambda - ab.theta/a)*(ab.theta < a*lambda)
    penalty_val <- b1#*sn.theta
  }
  if(penalty == "LASSO"){
    penalty_val <- lambda*rep(1, length(theta))
  }
  if(penalty == "Ridge"){
    penalty_val <- 2*lambda*ab.theta
  }
  if(penalty == "TLP"){  ###Truncated Lasso 
    penalty_val <- lambda *(ab.theta < a)  
  }
  return(penalty_val)
}



PHisCoM_POM <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, penalty, eps, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(nobs/(nobs-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  w_Kro_idx        = which(t(W2) == 99)
  #set.seed(2)
  W2[w_star_99] <- runif(num_windex)
  F_mat <- X_mat %*% W2
  
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  
  family <- make.link("logit") 
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  exclude <- seq(ncat, nobs * ncat, ncat)
  ncat1 <- ncat-1
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  
  ans <- diag(ncat1)
  ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
  
  while(iter < maxiter){
    W_old            <- W_new
    W2[w_star_99]    <- W_old
    beta_old         <- beta_new
    est_old          <- c(W_old, beta_old)
    F_mat            <- X_mat %*% W2
    eta              <- drop(F_mat %*% beta_old)   
    
    fitproball       <- mu_prob(linkinv(eta), nobs, ncat1)  ###linkinv(eta) produce cumulative probability
    fitprob          <- fitproball[-exclude]
    dummy            <- mu.eta(eta)
    pi               <- fitprob
    resids           <- y_mat - pi
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      
      Wj_old <- W2[k:kk,(ncat1+j)]
      
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector         <- Nsubs == i
        mueta_i1         <- dummy[selector]
        mueta_i          <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
        dpi_eta          <- ans*mueta_i
        B_mat_i          <- B_mat[selector,]
        Pi_vct           <- pi[selector]
        V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        V_inv            <- ginv(V_mat)  #inversematrix
        res_i            <- resids[selector]
        G_i              <- t(dpi_eta) %*% V_inv %*% dpi_eta
        Z_i              <- eta[selector] + (ginv(dpi_eta)%*%res_i)
        Z_1i             <- Z_i - z2[selector]
        Sw_300           <- (t(B_mat_i) %*% G_i %*% Z_1i)   
        Sw_mat           <- Sw_mat + Sw_300
        Hw_300           <- (t(B_mat_i) %*% G_i %*% B_mat_i)
        Hw_mat           <- Hw_mat + Hw_300
      }
      ##
      E_w   <-  diag(QPenalty(abs(as.vector(Wj_old)), lambda1, penalty)/(abs(as.vector(Wj_old)) + eps))
      ##
      w_j                <-   ginv(Hw_mat + E_w) %*% Sw_mat 
      w_j                <-   sqrt(nobs)*w_j/ norm(X_j%*%w_j) 
      if(any(w_j == 0)){
        w_j <- w_j + 10^-5
      }else{
        w_j <- w_j
      }
      W2[k:kk,(ncat1+j)] <- w_j
      F_mat[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    
    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      mueta_i1 <- dummy[selector]
      mueta_i <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
      dpi_eta <- ans*mueta_i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      G_i <- t(dpi_eta) %*% V_inv %*% dpi_eta
      Z_i <- eta[selector] + (ginv(dpi_eta)%*%res_i)   
      Sb_300 <- (t(F_mat_i) %*% G_i %*% Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% G_i %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    pindex <- c(seq(1:ncat1), indx)
    
    E2 <- diag(QPenalty(abs(as.vector(beta_old)), lambda2, penalty)/(abs(as.vector(beta_old)) + eps))
    
    if(is.null(pindex) == TRUE){
      E_beta <- E2
    } else
      if(is.null(pindex) !=TRUE){
        E2[,pindex] <- 0
        E_beta <- E2
      }
    
    
    #p_beta <- rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    #p_beta[c(1:ncat1, indx)] <- 0
    beta_new <- ginv(Hb_mat + E_beta) %*% Sb_mat
    
    if(any(beta_new == 0)){
      beta_new <- beta_new + 10^-5
    }else{
      beta_new <- beta_new
    }
    
    
    est_new <- c(W_new, beta_new)
    crit <- sum(abs(est_new - est_old))
    iter <- iter + 1
    if(iter%%5==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  
  ####Information Matrix######
  
  Phi_mat    <-  kronecker(X_mat, t(beta_coef))[,w_Kro_idx]
  
  nw <- ncol(Phi_mat)
  nb <- ncol(F_mat)
  
  I11_mat  <- matrix(0, nw, nw, FALSE)  
  I12_mat  <- matrix(0, nw, nb, FALSE) 
  I22_mat  <- matrix(0, nb, nb, FALSE) 
  for(i in 1:nobs){
    selector         <- Nsubs == i
    mueta_i1         <- dummy[selector]
    mueta_i          <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
    dpi_eta          <- ans*mueta_i
    G_i              <- t(dpi_eta) %*% V_inv %*% dpi_eta
    
    
    Pi_vct           <- pi[selector]
    V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
    V_inv            <- ginv(V_mat)  #inversematrix
    res_i            <- resids[selector]
    
    
    Phi_mat_i          <- Phi_mat[selector,]
    F_mat_i          <- F_mat[selector,]
    
    I11_300           <- (t(Phi_mat_i ) %*% G_i %*% Phi_mat_i )
    I11_mat           <- I11_mat + I11_300
    
    I22_300           <- (t(F_mat_i ) %*% G_i %*% F_mat_i )
    I22_mat           <- I22_mat + I22_300
    
    I12_300           <- (t(Phi_mat_i ) %*% G_i %*% F_mat_i )
    I12_mat           <- I12_mat + I12_300
  }
  
  
  I_F <- rbind(cbind(I11_mat , I12_mat), cbind(t(I12_mat) , I22_mat))
  I <- I_F + diag(c(rep(lambda1,nw), rep(0, ncat1), rep(lambda2,(nb-ncat1-length(indx))), rep(0, length(indx))))
  
  Vcov <- ginv(I) %*%I_F %*%ginv(I)
  beta.std <- sqrt(diag(Vcov))[-(1:nw)]
  ##
  rr <- rankMatrix(Vcov)
  
  z.beta <- beta_coef/beta.std
  p.beta <- 2*pnorm(abs(z.beta),lower.tail=FALSE)#, scientific=TRUE, digits = 5)
  #p.beta <- format(pchisq(z.beta**2, df = df1,lower.tail=FALSE), scientific=TRUE, digits = 5)
  beta_rslt <- data.frame(round(beta_coef,5), beta.std, 
                          round(z.beta,5), p.beta)
  colnames(beta_rslt) <- c( "Estimate", "Std.Error", "z value", "Pr(>|z|)")
  
  ## deviance calculate
  #eta_F <- drop(X_mat %*% W2 %*% beta_coef)
  #fitproball <- mu_prob(linkinv(eta_F), nobs, ncat1)  
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001)) 
  ##
  fit <- list()
  fit$xnames <- xnames
  fit$nvar <- nvar
  fit$beta_coef <- as.numeric(beta_coef)
  fit$weight_coef <- weight_coef
  #fit$sd <- sdP
  fit$W <- W2
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  fit$beta <- beta_rslt
  fit
}

