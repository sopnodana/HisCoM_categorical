
############
########################################################
###################Penalty Function#####################
########################################################
QPenalty <- function(theta, lambda, penalty){
  ab.theta <- abs(theta)
  sn.theta <- sign(theta)
  if(penalty == "SCAD"){
    a=3.7
    tem0.vec = (lambda) * (ab.theta < lambda)
    tem1.vec = ((a*lambda-ab.theta)/(a-1)) * (ab.theta>=lambda)*(ab.theta<a*lambda)
    penalty_val <- (tem0.vec + tem1.vec)#*sn.theta
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


########################################################
####################Internal Detrivation###############
########################################################

####Probability calculate
mu_prob <- function(cumprob, nobs, ncat1) {
  mat <- matrix(cumprob, nobs, ncat1, TRUE)
  mat <- rbind(mat[, 1], diff(t(mat)), 1 - mat[, ncat1])
  mat <- c(mat)
  return(mat)
}
