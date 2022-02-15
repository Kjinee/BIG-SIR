######################################################################
# Simulation : Comparison of SAVE and BIG-SAVE approaches
######################################################################

library(edrGraphicalTools)
library(doSNOW)
library(dplyr)
library(ggplot2)

######################################################################
# quality measure
######################################################################

matpower <- function(a, alpha){
  a <- (a + t(a))/2 # symmetrize a
  tmp <- eigen(a)
  power <- tmp$vectors %*% diag((tmp$values)^alpha) %*% t(tmp$vectors)
  return(power)
}

dist <- function(v1, v2) {
  v1 <- as.matrix(v1)
  v2 <- as.matrix(v2)
  if (dim(v1)[2] == 1) {
    p1 <- v1 %*% t(v1) / c(t(v1)%*%v1)
    p2 <- v2 %*% t(v2) / c(t(v2)%*%v2)
  } else {
    p1 <- v1 %*% matpower(t(v1)%*%v1, -1) %*% t(v1)
    p2 <- v2 %*% matpower(t(v2)%*%v2, -1) %*% t(v2)
  }
  d <- sqrt( sum((p1-p2)*(p1-p2)) )
  return(d)
}

######################################################################
# function to calculate BIG-SAVE estimator (using loop)
# x, y : data
# k : single index(k=1) or multiple indices(k>1)
# ng : number of blocks(chunks)
######################################################################

BIG_SAVE <- function(x,y,k,ng,p=10){
  
  matEDR.block <- function(x){
    bhat <- matrix(x / sqrt((sum(x**2))), ncol=k)
    return(bhat %*% t(bhat) / k)
  }
  
  Scalable.save <- function(ng, data, size.chunk){
    rows <- ((ng-1)*size.chunk+1):(ng*size.chunk)
    matEDR.block(edr(data[rows,1], data[rows,-1], H=8, K=k, method="SAVE")$matEDR[,1:k])
  }
  
  size.chunk <- nrow(x)/ng
  
  BIGsave <- matrix(0, nrow=p, ncol=p)
  for (g in 1:ng){
    temp <- Scalable.save(g, cbind(y,x), size.chunk)
    BIGsave <- BIGsave + temp
  }
  
  BIGsave.estimator <- eigen(BIGsave)$vectors[,1:k]
  return(BIGsave.estimator)
}

######################################################################
# function to generate data

# (model 1) y = 2 (X'b1)^2 + e
# (model 2) y = 0.4 (X'b1)^2 + (X'b2)^2 + e
# (model 3) y = sin(X'b1) + e
# (model 4) y = sin(X'b1) + (x'b2)^2 + e
# (model 5) y = 0.5 + X'b1^3 + e
# (model 6) y = X'b1 / (0.5 + X'b2^2) + e
######################################################################

gen.xy <- function(model, n, p=10){
  beta <- cbind(c(1,-1,2,-2,0,0,0,0,0,0),
                     c(0,0,0,0,0,0,1,-1,2,-2))/sqrt(10)
  
  sig <- sqrt(2)
  eps <- rnorm(n, 0, sig^2)
  A <- matrix(runif(p^2, -1, 1), nrow=p)
  SIG <- A %*% t(A) + diag(p)
  
  x <- rmvnorm(n, mean=rep(0,p), sigma=SIG)
  
  if (model==1){
    y <- 2 * (x %*% beta[,1])^2 + eps ; k <- 1
  } else if (model==2){
    y <- (x %*% beta[,1])^2 + (x %*% beta[,2])^2 + eps ; k <- 2
  } else if (model==3){
    y <- sin(x %*% beta[,1]) + eps ; k <- 1
  } else if (model==4){
    y <- sin(x %*% beta[,1]) + (x %*% beta[,2])^2 + eps ; k <- 2
  } else if (model==5){
    y <- 0.5 + (x %*% beta[,1])^3 + eps ; k <- 1
  } else if (model==6){
    y <- (x %*% beta[,1]) / (0.5 + (x %*% beta[,2])^2) + eps ; k <- 2
  }
  
  return(list(x=x, y=y, k=k, beta=beta[,1:k]))
}

######################################################################
# function to draw box plot to compare SAVE & BIG-SAVE for each model
######################################################################

draw.boxplot <- function(model){
  set.seed(1234)
  
  qual.list <- c()
  n.list <- c(10^3, 5*10^3, 10^4, 5*10^4)
  
  for (n in n.list){
    if (n == 10^3){
      ng.list <- c(1, 10, 20)
    } else {
      ng.list <- c(1, 10, 50, 100)
    }
    
    # 500 repetitions
    for (m in 1:500){
      data <- gen.xy(model, n)
      x <- data$x ; y <- data$y ; k <- data$k ; true.beta <- data$beta
      # calculate the estimator & measure the quality
      for (ng in ng.list){
        beta <- BIG_SAVE(x,y,k,ng)
        qual <- dist(beta, true.beta)
        qual.list <- rbind(qual.list, c(n,ng,qual))
      }
    }
  }
  
  colnames(qual.list) <- c("n", "ng", "quality")
  qual.list <- qual.list %>%
    as.data.frame() %>%
    mutate(ng = factor(ng, levels=c("1","10","20","50","100")),
           quality = as.numeric(quality))
  qual.list$ng <- recode_factor(qual.list$ng, "1"="SAVE")
  
  ggplot(qual.list, aes(x=ng, y=quality)) +
    geom_boxplot(aes(group=ng)) +
    facet_wrap(~ n, ncol=4, scales = "free_x") +
    labs(title=paste("BIG-SAVE: model",model), x="G", y="Quality measure")
}

draw.boxplot(1)
draw.boxplot(2)
draw.boxplot(3)
draw.boxplot(4)
draw.boxplot(5)
draw.boxplot(6)
