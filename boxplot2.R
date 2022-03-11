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

# (model 1) y = (x'b1) + e
# (model 2) y = 2 (X'b1)^2 + e
# (model 3) y = (X'b1)^2 + (X'b2)^2 + e
# (model 4) y = 0.5 + (x'b1)^3 + e
# (model 5) y = exp(x'b1) + e
# (model 6) y = X'b1 / (0.5 + X'b2^2) + e x'b1
######################################################################

gen.xy <- function(model, n, p){
  beta <- cbind(c(1,-1,2,-2,rep(0,p-4)),
                c(rep(0,p-4),1,-1,2,-2))/sqrt(10)
  
  sig <- sqrt(2)
  eps <- rnorm(n, 0, sig^2)
  A <- matrix(runif(p^2, -1, 1), nrow=p)
  SIG <- A %*% t(A) + diag(p)
  
  x <- rmvnorm(n, mean=rep(0,p), sigma=SIG)
  
  if (model==1){
    y <- x %*% beta[,1] + eps ; k <- 1
  } else if (model==2){
    y <- 2 * (x %*% beta[,1])^2 + eps ; k <- 1
  } else if (model==3){
    y <- (x %*% beta[,1])^2 + (x %*% beta[,2])^2 + eps ; k <- 2
  } else if (model==4){
    y <- 0.5 + (x %*% beta[,1])^3 + eps ; k <- 1
  } else if (model==5){
    y <- exp(x %*% beta[,1]) + eps ; k <- 1
  } else if (model==6){
    y <- x %*% beta[,1] / (0.5 + (x %*% beta[,2])^2) + eps * x %*% beta[,1] ; k <- 2
  }
  
  return(list(x=x, y=y, k=k, beta=beta[,1:k]))
}

######################################################################
# function to compare SAVE & BIG-SAVE for each model
######################################################################

# calculate quality measure for each n & ng
compare.qual <- function(model, p, n.list){
  set.seed(1234)
  ng.list <- c(1, 10, 50, 100)
  qual.list <- c()
  
  # 500 repetitions for each n
  for (n in n.list){
    for (m in 1:500){
      data <- gen.xy(model, n, p)
      x <- data$x ; y <- data$y ; k <- data$k ; true.beta <- data$beta
      # calculate the estimator & measure the quality
      for (ng in ng.list){
        beta <- BIG_SAVE(x,y,k,ng,p)
        qual <- dist(beta, true.beta)
        qual.list <- rbind(qual.list, c(n,ng,qual))
      }
    }
  }
  
  colnames(qual.list) <- c("n", "ng", "quality")
  qual.list <- qual.list %>%
    as.data.frame() %>%
    mutate(ng = factor(ng, levels=c("1","10","50","100")),
           quality = as.numeric(quality))
  qual.list$ng <- recode_factor(qual.list$ng, "1"="SAVE")
  
  return(qual.list)
}

# draw boxplot of quality measure
draw.boxplot <- function(model, qual.list) {
  ggplot(qual.list, aes(x=ng, y=quality)) +
    geom_boxplot(aes(group=ng)) +
    facet_wrap(~ n, ncol=4, scales = "free_x") +
    labs(title=paste("BIG-SAVE: model",model), x="G", y="Quality measure")
}


n.list <- c(10^4, 5*10^4, 10^5, 5*10^5)
model1 <- compare.qual(model=1, p=10, n.list)
draw.boxplot(1, model1)

n.list <- c(5*10^3, 10^4, 5*10^4, 10^5)
model2 <- compare.qual(model=2, p=20, n.list)
draw.boxplot(2, model2)

n.list <- c(5*10^3, 10^4, 5*10^4, 10^5)
model3 <- compare.qual(model=3, p=25, n.list)
draw.boxplot(3, model3)

n.list <- c(5*10^3, 10^4, 5*10^4, 10^5)
model4 <- compare.qual(model=4, p=15, n.list)
draw.boxplot(4, model4)

n.list <- c(5*10^3, 10^4, 5*10^4, 10^5)
model5 <- compare.qual(model=5, p=15, n.list)
draw.boxplot(5, model5)

n.list <- c(10^5, 5*10^5, 10^6, 3*10^6)
model6 <- compare.qual(model=6, p=10, n.list)
draw.boxplot(6, model6)
