######################################################################
# Simulation : Comparison of SAVE and BIG-SAVE approaches
######################################################################

library(edrGraphicalTools)
library(foreach)
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
# function to calculate BIG-SAVE estimator (using foreach)
# x, y : data
# k : single index(k=1) or multiple indices(k>1)
# ng : number of blocks(chunks)
######################################################################

BIG_SAVE <- function(x,y,k,ng){
  
  matEDR.block <- function(x){
    bhat <- matrix(x / sqrt((sum(x**2))), ncol=k)
    return(bhat %*% t(bhat) / k)
  }
  
  Scalable.save <- function(ng, data, size.chunk){
    rows <- ((ng-1)*size.chunk+1):(ng*size.chunk)
    matEDR.block(edr(data[rows,1], data[rows,-1], H=8, K=k,
                     method="SAVE")$matEDR[,1:k])
  }
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  
  size.chunk <- nrow(x)/ng
  
  BIGsave <- foreach(g=1:ng, .combine="+")%dopar%{
    require("edrGraphicalTools")
    Scalable.save(g, cbind(y,x), size.chunk)
  }
  
  stopCluster(cl)
  BIGsave.estimator <- eigen(BIGsave)$vectors[,1:k]
  return(BIGsave.estimator)
}

######################################################################
# model 1 : single index, asymmetric
######################################################################

set.seed(1234)
true.beta <- c(1,-1,2,-2,0,0,0,0,0,0)/sqrt(10)

# single-index model (w/ sample size n)
gen.xy <- function(n, p=10) {
  beta <- c(1,-1,2,-2,rep(0,p-4))/sqrt(10)
  
  sig <- sqrt(2)
  eps <- rnorm(n, 0, sig^2)
  A <- matrix(runif(p^2, -1, 1), nrow=p)
  SIG <- A %*% t(A) + diag(p)
  
  x <- rmvnorm(n, mean=rep(0,p), sigma=SIG)
  y <- (4/10) * (x %*% beta)^3 + eps
  
  return(list(x=x, y=y))
}

# calculate quality measure for each n, G
qual.list1 <- c()
n.list <- c(10^3, 5*10^3, 10^4, 5*10^4)

for (n in n.list){
  
  if (n == 10^3){
    g.list <- c(10, 20)
  } else {
    g.list <- c(10, 50, 100)
  }
  
  # 200 repetitions
  for (m in 1:200){
    xy <- gen.xy(n)
    x <- xy$x ; y <- xy$y
    
    # SAVE
    beta <- edr(y,x,H=8,K=1,method="SAVE")$matEDR[,1]
    qual <- dist(beta, true.beta)
    qual.list1 <- rbind(qual.list1, c(n,"SAVE",qual))
    
    # BIG-SAVE
    for (ng in g.list){
      beta <- BIG_SAVE(x,y,1,ng)
      qual <- dist(beta, true.beta)
      qual.list1 <- rbind(qual.list1, c(n,ng,qual))
    }
  }
}

colnames(qual.list1) <- c("n","g","quality")
qual.list1 <- qual.list1 %>%
  as.data.frame() %>%
  mutate(g = factor(g, levels=c("SAVE","10","20","50","100")),
         quality = as.numeric(quality))

# box-plot of quality measure
ggplot(data = qual.list1, aes(x=g, y=quality)) + 
  geom_boxplot(aes(group=g)) +
  facet_wrap(~ n, ncol=4, scales = "free_x") +
  labs(title="BIG-SAVE: Single-index Model", x="g", y="Quality measure")


######################################################################
# model 2 : multi index(K=2), symmetric
######################################################################

set.seed(4321)
true.beta <- cbind(c(1,-1,2,-2,0,0,0,0,0,0),
                   c(0,0,0,0,0,0,1,-1,2,-2))/sqrt(10)

# multi-index model (w/ sample size n)
gen.xy2 <- function(n, p=10) {
  beta1 <- c(1,-1,2,-2,0,0,0,0,0,0)/sqrt(10)
  beta2 <- c(0,0,0,0,0,0,1,-1,2,-2)/sqrt(10)
  
  sig <- sqrt(2)
  eps <- rnorm(n, 0, sig^2)
  A <- matrix(runif(p^2, -1, 1), nrow=p)
  SIG <- A %*% t(A) + diag(p)
  
  x <- rmvnorm(n, mean=rep(0,p), sigma=SIG)
  y <- (x%*%beta1)^2 + (x%*%beta2)^2 + eps
  
  return(list(x=x, y=y))
}

# calculate quality measure for each n, G
qual.list2 <- c()
n.list <- c(10^3, 5*10^3, 10^4, 5*10^4)

for (n in n.list){
  
  if (n == 10^3){
    g.list <- c(10, 20)
  } else {
    g.list <- c(10, 50, 100)
  }
  
  # 200 repetitions
  for (m in 1:200){
    xy <- gen.xy2(n)
    x <- xy$x ; y <- xy$y
    
    # SAVE
    beta <- edr(y,x,H=8,K=2,method="SAVE")$matEDR[,1:2]
    qual <- dist(beta, true.beta)
    qual.list2 <- rbind(qual.list2, c(n,"SAVE",qual))
    
    # BIG-SAVE
    for (ng in g.list){
      beta <- BIG_SAVE(x,y,2,ng)
      qual <- dist(beta, true.beta)
      qual.list2 <- rbind(qual.list2, c(n,ng,qual))
    }
  }
}

colnames(qual.list2) <- c("n","g","quality")
qual.list2 <- qual.list2 %>%
  as.data.frame() %>%
  mutate(g = factor(g, levels=c("SAVE","10","20","50","100")),
         quality = as.numeric(quality))

# box-plot of quality measure
ggplot(data = qual.list2, aes(x=g, y=quality)) + 
  geom_boxplot(aes(group=g)) +
  facet_wrap(~ n, ncol=4, scales = "free_x") +
  labs(title="BIG-SAVE: Multi-indices Model", x="g", y="Quality measure")
