#
# 4.2.1 Comparison of SIR and BIG-SIR approaches
# (single-index model)
#

library(edrGraphicalTools)
library(bigmemory)
library(foreach)
library(doSNOW)
library(ggplot2)

set.seed(1234)

# quality measure
dist <- function(v1, v2) {
  v1 <- as.matrix(v1)
  v2 <- as.matrix(v2)
  p1 <- v1 %*% t(v1) / c(t(v1)%*%v1)
  p2 <- v2 %*% t(v2) / c(t(v2)%*%v2)
  d <- sqrt( sum((p1-p2)*(p1-p2)) )
  return(d)
}


# single-index model (w/ sample size n)
true.beta <- c(1,-1,2,-2,0,0,0,0,0,0)/sqrt(10)
gen.xy <- function(n) {
  p <- 10
  beta <- c(1,-1,2,-2,0,0,0,0,0,0)/sqrt(10)
  
  sig <- sqrt(2)
  eps <- rnorm(n, 0, sig^2)
  A <- matrix(runif(p^2, -1, 1), nrow=p)
  SIG <- A %*% t(A) + diag(p)
  
  x <- rmvnorm(n, mean=rep(0,p), sigma=SIG)
  y <- (4/10) * (x %*% beta)^3 + eps
  
  return(list(x=x, y=y))
}


# function to calculate BIG-SIR estimator (using bigmemory + foreach)
BIG_SIR <- function(x,y,ng){
  dataYX <- as.big.matrix(cbind(y,x), type="double")
  BIGmatdes <- describe(dataYX)
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  
  x <- attach.big.matrix(BIGmatdes)
  
  matEDR.block <- function(x){
    bhat <- matrix(x / sqrt((sum(x**2))), ncol=1)
    return(bhat %*% t(bhat))
  }
  
  Scalable.sir <- function(ng, data, size.chunk){
    rows <- ((ng-1)*size.chunk+1):(ng*size.chunk)
    matEDR.block(edr(data[rows,1], data[rows,-1], H=8, K=1, method="SIR-I")$matEDR[,1])
  }
  
  size.chunk <- nrow(x)/ng
  
  BIGsir <- foreach(g=1:ng, .combine="+")%dopar%{
    require("edrGraphicalTools")
    require("bigmemory")
    x <- attach.big.matrix(BIGmatdes)
    Scalable.sir(g, x, size.chunk)
  }
  
  stopCluster(cl)
  BIGsir.estimator <- eigen(BIGsir)$vectors[,1]
  return(BIGsir.estimator)
}


# calculate quality measure for each n, G
n.list <- c(10^3, 5*10^3, 10^4, 5*10^4)
qual.list <- c()
for (n in n.list){
  
  if (n == 10^3){
    g.list <- c(10, 20)
  } else {
    g.list <- c(10, 50, 100)
  }
  
  # 500 repetitions
  for (m in 1:500){
    xy <- gen.xy(n)
    x <- xy$x ; y <- xy$y
    
    # SIR
    beta <- edr(y,x,H=8,K=1,method="SIR-I")$matEDR[,1]
    qual <- dist(beta, true.beta)
    qual.list <- rbind(qual.list, c(n,"SIR",qual))
    
    # BIG-SIR
    for (ng in g.list){
      beta <- BIG_SIR(x,y,ng)
      qual <- dist(beta, true.beta)
      qual.list <- rbind(qual.list, c(n,ng,qual))
    }
  }
}

qual.list <- as.data.frame(qual.list)
colnames(qual.list) <- c("n","g","quality")
qual.list$g <- as.factor(qual.list$g, levels=c("SIR","10","20","50","100"))
qual.list$quality <- as.numeric(qual.list$quality)

# box-plot of quality measure
ggplot(data = qual.list, aes(x=g, y=quality)) + 
  geom_boxplot(aes(group=g)) +
  facet_wrap(~ n, ncol=4, scales = "free_x") +
  labs(title="Single Index Model", x="g", y="Quality measure")
