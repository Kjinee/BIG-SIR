######################################################################
# Simulation : Computing strategies and running time
######################################################################

library(edrGraphicalTools)
library(bigmemory)
library(foreach)
library(doSNOW)
library(ggplot2)


######################################################################
# function to generate data
######################################################################

gen.xy <- function(n, p) {
  beta <- c(1,-1,2,-2,rep(0,p-4))/sqrt(10)
  
  sig <- sqrt(2)
  eps <- rnorm(n, 0, sig^2)
  A <- matrix(runif(p^2, -1, 1), nrow=p)
  SIG <- A %*% t(A) + diag(p)
  
  x <- rmvnorm(n, mean=rep(0,p), sigma=SIG)
  y <- 2 * (x %*% beta)^2 + eps
  
  return(list(x=x, y=y))
}


######################################################################
# function to calculate mean running time for each strategies
######################################################################

compare.time <- function(key){

  matEDR.block <- function(x){
    bhat <- matrix(x/sqrt((sum(x**2))), ncol=1)
    return(bhat%*%t(bhat))
  }
  Scalable.save <- function(g, data, size.chunk){
    rows <- ((g-1)*size.chunk+1):(g*size.chunk)
    matEDR.block(edr(data[rows,1],data[rows,-1],H=8,K=1,method="SAVE")$matEDR[,1])
  }
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  
  runtime <- matrix(0, nrow=20, ncol=4) # 4 strategies, 20 reps
  for (m in 1:20){
    xy <- gen.xy(n,p)
    x <- xy$x ; y <- xy$y
    size.chunk <- n/ng
    
    #SAVE
    runtime[m,1] <- system.time({
      est <- edr(y,x,H=8,K=1,method="SAVE")$matEDR[,1]
    })[3]
    
    # BIG-SAVE : loop
    runtime[m,2] <- system.time({
      BIGsave <- matrix(0, nrow=p, ncol=p)
      for (g in 1:ng){
        temp <- Scalable.save(g, cbind(y,x), size.chunk)
        BIGsave <- BIGsave + temp
      }
      est <- eigen(BIGsave)$vectors[,1]
    })[3]
    
    # BIG-SAVE : foreach
    runtime[m,3] <- system.time({
      BIGsave <- foreach(g=1:ng, .combine="+")%dopar%{
        require("edrGraphicalTools")
        Scalable.save(g, cbind(y,x), size.chunk)}
      est <- eigen(BIGsave)$vectors[,1]
    })[3]
    
    # BIG-SAVE : bigmemory + foreach
    dataYX <- as.big.matrix(cbind(y,x), backingpath="E:/temp",
                            backingfile=paste0(key+m,"test-save.bin"),
                            descriptorfile=paste0(key+m,"test-save.desc"), type="double")
    BIGmatdes <- describe(dataYX)
    x <- attach.big.matrix(BIGmatdes)
    
    runtime[m,4] <- system.time({
      BIGsave <- foreach(g=1:ng, .combine="+")%dopar%{
        require("edrGraphicalTools")
        require("bigmemory")
        x <- attach.big.matrix(BIGmatdes)
        Scalable.save(g, x, size.chunk)
      }
      est <- eigen(BIGsave)$vectors[,1]
    })[3]
  }
  
  stopCluster(cl)
  
  # return mean running time of 20 reps
  result <- colMeans(runtime)
  return(result)
}


######################################################################
# calculate running time for each n, p, ng
######################################################################

# for n
p <- 15 ; ng <- 10
n.list <- c(5*10^5, 10^6, 5*10^6, 10^7)
n.time <- c()
key <- 1000 # backingfile num
for (n in n.list){
  result <- compare.time(key)
  n.time <- rbind(n.time,
                  data.frame(n=n,
                             strategy=c("SAVE","BIG-SAVE:loop","BIG-SAVE:foreach","BIG-SAVE:bigmemory+foreach"),
                             runtime=result))
  key <- key + 100
}

n.time$strategy <- as.factor(n.time$strategy)
ggplot(n.time, aes(x=as.factor(n), y=runtime, col=strategy, group=strategy)) +
  geom_line() +
  geom_point(aes(shape=strategy)) +
  scale_color_manual(values=c(2,3,4,1)) +
  labs(title="running time vs n", x="sample size", y="running time")



# for p
n <- 10^6 ; ng <- 10
p.list <- c(10,30,50,80)
p.time <- c()
key <- 2000
for (p in p.list){
  result <- compare.time(key)
  p.time <- rbind(p.time,
                  data.frame(p=p,
                             strategy=c("SAVE","BIG-SAVE:loop","BIG-SAVE:foreach","BIG-SAVE:bigmemory+foreach"),
                             runtime=result))
  key <- key + 100
}

p.time$strategy <- as.factor(p.time$strategy)
ggplot(p.time, aes(x=as.factor(p), y=runtime, col=strategy, group=strategy)) +
  geom_line() +
  geom_point(aes(shape=strategy)) +
  scale_color_manual(values=c(2,3,4,1)) +
  labs(title="running time vs p", x="covariate dimension p", y="running time")



# for g
p <- 15 ; n <- 10^6
g.list <- c(10,50,100,500)
g.time <- c()
key <- 3000
for (ng in g.list){
  result <- compare.time(key)[-1] # SAVE 결과 제외
  g.time <- rbind(g.time,
                  data.frame(g=ng,
                             strategy=c("BIG-SAVE:loop","BIG-SAVE:foreach","BIG-SAVE:bigmemory+foreach"),
                             runtime=result))
  key <- key + 100
}

g.time$strategy <- as.factor(g.time$strategy)
ggplot(g.time, aes(x=as.factor(g), y=runtime, col=strategy, group=strategy)) +
  geom_line() +
  geom_point(aes(shape=strategy)) +
  scale_color_manual(values=c(2,3,4)) +
  labs(title="running time vs G", x="number of chunks", y="running time")
