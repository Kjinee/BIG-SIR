#
# 4.2.2 Computing strategies and running time
#

set.seed(4321)

# function to calculate mean running time for each strategies
compare.time <- function(key){
  
  matEDR.block <- function(x){
    bhat <- matrix(x/sqrt((sum(x**2))), ncol=1)
    return(bhat%*%t(bhat))
  }
  Scalable.sir <- function(g, data, size.chunk){
    rows <- ((g-1)*size.chunk+1):(g*size.chunk)
    matEDR.block(edr(data[rows,1],data[rows,-1],H=8,K=1,method="SIR-I")$matEDR[,1])
  }
  
  runtime <- matrix(0, nrow=50, ncol=4) # 4 strategies, 50 reps
  for (m in 1:50){
    xy <- gen.xy(n,p)
    x <- xy$x ; y <- xy$y
    
    #SIR
    start.time <- Sys.time()
    
    est <- edr(y,x,H=8,K=1,method="SIR-I")$matEDR[,1]
    
    end.time <- Sys.time()
    runtime[m,1] <- end.time - start.time
    
    # BIG-SIR : loop
    start.time <- Sys.time()
    
    size.chunk <- n/ng
    mg.hat <- matrix(0, nrow=p, ncol=p)
    for (g in 1:ng){
      temp <- Scalable.sir(g, cbind(y,x), size.chunk)
      mg.hat <- mg.hat + temp
    }
    est <- eigen(mg.hat)$vectors[,1]
    
    end.time <- Sys.time()
    runtime[m,2] <- end.time - start.time
    
    # BIG-SIR : foreach
    start.time <- Sys.time()
    
    cl <- makeCluster(4)
    registerDoSNOW(cl)
    
    size.chunk <- n/ng
    BIGsir <- foreach(g=1:ng, .combine="+")%dopar%{
      require("edrGraphicalTools")
      Scalable.sir(g, cbind(y,x), size.chunk)}
    
    stopCluster(cl)
    est <- eigen(BIGsir)$vectors[,1]
    
    end.time <- Sys.time()
    runtime[m,3] <- end.time - start.time
    
    # BIG-SIR : bigmemory + foreach
    start.time <- Sys.time()
    
    dataYX <- as.big.matrix(cbind(y,x), backingpath="E:/temp",
                            backingfile=paste0(key+m,"test-sir.bin"),
                            descriptorfile=paste0(key+m,"test-sir.desc"), type="double")
    BIGmatdes <- describe(dataYX)
    
    cl <- makeCluster(4)
    registerDoSNOW(cl)
    
    x <- attach.big.matrix(BIGmatdes)
    
    size.chunk <- n/ng
    
    BIGsir <- foreach(g=1:ng, .combine="+")%dopar%{
      require("edrGraphicalTools")
      require("bigmemory")
      x <- attach.big.matrix(BIGmatdes)
      Scalable.sir(g, x, size.chunk)
    }
    
    stopCluster(cl)
    est <- eigen(BIGsir)$vectors[,1]
    
    end.time <- Sys.time()
    runtime[m,4] <- end.time - start.time
  }
  # return mean running time of 50 reps
  return(colMeans(runtime))
}


# for n
p <- 10 ; ng <- 10
n.list <- c(5*10^5, 10^6, 5*10^6, 10^7)
n.time <- c()
key <- 1000 # code for backingfile name
for (n in n.list){
  result <- compare.time(key)
  n.time <- rbind(n.time,
                  data.frame(n=n,
                             strategy=c("SIR","BIG.SIR1","BIG.SIR2","BIG.SIR3"),
                             runtime=result))
  key <- key + 100
}

n.time$strategy <- as.factor(n.time$strategy)
ggplot(n.time, aes(x=n, y=runtime, col=strategy)) +
  geom_line() +
  geom_point(aes(shape=strategy)) +
  labs(title="running time vs n", x="sample size", y="running time")


# for p
n <- 10^6 ; ng <- 10
p.list <- c(5,10,20,30,50,80)
p.time <- c()
key <- 3000
for (p in p.list){
  result <- compare.time(key)
  p.time <- rbind(p.time,
                  data.frame(p=p,
                             strategy=c("SIR","BIG.SIR1","BIG.SIR2","BIG.SIR3"),
                             runtime=result))
  key <- key + 100
}

p.time$strategy <- as.factor(p.time$strategy)
ggplot(p.time, aes(x=p, y=runtime, col=strategy)) +
  geom_line() +
  geom_point(aes(shape=strategy)) +
  labs(title="running time vs p", x="covariate dimension p", y="running time")


# remove SIR from compare.time (comparison of g doesn't include SIR)
compare.time <- function(key){
  matEDR.block <- function(x){
    bhat <- matrix(x/sqrt((sum(x**2))), ncol=1)
    return(bhat%*%t(bhat))
  }
  Scalable.sir <- function(g, data, size.chunk){
    rows <- ((g-1)*size.chunk+1):(g*size.chunk)
    matEDR.block(edr(data[rows,1],data[rows,-1],H=8,K=1,method="SIR-I")$matEDR[,1])
  }
  runtime <- matrix(0, nrow=50, ncol=3) # 3 methods, 50 reps
  for (m in 1:50){
    xy <- gen.xy(n,p)
    x <- xy$x ; y <- xy$y
    
    # BIG-SIR : loop
    start.time <- Sys.time()
    size.chunk <- n/ng
    mg.hat <- matrix(0, nrow=p, ncol=p)
    for (g in 1:ng){
      temp <- Scalable.sir(g, cbind(y,x), size.chunk)
      mg.hat <- mg.hat + temp
    }
    est <- eigen(mg.hat)$vectors[,1]
    end.time <- Sys.time()
    runtime[m,1] <- end.time - start.time
    
    # BIG-SIR : foreach
    start.time <- Sys.time()
    cl <- makeCluster(4)
    registerDoSNOW(cl)
    size.chunk <- n/ng
    BIGsir <- foreach(g=1:ng, .combine="+")%dopar%{
      require("edrGraphicalTools")
      Scalable.sir(g, cbind(y,x), size.chunk)}
    stopCluster(cl)
    est <- eigen(BIGsir)$vectors[,1]
    end.time <- Sys.time()
    runtime[m,2] <- end.time - start.time
    
    # BIG-SIR : bigmemory + foreach
    start.time <- Sys.time()
    dataYX <- as.big.matrix(cbind(y,x), backingpath="E:/temp",
                            backingfile=paste0(key+m,"test-sir.bin"),
                            descriptorfile=paste0(key+m,"test-sir.desc"), type="double")
    BIGmatdes <- describe(dataYX)
    cl <- makeCluster(4)
    registerDoSNOW(cl)
    x <- attach.big.matrix(BIGmatdes)
    size.chunk <- n/ng
    BIGsir <- foreach(g=1:ng, .combine="+")%dopar%{
      require("edrGraphicalTools")
      require("bigmemory")
      x <- attach.big.matrix(BIGmatdes)
      Scalable.sir(g, x, size.chunk)
    }
    stopCluster(cl)
    est <- eigen(BIGsir)$vectors[,1]
    end.time <- Sys.time()
    runtime[m,3] <- end.time - start.time
  }
  return(colMeans(runtime))
}


# for g
p <- 10 ; n <- 10^7
g.list <- c(10,50,100,150,500)
g.time <- c()
key <- 2000
for (ng in g.list){
  result <- compare.time(key)
  g.time <- rbind(g.time,
                  data.frame(g=ng,
                             strategy=c("BIG.SIR1","BIG.SIR2","BIG.SIR3"),
                             runtime=result))
  key <- key + 100
}

g.time$strategy <- as.factor(g.time$strategy)
ggplot(g.time, aes(x=g, y=runtime, col=strategy)) +
  geom_line() +
  geom_point(aes(shape=strategy)) +
  labs(title="running time vs G", x="number of chunks", y="running time")
