# In the absence of an actual mesh of points covering region of interest, make one up:
nx = ny = 10
dx = dy = 1
xs = seq(from=1,to=nx,by=dx)
ys = seq(from=1,to=ny,by=dy)
mesh = data.frame(
  x=rep(xs,ny),
  y=rep(ys,rep(nx,ny))
)

# In the absence of actual detectors, make some up:
n = 50 # no. of detectors
set.seed(15) # for reproducibility
# generate random coordinates in region
xmax = max(xs)+dx/2
ymax = max(ys)+dy/2
xmin = dx/2
ymin = dy/2
detectors = data.frame(
  x=runif(n,min=xmin,max=xmax),
  y=runif(n,min=ymin,max=ymax)
)

plot(detectors,pch="*",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
points(mesh,pch="+",col="gray")

# Calculate number of detectors associated with each mesh point
# -------------------------------------------------------------
nmesh = dim(mesh)[1]
ndets = dim(detectors)[1]
# get distances from detectors to each mesh point:
det2meshdist = as.matrix(dist(rbind(detectors,mesh)))[1:ndets,(ndets+1):(ndets+nmesh)]
# find which mesh point is closest:
mesh4det = rep(NA,ndets)
for(i in 1: ndets) mesh4det[i] = which(det2meshdist[i,]==min(det2meshdist[i,]))
ndets.at.mesh = tabulate(mesh4det) # number of detectors at each mesh point
# fill in zeros for any meshes beyond the last with a detector
if(length(ndets.at.mesh)<nmesh) ndets.at.mesh = c(ndets.at.mesh,rep(0,nmesh-length(ndets.at.mesh)))

plot(mesh,pch="+",xlim=c(xmin,xmax), ylim=c(ymin,ymax),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh[i]),col="red")
}

#inverted for corrplot later on
plot(mesh,pch="+",xlim=c(xmin,xmax), xaxt = "n", ylim=rev(c(ymin,ymax)),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh[i]),col="red")
}
axis(3)

# Make a transition probability matrix 
# with value lambda at all transitions within distance dmax, then scale:
# ---------------------------------------------------------------------
# First find the mesh points that are within distance dmax of each other
dmax = sqrt(2)+0.0001
meshmat = as.matrix(mesh)
dists = as.matrix(dist(meshmat)) # distances between all points in mesh
close = which(dists<dmax,arr.ind=TRUE) # indices of close points 

lambda = 0.75 # starting value for transition probabilities
closemat = dists*0 # initialise a matrix to hold binary indicator of being close
for(i in 1:dim(close)[1]) closemat[close[i,1],close[i,2]] = 1 # mark the close 
tpm = closemat*lambda # make the transition probability matrix
# Need rows to add to 1:
tpmsum = apply(tpm,1,sum)
for(i in 1:dim(tpm)[1]){
  tpm[i,] = tpm[i,]/tpmsum[i]
}
tpm # look at it




require(corrplot)
corrplot(tpm,method="ellipse")
corrplot(tpm, method = "ellipse", cl.lim = c(0, 0.5), is.corr = F)

corrplot(closemat,method="number")

###################
#More detectors
##################
# In the absence of actual detectors, make some up:
n = 100 # no. of detectors
set.seed(15) # for reproducibility
# generate random coordinates in region
xmax = max(xs)+dx/2
ymax = max(ys)+dy/2
xmin = dx/2
ymin = dy/2
detectors = data.frame(
  x=runif(n,min=xmin,max=xmax),
  y=runif(n,min=ymin,max=ymax)
)

plot(detectors,pch="*",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
points(mesh,pch="+",col="gray")

# Calculate number of detectors associated with each mesh point
# -------------------------------------------------------------
nmesh = dim(mesh)[1]
ndets = dim(detectors)[1]
# get distances from detectors to each mesh point:
det2meshdist = as.matrix(dist(rbind(detectors,mesh)))[1:ndets,(ndets+1):(ndets+nmesh)]
# find which mesh point is closest:
mesh4det = rep(NA,ndets)
for(i in 1: ndets) mesh4det[i] = which(det2meshdist[i,]==min(det2meshdist[i,]))
ndets.at.mesh3 = tabulate(mesh4det) # number of detectors at each mesh point
# fill in zeros for any meshes beyond the last with a detector
if(length(ndets.at.mesh3)<nmesh) ndets.at.mesh = c(ndets.at.mesh,rep(0,nmesh-length(ndets.at.mesh)))


##########################
#Likelihood
#########################



theta.fn <- function(ndets.at.mesh, pdet){
  #creates det. prob. at each location (state)  
  #pdet is baseline prob of detection
  theta <-rep(NA, length(ndets.at.mesh))
  n <- length(ndets.at.mesh)
  
  for (i in 1:n){
    
    theta[i] <- 1-(1-pdet)^(ndets.at.mesh[i])
  }
  return(theta)
}

theta.fn(ndets.at.mesh = ndets.at.mesh, pdet = 0.5)

create.p.mat <- function(x, N, ndets.at.mesh, pdet){
  
  #creates prob matrix
  #theta is prob of detection at each grid square
  
  theta <- theta.fn(ndets.at.mesh = ndets.at.mesh, pdet = pdet)
  p.mat <- diag(nrow=N)
  
  for (i in 1:N){
    if (sum(x) == 0){
      if (ndets.at.mesh[i] > 0){
        p.mat[i,i] <- 1-theta[i]
      }
      if (ndets.at.mesh[i] == 0){
        p.mat[i,i] <- 1
      }
    }
    if (sum(x) > 0){
      if (ndets.at.mesh[i] > 0 && x[i] == 1){
        p.mat[i,i] <- theta[i]
      } else {
        p.mat[i,i] <- 0
      }
    }
  }
  
  return(p.mat)
}

param.transform <- function(pars, closemat){
  
  #transforms parameters for likelihood
  #and creates transition proability matrix
  
  lambda <- plogis(pars[1])
  pdet <- plogis(pars[2])
  tpm = closemat*lambda
  tpmsum = apply(tpm,1,sum)
  for(i in 1:dim(tpm)[1]){
    tpm[i,] = tpm[i,]/tpmsum
  }
  delta <- solve(t(diag(nrow(tpm)) - tpm + 1), rep(1, nrow(tpm)))
  
  return(list(lambda = lambda, pdet = pdet, tpm = tpm, delta = delta))
}

L <- function(par.vec, x, ndets.at.mesh, closemat){
  
  #likelihood function
  
  N <- ncol(x)
  pars <- param.transform(par.vec, closemat) #transforms parameters
  prob.mat <- create.p.mat(x[1,], N=N, ndets.at.mesh, pdet = pars$pdet) #creates P matrix for observation
  foo <- pars$delta %*% prob.mat
  l <- log(sum(foo)) #scaling to avoid underflow
  phi <- foo/sum(foo)
  
  for (i in 2:nrow(x)){
    prob.mat <- create.p.mat(x[i,], N=N, ndets.at.mesh, pdet = pars$pdet)
    foo <- phi%*%pars$tpm%*%prob.mat
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  
  return(-l)
}

pars <- c(qlogis(0.75), qlogis(0.75))
dat<-dat.sim(200, n = 100, pars, ndets.at.mesh, closemat)
mod <- nlm(L, p = pars, closemat = closemat, x = dat, N = 16, ndets.at.mesh = ndets.at.mesh)
mod <- optim(par = pars, L, closemat = closemat, x = dat$data, ndets.at.mesh = ndets.at.mesh, hessian = T)

mod.HMM <-function(x, pars, closemat, ndets.at.mesh){
  
  #returns estimated parameters, tpm, stationary dist, and 
  #number of detectors at each square
  
  N <- ncol(x$data)
  x <- x$data
  mod <- optim(par = pars, L, closemat = closemat, x = x, 
               ndets.at.mesh = ndets.at.mesh, hessian = T)
  par.vec <- param.transform(mod$par, closemat)
  
  return(list(lambda = par.vec$lambda, pdet = par.vec$pdet, 
              tpm = par.vec$tpm, delta = par.vec$delta, 
              ndets.at.mesh = ndets.at.mesh, hessian = mod$hessian))
}

#simulations

dat.sim <- function(t, n, pars, ndets.at.mesh, closemat){
  
  #simulates data and records true path
  
  pars <- plogis(pars)
  tpm = closemat*pars[1]
  tpmsum = apply(tpm,1,sum)
  for(i in 1:dim(tpm)[1]){
    tpm[i,] = tpm[i,]/tpmsum   #creates tpm
  }
  path <- rep(NA, t)
  state <- sample.int(n, 1) #random starting point
  prob.of.det <- theta.fn(ndets.at.mesh, pars[2])
  data <- matrix(nrow = t, ncol = n)
  data1 <- as.vector(rmultinom(1, 1, tpm[state,]))  #movement to next state
  state <- path[1] <- which(data1 == 1)
  detect <- rbinom(1, 1, prob = prob.of.det[state])  #bernoulli trial to see if detected or not
  data1[which(data1 == 1)] <- detect
  data[1,] <- data1
  
  for (i in 2:t){
    data.vec <- as.vector(rmultinom(1, size = 1, prob = tpm[state,]))
    state <- path[i] <- which(data.vec == 1)
    detect <- rbinom(1, 1, prob = prob.of.det[state])
    data.vec[which(data.vec == 1)] <- detect
    data[i,] <- data.vec
  }
  
  return(list(data = data, path = path))
}

set.seed(101976)
dat<-dat.sim(200, n= 100, pars = pars, ndets.at.mesh, closemat)
mod1 <- mod.HMM(dat, pars, closemat, ndets.at.mesh)
#lambda = 0.773355
#pdet = 0.76181998

bootstrap <- function(N, n, nobs, pars, ndets.at.mesh, closemat){
  
  #bootstraps model
  #generates data and fits model N times
  #creates CI through quantile
  #takes mean of bootstraps to help estimate bias
  #takes sd of parameter estimates
  #returns all parameter estimates from bootstraps
  
  pdet <- pars[2]
  lambda <- pars[1]
  pdet.ests <- rep(NA, N)
  lambda.ests <- rep(NA, N)
  
  for (i in 1:N){
    
    dat <- dat.sim(nobs, n = n, pars = pars, ndets.at.mesh = ndets.at.mesh,
                   closemat = closemat)
    mod <- mod.HMM(x = dat, pars = pars, closemat = closemat, 
                   ndets.at.mesh = ndets.at.mesh)
    pdet.ests[i] <- mod$pdet
    lambda.ests[i] <- mod$lambda
    
  }
  
  pdet.B <- quantile(pdet.ests, probs = c(0.025, 0.975))
  lambda.B <- quantile(lambda.ests, probs = c(0.025, 0.975))
  lambda.mean <- mean(lambda.ests)
  pdet.mean <- mean(pdet.ests)
  sd.lambda <- sd(lambda.ests)
  sd.pdet <- sd(pdet.ests)
  
  return(list(CI.B.pdet = pdet.B, CI.B.lambda = lambda.B, mean.lambda = lambda.mean,
              mean.pdet = pdet.mean, sd.lambda = sd.lambda,
              sd.pdet = sd.pdet, pdet.ests = pdet.ests, lambda.ests = lambda.ests))
}

dat3 <- dat.sim(1000, 100, pars, ndets.at.mesh, closemat)
mod3 <- mod.HMM(dat3, pars, closemat, ndets.at.mesh)

dat4 <- dat.sim(200, 100, pars, ndets.at.mesh3, closemat)
mod4 <- mod.HMM(dat4, pars, closemat, ndets.at.mesh3)

dat8 <- dat.sim(1000, 100, pars, ndets.at.mesh3, closemat)
mod8 <- mod.HMM(dat8, pars, closemat, ndets.at.mesh3)

dat5 <- dat.sim(1000, 100, pars, ndets.at.mesh2, closemat)
mod5 <- mod.HMM(dat5, pars, closemat, ndets.at.mesh2)

dat6 <- dat.sim(200, 100, pars, ndets.at.mesh4, closemat)
mod6 <- mod.HMM(dat6, pars, closemat, ndets.at.mesh4)

dat7 <- dat.sim(1000, 100, pars, ndets.at.mesh4, closemat)
mod7 <- mod.HMM(dat7, pars, closemat, ndets.at.mesh4)


set.seed(101976)
b1 <- bootstrap(100, 100, 200, pars, ndets.at.mesh, closemat)
set.seed(101986)
b2 <- bootstrap(200, 100, 1000, pars, ndets.at.mesh, closemat)
set.seed(101976)
b3 <- bootstrap(200, 100, 200, pars, ndets.at.mesh2, closemat)
set.seed(101986)
b4 <- bootstrap(200, 100, 1000, pars, ndets.at.mesh2, closemat)

hist(b1$lambda.ests, freq = F, ylim = c(0, 15), main ="", xlab = "lambda")
lines(density(b1$lambda.ests))
hist(b1$pdet.ests, freq = F, ylim = c(0, 10), main ="", xlab = "p")
lines(density(b1$pdet.ests))


set.seed(101976)
b1. <- bootstrap(500, 100, 200, pars, ndets.at.mesh, closemat)
set.seed(101986)
b2. <- bootstrap(500, 100, 1000, pars, ndets.at.mesh, closemat)
set.seed(101976)
b3. <- bootstrap(500, 100, 200, pars, ndets.at.mesh2, closemat)
set.seed(101986)
b4. <- bootstrap(500, 100, 1000, pars, ndets.at.mesh2, closemat)

set.seed(1017)
dat3<-dat.sim(200, n= 100, pars = pars, ndets.at.mesh3, closemat)
mod3 <- mod.HMM(dat3, pars, closemat, ndets.at.mesh3)

set.seed(10176)
b1_ <- bootstrap(500, 100, 200, pars, ndets.at.mesh3, closemat)
set.seed(10186)
b2_ <- bootstrap(500, 100, 1000, pars, ndets.at.mesh3, closemat)
set.seed(10176)
b3_ <- bootstrap(500, 100, 200, pars, ndets.at.mesh4, closemat)
set.seed(10186)
b4_ <- bootstrap(500, 100, 1000, pars, ndets.at.mesh4, closemat)


###################################
#Decoding and Prediction
###################################

HMM.forward <- function(x, mod){
  
  #computes log forward probabilities
  
  N <- ncol(x)
  t <- nrow(x)
  lalpha <- matrix(NA, N, t)
  P <- diag(create.p.mat(x[1,], N, mod$ndets.at.mesh, mod$pdet))
  foo <- mod$delta*P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- lscale + log(foo)
  
  for (i in 2:t){
    P <- diag(create.p.mat(x[i,], N, mod$ndets.at.mesh, mod$pdet))
    foo <- foo%*%mod$tpm*P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  
  return(lalpha)
}

a.f <- HMM.forward(dat1, mod1)


HMM.backward <- function(x, mod){
  
  #computes log backward probabilities
  
  N <- ncol(x)
  t <- nrow(x)
  lbeta <- matrix(NA, N, t)
  lbeta[,t] <- rep(0, N)
  foo <- rep(1/N, N)
  lscale <- log(N)
  
  for (i in (t-1):1){
    P <- diag(create.p.mat(x[i+1,], N, mod$ndets.at.mesh, mod$pdet))
    P.foo <- P*foo
    foo <- mod$tpm%*%P.foo
    lbeta[,i] <- log(foo) +lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  return(lbeta)
}
l.b <- HMM.backward(dat, mod1)


HMM.conditional <- function(x, mod){
  
  N <- ncol(x)
  t <- nrow(x)
  dxc <- matrix(NA, nrow = t, ncol = t)
  Px <- matrix(NA, nrow = N, ncol = t)

  
  for (j in 1:t){
    Px[,j] <- diag(create.p.mat(x[j,], N, mod$ndets.at.mesh, mod$pdet))
  }
  la <- HMM.forward(x, mod)
  lb <- HMM.backward(x, mod)
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  
  for (i in 1:t){
    foo <- (exp(la[,i]-lafact[i])%*%mod$tpm)*exp(lb[,i]-lbfact[i])
    foo <- foo/sum(foo)
    dxc[,i] <- foo%*%Px
  }
  
  return(dxc)
}


state.probs <- function(x, mod){
  
  #calculates probability of being in state i at time t
  #given our observations
  
  t <- nrow(x)
  N <- ncol(x)
  la <- HMM.forward(x, mod)
  lb <- HMM.backward(x, mod)
  c <- max(la[,N])
  llk <- c + log(sum(exp(la[,N] - c)))
  stateprobs <- matrix(NA, ncol = t, nrow = N)
  
  for (i in 1:t){
    stateprobs[,i] <- exp(la[,i] + lb[,i] - llk)
  }
  
  return(stateprobs)
}



local.decoding <- function(x, mod){
  
  #carries out local decoding on a set of observations given a model
  
  t <- nrow(x)
  stateprobs <- state.probs(x, mod)
  ild <- rep(NA, t)
  
  for (i in 1:t){
    ild[i] <- which.max(stateprobs[,i])
  }
  
  ild
}
ld <- local.decoding(dat, mod1)

state.pred <- function(h = 1, x, mod){
  
  #carries out state predications for time t+h given data and a model
  
  t <- nrow(x$data)
  N <- ncol(x$data)
  la <- HMM.forward(x$data, mod)
  c <- max(la[,t])
  llk <- c + log(sum(exp(la[,t] - c)))
  statepreds <- matrix(NA, ncol = h, nrow = N)
  foo <- exp(la[,t] - llk)
  
  for (i in 1:h){
    foo <- foo%*%mod$tpm
    statepreds[,i] <- foo
  }
  
  return(statepreds)
}

mean.group.pred <- function(S, h, x){
  
  #x is to be given as a list of the data for each individual
  #mod is a list of the fitted model for each individual
  #function calculates prediction for each individual h steps in future
  #then averages these for each state
  
  dat <- x[[1]]
  mod <- x[[2]]
  
  group <- matrix(NA, nrow = 100, ncol = S)
  
  for (i in 1:S){
    ind <-state.pred(h, dat[[i]], mod[[i]])
    group[,i] <- ind[,h]
  }
  
  mean.pred <- apply(group, 1, mean)
  
  return(mean.pred)
}

dat.gen <- function(N, t, n, pars, ndets.at.mesh, closemat){
  
  data <- list()
  mods <- list()
  for(i in 1:N){
    data[[i]] <- dat.sim(t, n, pars, ndets.at.mesh, closemat)
    mods[[i]] <- mod.HMM(data[[i]], pars, closemat, ndets.at.mesh)
  }
  
  return(list(data = data, mods = mods))
}
y <- dat.gen(100, 200, 100, pars, ndets.at.mesh, closemat)
h <- 4
mp <- mean.group.pred(100, 1, y)
mp. <- matrix(mp, nrow = 10, ncol = 10, byrow = T)
corrplot(mp., method = "color", cl.lim = c(0, max(mp.)), is.corr = F,
         xlab = "State", ylab = "probability", mar = c(1,0,2,0))

p <- state.pred(4, dat, mod1)
h <- c(1,2,3,4)

library(corrplot)

pdf("prednew.pdf",w =10, h = 9)
par(oma=c(1,0,2,0), mfrow = c(2,2),las = 1)
for (i in 1:4){
  ps <- p[,i]
  pmat <- matrix(ps, nrow = 10, ncol = 10, byrow = T)
  corrplot(pmat, method = "color", cl.lim = c(0, max(pmat)), is.corr = F, 
           main=paste("State predictions for h = ", h[i]),
           xlab = "State", ylab = "probability", mar = c(1,0,2,0))
}
dev.off()

pmat1 <- matrix(p[,1], nrow = 10, ncol = 10, byrow = T)
corrplot(pmat1, method = "color", cl.lim = c(0, 0.3), is.corr = F)

pdf("preds50new.pdf")
sp <- state.pred(500, dat, mod1)
h <- 50
ps50 <- sp[,50]
pmat50 <- matrix(ps50, nrow = 10, ncol = 10, byrow = T)
corrplot(pmat50, method = "color", cl.lim = c(0, max(pmat50)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()

sp2 <- state.pred(50, dat, mod1)
h <- 50
ps250 <- sp2[,50]
plot(1:100, ps250, type="h", main="State predictions for h = 50", 
     xlab = "State", ylab = "probability", lwd = 3)
pmat250 <- matrix(ps250, nrow = 10, ncol = 10, byrow = T)
corrplot(pmat250, method = "color", cl.lim = c(0, max(pmat250)), is.corr = F,
         xlab = "State", ylab = "probability")
title("State predictions for h=50", line=3)

pdf("stat_new2.pdf", w=8,h=8)
lambda <- mod1$lambda
pdet <- mod1$pdet
delta <- solve(t(diag(100)-mod1$tpm+1), rep(1, 100))
plot(1:100, delta, type="h", main="Stationary distribution",
     xlab = "State", ylab = "probability", lwd = 3)
pmatd <- matrix(delta, nrow = 10, ncol = 10, byrow = T)
corrplot(pmatd, method = "color", cl.lim = c(0, 0.02), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()


viterbi <- function(x, mod){
  
  #vitberi algorithm
  #carries out global decoding
  #returns vector of most likely state sequence
  
  t <- nrow(x)
  N <- ncol(x)
  xi <- matrix(0, t, N)
  P <- create.p.mat(x[1,], N, mod$ndets.at.mesh, mod$pdet)
  foo <- mod$delta*diag(P)
  xi[1,] <- foo/sum(foo)
  
  for (i in 2:t){
    P <- create.p.mat(x[i,], N, mod$ndets.at.mesh, mod$pdet)
    foo <- apply(xi[i-1,]*mod$tpm, 2, max)*diag(P)
    xi[i,] <- foo/sum(foo)
  }
  
  iv <- numeric(t)
  iv[t] <- which.max(xi[t,])
  
  for(i in (t-1):1){
    iv[i] <- which.max(mod$tpm[,iv[i+1]]*xi[i,])
  }
  
  return(iv)
}

g.transform <- function(v){
  
  #this function just transforms state sequence into a mapping
  #so that it can be plotted with the igraph package
  #igraph requires 
  #so a sequence of 1,2,3,4 becomes
  #1->2, 2->3, 3->4
  
  n <- length(v)
  v.star <- rep(NA, 2*(n-1))
  v.star[1] <- v[1]
  v.star[2*(n-1)] <- v[n]
  
  for (i in 2:(n-1)){
    x <- rep(v[i], 2)
    v.star[2*(i-1)] <- x[1]
    v.star[2*(i-1)+1] <- x[2]
  }
  
  return(v.star)
}

library(igraph)

pdf("comparison.pdf", w=8, h=4)
par(mfrow=c(1,2))
v <- viterbi2(dat$data, mod1)
v1 <- g.transform(v)
g <- graph(edges = v1, n = ncol(dat$data), directed = T)
a <- permute.vertices(g, c(1:100))
plot(a, layout= layout_on_grid(a), edge.arrow.size = 0.1, vertex.label = NA, 
     main = "Global Decoding")
#ld1 <- g.transform(ld)
#g1 <- graph(edges = ld1, n = ncol(dat$data), directed = T)
#a1 <- permute.vertices(g1, c(100:1))
#plot(a1, layout= layout_on_grid(a), edge.arrow.size = 0.1, vertex.label = NA)

v2 <- g.transform(dat$path)
g <- graph(edges = v2, n = ncol(dat$data), directed = T)
a <- permute.vertices(g, c(1:100))
plot(a, layout= layout_on_grid(a), edge.arrow.size = 0.1, vertex.label = NA,
     main = "True Path")
dev.off()

comparison <- function(vit, path){
  
  #looks at the proportion of time viterbi gets the most 
  #likely path right, by comparing to the true simulated path
  
  TF <- path == vit
  TS <- sum(TF, na.rm = T)
  prop <- TS/length(path)
  
  return(prop)
}

comp <- comparison(v, dat$path) #0.53

comp.sim <- function(N, mod, closemat){
  
  #looks at proportion of time global decoding is right compared
  #to true simulated path for multiple individuals
  
  c <- rep(NA, N)
  pars <- c(mod$lambda, mod$pdet)
  
  for (i in 1:N){
    dat <- dat.sim(200, 100, pars, mod$ndets.at.mesh, closemat)
    v <- viterbi(dat$data, mod)
    c[i] <- comparison(v, dat$path)
  }
  
  return(c)
}

c <- comp.sim(1000, mod1, closemat)
hist(c, col ="blue", xlab = "Proportion of correct decoded states", main = "")

comparison.r <- function(vit, path, radius = 1){
  
  #looks at proportion of time that global decoding was correct
  #when compared to true simulated path within radius of 
  # 1 to 4 states
  
  diff <- path - vit
  n <- length(path)
  close <- rep(NA, n)
  in.r <- c(9,10,11,-1,1,0,-9,-10,-11)
  rad <- c(radius*in.r, (radius-2)*in.r, (radius-1)*in.r, in.r)
  
  for (i in 1:n){
    if (diff[i] %in% rad){
      close[i] <- 1
    }
    else {
      close[i] <- 0
    }
  }
  
  close <- sum(close)
  prop <- close/n
  
  return(prop)
}

comp.r <- comparison.r(v, dat$path, 1) #0.69
comp.r2 <- comparison.r(v, dat$path, 2) #0.765
comp.r3 <- comparison.r(v, dat$path, 3)#0.79
comp.r4 <- comparison.r(v, dat$path, 4) #0.815

comp.bs <- function(N, r, mod, closemat){
  
  #looks at proportion of time that global decoding was correct
  #when compared to true simulated path within radius of 
  # 1 to 4 states for multiple individuals
  
  comp <- rep(NA, N)
  pars <- c(mod$lambda, mod$pdet)
  
  for (i in 1:N){
    data <- dat.sim(200, 100, pars, mod$ndets.at.mesh, closemat)
    v <- viterbi(data$data, mod)
    comp[i] <- comparison.r(v, data$path, r)
  }
  return(comp)
}

bscomp <- comp.bs(1000, 1, mod1, closemat)
hist(bscomp, col = "red", xlab = "Proportion of correct decoded states within a radius of one state",
     main = "")
bs.comp2 <- comp.bs(1000, 2, mod1, closemat)
hist(bs.comp2, col = "green", xlab = "Proportion of correct decoded states within a radius of two states",
     main = "")

usage <- function(N, pars, ndets.at.mesh, closemat){
  
  #looks at proportion of time spent in each state in the global
  #decoding sequence for one individual
  
  dat<-dat.sim(N, n= 100, pars = pars, ndets.at.mesh, closemat)
  mod <- mod.HMM(dat, pars, closemat, ndets.at.mesh)
  v <- viterbi(dat$data, mod)
  t <- table(factor(v, levels = c(1:100)))
  u <- t/N
  p <- table(factor(dat$path, levels = c(1:100))) 
  path.u <- p/N
  
  return(list(u=u, path.u=path.u))
}

par(mfrow=c(1,2))
u <- usage(2000, pars, ndets.at.mesh, closemat)
u. <- matrix(u[[1]], nrow = 10, ncol = 10, byrow = T)
corrplot(u., method = "color", cl.lim = c(0, max(u.)), is.corr = F,
         xlab = "State", ylab = "probability")
p. <- matrix(u[[2]], nrow = 10, ncol = 10, byrow = T)
corrplot(p., method = "color", cl.lim = c(0, max(p.)), is.corr = F,
         xlab = "State", ylab = "probability")

usage.group <- function(S, N, n, pars, ndets.at.mesh, closemat){
  
  #looks at proportion of time spent in each state in the global
  #decoding sequence and for the true path for several individuals
  
  
  u <- matrix(NA, nrow = S, ncol = n)
  p <- matrix(NA, nrow = S, ncol = n)
  mean.u <- rep(NA, n)
  mean.p <- rep(NA, n)
  
  for (i in 1:S){
    use <- usage(N, pars, ndets.at.mesh, closemat)
    u[i,] <- use[[1]]
    p[i,] <- use[[2]]
  }
  
  for (j in 1:n){
    mean.u[j] <- mean(u[,j])
    mean.p[j] <- mean(p[,j])
    
  }
  
  return(list(mean.use = mean.u, mean.path = mean.p))
}

usage.sim <- usage.group(50, 2000, 100, pars, ndets.at.mesh, closemat)
pdf("pathvsusegroup.pdf", w=10, h=5)
par(mfrow=c(1,2))
us. <- matrix(usage.sim[[1]], nrow = 10, ncol = 10, byrow = T)
corrplot(us., method = "color", cl.lim = c(0, max(us.)), is.corr = F,
         xlab = "State", ylab = "probability")
use. <- matrix(usage.sim[[2]], nrow = 10, ncol = 10, byrow = T)
corrplot(use., method = "color", cl.lim = c(0, max(use.)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()

pdf("avguse.pdf", w=4, h=4)
corrplot(us., method = "color", cl.lim = c(0, max(us.)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()

pdf("avgusedets.pdf", w=8, h=4)
par(mfrow=c(1,2))
corrplot(us., method = "color", cl.lim = c(0, max(us.)), is.corr = F,
         xlab = "State", ylab = "probability")
plot(mesh,pch="+",xlim=c(xmin,xmax), xaxt = "n", ylim=rev(c(ymin,ymax)),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh[i]),col="red")
}
axis(3)
dev.off()

d <- dat.sim(10000, 100, pars, ndets.at.mesh, closemat)
z <- table(factor(d$path, levels = c(1:100)))
use <- z/length(d$path)
use. <- matrix(use, nrow = 10, ncol = 10, byrow = T)
corrplot(use., method = "color", cl.lim = c(0, max(use.)), is.corr = F,
         xlab = "State", ylab = "probability")
use_ <- matrix(use.v, nrow = 10, ncol = 10, byrow = T)
use.v <- usage(10000, pars, ndets.at.mesh, closemat)
corrplot(use_, method = "color", cl.lim = c(0, max(use_)), is.corr = F,
         xlab = "State", ylab = "probability")


#------------------------EDGE--------------------------
layout <- c(99, 99, 96, 96, 92, 23, 21, 61, 40, 70, 82, 48, 59, 42, 40, 21, 62, 80, 95, 89,
            91,91,91,91,91,10,10,10,10,10,93,93,93,93,93,94,94,94,82,82,83,
            30,30,50,50,84,84,84,61,61) #50 detectors
layout2 <- c(layout, rep(61, 5), rep(98, 4), 97, rep(11, 3), rep(51,2), rep(41,5),
             rep(86, 3), rep(87, 3), rep(88,3), 1, rep(12,2), rep(100,2), rep(90,2),
             rep(60, 2), rep(20, 4), 72, 68, 49, 83, rep(42, 2), rep(39,2)) #100 detectors
ndets.at.mesh2 = tabulate(layout)
ndets.at.mesh2[100] <- 0  # number of detectors at each mesh point
plot(mesh,pch="+",xlim=c(xmin,xmax),ylim=c(ymin,ymax),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh2[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh2[i]),col="red")
}
ndets.at.mesh4 = tabulate(layout2)
ndets.at.mesh4[100] <- 0  # number of detectors at each mesh point
plot(mesh,pch="+",xlim=c(xmin,xmax),ylim=c(ymin,ymax),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh4[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh4[i]),col="red")
}


plot(mesh,pch="+",xlim=c(xmin,xmax), xaxt = "n", ylim=rev(c(ymin,ymax)),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh2[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh2[i]),col="red")
}
axis(3)

library(igraph)

pdf("compedge.pdf", w=8, h=4)
par(mfrow=c(1,2))
pars <- c(qlogis(0.75), qlogis(0.75))
set.seed(9283)
dat2<-dat.sim(200, n= 100, pars = pars, ndets.at.mesh2, closemat)
mod2 <- mod.HMM(dat2, pars, closemat, ndets.at.mesh2)
v <- viterbi(dat2$data, mod2)
v2 <- g.transform(v)
g2 <- graph(edges = v2, n = ncol(dat2$data), directed = T)
a2 <- permute.vertices(g2, c(1:100))
plot(a2, layout= layout_on_grid(a2), edge.arrow.size = 0.2, vertex.label = NA,
     main = "Global Decoding Sequence")

v2 <- g.transform(dat2$path)
g <- graph(edges = v2, n = ncol(dat2$data), directed = T)
a <- permute.vertices(g, c(1:100))
plot(a, layout= layout_on_grid(a), edge.arrow.size = 0.1, vertex.label = NA,
     main = "True Path Sequence")
dev.off()

pars <- c(qlogis(0.75), qlogis(0.75))
dat3<-dat.sim(2000, n= 100, pars = pars, ndets.at.mesh2, closemat)
mod3 <- mod.HMM(dat3, pars, closemat, ndets.at.mesh2)
v3 <- viterbi(dat3, mod3)
v3 <- g.transform(v3)
g3 <- graph(edges = v3, n = ncol(dat3), directed = T)
a3 <- permute.vertices(g3, c(1:100))
plot(a3, layout= layout_on_grid(a3), edge.arrow.size = 0.1, vertex.label = NA)

comp.r_ <- comparison.r(v, dat2$path, r =1) #0.475
comp.r2_ <- comparison.r(v, dat2$path, r =2) #0.535
comp.r3_ <- comparison.r(v, dat2$path, r=3)
comp.r4_ <- comparison.r(v, dat2$path, r=4)
compr <- comparison(v, dat2$path) #0.35


bscomp_ <- comp.bs(1000, 1, mod2, closemat)
bs.comp2_ <- comp.bs(1000, 2, mod2, closemat)
hist(bscomp_)
hist(bscomp_, col = "red", xlab = "Proportion of correct decoded states within a radius of two states",
     main = "")
hist(bs.comp2_)
hist(bs.comp2_, col = "green", xlab = "Proportion of correct decoded states within a radius of two states",
     main = "")

u2 <- usage(2000, pars, ndets.at.mesh2, closemat)
u2. <- matrix(u2[[1]], nrow = 10, ncol = 10, byrow = T)
pdf("usegraphdecode2", w=8, h=4)
par(mfrow=c(1,2))
corrplot(u2., method = "color", cl.lim = c(0, max(u2.)), is.corr = F,
         xlab = "State", ylab = "probability")
plot(a3, layout= layout_on_grid(a3), edge.arrow.size = 0.1, vertex.label = NA)
dev.off()

pdf("usegraph2", w=8, h=4)
par(mfrow=c(1,2))
corrplot(u2., method = "color", cl.lim = c(0, max(u2.)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()

pdf("usegraphpath2", w=8, h=4)
par(mfrow=c(1,2))
corrplot(u2., method = "color", cl.lim = c(0, max(u2.)), is.corr = F,
         xlab = "State", ylab = "probability")
u2_ <- matrix(u2[[2]], nrow = 10, ncol = 10, byrow = T)
corrplot(use., method = "color", cl.lim = c(0, max(u2_)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()


usage.sim2 <- usage.group(50, 2000, 100, pars, ndets.at.mesh2, closemat)
us2. <- matrix(usage.sim2[[1]], nrow = 10, ncol = 10, byrow = T)
pdf("avguse2.pdf", w=4, h=4)
corrplot(us2., method = "color", cl.lim = c(0, max(us2.)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()

usage.sim2 <- usage.group(50, 2000, 100, pars, ndets.at.mesh2, closemat)
pdf("pathvsusegroupavg2.pdf", w=10, h=5)
par(mfrow=c(1,2))
us2. <- matrix(usage.sim[[1]], nrow = 10, ncol = 10, byrow = T)
corrplot(us2., method = "color", cl.lim = c(0, max(us2.)), is.corr = F,
         xlab = "State", ylab = "probability")
use2. <- matrix(usage.sim2[[2]], nrow = 10, ncol = 10, byrow = T)
corrplot(use2., method = "color", cl.lim = c(0, max(use2.)), is.corr = F,
         xlab = "State", ylab = "probability")
dev.off()

y2 <- dat.gen(100, 200, 100, pars, ndets.at.mesh2, closemat)
h <- 1
mp2 <- mean.group.pred(100, 1, y2)
mp2. <- matrix(mp2, nrow = 10, ncol = 10, byrow = T)
corrplot(mp2., method = "color", cl.lim = c(0, max(mp2.)), is.corr = F,
         xlab = "State", ylab = "probability")


pdf("avguse2det.pdf", w=8, h=4)
par(mfrow=c(1,2))
corrplot(us2., method = "color", cl.lim = c(0, max(us2.)), is.corr = F,
         xlab = "State", ylab = "probability")
plot(mesh,pch="+",xlim=c(xmin,xmax), xaxt = "n", ylim=rev(c(ymin,ymax)),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh2[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh2[i]),col="red")
}
axis(3)
dev.off()

#####################################
#GOF
#####################################

matsplitter<-function(M, r, c){
  
  #divides a matrix up into smaller matrices with r rows and c columns
  
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  
  return(cv)
} 

exp.fn <- function(dat, mod){
  
  #creates expected matrix
  
  vit <- viterbi(dat, mod) #global decoding sequence for data and model
  n <- length(vit)
  exp.mat <- matrix(NA, nrow = n, ncol = ncol(dat))
  
  for (i in 1:n){
    
    state <- vit[i]
    pmat <- create.p.mat(dat[i,], N = nrow(mod$tpm), ndets.at.mesh = mod$ndets.at.mesh,
                         pdet = mod$pdet) 
    exp.mat[i,] <- pmat[state,]
    
  }
  
  return(exp.mat)
}

gof <- function(data, mod, B){
  
  #calculates gof statistic for data under model, mod, and splits it into B bins
  
  bin.length <- nrow(data$data)/B
  mats <- matsplitter(data$data, bin.length, ncol(data$data))
  exp.mat <- exp.fn(data$data, mod)
  exp.bins <- matsplitter(exp.mat, bin.length, ncol(data$data))
  stat.vec <- rep(NA, B)

  
  for (b in 1:B){
      obs <- colSums(mats[,,b])
      exp <- colSums(exp.bins[,,b])
      stat.vec[b] <- sum(((obs - exp)^2))
  }
  stat <- sum(stat.vec)
  
  return(stat)
}


sim.gof <- function(mod, B, N, n, closemat){
  
  #creates gof statistic for multiple individuals under same model
  
  sim <- rep(NA, N)
  pars <- c(mod$lambda, mod$pdet)

  for (i in 1:N){
    data <- dat.sim(200, n, pars, mod$ndets.at.mesh, closemat)
    sim[i] <- gof(data, mod, B)
  }
  return(sim)
}
sim<-sim.gof(mod1, 10, 500, 100, closemat)
sim1 <- sim.gof(mod2, 10, 500, 100, closemat)
sim2 <- sim.gof(mod3, 100, 500, 100, closemat)
sim3 <- sim.gof(mod4, 10, 500, 100, closemat)
sim4 <- sim.gof(mod5, 100, 500, 100, closemat)
sim5 <- sim.gof(mod6, 10, 500, 100, closemat)
sim6 <- sim.gof(mod7, 100, 500, 100, closemat)
sim7 <- sim.gof(mod8, 100, 500, 100, closemat)


hist(sim, main = "Distibution of test statistic A", xlab = "A", col = "blue")
hist(sim1, main = "Distibution of test statistic A", xlab = "A", col = "blue")

cr <- quantile(sim, 0.95)  #1153.116
cr1 <- quantile(sim1, 0.95) #2313.63
cr2 <- quantile(sim2, 0.95) #1186.661
cr3 <- quantile(sim3, 0.95) #702.0949
cr4 <- quantile(sim4, 0.95) #2277.788
cr5 <- quantile(sim5, 0.95) #1795.939
cr6 <- quantile(sim6, 0.95) #1708.28
cr7 <- quantile(sim7, 0.95) #713.1083


stat <- gof(dat, mod1, 10) #579.1252
stat2 <- gof(dat2, mod2, 10) #2059.049
stat3 <- gof(dat3, mod3, 10)
stat4 <- gof(dat4, mod4, 10) #604.5059
stat5 <- gof(dat5, mod5, 10) 
stat6 <- gof(dat6, mod6, 10) #124.0035
stat7 <- gof(dat7, mod7, 10)
stat8 <- gof(dat8, mod8, 10)


A <- gof(dat, mod1, 10)
sim<-sim.gof(mod1, 10, 1000, 100, closemat = closemat)
hist(sim7, main = "Distibution of test statistic A", xlab = "A", col = "blue")
cr <- quantile(sim, 0.95)

pars <-c(mod2$lambda, mod2$pdet)
sim<-sim.gof(mod2, 10, 500, 100, pars, closemat)
hist(sim, main = "Distibution of test statistic A", xlab = "A", col = "blue")
cr <- quantile(sim, 0.95)

#-------------------more detectors----------------
nx = ny = 10
dx = dy = 1
xs = seq(from=1,to=nx,by=dx)
ys = seq(from=1,to=ny,by=dy)
mesh = data.frame(
  x=rep(xs,ny),
  y=rep(ys,rep(nx,ny))
)

# In the absence of actual detectors, make some up:
n = 50 # no. of detectors
set.seed(15) # for reproducibility
# generate random coordinates in region
xmax = max(xs)+dx/2
ymax = max(ys)+dy/2
xmin = dx/2
ymin = dy/2
detectors = data.frame(
  x=runif(n,min=xmin,max=xmax),
  y=runif(n,min=ymin,max=ymax)
)

plot(detectors,pch="*",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
points(mesh,pch="+",col="gray")

# Calculate number of detectors associated with each mesh point
# -------------------------------------------------------------
nmesh = dim(mesh)[1]
ndets = dim(detectors)[1]
# get distances from detectors to each mesh point:
det2meshdist = as.matrix(dist(rbind(detectors,mesh)))[1:ndets,(ndets+1):(ndets+nmesh)]
# find which mesh point is closest:
mesh4det = rep(NA,ndets)
for(i in 1: ndets) mesh4det[i] = which(det2meshdist[i,]==min(det2meshdist[i,]))
ndets.at.mesh = tabulate(mesh4det) # number of detectors at each mesh point
# fill in zeros for any meshes beyond the last with a detector
if(length(ndets.at.mesh)<nmesh) ndets.at.mesh = c(ndets.at.mesh,rep(0,nmesh-length(ndets.at.mesh)))

plot(mesh,pch="+",xlim=c(xmin,xmax), ylim=c(ymin,ymax),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh[i]),col="red")
}

#inverted plot
plot(mesh,pch="+",xlim=c(xmin,xmax), xaxt = "n", ylim=rev(c(ymin,ymax)),col="gray")
for(i in 1:nmesh) {
  if(ndets.at.mesh[i]>0) 
    text(mesh$x[i],mesh$y[i],labels=as.character(ndets.at.mesh[i]),col="red")
}
axis(3)

