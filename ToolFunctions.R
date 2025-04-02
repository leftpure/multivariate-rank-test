require(clue, quietly=T)
require(energy, quietly = T)
require(randtoolbox, quietly = T)
require(pracma)
require(kernlab)

####------multivariate rank------#### 
compurank=function(x)
{
  m=dim(x)[1]
  dim=dim(x)[2]
  
  gridch=halton(m,dim)
  if(dim==1)
    gridch=matrix((1:m)/m)
  
  distmat=matrix(0,nrow=m,ncol=m)
  for(i in 1:(m))
    distmat[i,]=apply((x[i,]-t(gridch)),2,Norm)^2
  assignmentFUN=solve_LSAP(distmat)
  assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
  return(gridch[assignmentSOL[,2],])
  gridch[assignmentSOL[,2],]
}


####------Generating RMMD Universal Distribution------#### 
# M: X Sample Size
# N: Y Sample Size
# dime: dimension of X and Y
##
dime=2
M=100
N=100
generdist=function(M,N,dime)
{
  niter=5000
  tstat=matrix(0,niter,2)
  for(i in 1:niter)
  {
    x1 <- matrix(runif(dime*M),M)
    y1 <- matrix(runif(dime*N),N) 
    
    z1=as.matrix(compurank(rbind(x1,y1))) 
    
    print(i)
    tstat[i,]=kmmd(as.matrix(z1[1:M,]),as.matrix(z1[(M+1):(M+N),]),ntimes=1)@mmdstats
  }
  return(tstat)
}

####------Generating Rank Energy Universal Distribution-torus------#### 
# M: X Sample Size
# N: Y Sample Size
# dim: dimension of X and Y
# sequence:torus
gensamdist=function(M,N,dim,niter=5000,fixgrid=torus(M+N,dim))
{
  energytstat=numeric(niter)
  for(i in 1:niter)
  {
    ranper=sample(M+N)
    print(i)
    energytstat[i]=eqdist.etest(fixgrid[ranper,],sizes=c(M,N),R=1)$statistic
  }
  return(energytstat)
}

####------Generating Rank Energy Universal Distribution-halton-----#### 
# M: X Sample Size
# N: Y Sample Size
# dim: dimension of X and Y
# sequence:halton
gensamdist_h=function(M,N,dim,niter=5000,fixgrid=halton(M+N,dim))
{
  tstat=numeric(niter)
  for(i in 1:niter)
  {
    ranper=sample(M+N)
    print(i)
    tstat[i]=eqdist.etest(fixgrid[ranper],sizes=c(M,N),R=1)$statistic
  }
  return(tstat)
}


####------Generating MMD Universal Distribution-----#### 
# M: X Sample Size
# N: Y Sample Size
# dim: dimension of X and Y
# sequence:halton
generdist1=function(M,N,dim)
{
  niter=5000
  tstatmmd=matrix(0,niter,1)
  
  for(i in 1:niter)
  {
    x1 <- matrix(runif(dim*M),M)
    y1 <- matrix(runif(dim*N),N) 
    print(i)
    tstatmmd[i]=kmmd(x1,y1,ntimes=150)@mmdstats[2]
  }
  return(tstatmmd)
}

####------Generating Rank Energy Universal Distribution-halton-2-----#### 
# M: X Sample Size
# N: Y Sample Size
# dim: dimension of X and Y
generdistre=function(M,N,dim)
{
  niter=5000
  tstatre=matrix(0,niter,1)
  
  for(i in 1:niter)
  {
    x1 <- matrix(runif(dim*M),M)
    y1 <- matrix(runif(dim*N),N) 
    print(i)
    z1=as.matrix(compurank(rbind(x1,y1))) 
    tstatre[i]=eqdist.etest(as.matrix(z1),sizes=c(M,N),R=1)$statistic#eqdist.etest(as.matrix(z1[1:M,]),as.matrix(z1[(M+1):(M+N),]))
  }
  return(tstatre)
}

####------ Computing Rank Energy Statistic ------#### 
##
computestatistic=function(x,y,m=nrow(x),n=nrow(y),dim=ncol(x),gridch=halton(m+n,dim))
{
  comdata=rbind(x,y)
  distmat=matrix(0,nrow=m+n,ncol=m+n)
  for(i in 1:(m+n))
    distmat[i,]=apply((comdata[i,]-t(gridch)),2,Norm)^2
  assignmentFUN=solve_LSAP(distmat)
  assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
  randenergySTAT=eqdist.etest(gridch[assignmentSOL[,2],],sizes = c(m,n), R=1)
  return(randenergySTAT$statistic)
}

