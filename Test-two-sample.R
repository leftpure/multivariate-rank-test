##------------------------------------------------------
require(clue, quietly=T)
require(energy, quietly = T)
require(randtoolbox, quietly = T)
require(pracma, quietly = T)
require(kernlab, quietly = T)
require(crossmatch, quietly = T)
require(HHG, quietly = T)
require(gTests, quietly = T)
require(ramify)
require(ggplot2)
require(mlbench)
require(dprep)
library(mlbench)
library(export)
source("ToolFunctions.R")


####-----------------generating critical value-----------------###
# M: X Sample Size
# N: Y Sample Size
# dim: dimension of X and Y
M=100
N=100
dim=2

#RE
energytstat=gensamdist(M,N,dim)
par(mar = c(2, 2, 2, 2)) 
hist(energytstat,breaks = 35)
#RMMD
tstatest=generdist(M,N,3)
hist(tstatest[,2])
#MMD
tstatmmd2=generdist1(M,N,2)

#0.9、0.95、0.99 quantiles
quantile_value0 <- quantile(tstat[,2], 0.9)[[1]]
quantile_value1 <- quantile(tstat[,2], 0.95)[[1]]
quantile_value2 <- quantile(tstat[,2], 0.99)[[1]]

desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")

####-----------------simulation-----------------###
for(p in c(6))
{ 
  energytstat=gensamdist(100,100,p)
  tstat=generdist(100,100,p)
  tstatmmd2=generdist1(100,100,p)
  for(c in c(0.4))
  {
    perfvec=matrix(0,4,1003)
    for(j in 1:1000)
    {
      #setting
      p=25
      n_values=25
      beta_values <- seq(0, 1, by = 0.1)
      rho <- 0.5
      
      Sigma <- matrix(rho^abs(outer(1:p, 1:p, "-")), nrow = p, ncol = p)

      ZX <- matrix(rnorm(p * max(n_values)), nrow = max(n_values), ncol = p)
      ZY <- matrix(rnorm(p * max(n_values)), nrow = max(n_values), ncol = p)
      
      for (n in n_values) {
        for (m in n_values) {
          for (beta in beta_values) {
            X <- as.matrix(sqrtm(Sigma))%*% ZX
            Y <- sqrtm(Sigma) %*% ZY[1:m, , drop = FALSE] + 0.15 * beta * rep(1, p)
            print(dim(X))
            print(dim(Y))
          }
        }
      }
      ###--local alternative--###
      #uniform
      # setdata1 <- matrix(runif((p*100), min = -1, max = 1),100)
      # setdata2 <- matrix(runif((p*100), min = -1, max = 1+c),100) #H0
      #Cauchy setup independent
      # setdata1=array(0,dim=c(100,p))
      # setdata2=array(0,dim=c(100,p))
      # for(i in 1:p)
      # {
      #   #prep=rnorm(n,0,1)
      #   setdata1[,i]=rcauchy(100,0,0.5)
      #   setdata2[,i]=rcauchy(100,0+c,0.5)
      # }

      mat3=matrix(setdata1[,], nrow = 100, ncol = p)
      mat4=matrix(setdata2[,], nrow = 100, ncol = p)
      
      #Rank Energy
      T1re2=computestatistic(mat3,mat4)
      rsetdata=as.matrix(compurank(rbind(mat3,mat4))) 
      T1re=eqdist.etest(as.matrix(rsetdata),sizes=c(100,100),R=1)$statistic
      #T1re=computestatistic(mat3,mat4)
      perfvec[1,j]=length(which(energytstat>T1re))/5000
      
      #Energy
      m1=nrow(mat3)
      n1=nrow(mat4)
      comdata=dist(rbind(mat3,mat4))
      T2energy=eqdist.etest(comdata,sizes = c(m1,n1), R=199)$p.value 
      perfvec[2,j]=T2energy
        
      #RMMD
      rsetdata=as.matrix(compurank(rbind(mat3,mat4))) 
      T3Rmmd=kmmd(rsetdata[1:100,],rsetdata[(100+1):(2*100),],ntimes=1)@mmdstats[2]
      perfvec[3,j]=length(which(tstat[,2]>T3Rmmd))/5000 #T3Rmmd
     
      #MMD
      T4mmd=kmmd(mat3,mat4,ntimes=150)@mmdstats[2]
      perfvec[4,j]=length(which(tstatmmd2>T4mmd))/5000
  
      print(j)
    }
    for(i in 1:4)
    {
      for(j in 1:1000)
      {
        if(perfvec[i,j]<0.01)
        {perfvec[i,1001]=1+perfvec[i,1001]}
        if(perfvec[i,j]<0.05)
        {perfvec[i,1002]=1+perfvec[i,1002]}
        if(perfvec[i,j]<0.1)
        {perfvec[i,1003]=1+perfvec[i,1003]}
        
      }
      #power
      if(j==1000)
      {
        perfvec[i,1001]=perfvec[i,1001]/1000
        perfvec[i,1002]=perfvec[i,1002]/1000
        perfvec[i,1003]=perfvec[i,1003]/1000
      }
    }
    file_path <- file.path(desktop_path, paste0('Rtest','100_uniform_p', p, "_c", c,"_homo.csv"))
    write.csv(perfvec, file = file_path)
    }
}

####-----------------Real data analysis-----------------###
data("Sonar")
#let n=40
random_numbers1 <- sample(1:97, 97)
random_numbers2 <- sample(98:208, 111)
class1=Sonar[random_numbers1,1:60]
class2=Sonar[random_numbers2,1:60]
class1=as.matrix(class1)
class2=as.matrix(class2)
sonartest=matrix(0,4,1)
#Rank Energy
T1re=computestatistic(class1,class2)
sonartest[1]=length(which(energytstats>T1re))/30000
#Energy
m1=nrow(class1)
n1=nrow(class2)
comdata=dist(rbind(class1,class2))
T2senergy=eqdist.etest(comdata,sizes = c(m1,n1), R=199)$p.value
sonartest[2]=T2senergy
#RMMD
rsetdata=as.matrix(compurank(rbind(class1,class2))) 
T3Rmmd=kmmd(rsetdata[1:97,],rsetdata[(97):(2*97),],ntimes=1)@mmdstats[2]
sonartest[3]=length(which(tstats[,2]>T3Rmmd))/5000 #T3Rmmd
#MMD
T4mmd=kmmd(class1,class2,ntimes=150)@mmdstats[2]
sonartest[4]=length(which(tstatmmd2s>T4mmd))/5000
sonartest


data("Ionosphere")
good_reflection <- Ionosphere[Ionosphere$Class == "good", ]
bad_reflection <- Ionosphere[Ionosphere$Class == "bad", ]
random_numbers1 <- sample(1:225, 225)
random_numbers2 <- sample(1:126, 126)
rclass1=good_reflection[random_numbers1,1:34]
rclass2=bad_reflection[random_numbers2,1:34]
class1=as.matrix(rclass1)
class2=as.matrix(rclass2)
class1=apply(class1, c(1, 2), as.numeric)
class2=apply(class2, c(1, 2), as.numeric)
energytstats=gensamdist(97,111,60)
tstats=generdist(97,111,60)
tstatmmd2s=generdist1(97,111,60)

energytstationo=gensamdist(225,126,34)
tstationo=generdist(225,126,34)
tstatmmd2iono=generdist1(225,126,34)
ionospheretest=matrix(0,4,1)

  
 