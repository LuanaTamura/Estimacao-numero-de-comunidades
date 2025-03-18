library(RMTstat)
library(mgcv)

library(randnet)

source("D:/UFSCar - USP/Dissertação - Andressa Cerqueira/Códigos do R/2022.2/2. Bernoulli_TH.R", local=TRUE)

q <- c(1/3,1/3,1/3)

a <- 0.3#prob dentro da comunidade
b <- 0.2 #prob fora das comunidades


#testar para a=0.7 b=0.2 e a diferen?a ser? maior 

P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))

#vetor de comunidades verdadeiro
z0 <- sample.z(n, q)

#gerar matriz de adjac?ncia
X <- sample.g(z0,P)
#res <- SeqTestSBM(X)
S <- 100 #numero de vezes para realizar a simulação


# Teste novo 

#Para a=8 e b=3
##### 
#n=100
#####
mc2_n100 <- for(i in 1:S){
  n <- 100
  z0 <- sample.z(n,q)
  X <- sample.g(z0,P)
  test.hip(n)
  print(i)
  if(is.na(test.hip[i])) break
}

#Antigo
  
  q <- c(1/3,1/3,1/3)
  
  a <- 0.7 #prob dentro da comunidade
  b <- 0.2 #prob fora das comunidades
  
  
  #testar para a=0.7 b=0.2 e a diferença será maior 
  
  P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
  
  #vetor de comunidades verdadeiro
  z0 <- sample.z(n, q)
  
  #gerar matriz de adjacência
  
  X <- sample.g(z0,P)
  
  res <- SeqTestSBM(X)
  return(res)
}

#Para a=0.3 e b=0.2
##### 
#n=100
#####
mc_n100 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 100)}
))
write.table(mc_n100, "mc_n100.txt")
mc_n100 <- read.table("mc_n100.txt")

table(
  unlist(mc_n100[,1])
)

table(
  unlist(mc_n100[,2])
)

table(
  unlist(mc_n100[,1]),
  unlist(mc_n100[,2])
)

#####
#n=200
#####
mc_n200 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 200)}
))
write.table(mc_n200, "mc_n200.txt")
mc_n200 <- read.table("mc_n200.txt")
table(
  unlist(mc_n200[,1])
)
table(
  unlist(mc_n200[,2])
)
table(
  unlist(mc_n200[,1]),
  unlist(mc_n200[,2])
)

#####
#n=300
#####
mc_n300 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 300)}
))
write.table(mc_n300, "mc_n300.txt")
mc_n300 <- read.table("mc_n300.txt")
table(
  unlist(mc_n300[,1])
)
table(
  unlist(mc_n300[,2])
)
table(
  unlist(mc_n300[,1]),
  unlist(mc_n300[,2])
)
#####
#n=400
#####
mc_n400 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 400)}
))
write.table(mc_n400, "mc_n400.txt")
mc_n400 <- read.table("mc_n400.txt")
table(
  unlist(mc_n400[,1])
)
table(
  unlist(mc_n400[,2])
)
table(
  unlist(mc_n400[,1]),
  unlist(mc_n400[,2])
)

#####
#n=500
#####
mc_n500 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 500)}
))
write.table(mc_n500, "mc_n500.txt")
mc_n500 <- read.table("mc_n500.txt")
table(
  unlist(mc_n500[,1])
)
table(
  unlist(mc_n500[,2])
)
table(
  unlist(mc_n500[,1]),
  unlist(mc_n500[,2])
)





#Para a=0.7 e b=0.2
##### 
#n=40
#####
mc2_n40 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 40)}
))
write.table(mc2_n40, "mc2_n40.txt")
mc2_n40 <- read.table("mc2_n40.txt")

table(
  unlist(mc2_n40[,1])
)

table(
  unlist(mc2_n40[,2])
)

table(
  unlist(mc2_n40[,1]),
  unlist(mc2_n40[,2])
)

##### 
#n=50
#####
mc2_n50 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 50)}
))
write.table(mc2_n50, "mc2_n50.txt")
mc2_n50 <- read.table("mc2_n50.txt")

table(
  unlist(mc2_n50[,1])
)

table(
  unlist(mc2_n50[,2])
)

table(
  unlist(mc2_n50[,1]),
  unlist(mc2_n50[,2])
)

##### 
#n=60
#####
mc2_n60 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 60)}
))
write.table(mc2_n60, "mc2_n60.txt")
mc2_n60 <- read.table("mc2_n60.txt")

table(
  unlist(mc2_n60[,1])
)

table(
  unlist(mc2_n60[,2])
)

table(
  unlist(mc2_n60[,1]),
  unlist(mc2_n60[,2])
)
##### 
#n=70
#####
mc2_n70 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 70)}
))
write.table(mc2_n70, "mc2_n70.txt")
mc2_n70 <- read.table("mc2_n70.txt")

table(
  unlist(mc2_n70[,1])
)

table(
  unlist(mc2_n70[,2])
)

table(
  unlist(mc2_n70[,1]),
  unlist(mc2_n70[,2])
)

##### 
#n=80
#####
mc2_n80 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 80)}
))
write.table(mc2_n80, "mc2_n80.txt")
mc2_n80 <- read.table("mc2_n80.txt")

table(
  unlist(mc2_n80[,1])
)

table(
  unlist(mc2_n80[,2])
)

table(
  unlist(mc2_n80[,1]),
  unlist(mc2_n80[,2])
)

##### 
#n=100
#####
mc2_n100 <- t(replicate(
  n = 100,
  expr = {test.hip(n = 100)}
))
write.table(mc2_n100, "mc2_n100.txt")
mc2_n100 <- read.table("mc2_n100.txt")

table(
  unlist(mc2_n100[,1])
)

table(
  unlist(mc2_n100[,2])
)

table(
  unlist(mc2_n100[,1]),
  unlist(mc2_n100[,2])
)






