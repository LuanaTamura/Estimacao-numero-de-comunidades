install.packages("RMTstat")
install.packages("mgcv")
library(RMTstat)
library(mgcv)

p = 1
beta1 <- 1 
beta2 <- 2 
beta4 <- 4
x <-  seq(-6,4, length = 200)
y1 <- dtw(x, beta=beta1)
y2 <- dtw(x, beta=beta2)
y3 <- dtw(x, beta=beta4)

plot(x,y1,type="l",xlab="x", ylab="densidade", ylim=c(0,max(y1,y2,y3)))
lines(x,y2,col="red")
lines(x,y3,col="blue")
legend(2,0.35, legend=c(expression(beta~" = 1"),expression(beta~" = 2"),expression(beta~" = 4")),
       col=c("black","red","blue"), cex=0.8, lty = c(1,1,1),
       box.lty=0)

#Teste de Hipoteses

sample.z <- function(n,q){
  k <- length(q)
  return(sample(x = seq(1:k), n, replace = T, prob = q))
}

#Matriz de adjacencia simetrica do grafo 
sample.g <- function(P,n){
  mat.adj <- matrix(2,n,n)
  for(i in 1:n){
    for(j in i:n){
      mat.adj[i,j] <- mat.adj[j,i] <- rbinom(1,1,P)
    }
  }
  # diag(mat.adj) <- 0 #A diagonal ? 0 
  return(mat.adj)
}

P <- 0.5 # probabilidade de conexão
n1 <- 100
n2 <- 300 
n3 <- 1000
S <- 1000 # número de simulações que resulta em autovalores

######
#Para n=100
######

#Maior e menor autovalor, convergindo a Tracy Widom
L_max1 <- vector()
L_min1 <- vector()
  for(i in 1:S){
  A1 <- sample.g(P, n1)
  A.star1 <- (A1-P)/sqrt((n1-1)*P*(1-P))
  L_max1[i] <- max(eigen(A.star1)$values)
  L_min1[i] <- min(eigen(A.star1)$values)
  print(i)
  }


######
#Para n=300
######

#Maior e menor autovalor, convergindo a Tracy Widom
L_max2 <- vector()
L_min2 <- vector()
for(i in 1:S){
  A2 <- sample.g(P, n2)
  A.star2 <- (A2-P)/sqrt((n2-1)*P*(1-P))
  L_max2[i] <- max(eigen(A.star2)$values)
  L_min2[i] <- min(eigen(A.star2)$values)
}


######
#Para n=1000
######

#Maior e menor autovalor, convergindo a Tracy Widom
L_max3 <- vector()
L_min3 <- vector()
for(i in 1:S){
  A3 <- sample.g(P, n3)
  A.star3 <- (A3-P)/sqrt((n3-1)*P*(1-P))
  L_max3[i] <- max(eigen(A.star3)$values)
  L_min3[i] <- min(eigen(A.star3)$values)
}

write.table(L_min1,"L_min n=100.txt")
write.table(L_max1,"L_max n=100.txt")
write.table(L_min2,"L_min n=300.txt")
write.table(L_max2,"L_max n=300.txt")
write.table(L_min3,"L_min n=1000.txt")
write.table(L_max3,"L_max n=1000.txt")

L_min1 <- read.table("L_min n=100.txt")
L_max1 <- read.table("L_max n=100.txt")
L_min2 <- read.table("L_min n=300.txt")
L_max2 <- read.table("L_max n=300.txt")
L_min3 <- read.table("L_min n=1000.txt")
L_max3 <- read.table("L_max n=1000.txt")


L_min1 <- L_min1$x
L_max1 <- L_max1$x
L_min2 <- L_min2$x
L_max2 <- L_max2$x
L_min3 <- L_min3$x
L_max3 <- L_max3$x


######
#Histogramas dos menores autovalores
######

par(mfrow=c(1,3))

hist((-L_min1-2)*n1^(2/3), freq = FALSE, xlab=expression(lambda), ylab="Densidade",
     main=bquote(lambda[1] ~ "com n = 100"), xlim = c(-6, 4), ylim = c(0, 0.4), col="grey")
x <- seq(min((-L_min1-2)*n1^(2/3)),max((-L_min1-2)*n1^(2/3)), 0.005)
y <- dtw(x, beta=1)
lines(x,y,col="blue")

hist((-L_min2-2)*n2^(2/3), freq = FALSE, xlab=expression(lambda), ylab="Densidade",
     main=bquote(lambda[1] ~ "com n = 300"), xlim = c(-6, 4), ylim = c(0, 0.4), col='grey')
x <- seq(min((-L_min2-2)*n2^(2/3)),max((-L_min2-2)*n2^(2/3)), 0.005)
y <- dtw(x, beta=1)
lines(x,y,col="blue")

hist((-L_min3-2)*n3^(2/3), freq = FALSE, xlab=expression(lambda), ylab="Densidade", 
     main=bquote(lambda[1] ~ "com n = 1000"), xlim = c(-6, 4), ylim = c(0, 0.4), col='grey')
x <- seq(min((-L_min3-2)*n3^(2/3)),max((-L_min3-2)*n3^(2/3)), 0.005)
y <- dtw(x, beta=1)
lines(x,y,col="blue")

legend(-1,0.4, legend=c(expression(Simulado), expression(TW[1])),
       col=c("grey", "blue"), cex=0.8, lty = c(1,1),
       box.lty=0)

#Histogramas dos maiores autovalores
par(mfrow=c(1,3))
hist((L_max1-2)*n1^(2/3), freq = FALSE, xlab=expression(lambda), ylab="Densidade", 
     main=bquote(lambda[n] ~ "com n = 100"), xlim = c(-6, 4), ylim = c(0, 0.4), col='grey') 
x <- seq(min((L_max1-2)*n1^(2/3)),max((L_max1-2)*n1^(2/3)), 0.005)
y <- dtw(x, beta=1)
lines(x,y, col="blue")

hist((L_max2-2)*n2^(2/3), freq = FALSE, xlab=expression(lambda), ylab="Densidade", 
     main=bquote(lambda[n] ~ "com n = 300"), xlim = c(-6, 4), ylim = c(0, 0.4), col='grey') 
x <- seq(min((L_max2-2)*n2^(2/3)),max((L_max2-2)*n2^(2/3)), 0.005)
y <- dtw(x, beta=1)
lines(x,y, col="blue")

hist((L_max3-2)*n3^(2/3), freq = FALSE, xlab=expression(lambda), ylab="Densidade", 
     main=bquote(lambda[n] ~ "com n = 1000"), xlim = c(-6, 4), ylim = c(0, 0.4), col='grey')
     x <- seq(min((L_max3-2)*n3^(2/3)),max((L_max3-2)*n3^(2/3)), 0.005)
y <- dtw(x, beta=1)
lines(x,y, col="blue")

legend(-1,0.4, legend=c(expression(Simulado), expression(TW[1])),
       col=c("grey", "blue"), cex=0.8, lty = c(1,1),
       box.lty=0)







