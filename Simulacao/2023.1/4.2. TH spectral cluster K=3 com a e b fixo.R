set.seed(030699)
library(nlme)
library(mgcv)
install.packages("randnet")
library(randnet)
install.packages("entropy")
install.packages("AUC")
install.packages("poweRlaw")
install.packages("RSpectra")
install.packages("irlba")

source("D:/UFSCar - USP/Dissertação - Andressa Cerqueira/Códigos do R/2023.1/2. Teste de Hipóteses Poisson.R", local=TRUE)

#Variando n com |a-b| = 1, para caso balanceado e desbalanceado

#####
#1) Caso balanceado K=3
#Teste para q=(1/3,1/3,1/3), a e b próximos: a=4 e b=3 
#####

  q <- c(1/3,1/3,1/3)
  a <- 4 #prob dentro da comunidade
  b <- 3 #prob fora das comunidades
  P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
  S <- 100 #numero de vezes para realizar a simulação

N <- seq(40, 165, 5)
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaobootna <- c()
proporcaosd <- c()
proporcaobootsd <- c()
media <- c()
mediaboot <- c()
contador <- 1
for(n in N){
for(i in 1:S){
  z0 <- sample.z(n,q)
  X <- sample.g(z0,P)
  result <- SeqTestSBM(X, K.min = 1, K.max = 5)
  z_hat[i] <- result$K.est
  z_hatboot[i] <- result$K.est.boot
  cat(n, ':', paste0(round(100*i/S), '%'), '\r')
  #if(is.na(test.hip[i])) break
}
  proporcaona[contador] <- sum(is.na(z_hat))/S
  proporcaobootna[contador] <- sum(is.na(z_hatboot))/S
  z_hat[is.na(z_hat)] <- 0
  z_hatboot[is.na(z_hatboot)] <- 0
  proporcao[contador] <- mean(z_hat==3)
  proporcaoboot[contador] <- mean(z_hatboot==3)
  media[contador] <- mean(z_hat)
  mediaboot[contador] <- mean(z_hatboot)
  proporcaosd[contador] <- sd(z_hat)
  proporcaobootsd[contador] <- sd(z_hatboot)
  contador <- contador + 1
}

proporcao 
proporcaoboot
proporcaona 
proporcaobootna 
proporcaosd 
proporcaobootsd 
media 
mediaboot 

#####
#Dados 1
#####
N <- seq(50,165,5)
valorx <- c(50,165,5)
proporcaoboot <- read.table("3.1 proporcaoboot1.txt")
proporcao <- read.table("3.1 proporcao1.txt")

x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)


plot(x,y1,type="l",xlab="Número de vértices no grafo", ylab="Proporção", xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(130,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
N <- seq(40,165,5)
valorx <- c(40,165,5)
proporcao <- read.table("3.1 proporcao4.txt")
proporcaoboot <- read.table("3.1 proporcaoboot4.txt")
proporcaona <- read.table("3.1 proporcao4na.txt")
proporcaobootna <- read.table("3.1 proporcaoboot4na.txt")
media <- read.table("3.1 media4.txt")
mediaboot <- read.table("3.1 mediaboot4.txt")
sd <- read.table("3.1 proporcao4sd.txt")
sdboot <- read.table("3.1 proporcaoboot4sd.txt")

tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Média" = round(media$x,3),
                      "SD" = round(sd$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      "SD boot" = round(sdboot$x,3),
                      "Média boot" = round(mediaboot$x,3),
                      row.names = N)
print(tabela)

x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
media <- (media$x)
mediaboot <- (mediaboot$x) 
sd <- (sd$x)
sdboot <- (sdboot$x)
y <- expression(\hat(K))
 
plot(x,y1,type="l",xlab="n", ylab=expression("Taxa de acertos para" ~ hat(K)), xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(125,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab="n", ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(125,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=c(0.5,4.5),
     pch=16, xlab="n", ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(125,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

#####
#2) Caso desbalanceado K=3
#Testar para q=(0.4,0.3,0.3) a e b próximos: a=4 e b=3 
#####

q <- c(0.4,0.3,0.3)
a <- 4 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
S <- 100 #numero de vezes para realizar a simulação

N <- seq(40,194,7)
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaobootna <- c()
proporcaosd <- c()
proporcaobootsd <- c()
media <- c()
mediaboot <- c()
contador <- 1
for(n in N){
  for(i in 1:S){
    z0 <- sample.z(n,q)
    X <- sample.g(z0,P)
    result <- SeqTestSBM(X, K.min = 1, K.max = 10)
    z_hat[i] <- result$K.est
    z_hatboot[i] <- result$K.est.boot
    cat(n, ':', paste0(round(100*i/S), '%'), '\r')
    #if(is.na(test.hip[i])) break
  }
  proporcaona[contador] <- sum(is.na(z_hat))/S
  proporcaobootna[contador] <- sum(is.na(z_hatboot))/S
  z_hat[is.na(z_hat)] <- 0
  z_hatboot[is.na(z_hatboot)] <- 0
  proporcao[contador] <- mean(z_hat==3)
  proporcaoboot[contador] <- mean(z_hatboot==3)
  media[contador] <- mean(z_hat)
  mediaboot[contador] <- mean(z_hatboot)
  proporcaosd[contador] <- sd(z_hat)
  proporcaobootsd[contador] <- sd(z_hatboot)
  contador <- contador + 1
}

#####
#Dados 1
#####
N <- seq(40,194,7)
valorx <- c(40,194,7)
proporcao <- read.table("3.1 proporcao2.txt")
proporcaona <- read.table("3.1 proporcao2na.txt")
proporcaoboot <- read.table("3.1 proporcaoboot2.txt")
proporcaobootna <- read.table("3.1 proporcaoboot2na.txt")
tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      row.names = N)
x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

plot(x,y1,type="l",xlab="Número de vértices no grafo", ylab="Proporção", xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(150,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
N <- seq(40,194,7)
valorx <- c(40,195,5)
proporcao <- read.table("3.1 proporcao5.txt")
proporcaoboot <- read.table("3.1 proporcaoboot5.txt")
proporcaona <- read.table("3.1 proporcao5na.txt")
proporcaobootna <- read.table("3.1 proporcaoboot5na.txt")
media <- read.table("3.1 media5.txt")
mediaboot <- read.table("3.1 mediaboot5.txt")
sd <- read.table("3.1 proporcao5sd.txt")
sdboot <- read.table("3.1 proporcaoboot5sd.txt")

tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Média" = round(media$x,3),
                      "SD" = round(sd$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      "SD boot" = round(sdboot$x,3),
                      "Média boot" = round(mediaboot$x,3),
                      row.names = N)
print(tabela)

x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
media <- (media$x)
mediaboot <- (mediaboot$x) 
sd <- (sd$x)
sdboot <- (sdboot$x)

plot(x,y1,type="l",xlab="n", ylab=expression("Taxa de acertos para" ~ hat(K)), xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(145,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B", xaxp=valorx,
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab="n", ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(145,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B", xaxp=valorx,
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab="n", ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(145,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

#####
#3) Caso desbalanceado K=3
#Testar para q=(0.5,0.3,0.2) a e b próximos: a=4 e b=3 
#####

q <- c(0.5,0.3,0.3)
a <- 4 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
S <- 100 #numero de vezes para realizar a simulação

N <- seq(40,280,10)
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaobootna <- c()
proporcaosd <- c()
proporcaobootsd <- c()
media <- c()
mediaboot <- c()
contador <- 1
for(n in N){
  for(i in 1:S){
    z0 <- sample.z(n,q)
    X <- sample.g(z0,P)
    result <- SeqTestSBM(X, K.min = 1, K.max = 10)
    z_hat[i] <- result$K.est
    z_hatboot[i] <- result$K.est.boot
    cat(n, ':', paste0(round(100*i/S), '%'), '\r')
    #if(is.na(test.hip[i])) break
  }
  proporcaona[contador] <- sum(is.na(z_hat))/S
  proporcaobootna[contador] <- sum(is.na(z_hatboot))/S
  z_hat[is.na(z_hat)] <- 0
  z_hatboot[is.na(z_hatboot)] <- 0
  proporcao[contador] <- mean(z_hat==3)
  proporcaoboot[contador] <- mean(z_hatboot==3)
  media[contador] <- mean(z_hat)
  mediaboot[contador] <- mean(z_hatboot)
  proporcaosd[contador] <- sd(z_hat)
  proporcaobootsd[contador] <- sd(z_hatboot)
  contador <- contador + 1
}

#####
#Dados 1  
#####
N <- seq(40,280,10)
valorx <- c(40,280,10)
proporcao <- read.table("3.1 proporcao3.txt")
proporcaona <- read.table("3.1 proporcao3na.txt")
proporcaoboot <- read.table("3.1 proporcaoboot3.txt")
proporcaobootna <- read.table("3.1 proporcaoboot3na.txt")
tabela <- data.frame ("Prop" = round(proporcao$x,3), 
            "Prop NA" = round(proporcaona$x,3),
            "Prop boot" = round(proporcaoboot$x,3),
            "Prop boot NA" = round(proporcaobootna$x,3),
            row.names = N)

x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

plot(x,y1,type="l",xlab="Número de vértices no grafo", ylab="Proporção", xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(220,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
N <- seq(40,230,10)
valorx <- c(40,230,5)
proporcao <- read.table("3.1 proporcao6.txt")
proporcaoboot <- read.table("3.1 proporcaoboot6.txt")
proporcaona <- read.table("3.1 proporcao6na.txt")
proporcaobootna <- read.table("3.1 proporcaoboot6na.txt")
media <- read.table("3.1 media6.txt")
mediaboot <- read.table("3.1 mediaboot6.txt")
sd <- read.table("3.1 proporcao6sd.txt")
sdboot <- read.table("3.1 proporcaoboot6sd.txt")

tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Média" = round(media$x,3),
                      "SD" = round(sd$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      "SD boot" = round(sdboot$x,3),
                      "Média boot" = round(mediaboot$x,3),
                      row.names = N)
print(tabela)

x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
media <- (media$x)
mediaboot <- (mediaboot$x) 
sd <- (sd$x)
sdboot <- (sdboot$x)

plot(x,y1,type="l",xlab="n", ylab=expression("Taxa de acertos para" ~ hat(K)), xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(170,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B", xaxp=valorx,
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab="n", ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(170,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B", xaxp=valorx,
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab="n", ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(170,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

