set.seed(030699)
library(nlme)
library(mgcv)

source("D:/UFSCar - USP/Dissertação - Andressa Cerqueira/Códigos do R/2023.1/2. Teste de Hipóteses Poisson.R", local=TRUE)

#Caso balanceado para K=3, q=(1/3,1/3,1/3), n=500 e variando |a-b|=(0,1)

#####
#Variando a diferença entre a e b, com b=3 (0.25; 0.45) e com n = 500.
#####

q <- c(1/3,1/3,1/3)
a <- 3 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
S <- 100 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
d <- seq(0.25,0.45,0.01)
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaonaboot <- c()
proporcaosd <- c()
proporcaobootsd <- c()
media <- c()
mediaboot <- c()
contador <- 1
for(j in 1:length(d)){
  a <- b+d[j]
  P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
  {
    for(i in 1:S){
      z0 <- sample.z(n,q)
      X <- sample.g(z0,P)
      result <- SeqTestSBM(X, K.min = 1, K.max = 10)
      z_hat[i] <- result$K.est
      z_hatboot[i] <- result$K.est.boot
      cat(j, ':', paste0(round(100*i/S), '%'), '\r')
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
  }}

#####
#Dados 1
#####
proporcaoboot <- read.table("3.2 proporcaoboot2.txt")
proporcao <- read.table("3.2 proporcao2.txt")
d <- seq(0.25,0.45,0.02)

x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

plot(x,y1,type="l",xlab="Diferença das probabilidades entre comunidades", ylab="Proporção", ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(0.39,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
d <- seq(0.25,0.45,0.01) 
proporcao <- read.table("3.2 proporcao3.txt")
proporcaoboot <- read.table("3.2 proporcaoboot3.txt")
proporcaona <- read.table("3.2 proporcao3na.txt")
proporcaobootna <- read.table("3.2 proporcaoboot3na.txt")
media <- read.table("3.2 media3.txt")
mediaboot <- read.table("3.2 mediaboot3.txt")
sd <- read.table("3.2 proporcao3sd.txt")
sdboot <- read.table("3.2 proporcaoboot3sd.txt")

tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Média" = round(media$x,3),
                      "SD" = round(sd$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      "SD boot" = round(sdboot$x,3),
                      "Média boot" = round(mediaboot$x,3),
                      row.names = d)

#Gráfico
x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
media <- (media$x)
mediaboot <- (mediaboot$x) 
sd <- (sd$x)
sdboot <- (sdboot$x)

plot(x,y1,type="l",xlab=expression(a-b), ylab=expression("Hit rate of" ~ hat(K)), ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(0.38,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(0.38,1.5, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(0.38,1.5, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

#####
#Variando a diferença entre a e b, com b=5 (0.35; 0.55) e com n = 500.
#####

q <- c(1/3,1/3,1/3)
a <- 5 #prob dentro da comunidade
b <- 5 #prob fora das comunidades
S <- 100 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
d <- seq(0.35,0.57,0.01)
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaonaboot <- c()
proporcaosd <- c()
proporcaobootsd <- c()
media <- c()
mediaboot <- c()
contador <- 1
for(j in 1:length(d)){
  a <- b+d[j]
  P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))
  {
    for(i in 1:S){
      z0 <- sample.z(n,q)
      X <- sample.g(z0,P)
      result <- SeqTestSBM(X, K.min = 1, K.max = 10)
      z_hat[i] <- result$K.est
      z_hatboot[i] <- result$K.est.boot
      cat(j, ':', paste0(round(100*i/S), '%'), '\r')
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
  }}

#####
#Dados 1
#####
proporcaoboot <- read.table("3.2 proporcaoboot1.txt")
proporcao <- read.table("3.2 proporcao1.txt")
d <- seq(0.35,0.57,0.02)

x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

plot(x,y1,type="l",xlab="Diferença das probabilidades entre comunidades", ylab="Proporção", ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(0.51,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
d <- seq(0.35,0.57,0.01) 
proporcao <- read.table("3.2 proporcao4.txt")
proporcaoboot <- read.table("3.2 proporcaoboot4.txt")
proporcaona <- read.table("3.2 proporcao4na.txt")
proporcaobootna <- read.table("3.2 proporcaoboot4na.txt")
media <- read.table("3.2 media4.txt")
mediaboot <- read.table("3.2 mediaboot4.txt")
sd <- read.table("3.2 proporcao4sd.txt")
sdboot <- read.table("3.2 proporcaoboot4sd.txt")

tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Média" = round(media$x,3),
                      "SD" = round(sd$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      "SD boot" = round(sdboot$x,3),
                      "Média boot" = round(mediaboot$x,3),
                      row.names = d)

#Gráfico
x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
media <- (media$x)
mediaboot <- (mediaboot$x) 
sd <- (sd$x)
sdboot <- (sdboot$x)

plot(x,y1,type="l",xlab=expression(a-b), ylab=expression("Taxa de acertos para" ~ hat(K)), ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(0.5,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(0.5,1.5, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(0.5,1.5, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)


