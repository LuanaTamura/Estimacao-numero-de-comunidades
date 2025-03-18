set.seed(030699)
library(nlme)
library(mgcv)

source("D:/UFSCar - USP/Dissertação - Andressa Cerqueira/Códigos do R/2023.1/2. Teste de Hipóteses Poisson.R", local=TRUE)


#Caso desbalanceado para K=3, q=(0.5,0.3,0.2), n=500 e variando |a-b|=(0,1)

#####
#Variando a diferença entre a e b, com b=3 (0.2; 0.55) e com n = 500.
#####

q <- c(0.5,0.3,0.2)
a <- 3 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
S <- 10 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaonaboot <- c()
contador <- 1
d <- 0.5
d <- seq(0.2,0.54,0.02) 
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
    proporcaonaboot[contador] <- sum(is.na(z_hatboot))/S
    proporcao[contador] <- mean(z_hat==3, na.rm = TRUE)
    proporcaoboot[contador] <- mean(z_hatboot==3, na.rm = TRUE)
    contador <- contador + 1
  }}

proporcao
proporcaoboot
proporcaona
proporcaonaboot


#####
#Dados 1
#####
d <- seq(0.2,0.54,0.02) 
proporcao <- read.table("3.4 proporcao1.txt")
proporcaoboot <- read.table("3.4 proporcaoboot1.txt")
proporcaona <- read.table("3.4 proporcao1na.txt")
proporcaobootna <- read.table("3.4 proporcaoboot1na.txt")
tabela7 <- data.frame ("Prop" = round(proporcao$x,3), 
                       "Prop NA" = round(proporcaona$x,3),
                       "Prop boot" = round(proporcaoboot$x,3),
                       "Prop boot NA" = round(proporcaobootna$x,3),
                       row.names = d)
print(tabela7)

#Gráfico
d <- seq(0.2,0.7,0.02) 
x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

plot(x,y1,type="l",xlab="Diferença das probabilidades entre comunidades", ylab="Proporção", ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(0.55,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
d <- seq(0.2,0.68,0.02)
proporcao <- read.table("3.4 proporcao3.txt")
proporcaoboot <- read.table("3.4 proporcaoboot3.txt")
proporcaona <- read.table("3.4 proporcao3na.txt")
proporcaobootna <- read.table("3.4 proporcaoboot3na.txt")
media <- read.table("3.4 media3.txt")
mediaboot <- read.table("3.4 mediaboot3.txt")
sd <- read.table("3.4 proporcao3sd.txt")
sdboot <- read.table("3.4 proporcaoboot3sd.txt")

tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop NA" = round(proporcaona$x,3),
                      "Média" = round(media$x,3),
                      "SD" = round(sd$x,3),
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop boot NA" = round(proporcaobootna$x,3),
                      "SD boot" = round(sdboot$x,3),
                      "Média boot" = round(mediaboot$x,3),
                      row.names = d)
print(tabela)

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
legend(0.52,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(0.52,1.5, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(0.52,1.5, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

#####
#Caso desbalanceado para K=3, q=(0.5,0.3,0.2), n=500 e variando |a-b|=(0,1)

#####
#Variando a diferença entre a e b, com b=5 (0.25; 0.45) e com n = 500.
#####

q <- c(0.5,0.3,0.2)
a <- 5 #prob dentro da comunidade
b <- 5 #prob fora das comunidades
S <- 10 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
z_hat <- vector()
z_hatboot <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcaona <- c()
proporcaonaboot <- c()
contador <- 1
d <- seq(0.32,0.84,0.03) 
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
    proporcaonaboot[contador] <- sum(is.na(z_hatboot))/S
    proporcao[contador] <- mean(z_hat==3, na.rm = TRUE)
    proporcaoboot[contador] <- mean(z_hatboot==3, na.rm = TRUE)
    contador <- contador + 1
  }}

proporcao
proporcaoboot
proporcaona
proporcaonaboot

#####
#Dados 1
#####

d <- seq(0.32,0.84,0.03) 
proporcao <- read.table("3.4 proporcao2.txt")
proporcaoboot <- read.table("3.4 proporcaoboot2.txt")
proporcaona <- read.table("3.4 proporcao2na.txt")
proporcaobootna <- read.table("3.4 proporcaoboot2na.txt")
tabela8 <- data.frame ("Prop" = round(proporcao$x,3), 
                       "Prop NA" = round(proporcaona$x,3),
                       "Prop boot" = round(proporcaoboot$x,3),
                       "Prop boot NA" = round(proporcaobootna$x,3),
                       row.names = d)
print(tabela8)

#Gráfico
d <- seq(0.32,0.84,0.03) 
x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

plot(x,y1,type="l",xlab="Diferença das probabilidades entre comunidades", ylab="Proporção", ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(0.7,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#Dados 2
#####
d <- seq(0.32,0.87,0.03) 
proporcao <- read.table("3.4 proporcao4.txt")
proporcaoboot <- read.table("3.4 proporcaoboot4.txt")
proporcaona <- read.table("3.4 proporcao4na.txt")
proporcaobootna <- read.table("3.4 proporcaoboot4na.txt")
media <- read.table("3.4 media4.txt")
mediaboot <- read.table("3.4 mediaboot4.txt")
sd <- read.table("3.4 proporcao4sd.txt")
sdboot <- read.table("3.4 proporcaoboot4sd.txt")

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
legend(0.69,0.4, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.7, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(0.69,1.7, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab=expression(a-b), ylab=expression(hat(K)),
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(0.69,1.7, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.7, lty = c(0,1), pch=c(16,NA),
       box.lty=0)
