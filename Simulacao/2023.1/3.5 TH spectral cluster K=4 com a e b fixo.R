set.seed(030699)
library(nlme)
library(mgcv)

source("D:/UFSCar - USP/Dissertação - Andressa Cerqueira/Códigos do R/2023.1/2. Teste de Hipóteses Poisson.R", local=TRUE)

#Variando n com |a-b| = 1, para caso balanceado e desbalanceado

#####
#1) Caso balanceado K=4
#Teste para q=(1/4,1/4,1/4,1/4), a e b próximos: a=4 e b=3 
#####

q <- c(1/4,1/4,1/4,1/4)
a <- 4 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
P <- rbind(c(a,b,b,b),c(b,a,b,b),c(b,b,a,b),c(b,b,b,a))
S <- 100 #numero de vezes para realizar a simulação

N <- seq(70, 260, 10)
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
  proporcao[contador] <- mean(z_hat==4)
  proporcaoboot[contador] <- mean(z_hatboot==4)
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

write.table(proporcao,"3.5 proporcao1.txt")
write.table(proporcaoboot,"3.5 proporcaoboot1.txt")
write.table(proporcaona,"3.5 proporcao1na.txt")
write.table(proporcaobootna,"3.5 proporcaoboot1na.txt")
write.table(proporcaosd,"3.5 proporcao1sd.txt")
write.table(proporcaobootsd,"3.5 proporcaoboot1sd.txt")
write.table(media, "3.5 media1.txt")
write.table(mediaboot, "3.5 mediaboot1.txt")

#Tabela
N <- seq(70, 260 , 10)
valorx <- c(70, 260 , 10)
proporcao <- read.table("3.5 proporcao1.txt")
proporcaoboot <- read.table("3.5 proporcaoboot1.txt")
proporcaona <- read.table("3.5 proporcao1na.txt")
proporcaobootna <- read.table("3.5 proporcaoboot1na.txt")
media <- read.table("3.5 media1.txt")
mediaboot <- read.table("3.5 mediaboot1.txt")
sd <- read.table("3.5 proporcao1sd.txt")
sdboot <- read.table("3.5 proporcaoboot1sd.txt")

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

plot(x,y1,type="l",xlab="Número de vértices no grafo", ylab="Proporção", xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(215,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab="Número de vértices no grafo", ylab="N° de comunidades",
     abline(h=4, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(215,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.8, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab="Número de vértices no grafo", ylab="N° de comunidades",
     abline(h=4, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(215,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.8, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

#####
#Gráfico
#####
x <- N
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)

#Gráfico de proporções para balanceado K=4
plot(x,y1,type="l",xlab="Número de vértices no grafo", ylab="Proporção", xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(210,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

#####
#2) Caso balanceado K=4
#Teste para q=(1/4,1/4,1/4,1/4), a e b próximos: a=5 e b=4 
#####

q <- c(1/4,1/4,1/4,1/4)
a <- 6 #prob dentro da comunidade
b <- 5 #prob fora das comunidades
P <- rbind(c(a,b,b,b),c(b,a,b,b),c(b,b,a,b),c(b,b,b,a))
S <- 100 #numero de vezes para realizar a simulação

N <- seq(100, 310 , 5)
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
  proporcao[contador] <- mean(z_hat==4)
  proporcaoboot[contador] <- mean(z_hatboot==4)
  media[contador] <- mean(z_hat)
  mediaboot[contador] <- mean(z_hatboot)
  proporcaosd[contador] <- sd(z_hat)
  proporcaobootsd[contador] <- sd(z_hatboot)
  contador <- contador + 1
}
write.table(proporcao,"3.5 proporcao2.txt")
write.table(proporcaoboot,"3.5 proporcaoboot2.txt")
write.table(proporcaona,"3.5 proporcao2na.txt")
write.table(proporcaobootna,"3.5 proporcaoboot2na.txt")
write.table(proporcaosd,"3.5 proporcao2sd.txt")
write.table(proporcaobootsd,"3.5 proporcaoboot2sd.txt")
write.table(media, "3.5 media2.txt")
write.table(mediaboot, "3.5 mediaboot2.txt")

#Tabela
N <- seq(100, 310 , 5)
valorx <- c(100, 310 , 5)
proporcao <- read.table("3.5 proporcao2.txt")
proporcaoboot <- read.table("3.5 proporcaoboot2.txt")
proporcaona <- read.table("3.5 proporcao2na.txt")
proporcaobootna <- read.table("3.5 proporcaoboot2na.txt")
media <- read.table("3.5 media2.txt")
mediaboot <- read.table("3.5 mediaboot2.txt")
sd <- read.table("3.5 proporcao2sd.txt")
sdboot <- read.table("3.5 proporcaoboot2sd.txt")
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

plot(x,y1,type="l",xlab="Número de vértices no grafo", ylab="Proporção", xaxp=valorx, ylim=c(0,1), col="#00008B")
lines(x,y2,col="#FF007F")
legend(130,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap")),
       col=c("#00008B","#FF007F"), cex=0.8, lty = c(1,1),
       box.lty=0)

plot(x, media, col="#00008B",
     ylim=range(c(media-sd, media+sd)),
     pch=16, xlab="Número de vértices no grafo", ylab="N° de comunidades",
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, media-sd, x, media+sd, length=0.05, angle=90, code=3)
legend(135,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.8, lty = c(0,1), pch=c(16,NA),
       box.lty=0)

plot(x, mediaboot, col="#00008B",
     ylim=range(c(mediaboot-sdboot, mediaboot+sdboot)),
     pch=16, xlab="Número de vértices no grafo", ylab="N° de comunidades",
     abline(h=3, col="#FF007F", lty=2)
)
arrows(x, mediaboot-sdboot, x, mediaboot+sdboot, length=0.05, angle=90, code=3)
legend(135,2, legend=c(expression("Média"),expression("Desvio Padrão")),
       col=c("#00008B","black"), cex=0.8, lty = c(0,1), pch=c(16,NA),
       box.lty=0)