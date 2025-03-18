set.seed(030699)
install.packages("RMTstat")
install.packages("randnet")
install.packages("entropy")
install.packages("AUC")
install.packages("Matrix")
install.packages("sparseFLMM")
install.packages("refund")
install.packages("mgcv")
install.packages("nlme")
install.packages("igraph")

library(nlme)
library(mgcv)
library(Matrix)
library(entropy)
library(AUC)
library(sparseFLMM)
library(refund)
library(randnet)
library(igraph)

source("D:/UFSCar - USP/Dissertação - Andressa Cerqueira/Códigos do R/2023.1/2. Teste de Hipóteses Poisson.R", local=TRUE)

#Caso balanceado para K=3, q=(1/3,1/3,1/3), n=500 e variando |a-b|=(0,1)

#####
#Variando a diferença entre a e b, com b=3 (0.25; 0.45) e com n = 500.
#####

q <- c(1/3,1/3,1/3)
a <- 3 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
S <- 50 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
d <- seq(0.25,0.45,0.01)
z_hat <- vector()
z_hatboot <- vector()
fastgreedy <- vector()
louvain <- vector()
walktrap <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcao_fastgreedy <- c()
proporcao_louvain <- c()
proporcao_walk <- c()
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
    G <- graph_from_adjacency_matrix(X, weighted = TRUE, mode="undirected")
    fg <- cluster_fast_greedy(G)
    fastgreedy[i] <- max(fg$membership)
    louv <- cluster_louvain(G, weights = NULL)
    louvain[i] <- max(louv$membership)
    walk <- cluster_walktrap(G, weights = NULL, steps = 4,
                             merges = TRUE, modularity = TRUE, membership = TRUE)
    walktrap[i] <- max(walk$membership)
    cat(j, ':', paste0(round(100*i/S), '%'), '\r')
    #if(is.na(test.hip[i])) break
  }
  z_hat[is.na(z_hat)] <- 0
  z_hatboot[is.na(z_hatboot)] <- 0
  fastgreedy[is.na(fastgreedy)] <- 0
  louvain[is.na(louvain)] <- 0
  walktrap[is.na(walktrap)] <- 0
  proporcao[contador] <- mean(z_hat==3)
  proporcaoboot[contador] <- mean(z_hatboot==3)
  proporcao_fastgreedy[contador] <- mean(fastgreedy==3)
  proporcao_louvain[contador] <- mean(louvain==3)
  proporcao_walk[contador] <- mean(walktrap==3)
  contador <- contador + 1
}}

write.table(proporcao,"4.2 proporcao1.txt")
write.table(proporcaoboot,"4.2 proporcaoboot1.txt")
write.table(proporcao_fastgreedy,"4.2 proporcao_fastgreedy1.txt")
write.table(proporcao_louvain,"4.2 proporcao_louvain1.txt")
write.table(proporcao_walk,"4.2 proporcao_walk1.txt")

proporcao <- read.table("4.2 proporcao1.txt")
proporcaoboot <- read.table("4.2 proporcaoboot1.txt")
proporcao_fastgreedy <- read.table("4.2 proporcao_fastgreedy1.txt")
proporcao_louvain <- read.table("4.2 proporcao_louvain1.txt")
proporcao_walk <- read.table("4.2 proporcao_walk1.txt")

d <- seq(0.25,0.45,0.01)
tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop fast greedy" = round(proporcao_fastgreedy$x,3),
                      "Prop louvain" = round(proporcao_louvain$x,3),
                      "Prop walktrap" = round(proporcao_walk$x,3),
                      row.names = d)
print(tabela)

x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
y3 <- (proporcao_fastgreedy$x)
y4 <- (proporcao_louvain$x)
y5 <- (proporcao_walk$x)
y <- expression(\hat(K))

plot(x,y1,type="b",xlab="a-b", ylab=expression("Taxa de acertos para" ~ hat(K)), ylim=c(0,1), col="#00008B", pch=15)
lines(x,y2,col="#FF007F", type = "b", pch=16)
lines(x,y3,col="#20B2AA", type = "b", pch=17)
lines(x,y4,col="#D2691E", type = "b", pch=18)
lines(x,y5,col="#228B22", type = "b", pch=4)
legend(0.395,0.3, legend=c(expression("Sem bootstrap"),expression("Com bootstrap"), expression("Fast greedy"), expression("Louvain"), expression("Walktrap")),
       col=c("#00008B","#FF007F","#20B2AA","#D2691E","#228B22"), cex=0.7, lty = c(1,1), pch=c(15,16,17,18,4),
       box.lty=0)


#####
# Caso 2
#####

q <- c(0.4,0.3,0.3)
a <- 3 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
S <- 50 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
d <- seq(0.25,0.55,0.01) 
z_hat <- vector()
z_hatboot <- vector()
fastgreedy <- vector()
louvain <- vector()
walktrap <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcao_fastgreedy <- c()
proporcao_louvain <- c()
proporcao_walk <- c()
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
      G <- graph_from_adjacency_matrix(X, weighted = TRUE, mode="undirected")
      fg <- cluster_fast_greedy(G)
      fastgreedy[i] <- max(fg$membership)
      louv <- cluster_louvain(G, weights = NULL)
      louvain[i] <- max(louv$membership)
      walk <- cluster_walktrap(G, weights = NULL, steps = 4,
                               merges = TRUE, modularity = TRUE, membership = TRUE)
      walktrap[i] <- max(walk$membership)
      cat(j, ':', paste0(round(100*i/S), '%'), '\r')
      #if(is.na(test.hip[i])) break
    }
    z_hat[is.na(z_hat)] <- 0
    z_hatboot[is.na(z_hatboot)] <- 0
    fastgreedy[is.na(fastgreedy)] <- 0
    louvain[is.na(louvain)] <- 0
    walktrap[is.na(walktrap)] <- 0
    proporcao[contador] <- mean(z_hat==3)
    proporcaoboot[contador] <- mean(z_hatboot==3)
    proporcao_fastgreedy[contador] <- mean(fastgreedy==3)
    proporcao_louvain[contador] <- mean(louvain==3)
    proporcao_walk[contador] <- mean(walktrap==3)
    contador <- contador + 1
  }}

write.table(proporcao,"4.2 proporcao2.txt")
write.table(proporcaoboot,"4.2 proporcaoboot2.txt")
write.table(proporcao_fastgreedy,"4.2 proporcao_fastgreedy2.txt")
write.table(proporcao_louvain,"4.2 proporcao_louvain2.txt")
write.table(proporcao_walk,"4.2 proporcao_walk2.txt")

proporcao <- read.table("4.2 proporcao2.txt")
proporcaoboot <- read.table("4.2 proporcaoboot2.txt")
proporcao_fastgreedy <- read.table("4.2 proporcao_fastgreedy2.txt")
proporcao_louvain <- read.table("4.2 proporcao_louvain2.txt")
proporcao_walk <- read.table("4.2 proporcao_walk2.txt")

d <- seq(0.25,0.55,0.01)
tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop fast greedy" = round(proporcao_fastgreedy$x,3),
                      "Prop louvain" = round(proporcao_louvain$x,3),
                      "Prop walktrap" = round(proporcao_walk$x,3),
                      row.names = d)

x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
y3 <- (proporcao_fastgreedy$x)
y4 <- (proporcao_louvain$x)
y5 <- (proporcao_walk$x)
y <- expression(\hat(K))

plot(x,y1,type="b",xlab="a-b", ylab=expression("Taxa de acertos para" ~ hat(K)), ylim=c(0,1), col="#00008B", pch=15)
lines(x,y2,col="#FF007F", type = "b", pch=16)
lines(x,y3,col="#20B2AA", type = "b", pch=17)
lines(x,y4,col="#D2691E", type = "b", pch=18)
lines(x,y5,col="#228B22", type = "b", pch=4)
legend(0.45,0.5, legend=c(expression("Sem bootstrap"),expression("Com bootstrap"), expression("Fast greedy"), expression("Louvain"), expression("Walktrap")),
       col=c("#00008B","#FF007F","#20B2AA","#D2691E","#228B22"), cex=0.7, lty = c(1,1), pch=c(15,16,17,18,4),
       box.lty=0)

#####
# Caso 3
#####

q <- c(0.5,0.3,0.2)
a <- 3 #prob dentro da comunidade
b <- 3 #prob fora das comunidades
S <- 50 #numero de vezes para realizar a simulação
d <- vector()
n <- 500
d <- seq(0.2,0.68,0.02) 
z_hat <- vector()
z_hatboot <- vector()
fastgreedy <- vector()
louvain <- vector()
walktrap <- vector()
proporcao <- c()
proporcaoboot <- c()
proporcao_fastgreedy <- c()
proporcao_louvain <- c()
proporcao_walk <- c()
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
      G <- graph_from_adjacency_matrix(X, weighted = TRUE, mode="undirected")
      fg <- cluster_fast_greedy(G)
      fastgreedy[i] <- max(fg$membership)
      louv <- cluster_louvain(G, weights = NULL)
      louvain[i] <- max(louv$membership)
      walk <- cluster_walktrap(G, weights = NULL, steps = 4,
                               merges = TRUE, modularity = TRUE, membership = TRUE)
      walktrap[i] <- max(walk$membership)
      cat(j, ':', paste0(round(100*i/S), '%'), '\r')
      #if(is.na(test.hip[i])) break
    }
    z_hat[is.na(z_hat)] <- 0
    z_hatboot[is.na(z_hatboot)] <- 0
    fastgreedy[is.na(fastgreedy)] <- 0
    louvain[is.na(louvain)] <- 0
    walktrap[is.na(walktrap)] <- 0
    proporcao[contador] <- mean(z_hat==3)
    proporcaoboot[contador] <- mean(z_hatboot==3)
    proporcao_fastgreedy[contador] <- mean(fastgreedy==3)
    proporcao_louvain[contador] <- mean(louvain==3)
    proporcao_walk[contador] <- mean(walktrap==3)
    contador <- contador + 1
  }}

write.table(proporcao,"4.2 proporcao4.txt")
write.table(proporcaoboot,"4.2 proporcaoboot4.txt")
write.table(proporcao_fastgreedy,"4.2 proporcao_fastgreedy4.txt")
write.table(proporcao_louvain,"4.2 proporcao_louvain4.txt")
write.table(proporcao_walk,"4.2 proporcao_walk4.txt")

proporcao <- read.table("4.2 proporcao3.txt")
proporcaoboot <- read.table("4.2 proporcaoboot3.txt")
proporcao_fastgreedy <- read.table("4.2 proporcao_fastgreedy3.txt")
proporcao_louvain <- read.table("4.2 proporcao_louvain3.txt")
proporcao_walk <- read.table("4.2 proporcao_walk3.txt")

proporcao <- read.table("4.2 proporcao4.txt")
proporcaoboot <- read.table("4.2 proporcaoboot4.txt")
proporcao_fastgreedy <- read.table("4.2 proporcao_fastgreedy4.txt")
proporcao_louvain <- read.table("4.2 proporcao_louvain4.txt")
proporcao_walk <- read.table("4.2 proporcao_walk4.txt")

d <- seq(0.2,0.68,0.02)
tabela <- data.frame ("Prop" = round(proporcao$x,3), 
                      "Prop boot" = round(proporcaoboot$x,3),
                      "Prop fast greedy" = round(proporcao_fastgreedy$x,3),
                      "Prop louvain" = round(proporcao_louvain$x,3),
                      "Prop walktrap" = round(proporcao_walk$x,3),
                      row.names = d)
print(tabela)
x <- d
y1 <- (proporcao$x)
y2 <- (proporcaoboot$x)
y3 <- (proporcao_fastgreedy$x)
y4 <- (proporcao_louvain$x)
y5 <- (proporcao_walk$x)
y <- expression(\hat(K))

plot(x,y1,type="b",xlab="a-b", ylab=expression("Taxa de acertos para" ~ hat(K)), ylim=c(0,1), col="#00008B", pch=15)
lines(x,y2,col="#FF007F", type = "b", pch=16)
lines(x,y3,col="#20B2AA", type = "b", pch=17)
lines(x,y4,col="#D2691E", type = "b", pch=18)
lines(x,y5,col="#228B22", type = "b", pch=4)
legend(0.56,0.63, legend=c(expression("Sem bootstrap"),expression("Com bootstrap"), expression("Fast greedy"), expression("Louvain"), expression("Walktrap")),
       col=c("#00008B","#FF007F","#20B2AA","#D2691E","#228B22"), cex=0.7, lty = c(1,1), pch=c(15,16,17,18,4),
       box.lty=0)
