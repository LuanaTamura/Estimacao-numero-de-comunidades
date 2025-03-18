set.seed(030699)

data <- read.csv("date-hour-soo-dest-2022.csv")

names(data)

nomes.origem <- names(table(data$Origin.Station))
length(nomes.origem)
nomes.destino <- names(table(data$Destination.Station))
length(nomes.destino)
# Temos os mesmos nomes de origem e destino - ok
length(intersect(nomes.destino,nomes.origem))


tabela <- cbind(seq(1,50),nomes.origem)
colnames(tabela) <- c("vertice","nome_estacao")
head(tabela)
write.table(tabela,"ID_vertices_2022.txt",row.names = FALSE)

nomes <- nomes.origem

#matriz de adjacencia com os pesos
Adj <- matrix(0,nrow=50,ncol=50)

for(k in 1:dim(data)[1]){
  origem <- data$Origin.Station[k]
  i <- which(nomes==origem)
  destino <- data$Destination.Station[k]
  j <- which(nomes==destino)
  Adj[i,j] <- Adj[i,j] + data$Trip.Count[k]
  print(paste("linha ", k, " de ", dim(data)[1]))
}

write.csv(Adj,"Adj_sym_2022.csv")

write.table(Adj, file = "Adj_sym_2022.txt",
            row.names = FALSE, col.names = FALSE)

Adj_teste = read.table("Adj_sym_2022.txt")
matriz = as.matrix(Adj_teste)
isSymmetric(matriz)
matriz_sym <- matriz + t(matriz)
isSymmetric(matriz_sym)


isSymmetric(Adj)
#transformar em simétrica
Adj_sym <- Adj + t(Adj)
isSymmetric(Adj_sym)

#colocar a diagonal 0
diag(Adj_sym) <- 0
Adj_transf <- ceiling(log(Adj_sym+1))
Adj_zero <- floor(log(Adj_sym+1))

#peso da aresta e histograma
hist(Adj_sym, 
     col="#20B2AA",
     ylab="Frequência",
     xlab="pesos de cada aresta",
     ylim=c(0,1700),
     main=" ")
hist(Adj_zero, 
     col="#20B2AA",
     ylab="Frequência",
     xlab="pesos de cada aresta",
     ylim=c(0,1700),
     main=" ")

#peso do vértice e histograma
grau <- apply(Adj_sym,1,mean)
grau_transf <- apply(Adj_transf,1,mean)
grau_zero <- apply(Adj_zero,1,mean)
#par(mfrow = c(1,2))
hist(grau, col="#20B2AA",
     ylab="Frequência",
     xlab="pesos de cada vértice",
     ylim=c(0,20),
     main=" ")
hist(grau_zero, 
     col="#20B2AA",
     ylab="Frequência",
     xlab="pesos de cada vértice",
     ylim=c(0,20),
     main=" ")


#testando alguns métodos para estimar k
require("igraph")
library(randnet)

g <- graph_from_adjacency_matrix(Adj_zero, mode = "undirected", weighted = TRUE)

# Fast Greedy method (aplica-se apenas a grafos não direcionados)
fg_communities <- cluster_fast_greedy(g)
num_fg_communities <- length(unique(membership(fg_communities)))
num_fg_communities #2
table(membership(fg_communities))

# Louvain method (funciona com grafos não direcionados)
louvain_communities <- cluster_louvain(g)
num_louvain_communities <- length(unique(membership(louvain_communities)))
num_louvain_communities #2
table(membership(louvain_communities))

# Walktrap method (funciona com grafos direcionados)
walktrap_communities <- cluster_walktrap(g)
num_walktrap_communities <- length(unique(membership(walktrap_communities)))
num_walktrap_communities #2
table(membership(walktrap_communities))


#randnet normal
z.hat <- reg.SP(Adj_zero,K=2)$cluster
table(z.hat)

#randnet ESFERICO
z.hat.SC <- reg.SSP(Adj_zero,K=2)$cluster
table(z.hat.SC)




result <- SeqTestSBM(Adj_zero,K.min=1,K.max=10,n.b=10000)
result #1

result <- SeqTestSBM(Adj_zero,K.min=1,K.max=10,n.b=1000)
result #1

result <- SeqTestSBM(Adj_zero,K.min=1,K.max=12,n.b=100)
result #2

vert_comunidades <- SpecClust(Adj_zero, 2)
table(vert_comunidades)


edge.cor <- rep("grey",gsize(g))
#lay <- layout_nicely(g) #rodar uma vez só
col <- vector()
col[which(vert_comunidades==2)] <- "lightsalmon1"
col[which(vert_comunidades==1)] <- "mediumpurple1"
plot(g,
     layout=lay,
     vertex.size=15,
     vertex.color = col,
     vertex.label.cex = 1.2,
     edge.color = edge.cor,
     vertex.label.color = "black")

#sera q as arestas estao com os pesos?

pesos = E(g)$weight/max(E(g)$weight)
plot(g,
     layout=lay,
     vertex.size=15,
     vertex.color = col,
     vertex.label.cex = 0.55,
     edge.color = edge.cor,
     vertex.label.color = "black",
     edge.arrow.size=0.3,
     edge.width=pesos)


install.packages("igraphdata")
library(igraphdata)
