
#### bike
#### https://s3.amazonaws.com/capitalbikeshare-data/index.html

# Pode ser feito mês e mês para comparar se os grupos mudam ao longo do ano
# selecionar um ano
# Mes 1
data <- read.csv("202401-capitalbikeshare-tripdata.csv")

names(data)

dim(table(data$start_station_name))

nomes.origem <- names(table(data$start_station_name))
length(nomes.origem)
nomes.destino <- names(table(data$end_station_name))
length(nomes.destino)
# obter a interseção de origem e destino
length(intersect(nomes.destino,nomes.origem))

nomes <- intersect(nomes.destino,nomes.origem)

length(nomes)
#matriz de adjacencia com os pesos
Adj <- matrix(0,nrow=length(nomes),ncol=length(nomes))

for(k in 1:dim(data)[1]){
  origem <- data$start_station_name[k]
  i <- which(nomes==origem)
  destino <- data$end_station_name[k]
  j <- which(nomes==destino)
  Adj[i,j] <- Adj[i,j] + 1
  print(paste("linha ", k, " de ", dim(data)[1]))
}

dim(Adj)

isSymmetric(Adj)
#transformar em simétrica
Adj_sym <- Adj + t(Adj)
#colocar a diagonal 0
diag(Adj_sym) <- 0
Adj_sym <- log(Adj_sym+1)

isSymmetric(Adj_sym)

#testando alguns métodos para estimar k
require("igraph")

g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)

# Fast Greedy method (aplica-se apenas a grafos não direcionados)
fg_communities <- cluster_fast_greedy(g)
num_fg_communities <- length(unique(membership(fg_communities)))
num_fg_communities 

# Louvain method (funciona com grafos não direcionados)
louvain_communities <- cluster_louvain(g)
num_louvain_communities <- length(unique(membership(louvain_communities)))
num_louvain_communities 

# Walktrap method (funciona com grafos direcionados)
walktrap_communities <- cluster_walktrap(g)
num_walktrap_communities <- length(unique(membership(walktrap_communities)))
num_walktrap_communities


a <- max(Adj_sym)/10
A <- ceiling(Adj_sym/a)
result <- SeqTestSBM(A,K.min=40,K.max=70, n.b=1000)
result





############ BART -- Kaggle
############ https://www.kaggle.com/datasets/mrgeislinger/bart-ridership

data <- read.csv("date-hour-soo-dest-2023.csv")

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
write.table(tabela,"ID_vertices_2023.txt",row.names = FALSE)

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

write.table(Adj, file = "Adj_sym_2023.txt", sep=";",
            row.names = FALSE, col.names = FALSE)

Adj_teste = read.table("Adj_sym_2023.txt")
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

table(Adj_zero==0)

#testando alguns métodos para estimar k
require("igraph")

g <- graph_from_adjacency_matrix(Adj_zero, mode = "undirected", weighted = TRUE)

# Fast Greedy method (aplica-se apenas a grafos não direcionados)
fg_communities <- cluster_fast_greedy(g)
num_fg_communities <- length(unique(membership(fg_communities)))
num_fg_communities

# Louvain method (funciona com grafos não direcionados)
louvain_communities <- cluster_louvain(g)
num_louvain_communities <- length(unique(membership(louvain_communities)))
num_louvain_communities 

# Walktrap method (funciona com grafos direcionados)
walktrap_communities <- cluster_walktrap(g)
num_walktrap_communities <- length(unique(membership(walktrap_communities)))
num_walktrap_communities



result <- SeqTestSBM(Adj_zero,K.min=1,K.max=10,n.b=10000)
result

result <- SeqTestSBM(Adj_zero,K.min=1,K.max=15,n.b=1000)
result #2

result <- SeqTestSBM(Adj_zero,K.min=1,K.max=30,n.b=100)
result #2


vert_comunidades <- SpecClust(Adj_zero, 2)
table(vert_comunidades)


edge.cor <- rep("grey",gsize(g))
#lay <- as.matrix(read.table("layout_UKfaculty"))
#lay <- layout_nicely(g)
col <- vector()
col[which(vert_comunidades==1)] <- "lightsalmon1"
col[which(vert_comunidades==2)] <- "mediumpurple1"
plot(g,
     layout=lay,
     vertex.size=15,
     vertex.color = col,
     vertex.label.cex = 1.2,
     edge.color = edge.cor,
     vertex.label.color = "black")
