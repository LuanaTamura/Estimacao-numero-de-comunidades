set.seed(030699)

data <- read.csv("date-hour-soo-dest-2020.csv")
data <- read.csv("date-hour-soo-dest-2021.csv")
data <- read.csv("date-hour-soo-dest-2022.csv")
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

isSymmetric(Adj)
#transformar em simétrica
Adj_sym <- Adj + t(Adj)
isSymmetric(Adj_sym)

#colocar a diagonal 0
diag(Adj_sym) <- 0
Adj_transf <- ceiling(log(Adj_sym+1))
Adj_zero <- floor(log(Adj_sym+1))

#gerando imagem do grafo
vert_comunidades <- SpecClust(Adj_zero, 2)
table(vert_comunidades)

require("igraph")
g <- graph_from_adjacency_matrix(Adj_zero, mode = "undirected", weighted = TRUE)

edge.cor <- rep("grey",gsize(g))
#lay <- as.matrix(read.table("layout_UKfaculty"))
#lay <- layout_nicely(g)
col <- vector()
col[which(vert_comunidades==1)] <- "mediumpurple1"
col[which(vert_comunidades==2)] <- "lightsalmon1"
#col[which(vert_comunidades==1)] <- "lightsalmon1"
#col[which(vert_comunidades==2)] <- "mediumpurple1"
plot(g,
     layout=lay,
     vertex.size=15,
     vertex.color = col,
     vertex.label.cex = 1.2,
     edge.color = edge.cor,
     vertex.label.color = "black")
