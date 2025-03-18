

#pacote para manipular grafos/redes
library(igraph)
#pacote com alguns dados de redes
library(igraphdata)

#dados dispon?veis no pacote
data(package="igraphdata")

#vamos usar os dados da rede UKfaculty
data("UKfaculty")


pesos <- E(UKfaculty)$weight/max(E(UKfaculty)$weight)*5
edge.cor <- rep("grey",gsize(UKfaculty))
#lay <- as.matrix(read.table("layout_UKfaculty"))
lay <- layout_nicely(UKfaculty)
plot(UKfaculty,
     layout=lay,
     vertex.size=12,
     vertex.color = "lightblue",
     vertex.label.cex = 0.5,
     edge.color = edge.cor,
     vertex.label.color = "black",
     edge.arrow.size=0.3,
     edge.width=pesos,
     edge.curved = TRUE)

degree_out = degree(UKfaculty, mode = 'out') 
degree_in = degree(UKfaculty, mode = 'in')

par(mfrow=c(1,2))

xmax <- max(hist(degree_in)$breaks,hist(degree_out)$breaks)
ymax <- max(hist(degree_in)$density,hist(degree_out)$density)
hist(degree_out,col="lightblue", 
     main = "", freq=FALSE, 
     xlab="grau de sa?da",ylab="densidade",
     xlim=c(0,xmax), ylim=c(0,ymax))
hist(degree_in,col="lightblue", 
     main = "", freq=FALSE,
     xlab="grau de entrada",ylab="densidade",
     xlim=c(0,xmax),ylim=c(0,ymax))


#Vamos trabalhar com essa rede sem dire??o e sem pesos
un.UKfaculty <- as.undirected(UKfaculty, mode = c("collapse"))

#numero de arestas
gsize(un.UKfaculty)


hist(degree(un.UKfaculty),col="lightblue", 
     main = "Rede da Universidade Inglesa", freq=FALSE,
     xlab="grau",ylab="densidade",ylim=c(0,0.07))

#comunidades verdadeiras
z <- vertex_attr(UKfaculty, index = V(UKfaculty))$Group
col <- vector()
col[which(z==1)] <- "lightblue"
col[which(z==2)] <- "lightsalmon1"
col[which(z==3)] <- "mediumpurple1"
col[which(z==4)] <- "olivedrab3"
#grafo sem dire??o e com as comunidades
plot(un.UKfaculty,
     layout=lay,
     vertex.size=10,
     vertex.color = col,
     vertex.label.cex = 0.55,
     edge.color = edge.cor,
     vertex.label.color = "black")

#pacote para estima??o
library(randnet)

#transformando o objeto igraph em uma matriz de adjac?ncia
A <- as_adjacency_matrix(un.UKfaculty, type = c("both"),
attr = NULL, edges = FALSE, names = TRUE, sparse=0)

#Estimar as comunidades usando Spectral Clustering na rede un.UKfaculty
z.hat.SC <- reg.SP(A,K=4)$cluster

col <- vector()
col[which(z.hat.SC==1)] <- "olivedrab3"
col[which(z.hat.SC==2)] <- "lightsalmon1"
col[which(z.hat.SC==3)] <- "lightblue"
col[which(z.hat.SC==4)] <- "mediumpurple1"
plot(un.UKfaculty,
     layout=lay,
     vertex.size=7,
     vertex.color = col,
     vertex.label.cex = 0.55,
     edge.color = edge.cor,
     vertex.label.color = "black")

#Estimar as comunidades usando Spectral Clustering Regularizado esf?rico na rede un.UKfaculty
z.hat.SC <- reg.SSP(A,K=4)$cluster

col <- vector()
col[which(z.hat.SC==1)] <- "lightsalmon1"
col[which(z.hat.SC==2)] <- "mediumpurple1"
col[which(z.hat.SC==3)] <- "lightblue"
col[which(z.hat.SC==4)] <- "olivedrab3"
plot(un.UKfaculty,
     layout=lay,
     vertex.size=7,
     vertex.color = col,
     vertex.label.cex = 0.55,
     edge.color = edge.cor,
     vertex.label.color = "black")
