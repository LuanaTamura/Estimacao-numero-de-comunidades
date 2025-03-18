

#Funcao que gera as comunidades
#Entrada: n: n. de vertices
# q (pi): prob. de cada comunidade (vetor de dimensao k)
sample.z <- function(n,q){
  k <- length(q)
  return(sample(x = seq(1:k), n, replace = T, prob = q))
}

# Gera a matriz de ajd
# Entradas: z vetor de comunidades
# P: matriz k x k com as prob. de conexao dadas as comunidades
# Pab = prob de conexao entre vertices da comun a e b
sample.SBM <- function(z,P){
  #number of vertices
  n <- length(z)
  mat.adj <- matrix(0,n,n)
  k <- dim(P)[1]
  for(i in 1:n){
    for(j in (i):n){
      mat.adj[i,j] <- rbinom(1,1,P[z[i],z[j]])
      mat.adj[j,i] <- mat.adj[i,j]
    }
  }
  #no self loops
  diag(mat.adj) <- 0
  return(mat.adj)
}

require(igraph)
require(randnet)
#### Exemplo
q <- c(1/3,1/3,1/3) #balanceada
n <- 100
a <- 0.7
b <- 0.3
P <- rbind(c(a,b,b),c(b,a,b),c(b,b,a))

#gerar as comunidades verdadeiras
z <- sample.z(n,q)
A <- sample.SBM(z,P)

S <- 100 #numero de vezes para realizar a simulação

#k.est <- c(1,2,3,4,NA,5,6,7)
#S <- 7
#for(i in 1:S){
#z <- sample.z(n,q)
#A <- sample.SBM(z,P)
#teste k.est[i] <- teste..(A)
#print(i)
#if(is.na(k.est[i])) break
#}

getwd()
