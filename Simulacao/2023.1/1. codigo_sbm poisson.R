set.seed(030699)

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
sample.g <- function(z,P){
  #number of vertices
  n <- length(z)
  mat.adj <- matrix(0,n,n)
  k <- dim(P)[1]
  for(i in 1:n){
    for(j in (i):n){
      mat.adj[i,j] <- rpois(1,P[z[i],z[j]])
      mat.adj[j,i] <- mat.adj[i,j]
    }
  }
  #no self loops
  diag(mat.adj) <- 0
  return(mat.adj)
}



