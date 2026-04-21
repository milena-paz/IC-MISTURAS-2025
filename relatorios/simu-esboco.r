rMisturaNorm <- function(n,mi,dp,probs){
  k<- length(mi)
  X<- numeric(n)
  U <- runif(n)
  p <- cumsum(c(0,probs))
  for(i in 1:k){
    caso <- U <= p[i+1] & U> p[i]
    tam <- sum(caso)
    X[caso] <- rnorm(tam,mi[i],dp[i])
  }
  return(X)
}#Fim função geradora
p <- c(0.64,0.36)
mu <- c(160,173)
sigma <- c(15,12)
y<- rMisturaNorm(300,mu,sqrt(sigma),p)

hist(y,freq=F,breaks=20)

#. Escolha do Valor Inicial

#Muitas vezes o algoritmo iterativo será muito sensível aos valores iniciais
#escolhidos para realizar a aproximação. Isso se deve à sua convergência tender
# a máximos locais, que podem não ser os máximos globais, que pode levar a um erro.
# Uma forma de obter "chutes iniciais" mais educados é por meio de uma análise
# exploratória dos dados antes mesmo de codificar o algoritmo. Entre as possíveis
# técnicas a serem usadas, temos densidades kernel, gráficos de dispersão, quantis,
# e histogramas. Por exemplo, em 2.1,2.2 e 2.3, o histograma das amostras foi usado
# para estimar valores iniciais para os parâmetros de locação.

#Aqui iremos avaliar a sensibilidade à inicialização do modelo do exemplo 2.2
#(normal, g=2) observando a convergência para #$N=100$ repetições do algoritmo EM com diferentes valores iniciais.

#y <- r



p <- c(0.64,0.36)
mu <- c(160,173)
sigma <- c(15,12)
#Usando a amostra de 2.1

component <- function(x,mu,sig,g){
  dnorm(x,mu[g],sqrt(sig[g]))
}
l <- function(p,mu,sigma){
  sum(log(p[1]*dnorm(y,mu[1],sigma[1])+p[2]*dnorm(y,mu[2],sigma[2])))
}
crit <- function(ll,p,mu,sig){
  ll2 <- l(p,mu,sig)
  val <-abs(ll2-ll)/(abs(ll) + 1e-6) > 1e-6
  if(is.na(val))
    return(F)
  return(val)
}

EM <- function(p,mu,sig){
  tempo <- Sys.time()
  l1<-100
  z <- 1
  cont <- 1
   #adicionando um limite de iteracoes porque vai divergir provavelmente
  while(crit(l1,p,mu,sig)& cont < 50){
    z <- sapply(y,function(y){
      sapply(1:2,function(j) p[j]*component(y,mu,sig,j)/sum(p[1]*component(y,mu,sig,1),
                                                            p[2]*component(y,mu,sig,2)))
    })
    l1<- l(p,mu,sig)
    p <- apply(z,1,mean)
    mu <- apply(z,1,function(z) sum(z*y)/sum(z))
    sig <- sapply(1:2, function(j)
      sum(z[j,]*(y-mu[j])^2)/sum(z[j,]))
    cont <- cont +1
  }
  tempo <- Sys.time() - tempo
  return(list(cbind(p=p,mu=mu,sig=sig),tempo,cont==50))
}

N<-10
prmt <- matrix(c(p,mu,sigma),ncol=3)
#cores <- RColorBrewer:::brewer.pal(6,"Set3")
n<-300
simula <- function(p0=p,mu0=mu,sig0=sigma){
  saida <- replicate(N, EM(p0,mu0,sig0))
  diffs <- sapply(saida[1,], "-", prmt)
  eqm <- matrix(apply(diffs,1,function(d) sum(d^2)/n),nrow=2)
  colnames(eqm) <- c("p","mu","sigma2")
  cat("Valores iniciais:\n p0: ",p0," mu0: ",mu0," sig0: ",sig0)
  cat("\n Tempo de execução médio(seg.): ",mean(unlist(saida[2,])),
      "\n Número de Divergências: ",sum(unlist(saida[3,])),
      "\n Erro Quadrático Médio:\n")
  print(eqm)
}
#para mu0s aleatorios:
mu0<-matrix(runif(10,100,200),nrow=2)
res <-apply(mu0,2,function(m) simula(mu0=m))
#para sigma0s aleatorios:
mu0<-matrix(runif(10,100,200),nrow=2)
apply(mu0,2,function(m) simula(sig0=m))
#para p0s aleatorios:
mu0<-matrix(runif(10,100,200),nrow=2)
apply(mu0,2,function(m) simula(mu=m))