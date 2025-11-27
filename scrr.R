#fixando os parametros e tamanho amostral
n <- 30
ni <- 9
sigma2<- 3
mi <- 20

#gerando amostra t locacao-escala
y <- mi + rt(n, df=ni)*sqrt(sigma2)

l<- function(nu,mu,sig,u,logu){
  nu/2*(sum(logu-u)+n*log(nu/2))-n*lgamma(nu/2)-n/2*log(sig) - sum(u*(y-mu)^2)/(2*sig)
}#fim funcao logverossimilhanca (somente parte dependente de theta)

#funcao de checagem de criterio de parada
crit <- function(l1,nu,mu,sig,u,logu){
  abs(l1-l(nu,mu,sig,u,logu)) > 1e-9
}

CM2 <- function(nu,mu,sig){
  sum(digamma((nu+1)/2) - log((nu+(y-mu)^2/sig)/2)-
        (nu+1)/(nu+(y-mu)^2/sig))/(2*n)+log(nu/2)/2-digamma(nu/2)+1
}

ECM <- function(nu,mu,sig){
  cont<-1
  l1<-100
  u <- 1
  logu <- 1
  d2 <- (y-mu)^2/sig
  while(crit(l1,nu,mu,sig,u,logu)){
    #--------PASSO E--------#
    u <- (nu+1)/(nu+(y-mu)^2/sig)
    logu <- digamma((nu+1)/2) - log((nu+d2)/2)
    #-----------------------#
    l1<- l(nu,mu,sig,u,logu)
    #------PASSO CM-1-------#
    mu <- sum(u*y)/sum(u)
    sig <- sum(u*(y-mu)^2)/n
    #-----------------------#
    #------PASSO CM-2-------#
    nu <- uniroot(CM2,
          interval=c(0.01,10),
          mu=mu,
          sig=sig)$root
    #-----------------------#
    cat("Iteração ",cont,": muhat=",mu," sigma2hat=",sig," nuhat=",nu,"\n")
    cont<-cont+1
  }
  return(c(nuhat=nu,muhat=mu,sigmahat=sig))
}#fim funcao algoritmo EM

# chute inicial usando media e variancia amostrais
ECM(5,mean(y),var(y))
################################################