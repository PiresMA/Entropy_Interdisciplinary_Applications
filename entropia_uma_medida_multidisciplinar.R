

#-------------------------------------------------#
#-------------------------------------------------#
##### Entropia: uma medida multidisciplinar  ######
#-------------------------------------------------#
#-------------------------------------------------#

#-------------------------------------------------#
############ Entropia de Boltzmann-Gibbs  #########
#-------------------------------------------------#
# preparar grafico
par(mfrow=c(1,1))
Nmax=10; Wmax=fac(Nmax)/(fac(Nmax/2)*fac(Nmax/2))
maxY=ceiling(log(Wmax,  base = exp(1)))
plot((0:maxY)/maxY,0:maxY , type="n", xlab="ne/N", ylab="H")

W=ne=nd=p=c(NULL); q=1
for(N in seq(2,Nmax,2)){
  # inicializacao para cada iteracao
  ne[1]=N; nd[1]=0
  totEst=2^(N)
  W[1]=factorial(N)/(factorial(ne[1])*factorial(nd[1]))
  p[1]=W[1]/totEst
  # expansao do gas: 
  # uma particula passa da metade esquerda para a direira
  for(c in 1:N){
    ne[c+1]=N-c
    nd[c+1]=c
    W[c+1]=fac(N)/(fac(ne[c+1])*fac(nd[c+1]))
    p[c+1]=W[c+1]/totEst    
  }
  # calculo da entropia  
  H=log(W, base = exp(1))
  #plot
  lines(ne/N, H, type="o", col=q)  
  q=q+1
}
legend( "topleft", legend=sprintf("%s%d","N=",seq(2,Nmax,2)),
        col=1:(q-1) , pch=1)


#-------------------------------------------------#
############ Entropia na estatística   ############
#-------------------------------------------------#
x=rbinom(n=1000, size=50, prob=0.5)
barplot(table(x))

# dados amostrais
N=10^(4)
x1=rnorm(N,0,1)
x2=rnorm(N,0,0.5)

#calculo das entropias
H1=round(getStandartEntropy(x1, bin=0.5),2)
H2=round(getStandartEntropy(x2, bin=0.5),2)

#visualizar
par(mfrow=c(1,2))
hist(x1, xlim=c(-5,5), ylab="f", xlab="x",
     main=sprintf("%s%.2f","H1=",H1) )
hist(x2, xlim=c(-5,5), ylab="f", xlab="x",
     main=sprintf("%s%.2f","H2=",H2) )

#-------------------------------------------------#
############ Entropia e diversidade   ############
#-------------------------------------------------#

# Geração dos sistemas I e II
sis_I=c(3,3)
sis_II=c(2,2,2)

p=sis_I/sum(sis_I)
H1=-sum(p*log(p))
p=sis_II/sum(sis_II)
H2=-sum(p*log(p))

par(mfrow=c(1,2))
barplot(sis_I,ylim=c(0,3), names.arg=c("A","B"),
        xlab="Sistema I", ylab="número de indivíduos", 
        main=sprintf("%s%.2f","H1=",H1) )
barplot(sis_II,ylim=c(0,3),names.arg=c("A","B","C"),
        xlab="Sistema II", ylab="número de indivíduos",
        main=sprintf("%s%.2f","H2=",H2) )


#-------------------------------------------------#
############ Entropia e Informação     ############
#-------------------------------------------------#
# jogo de moeda
# declarações e inicializações
N=100
p=H=c(NULL)
H[1]=0; H[N]=0
p[1]=0; p[N]=1
#iteração
for( k in 2:(N-1)){
  #probabilidade de cara
  p[k]=k/N
  #probabilidade de coroa
  q=1-p[k]
  #entropia
  H[k]=-p[k]*log(p[k],base=2)-q*log(q,base=2)
}
#gráfico
plot(p,H, col="red")


#-------------------------------------------------#
############ Entropia em Nanociencia   ############
#-------------------------------------------------#

# perfil 1
perfil1=c(1,2,1); x=perfil1
namesBar=1:length(x);
barplot(x, names.arg=namesBar)
# calculo de probabilidades
p=x/sum(x)
#calculo da entropia normalizada
H1=-sum(p*log(p))/log(length(p))

# perfil 2
perfil2=c(1,2,2); x=perfil2
barplot(x, names.arg=namesBar)
# calculo de probabilidades
p=x/sum(x)
#calculo da entropia normalizada
H2=-sum(p*log(p))/log(length(p))

# perfil 3
perfil3=c(2,2,2); x=perfil3
barplot(x, names.arg=namesBar)
# calculo de probabilidades
p=x/sum(x)
#calculo da entropia normalizada
H3=-sum(p*log(p))/log(length(p))

par(mfrow=c(1,3))
barplot(perfil1, names.arg=namesBar, ylab="frequência", xlab="posições")
barplot(perfil2, names.arg=namesBar, ylab="frequência", xlab="posições")
barplot(perfil3, names.arg=namesBar, ylab="frequência", xlab="posições")

# mostrar valores
round(c(H1,H2,H3),3)

#-------------------------------------------------#
############ Entropia em Redes Sociais ############
#-------------------------------------------------#
#Inicialização
num_cell<-4
mat_adj<-matrix(0 ,nrow=num_cell, ncol=num_cell)

# Criação da rede
mat_adj[1,2]<-1
mat_adj[1,3]<-1
mat_adj[2,3]<-1
mat_adj[2,4]<-1
mat_adj[3,4]<-1
g <- graph.adjacency( mat_adj, mode=c("undirected") )
E(g)$weight=c(0.3,0.3,0.6,0.01,0.99)

# Cálculo da entropia
round(graph.diversity(g), 2)

# Visualizar
par(mfrow=c(1,1))
plot(g, layout=layout.circle,
     edge.label=E(g)$weight)

###-----------------------------#
