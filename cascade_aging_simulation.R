
## note: figure 7 is skippe din the present versiomn of the paper
# thus figure 8 here is figure 7 in the paper
# thus figure 9 here is figure 8 in the paper



####################
##  Figure 4: HSBM
####################

library(igraph)
library(colorspace)
library(qgraph)

hsbm <- function (cluster_sizes,strengths=NULL,plot=T)
{
  k=cluster_sizes # cluster sizes
  levels=length(k) # hierarchic levels
  n=prod(k) # nodes
  if(n > 2000) print('Warning: ', n, 'nodes', 'might be too much for qgraph')
  if(length(k)!=length(strengths)) { print ('cluster_sizes and strengths dont match');stop()}
  group=rep(0,n) # node group
  m=matrix(levels,n,n) # connection matrix
  for(l in (levels-1):1)
    for(i in 1:n)
    {for(j in 1:n)
      if((i-1)%/%prod(k[1:l])==(j-1)%/%prod(k[1:l])) m[i,j]=l
    group[i]=(i-1)%/%prod(k[1:l])
    }
  
  if(plot) image(1/(m+.5))
  
  if(length(strengths)<1) m <- 1/m^3 else
  {
    for(i in 1:levels)
      m[m==i]=strengths[i]
  }
  return(list(m,group))
}

#clusters=c(4,4,4,4,4);strengths = c(1,.3,.1,.02,.005)
clusters=c(4,4,4,4);strengths = c(.8,.1,.04,.008)
clusters=c(4,3,6,5);strengths = c(.8,.1,.04,.02)


h=hsbm(clusters,strengths)
m <- h[[1]]
group=h[[2]]

name2=paste('network',paste0(clusters,collapse=''),'.jpeg',sep='',col='')
jpeg(name2,width = 5, height = 3, units = 'in', res = 300)
layout(t(1))
#qgraph(m,groups=group, layout="spring",labels=F,vsize=1,edge.width=.4, bg=adjustcolor("white", alpha.f=0),border.width=.1) # collor version
qgraph(m,groups = group,layout = "spring",labels = FALSE,vsize = 1,edge.width = 0.3,bg = "white",palette='grey',border.color = "black",edge.color = "gray",label.color = "black",edge.width=.6) 
dev.off()




####################
## figure 6: An illustrative simulation of the cascading transitions model for aging. 
####################

# packages required
library(deSolve)
library(rootSolve)
library(FME)
#install.packages("remotes") # if grind is not installed
#remotes::install_github("hansschepers/grindr") # if grind is not installed
library(Grind) 
library(igraph)
library(qgraph)
library(colorspace)


set.seed(1)
N <- 90 # dimensions (nr of nodes) of aging
clusters <- 3 # in 3 clusters (physical, psychological, social)
p_within <- 20/N;p_between <- .003;rewiring <- .02 # parameters for Stochastic block model (SBM)
pm <- matrix(p_between,clusters,clusters);
diag(pm) <- p_within

# create stochastic block model using sample_sbm from Igraph
g <- sample_sbm(N, pref.matrix=pm, block.sizes=rep(N/clusters,clusters))
l <-layout_nicely(g)
adj <- as_adjacency_matrix(g,sparse=F)
adj <- as.matrix(adj)
clusters <- cluster_fast_greedy(g) # detect clusters in network
member <- membership(clusters) # assign memberships 


# cascading transition model (for odesolver in grind)
model <- function(t, state, parms){
  with(as.list(c(state,parms)),{
    X <- state[1:N]
    a0 <- parms[1:N]
    dX <-  a0 + (a* adj) %*% X + b*X -X^3 # connected cusps (equation 7 in paper)
    return(list(dX))
  })
}

# parameter settings
a0=rep(0,N) # intercepts of normal
a=matrix(rnorm(N*N,.4,0),N,N) # couplings effects on normal
b=rnorm(N,1,1) # random splitting values
b[1]=-1 # to illustrate continuous decline by setitng the splitting value  low

# 20% of the variables will decrease slowly after t = 20 (see run())
degrading=-.2 * sample(0:1,N,T,prob=c(.8,.2))

s.m=2 # for plotting y-axis
tmax=300 # time duration
X <- rep(1,N) # initials state of X
s <- c(X);p <- c(a0) # required for grind
# run the model:
Xm=as.matrix(run(tmax=tmax,ymin=-s.m,ymax=s.m,
                 after="if(t>20)parms<-parms+degrading;
                 state<-state+rnorm(N,mean=0,sd=0.05)",timeplot=FALSE,table=TRUE))
# Xm contains the values of all nodes over all time points

# plot settings
pdf('figure6.pdf',h=6,w=9)
widths=rep(1,12) 
widths[1]=.3
height=c(rep(.5,6),rep(1,6))
layout(matrix(c(1,2,3,4,5,6,7,7,7,7,7,7),2,6,byrow=T),widths=widths)
par(mar=c(0,0,0,0))
min.o=-3;max.o=3
t_maxvar=which.max(apply(Xm[,-1],1,var))

# create empty plot
plot(1, type="n", axes=FALSE, ann=FALSE)

# plot network at 5 points in time
for (t in c(10,tmax/4,2*tmax/4, 3*tmax/4,tmax))
{
  green_to_red <- colorRampPalette(c("green", "red"))
  color_palette <- green_to_red(256)
  colored_values <- color_palette[as.numeric(cut(-1 * Xm[t, -1], 
                                                 seq(min.o, max.o, length = 255),
                                                 include.lowest = TRUE))]
  qgraph(adj, layout = l, labels = FALSE, color=colored_values,vsize=3.5,title=paste('    time =',t),border.width=.1)
}

#plot the time series for all nodes
par(cex=1.1,mar=c(3,4,0,1),las=1,mgp=c(1.5, 1, 0))
matplot(Xm[2:tmax,-1],type='l',lty=1,lwd=.5,col=3*(degrading!=0)+1, ylab='',xlab='time',bty='n',axes=F)
lines(Xm[2:tmax,2],col='purple',lwd=1.5)
axis(1,at=c(0,100,200,300))
axis(2,at=c(2,0,-2))
mtext('X',2,1.5,at=1)

dev.off() # close plot file

####################
## figure 7
####################
library(igraph)
library(colorspace)
library(qgraph)
library(png)
library(ggplot2)
library(grid)
library(gridGraphics)

# Energy of a particular state of the network (equation 8)
hamiltonian=function(x,t,w) -sum(t*x)-sum(w*x%*%t(x)/2)

# Glauber dynamics (follows from equation 9)
glauber_step = function(x,n,t,w,beta)
{
  i = sample(1:n,size=1) # take a random node
  x_new=x;x_new[i]=x_new[i]*-1 # construct new state with flipped node
  p=1/(1+exp(beta*(hamiltonian(x_new,t,w)-hamiltonian(x,t,w))))  # update probability
  if(runif(1)<p) x=x_new # update state
  return(x)
}

# create (hierarchical) SBM
hsbm <- function (cluster_sizes,strengths=NULL,plot=T)
{
  k=cluster_sizes # cluster sizes
  levels=length(k) # hierarchic levels
  n=prod(k) # nodes
  if(n > 2000) print('Warning: ', n, 'nodes', 'might be too much for qgraph')
  if(length(k)!=length(strengths)) { print ('cluster_sizes and strengths dont match');stop()}
  group=rep(0,n) # node group
  m=matrix(levels,n,n) # connection matrix
  for(l in (levels-1):1)
    for(i in 1:n)
    {for(j in 1:n)
      if((i-1)%/%prod(k[1:l])==(j-1)%/%prod(k[1:l])) m[i,j]=l
    group[i]=(i-1)%/%prod(k[1:l])
    }
  
  if(plot) image(1/(m+.5))
  
  if(length(strengths)<1) m <- 1/m^3 else
  {
    for(i in 1:levels)
      m[m==i]=strengths[i]
  }
  return(list(m,group))
}

iterations=30000

    set.seed(1)
    clusters=c(30,6);strengths = c(.3,.0015) # creates the SBM
    t=seq(-1.5,-5,length=iterations)
    beta=2
    name1=''
    name2='Stochastic block\nmodel (SBM)'
  
  #build network  
  h=hsbm(clusters,strengths)
  m <- h[[1]]
  group=h[[2]]
  
  # delete 20% of nodes
  del=sample(1:nrow(m),round(nrow(m)/5))
  m=m[-del,-del]
  group=group[-del]
  
  
  n <- nrow(m) # nr of nodes
  x=rep(1,n) # start values positive
  
  m=ifelse(m>matrix(runif(n*n,0,1),n,n),1,0)
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m)=0
  group <- as_membership(cluster_louvain(
    graph_from_adjacency_matrix(m, mode = "undirected")))$membership
  
  # run simulation
  s=rep(0,iterations)
  
  ii=iterations/5
  xi=matrix(NA,floor(iterations/ii),n) # collect x values 
  j=0
  for(i in 1:iterations)
  {  x<-glauber_step(x,n,t=t[i],m,beta)
  s[i]=mean(x)
  if(i%%ii==0)
  {
    j=j+1
    print(i)
    xi[j,]=x
  }
  }
  
  ## creating figure requires advances plotting tricks
  
  capturePlot <- function() {
    grid.echo()
    grid.grab()
  }
  plot.new()
  plotGrob <- vector("list", length = j)
  for(jj in 1:j)
  {
    qgraph(m,color=2+(xi[jj,]+1)/2, layout="spring",labels=F,
           vsize=.4,edge.width=.3,bg=adjustcolor("white", alpha.f=0),border.width=.1)
    plotGrob[[jj]] <- capturePlot()
            print(jj)
  }
  qgraph(m,groups=group, layout="spring",labels=F,
         vsize=.6,edge.width=.3, bg=adjustcolor("white", alpha.f=0),border.width=.1) 
  plotGrob1 <- capturePlot()
  
  png("temp_plot.png",w=800,h=500,res=200)
  par(mar=c(5,5,3,3),mgp=c(1.5, 1, 0))
  plot(t,type='l',xlab='time',ylab=expression(tau),bty='o',
       lwd=1.5,col='black',cex.axis=1.5,cex.lab=1.2,
       axes=F,cex.main=1.2,
       main=paste('Decrease in external field'),
       bg=adjustcolor("white", alpha.f=0),border.width=.1)
  axis(1,labels=FALSE);axis(2,at=range(t))
  dev.off()
  img <- readPNG("temp_plot.png")
  plotGrob2 <- rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))
  
  
  pdf(paste('figure7.pdf',sep='',collapse=''),h=8,w=12)
  par(mar=c(5,6,3,6))
  plot(s,bty='n',type='l',xlab='time',ylab=expression(bar(x)),
       ylim=c(-1.2,1.1),lwd=5,col='blue',cex.axis=1.5,cex.lab=1.5,las=1,
       axes=F,cex.main=1.5,main=name1)
  axis(1,labels=FALSE);axis(2)
  text(iterations/6,-1.1, name2)
  
  for(jj in 1:j)
  {
    vp <- viewport(x=.9*(jj*ii/iterations), y=-.1+(.6*s[ii+(jj-1)*ii]+1)/2, 
                   width=0.2, height=0.35, just=c("center","bottom"))
    pushViewport(vp)
    grid.draw(plotGrob[[jj]])
    popViewport()
  }
  
  vp <- viewport(x=.25, y=.18, 
                 width=0.25, height=.35, just=c("center","bottom"))
  pushViewport(vp)
  grid.draw(plotGrob1)
  popViewport()
  
  vp <- viewport(x=.83, y=.65,
                 width=0.28, height=.28, just=c("center","bottom"))
  pushViewport(vp)
  grid.draw(plotGrob2)
  popViewport()
  
  dev.off()


####################
## figure 8
###################

library(igraph)
library(colorspace)
library(qgraph)
library(png)
library(ggplot2)
library(grid)
library(gridGraphics)

# Returns the indices in the connection matrix of the node
# that should be removed 
getMatrixIndices <- function(matrix, index) {
  nrows <- nrow(matrix)
  row_index <- ((index - 1) %% nrows) + 1
  col_index <- ((index - 1) %/% nrows) + 1
  return(c(row_index,col_index))
}


## hamiltonian: Computes the Ising model’s potential function 
# Input:
## x: the state vector with the network of nodes
## n:
## t: the tau parameter (external magnetic field, can be positive or negative)
## w: TODO: explanation about function of w (kind of interaction matrix to specify
#    the neighboring nodes/states in X to sum over. In formula division by 2 because
#    the network is undirected and each connection needs to counted once)
# Returns: the energy of X (the network of nodes)
hamiltonian=function(x,t,w) -sum(t*x)-sum(w*x%*%t(x)/2)


# Update of the nodes in the network with the probability distribution
# that is based on the potential function (Hamiltonian in the model)
glauber_step = function(x,n,t,w,beta)
{
  i = sample(1:n,size=1) # take a random node
  x_new=x;x_new[i]=x_new[i]*-1 # construct new state with flipped node
  p=1/(1+exp(beta*(hamiltonian(x_new,t,w)-hamiltonian(x,t,w))))  # update probability
  if(runif(1)<p) x=x_new # update state
  return(x)
}

# create Stochastic Block Model (hsbm). 
# Input:
## cluster_sizes: A vector with as length the number of hierarchical levels/clusters
##                and the values specify the number of nodes in each cluster. 
##                *goes from higher levels to lower levels?*
## strengths: The strength for each cluster. Default is NULL. 
## plot: Default is TRUE
# Returns: 
## A list of length 2. First element is m (connection matrix of the network) and 
## the second element is group (for each node n, to which group it belongs). m is 
## the size of the product of cluster sizes (e.g., cluster_sizes = c(4,2), nrow(m) = ncol(m) = 8)
# and elements of m have the strenght values. group is a vector of length product cluster sizes

hsbm <- function (cluster_sizes,strengths=NULL,plot=T)
{
  k=cluster_sizes # cluster sizes
  levels=length(k) # hierarchic levels
  n=prod(k) # nodes
  
  # only possible to *compute* hsbm over aging networks smaller than 2000 nodes
  if(n > 2000) print('Warning: ', n, 'nodes', 'might be too much for qgraph')
  if(length(k)!=length(strengths)) { print ('cluster_sizes and strengths dont match');stop()}
  
  # all nodes initialized to belong to group 0
  group=rep(0,n) # node group
  m=matrix(levels,n,n) # connection matrix
  
  # loop over all hierarchical levels, from the highest level to
  # the lowest level (i.e., the smallest sub cluster)
  for(l in (levels-1):1)
    
    # there is looped over each pair of nodes i,j 
    for(i in 1:n)
    {for(j in 1:n)
      if((i-1)%/%prod(k[1:l])==(j-1)%/%prod(k[1:l])) m[i,j]=l
    group[i]=(i-1)%/%prod(k[1:l])
    }
  
  if(plot) image(1/(m+.5))
  
  if(length(strengths)<1) m <- 1/m^3 else
  {
    for(i in 1:levels)
      m[m==i]=strengths[i]
  }
  return(list(m,group))
}


## Simulation of the aging model 

nr_of_persons = 10 # the number of persons for figure 7B 
iterations = 50000 # the number of iterations we want to use in our network 

clusters=c(30,6); strengths = c(.3,.0015)
t=seq(-1.5,-5,length=iterations)
beta=2
name1=''
name2='Stochastic block\nmodel (SBM)'

# initialize everything to build the stochastic block model network for each person

mis <- list(nr_of_persons)
linksis <- list(nr_of_persons)
xis <- list(nr_of_persons)
qs <- list(nr_of_persons)
m_origs <- list(nr_of_persons)
ss <- list(nr_of_persons)

# build for each person the network 
for (person in 1:nr_of_persons) {
  
  #build network  
  h=hsbm(clusters,strengths)
  m <- h[[1]]
  group=h[[2]]
  
  # delete 20% of nodes (divided through 5)
  del=sample(1:nrow(m),round(nrow(m)/5))
  m=m[-del,-del]
  group=group[-del] # remove the deleted nodes from the group vector
  
  n <- nrow(m) # nr of nodes
  x=rep(1,n) # each node a positive start value of 1 
  
  
  # put everywhere 0 and 1 in the connection matrix, in a random way,
  # but with more chance of having a 1 if the strength at a specific 
  # position was higher specified
  m=ifelse(m>matrix(runif(n*n,0,1),n,n),1,0)
  m[lower.tri(m)] <- t(m)[lower.tri(m)] # symmetric (e.g., from node 1 to 5 and 5 to 1 a connection of 1)
  diag(m)=0 # elements at diagonal 0
  
   group <- as_membership(cluster_louvain(
   graph_from_adjacency_matrix(m, mode = "undirected")))$membership
  
  # create the graph of the network 
  q=qgraph(m,groups=group, layout="spring",labels=F) 
  m_orig=m
  
  # run simulation
  s=rep(0,iterations)
  
  ii=iterations/5
  xi=matrix(NA,floor(iterations/ii),n) # collect x values 
  mi=array(NA,c(floor(iterations/ii),n,n)) # collect x values 
  linksi=rep(n,iterations)
  
  j=0
  for(i in 1:iterations)
  {  x<-glauber_step(x,n,t=t[1],m,beta)
  s[i]=mean(x)
  linksi[i]=sum(m)
  if(i%%(iterations/400)==0)
  {
    links=which(m==1)
    if(length(links>0))
    {
      link=sample(links,1)
      m[link]=0 # TODO: waarom hier en 2 regels verder gelijk stellen aan 0?
      link_index=getMatrixIndices(m,link)
      m[link_index[2],link_index[1]]=0
    } # if(length(links>0))
  } #  if(i%%(iterations/400)==0)
  
  if(i%%ii==0)
  {
    j=j+1
    print(i);print(sum(m))
    xi[j,]=x
    mi[j,,]=m
  } # if(i%%ii==0)
  } # for(i in 1:iterations)
  
  xis[[person]] <- xi
  mis[[person]] <- mi
  linksis[[person]] <- linksi
  qs[[person]] <- q
  m_origs[[person]] <- m_orig
  ss[[person]] <- s
} # persons loop 

pdf('figure 8gr.pdf',w=5,h=6)
# for all persons
par(mar=c(5,5,3,6))
for(i in 1:nr_of_persons){
  if(i==1){
    plot(unlist(ss[[i]]),bty='n',type='l',xlab='time',ylab=expression(bar(x)),
         ylim=c(-1.2,1.1),lwd=1,col='black',cex.axis=1.5,cex.lab=1.5,las=1,
         axes=F,cex.main=1.5,main=name1)
    axis(1,labels=FALSE);axis(2)
    # text(iterations/6,-1.1, name2)
  } else {
    lines(unlist(ss[[i]]))
  }
  
  
} 
dev.off()
#  advanced plot stuff (figure 7A)
capturePlot <- function() {
  grid.echo()
  grid.grab()
}
plot.new()
#stop()
plotGrob <- vector("list", length = j)
for(jj in 1:j)
{
  qgraph(mi[jj,,],color=2+(xi[jj,]+1)/2, layout=q$layout,labels=F,
         vsize=.4,edge.width=.3,bg=adjustcolor("white", alpha.f=0),border.width=.1)
  plotGrob[[jj]] <- capturePlot()
  print(jj)
}

qgraph(m_orig,groups=group, layout=q$layout,labels=F,
       vsize=.6,edge.width=.3, bg=adjustcolor("white", alpha.f=0),border.width=.1)
plotGrob1 <- capturePlot()

png("temp_plot.png",w=800,h=500,res=200)
par(mar=c(5,5,3,3),mgp=c(2, 1, 0))
plot(linksi,type='l',xlab='time',ylab="# connections",bty='o',
     lwd=1.5,col='black',cex.axis=1.5,cex.lab=1.2,
     axes=F,cex.main=1.2,
     main=paste('Decrease in connections'),
     bg=adjustcolor("white", alpha.f=0),border.width=.1)
axis(1,labels=FALSE);axis(2,at=range(linksi))
dev.off()
img <- readPNG("temp_plot.png")
plotGrob2 <- rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))



pdf(paste('figure8.pdf',sep='',collapse=''),h=8,w=12)
par(mar=c(5,5,3,6))
plot(s,bty='n',type='l',xlab='time',ylab=expression(bar(x)),
     ylim=c(-1.2,1.1),lwd=5,col='blue',cex.axis=1.5,cex.lab=1.5,las=1,
     axes=F,cex.main=1.5,main=name1)
axis(1,labels=FALSE);axis(2)
text(iterations/6,-1.1, name2)
# axis(side = 4, at = pretty(range(t)), las = 1)
# mtext("Right Y-Axis", side = 4, line = 3)
# par(new = TRUE);par(mar=c(5,5,3,6))
# plot(t, type = "l", col = "black", ylim = range(t), axes = FALSE, xlab = "", ylab = "")
# legend("topright", legend = c("X", "tau"), col = c("blue", "black"), lty = 1)

for(jj in 1:j)
{
  vp <- viewport(x=.9*(jj*ii/iterations), y=-.1+(.6*s[ii+(jj-1)*ii]+1)/2,
                 width=0.2, height=0.35, just=c("center","bottom"))
  pushViewport(vp)
  grid.draw(plotGrob[[jj]])
  popViewport()
}

vp <- viewport(x=.25, y=.18,
               width=0.25, height=.35, just=c("center","bottom"))
pushViewport(vp)
grid.draw(plotGrob1)
popViewport()

vp <- viewport(x=.83, y=.65,
               width=0.28, height=.28, just=c("center","bottom"))
pushViewport(vp)
grid.draw(plotGrob2)
popViewport()

dev.off()

####################
## figure 9
###################

## same as figure 8 but now with a hierarchical stochastic block

library(igraph)
library(colorspace)
library(qgraph)

getMatrixIndices <- function(matrix, index) {
  nrows <- nrow(matrix)
  row_index <- ((index - 1) %% nrows) + 1
  col_index <- ((index - 1) %/% nrows) + 1
  return(c(row_index,col_index))
}

hamiltonian=function(x,n,t,w) -sum(t*x)-sum(w*x%*%t(x)/2)
glauber_step = function(x,n,t,w,beta)
{
  i = sample(1:n,size=1) # take a random node
  x_new=x;x_new[i]=x_new[i]*-1 # construct new state with flipped node
  p=1/(1+exp(beta*(hamiltonian(x_new,n,t,w)-hamiltonian(x,n,t,w))))  # update probability
  if(runif(1)<p) x=x_new # update state
  return(x)
}

# create (hierarchical) SBM
hsbm <- function (cluster_sizes,strengths=NULL,plot=T)
{
  k=cluster_sizes # cluster sizes
  levels=length(k) # hierarchic levels
  n=prod(k) # nodes
  if(n > 2000) print('Warning: ', n, 'nodes', 'might be too much for qgraph')
  if(length(k)!=length(strengths)) { print ('cluster_sizes and strengths dont match');stop()}
  group=rep(0,n) # node group
  m=matrix(levels,n,n) # connection matrix
  for(l in (levels-1):1)
    for(i in 1:n)
    {for(j in 1:n)
      if((i-1)%/%prod(k[1:l])==(j-1)%/%prod(k[1:l])) m[i,j]=l
    group[i]=(i-1)%/%prod(k[1:l])
    }
  
  if(plot) image(1/(m+.5))
  
  if(length(strengths)<1) m <- 1/m^3 else
  {
    for(i in 1:levels)
      m[m==i]=strengths[i]
  }
  return(list(m,group))
}

iterations=50000

    set.seed(1)
    clusters=c(8,3,3,4);strengths = c(1,.1,.02,.0005)
    t=seq(-2.5,-5,length=iterations)
    beta=2
    name1=''
    name2='Hierarchical stochastic\nblock model (HSBM)'
  
  
  #build network  
  h=hsbm(clusters,strengths)
  m <- h[[1]]
  group=h[[2]]
  
  # delete 20% of nodes
  del=sample(1:nrow(m),round(nrow(m)/5))
  m=m[-del,-del]
  group=group[-del]
  
  n <- nrow(m) # nr of nodes
  x=rep(1,n) # start values positive
  
  m=ifelse(m>matrix(runif(n*n,0,1),n,n),1,0)
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m)=0
  
 group <- as_membership(cluster_louvain(
    graph_from_adjacency_matrix(m, mode = "undirected")))$membership
  
  q=qgraph(m,groups=group, layout="spring",labels=F,vsize=1) 
  m_orig=m
  
  # run simulation
  s=rep(0,iterations)
  
  ii=iterations/5
  xi=matrix(NA,floor(iterations/ii),n) # collect x values 
  mi=array(NA,c(floor(iterations/ii),n,n)) # collect x values 
  linksi=rep(n,iterations)
  
  j=0
  for(i in 1:iterations)
  {  x<-glauber_step(x,n,t=t[1],m,beta)
  s[i]=mean(x)
  linksi[i]=sum(m)
  if(i%%(iterations/400)==0)
  {
    links=which(m==1)
    if(length(links>0))
    {
      link=sample(links,1)
      m[link]=0
      link_index=getMatrixIndices(m,link)
      m[link_index[2],link_index[1]]=0
    }
  }
  if(i%%ii==0)
  {
    j=j+1
    print(i);print(sum(m))
    xi[j,]=x
    mi[j,,]=m
  }
  }
  
  #  advanced plot stuff (figure 8)
  capturePlot <- function() {
    grid.echo()
    grid.grab()
  }
  plot.new()
  #stop()
  plotGrob <- vector("list", length = j)
  for(jj in 1:j)
  {
    qgraph(mi[jj,,],color=2+(xi[jj,]+1)/2, layout=q$layout,labels=F,
           vsize=.4,edge.width=.3,bg=adjustcolor("white", alpha.f=0),border.width=.1)
    plotGrob[[jj]] <- capturePlot()
    print(jj)
  }
  
  qgraph(m_orig,groups=group, layout=q$layout,labels=F,
         vsize=.6,edge.width=.3, bg=adjustcolor("white", alpha.f=0),border.width=.1)
  plotGrob1 <- capturePlot()
  
  png("temp_plot.png",w=800,h=500,res=200)
  par(mar=c(5,5,3,3),mgp=c(2, 1, 0))
  plot(linksi,type='l',xlab='time',ylab="# connections",bty='o',
       lwd=1.5,col='black',cex.axis=1.5,cex.lab=1.2,
       axes=F,cex.main=1.2,
       main=paste('Decrease in connections'),
       bg=adjustcolor("white", alpha.f=0),border.width=.1)
  axis(1,labels=FALSE);axis(2,at=range(linksi))
  dev.off()
  img <- readPNG("temp_plot.png")
  plotGrob2 <- rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))
  
  
  
  pdf(paste('figure9.pdf',sep='',collapse=''),h=8,w=12)
  par(mar=c(5,5,3,6))
  plot(s,bty='n',type='l',xlab='time',ylab=expression(bar(x)),
       ylim=c(-1.2,1.1),lwd=5,col='blue',cex.axis=1.5,cex.lab=1.5,las=1,
       axes=F,cex.main=1.5,main=name1)
  axis(1,labels=FALSE);axis(2)
  text(iterations/6,-1.1, name2)
  # axis(side = 4, at = pretty(range(t)), las = 1)
  # mtext("Right Y-Axis", side = 4, line = 3)
  # par(new = TRUE);par(mar=c(5,5,3,6))
  # plot(t, type = "l", col = "black", ylim = range(t), axes = FALSE, xlab = "", ylab = "")
  # legend("topright", legend = c("X", "tau"), col = c("blue", "black"), lty = 1)
  
  for(jj in 1:j)
  {
    vp <- viewport(x=.9*(jj*ii/iterations), y=-.1+(.6*s[ii+(jj-1)*ii]+1)/2,
                   width=0.2, height=0.35, just=c("center","bottom"))
    pushViewport(vp)
    grid.draw(plotGrob[[jj]])
    popViewport()
  }
  
  vp <- viewport(x=.25, y=.18,
                 width=0.25, height=.35, just=c("center","bottom"))
  pushViewport(vp)
  grid.draw(plotGrob1)
  popViewport()
  
  vp <- viewport(x=.83, y=.65,
                 width=0.28, height=.28, just=c("center","bottom"))
  pushViewport(vp)
  grid.draw(plotGrob2)
  popViewport()
  
  dev.off()
  
  
  
  
  