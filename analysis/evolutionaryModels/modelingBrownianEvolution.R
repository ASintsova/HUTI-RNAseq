#Into to Phylogenetics in R


#Install Dependencies:

#install.packages("ape")
#install.packages("geiger")
#install.packages("phytools")

library(ape)
library(geiger)
library(phytools)


## read tree from string
tt="(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"
vert.tree<-read.tree(text=tt)
plot(vert.tree)
#cladogram
plot(vert.tree,type="cladogram")
#unrooted
plot(unroot(vert.tree),type="unrooted")
#fan
plot(vert.tree,type="fan")


tree<-read.tree(text="(((A,B),(C,D)),E);")
plot(tree,type="cladogram",edge.width=2)

tree$edge
tree$tip.label
tree$Nnode

#from files
write.tree(tree, "example.tre")
cat(readLines("example.tre"))

#from phytools
writeNexus(tree, "example.nex")
cat(readLines("example.nex"), sep="\n")

#Simulating, plotting, extracting clades, and dropping tips
##for repeatability

set.seed(1)
#simulate a birth-death tree(?)
tree <- pbtree(b=1, d=0.2, n=40)
plotTree(tree)
nodelabels()
#extract specific clade

tt62 <- extract.clade(tree, 62)
plotTree(tt62)
#drop 10 random tips
dtips <- sample(tree$tip.label, 10)
dt <- drop.tip(tree, dtips)
plotTree(dt)


## we could also, say, drop all tips that go extinct before the present
## this is a fun way, but not the only way to do this:
plotTree(et<-drop.tip(tree,getExtinct(tree)),cex=0.7)


#Simulating Browninan Motion

t <- 0:100
sig2 <- 0.01
#set of random derivatives
x <- rnorm(n=length(t)-1, sd  = sqrt(sig2))
#cumulative sum
x <- c(0, cumsum(x))
plot(t, x, type ="l", ylim = c(-2, 2))


nsim<-100
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2)),nsim,length(t)-1)
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
plot(t,X[1,],xlab="time",ylab="phenotype",ylim=c(-2,2),type="l")
apply(X[2:nsim,],1,function(x,t) lines(t,x),t=t)


t<- 100
n <- 30 # taxa
b <- (log(n) - log(2))/t
tree <- pbtree(b=b, n=n, t=t, type= "discrete")

plotTree(tree)



## simulate evolution along each edge
X<-lapply(tree$edge.length,function(x) c(0,cumsum(rnorm(n=x,sd=sqrt(sig2)))))
## reorder the edges of the tree for pre-order traversal
cw<-reorder(tree)
## now simulate on the tree
ll<-tree$edge.length+1
for(i in 1:nrow(cw$edge)){
        pp<-which(cw$edge[,2]==cw$edge[i,1])
        if(length(pp)>0) X[[i]]<-X[[i]]+X[[pp]][ll[pp]]
        else X[[i]]<-X[[i]]+X[[1]][1]
}
## get the starting and ending points of each edge for plotting
H<-nodeHeights(tree)
## plot the simulation
plot(H[1,1],X[[1]][1],ylim=range(X),xlim=range(H),xlab="time",ylab="phenotype")
for(i in 1:length(X)) lines(H[i,1]:H[i,2],X[[i]])
## add tip labels if desired
yy<-sapply(1:length(tree$tip.label),function(x,y) which(x==y),y=tree$edge[,2])
yy<-sapply(yy,function(x,y) y[[x]][length(y[[x]])],y=X)
text(x=max(H),y=yy,tree$tip.label)


## simulate Brownian evolution on a tree with fastBM
x<-fastBM(tree,sig2=sig2,internal=TRUE)
## visualize Brownian evolution on a tree
phenogram(tree,x,spread.labels=TRUE,spread.cost=c(1,0))

## linear regression model is y = b0 + b*x + e
b0<-5.3
b<-0.75
x<-rnorm(n=100,sd=1)
y<-b0+b*x+rnorm(n=100,sd=0.2)
plot(x,y)
fit<-lm(y~x)
abline(fit)