install.packages("splines2") 
install.packages("grpreg")
#install.packages("Matrix")
library(splines2)
library(Matrix)
library(grpreg)
library(igraph)
library(np)

#some useful fct
step_fct<-function(x,step)
{
  t<-length(x)
  x_new<-x[1]
  for (i in 1:(t-1))
  {
    x_int<-seq(x[i],x[i+1],by=step)[-1]
    x_new<-c(x_new,x_int)
  }
  x_new
}

#generate B-spline and its integration
X_big<-rep()
X_big_int<-rep()
degree<-3
step<-0.01
exp_index<-round(exp_index, 2)
for (i in 1:cluster)
{
  knots <- as.vector(summary(data_fitted[,i])[c(2,5)])
  bsMat <- bSpline(data_fitted[,i], knots = knots, degree = degree, intercept=F)
  x_int<-step_fct(as.vector(exp_index),step)
  
  d<-data.frame(data_clustered[,i],exp_index)
  colnames(d)<-c("y","x")
  bw <- npregbw(data_clustered[,i]~exp_index,regtype="ll")
  ghat <- npreg(bw,exdat=data.frame(exp_index=x_int),regtype="ll")
  #bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll",bwtype="generalized_nn")
  #bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=x_cov),ezdat=data.frame(t=t),regtype="ll")
  x_hat<-fitted(ghat)
  #bw<-regCVBwSelC(d$x,d$y, deg=1, kernel=EpaK, interval = c(0, 20))
  #bw<-regCVBwSelC(exp_index, data_clustered[,i], deg=1, kernel=EpaK, interval = c(2000, 11000))
  #x_hat <- locCuadSmootherC(d$x, d$y ,x_int, bw=bw, EpaK)$beta0
  
  #x_hat<-exp(cbind(1,log(x_int)) %*% beta[,i])
  basis_int<-bSpline(x_hat, knots = knots, degree = degree, intercept=F)
  base_int<-0
  for (j in 1:(length(exp_index)-1))
  {
    int_row<-apply(basis_int[x_int>=exp_index[j] & x_int<=exp_index[j+1],][-1,],2,function (x) sum(x)*step)
    base_int<-rbind(base_int,int_row)
  }
  base_int<-apply(base_int,2,cumsum)
  #base_int<-ibs(data_fitted[,i], knots = knots, degree = degree, intercept=F)*(max(exp_index)-min(exp_index)) #not correct
  X_big <- cbind(X_big,bsMat)
  X_big_int<-cbind(X_big_int,base_int)
}
dim(X_big)
dim(X_big_int)
num_cov<-degree+length(knots)
X_big_int_exp_intcep<-cbind(1,exp_index,X_big_int)


##################
##use idea in 2014
##################
#by EBIC
exp_index_scale<-(exp_index-min(exp_index))/(max(exp_index)-min(exp_index))
group<-rep(1:50,each=num_cov)
fit<-grpreg(X_big, data_derivative[,50], group, penalty="grLasso",nlambda=100)
beta_select<-select(fit, "EBIC")$beta
length(beta_select[beta_select!=0])
beta_select
fit$lambda

#by cv
cvfit <- cv.grpreg(X_big, data_derivative[,50], penalty="grLasso",group,nfolds=4,lambda.min=0.05)
plot(cvfit)
cvfit$lambda.min
summary(cvfit)
beta_select<-coef(cvfit)
index_nonzero<-which(beta_select!=0)
beta_select[index_nonzero]

X_big_int_exp_intcep<-cbind(1,exp_index,X_big_int)
X_big_int_exp<-cbind(exp_index,X_big_int)
dim(X_big_int_exp)
fitted<-X_big_int_exp[,c(index_nonzero)]%*%beta_select[index_nonzero]
fitted_self<-X_big_int_exp[,1]*beta_select[1]
fitted_1<-X_big_int_exp[,index_nonzero[2:(2+num_cov-1)]]%*%beta_select[index_nonzero[2:(2+num_cov-1)]]
fitted_2<-X_big_int_exp[,index_nonzero[(2+num_cov):(2+num_cov+num_cov-1)]]%*%beta_select[index_nonzero[(2+num_cov):(2+num_cov+num_cov-1)]]
fitted_3<-X_big_int_exp[,index_nonzero[(2+2*num_cov):(2+2*num_cov+num_cov-1)]]%*%beta_select[index_nonzero[(2+2*num_cov):(2+2*num_cov+num_cov-1)]]
fitted_4<-X_big_int_exp[,index_nonzero[(2+3*num_cov):(2+3*num_cov+num_cov-1)]]%*%beta_select[index_nonzero[(2+3*num_cov):(2+3*num_cov+num_cov-1)]]

par(mfrow = c(1, 1))
plot(x=exp_index,y=data_clustered[,50],ylim=c(-1,3))
lines(x=exp_index,y=fitted_self+mean(data_clustered[,50]-fitted),col="grey")
lines(x=exp_index,y=fitted+mean(data_clustered[,50]-fitted),col="yellow")
lines(x=exp_index,y=fitted_1,ylab="y",col="red")
lines(x=exp_index,y=fitted_2,ylab="y",col="blue")
lines(x=exp_index,y=fitted_3,ylab="y",col="green")
lines(x=exp_index,y=fitted_4,ylab="y",col="orange")

#################################################################

##################
##use idea in 2017
##################
#by EBIC
group<-c(0,rep(1:cluster,each=num_cov))
X_big_int_exp<-cbind(exp_index,X_big_int)
#lambda<-seq(0.001,0.1,length=100)
fit<-grpreg(X_big_int_exp, data_clustered[,50], group, penalty="grLasso",max.iter=100000)
beta_select<-select(fit, "BIC")$beta
length(beta_select[beta_select!=0])
beta_select[beta_select!=0]


##########by cv
alpha<-1
gene_whole<-rep()
self_size<-rep()
X_big_int_exp<-cbind(exp_index,X_big_int)
X_big_int_exp_intcep<-cbind(1,exp_index,X_big_int)
#268524
#seed=2367
#set.seed(23556)
#set.seed(235672)
#lambda<-seq(0.001,0.3,length=100)
group<-c(0,rep(1:cluster,each=num_cov))
lambda.min<-0.05
for (j in 1:cluster)
{
  gene_index<-j
  #cvfit <- cv.grpreg(X_big_int_exp, data_clustered[,gene_index], group=group,penalty="grLasso",nfolds=nrow(data_clustered),seed=23556,alpha=alpha)
  #beta_ini<-grpreg(X_big_int_exp, data_clustered[,gene_index], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha)$beta
  #group_id<-c(1,2,rep(3:(2+cluster),each=num_cov))
  #beta_each_g<-unlist(tapply(beta_ini,group_id,function(x) x))
  #weight_ini<-rep()
  #for(i in 1:ncol(data_clustered))
  #{
  #  m<-sum((beta_each_g[(3+(i-1)*num_cov):(3+(i-1)*num_cov+num_cov-1)])^2)
  #  weight_ini<-c(weight_ini,m)
  #}
  #weight<-(1/(sqrt(weight_ini)+0.0001))^0.3
  cvfit <- cv.grpreg(X_big_int_exp, data_clustered[,gene_index], group=group,penalty="grLasso",nfolds=10,seed=23556,alpha=alpha,lambda.min=lambda.min)
  fit<-grpreg(X_big_int_exp, data_clustered[,gene_index], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha)
  beta_select<-fit$beta
  index_nonzero<-which(beta_select!=0)
  beta_select[index_nonzero]
  
  find_index<-index_nonzero[-(1:2)]-2
  #influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
  #m<-which(influ==j)
  #if (length(m)==0)
  #{
  fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
  fitted_self<-X_big_int_exp_intcep[,1:2]%*%beta_select[1:2]
  fitted_self<-t(fitted_self)
  self_size<-rbind(self_size,fitted_self)
  #}
  #else
  #{
  #  fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
  #  fitted_self<-X_big_int_exp_intcep[,c(1:2,(2+(j-1)*num_cov+1):(2+(j-1)*num_cov+num_cov))]%*%beta_select[c(1:2,(2+(j-1)*num_cov+1):(2+(j-1)*num_cov+num_cov))]
  #  fitted_self<-t(fitted_self)
  #  self_size<-rbind(self_size,fitted_self)
  #  find_self<-which(influ==j)
  #  find_index<-find_index[-c((num_cov*(find_self-1)+1):(num_cov*(find_self-1)+num_cov))]
  #}
  par(mfrow = c(1, 1))
  plot(x=exp_index,y=data_clustered[,gene_index],ylim=c(-1,max(data_clustered[,gene_index])+0.2))
  lines(x=exp_index,y=fitted_self,col="black")
  lines(x=exp_index,y=fitted,col="red")
  
  num_gene<-length(find_index)/num_cov
  gene_cluster<-rep()
  fitted_gene_matrix<-rep()
  if (num_gene>0)
  {
    for (i in 1:num_gene)
    {
      gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
      gene_cluster<-c(gene_cluster,gene_relate)
      fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
      fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
      lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col=i+2)
    }
    gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
    gene_whole<-rbind(gene_whole,gene_cluster)
  }
  print(j)
}
#self fit
self_size<-cbind(rep(1:cluster),self_size)
dim(self_size)

#################################
##subject40: failure
#################################

###order phenotype
dim(data_phenotype)
total<-apply(data_gene_exp,2,sum)
data_phenotype_order<-data_phenotype[,order(total)]
head(data_phenotype_order)
i_40_0<-which(data_phenotype_order[1,]==40 & data_phenotype_order[2,]==0)
i_40_1<-which(data_phenotype_order[1,]==40 & data_phenotype_order[2,]==1)
i_40_7<-which(data_phenotype_order[1,]==40 & data_phenotype_order[2,]==7)
i_40_28<-which(data_phenotype_order[1,]==40 & data_phenotype_order[2,]==28)

####################################################
####construct the data we need to plot for subject 40: 0->1  draw graph!!!!!!
####################################################

###try not circle

####################################################
####construct the data we need to plot 0->1
####################################################
#lines(x=exp_index_scale,y=fitted,ylab="y")
par(mfrow = c(1, 3))
par(cex=0.8)
#par(mar = c(0, 0, 0, 0), oma = c(6, 6, 3, 3))
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
dim(gene_whole)
unique(gene_whole[,1])
table(gene_whole[,1])
gene_whole_threhold<-gene_whole
gene_whole_threhold[,i_40_0+2]<-gene_whole[,i_40_1+2]-gene_whole[,i_40_0+2]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,i_40_0+2])>0.005,]
dim(gene_whole_threhold)
for (i in i_40_0)
{
  #geneate link data with 25th tissue
  links_gene<-gene_whole_threhold[,c(1,2,i+2)]
  colnames(links_gene)<-c("from","to","tissue_name")
  links_gene<-as.data.frame(links_gene)
  links_gene$edge_type<-ifelse(gene_whole_threhold[,i+2]>0,1,2)
  rownames(links_gene)<-rep(1:nrow(links_gene))
  #generate nodes with 25th tissue
  self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
  #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
  nodes_gene<-data.frame(self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(53),1,ifelse(nodes_gene$id %in% (setdiff(self_size[,1],self_index)),2,3))
  net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
  class(net)
  net 
  
  # It's easy to access nodes, edges, and their attributes:
  E(net)
  V(net)$self_weight
  E(net)$tissue_name
  V(net) #total 34 clusters included in the network
  
  # Now we should be able to do this:
  #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
  
  #adjust size of nodes and width of edges
  colrs_nodes <- c("red","gray", "orange")
  V(net)$color <- colrs_nodes[V(net)$node_type]
  colrs_edge<-c("tomato","black")
  E(net)$color <- colrs_edge[E(net)$edge_type]
  V(net)$size <- (V(net)$self_weight)^0.4*3.5
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])*30
  
  #change coordinate to separate the nodes
  #set.seed(4343)
  net.bg<-net
  #l <- layout.fruchterman.reingold(net.bg)
  #l<-layout_with_kk(net.bg)
  #l<-layout_as_tree(net.bg)
  #l<-layout_nicely(net.bg)
  #l<-layout_in_circle(net.bg)
  #l<-layout_as_star(net.bg)
  #l<-layout_as_bipartite(net.bg)
  #l<-layout_on_grid(net.bg)
  #l<-layout_on_sphere(net.bg)
  #l<-layout_with_dh(net.bg)
  
  # Normalize them so that they are in the -1, 1 interval:
  #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #V(net) #total 34 clusters included in the network
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("lightcyan"),mark.border=NA)
}


####################################################
####construct the data we need to plot for subject 40: 0->7
####################################################

###try not circle

####################################################
####construct the data we need to plot 0->7
####################################################
#lines(x=exp_index_scale,y=fitted,ylab="y")

dim(gene_whole)
unique(gene_whole[,1])
table(gene_whole[,1])
gene_whole_threhold<-gene_whole
gene_whole_threhold[,i_40_0+2]<-gene_whole[,i_40_7+2]-gene_whole[,i_40_0+2]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,i_40_0+2])>0.005,]
dim(gene_whole_threhold)
for (i in i_40_0)
{
  #geneate link data with 25th tissue
  links_gene<-gene_whole_threhold[,c(1,2,i+2)]
  colnames(links_gene)<-c("from","to","tissue_name")
  links_gene<-as.data.frame(links_gene)
  links_gene$edge_type<-ifelse(gene_whole_threhold[,i+2]>0,1,2)
  rownames(links_gene)<-rep(1:nrow(links_gene))
  #generate nodes with 25th tissue
  self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
  #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
  nodes_gene<-data.frame(self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(53),1,ifelse(nodes_gene$id %in% (setdiff(self_size[,1],self_index)),2,3))
  net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
  class(net)
  net 
  
  # It's easy to access nodes, edges, and their attributes:
  E(net)
  V(net)$self_weight
  E(net)$tissue_name
  V(net) #total 34 clusters included in the network
  
  # Now we should be able to do this:
  #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
  
  #adjust size of nodes and width of edges
  colrs_nodes <- c("red","gray", "orange")
  V(net)$color <- colrs_nodes[V(net)$node_type]
  colrs_edge<-c("tomato","black")
  E(net)$color <- colrs_edge[E(net)$edge_type]
  V(net)$size <- (V(net)$self_weight)^0.4*3.5
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])*30
  
  #change coordinate to separate the nodes
  #set.seed(4343)
  net.bg<-net
  #l <- layout.fruchterman.reingold(net.bg)
  #l<-layout_with_kk(net.bg)
  #l<-layout_as_tree(net.bg)
  #l<-layout_nicely(net.bg)
  #l<-layout_in_circle(net.bg)
  #l<-layout_as_star(net.bg)
  #l<-layout_as_bipartite(net.bg)
  #l<-layout_on_grid(net.bg)
  #l<-layout_on_sphere(net.bg)
  #l<-layout_with_dh(net.bg)
  
  # Normalize them so that they are in the -1, 1 interval:
  #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #V(net) #total 34 clusters included in the network
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("moccasin"),mark.border=NA)
}



###try not circle

####################################################
####construct the data we need to plot 0->28
####################################################
#lines(x=exp_index_scale,y=fitted,ylab="y")
dim(gene_whole)
unique(gene_whole[,1])
table(gene_whole[,1])
gene_whole_threhold<-gene_whole
gene_whole_threhold[,i_40_0+2]<-gene_whole[,i_40_28+2]-gene_whole[,i_40_0+2]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,i_40_0+2])>0.005,]
dim(gene_whole_threhold)
for (i in i_40_0)
{
  #geneate link data with 25th tissue
  links_gene<-gene_whole_threhold[,c(1,2,i+2)]
  colnames(links_gene)<-c("from","to","tissue_name")
  links_gene<-as.data.frame(links_gene)
  links_gene$edge_type<-ifelse(gene_whole_threhold[,i+2]>0,1,2)
  rownames(links_gene)<-rep(1:nrow(links_gene))
  #generate nodes with 25th tissue
  self_index<-unique(c(gene_whole_threhold[,1],gene_whole_threhold[,2]))
  #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
  nodes_gene<-data.frame(self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(53,59,98,109),1,ifelse(nodes_gene$id %in% (setdiff(self_size[,1],self_index)),2,3))
  net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
  class(net)
  net 
  
  # It's easy to access nodes, edges, and their attributes:
  E(net)
  V(net)$self_weight
  E(net)$tissue_name
  V(net) #total 34 clusters included in the network
  
  # Now we should be able to do this:
  #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
  
  #adjust size of nodes and width of edges
  colrs_nodes <- c("red","gray", "orange")
  V(net)$color <- colrs_nodes[V(net)$node_type]
  colrs_edge<-c("tomato","black")
  E(net)$color <- colrs_edge[E(net)$edge_type]
  V(net)$size <- (V(net)$self_weight)^0.4*3.5
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])*30
  
  #change coordinate to separate the nodes
  set.seed(4343)
  net.bg<-net
  l <- layout.fruchterman.reingold(net.bg)
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #l<-layout_with_kk(net.bg)
  #l<-layout_as_tree(net.bg)
  #l<-layout_nicely(net.bg)
  #l<-layout_in_circle(net.bg)
  #l<-layout_as_star(net.bg)
  #l<-layout_as_bipartite(net.bg)
  #l<-layout_on_grid(net.bg)
  #l<-layout_on_sphere(net.bg)
  #l<-layout_with_dh(net.bg)
  # Normalize them so that they are in the -1, 1 interval:
  l[1,1]<--0.6
  l[109,]<-c(-0.45,-0.06)
  l[98,]<-c(-0.15,-0.1)
  l[26,]<-c(0.15,-0.2)
  l[122,2]<--0.19
  l[112,1]<--0.38
  l[25,2]<--0.35
  l[104,2]<--0.25
  l[135,]<-c(-0.23,-0.38)
  l[120,]<-c(-0.1,0)
  l[49,]<-c(-0.3,0.7)
  l[108,]<-c(-0.35,0.83)
  l[41,2]<-l[123,2]<-0.5
  l[90,2]<-0.3
  l[33,]<-c(0,0.6)
  l[99,]<-c(0.3,0.82)
  l[138,1]<-0.6
  l[45,1]<-0.55
  l[38,1]<-0.6
  l[50,1]<-0.62
  l[37,2]<-0.73
  l[77,2]<--0.9
  l[63,1]<-l[140,1]<-l[52,1]<--0.85
  l[117,2]<-l[124,2]<-0.8
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #V(net) #total 34 clusters included in the network
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("darkolivegreen1"),mark.border=NA)
}




###order phenotype
total<-apply(data_gene_exp,2,sum)
data_phenotype_order<-data_phenotype[,order(total)]
head(data_phenotype_order)
high<-data_phenotype_order[2,]==1
low<-data_phenotype_order[2,]==0

##############################
#plot four fitted values
##############################
gene_whole<-rep()
self_size<-rep()
group<-c(0,rep(1:cluster,each=num_cov))
cluster_select<-c(3,20,45,48)
par(mfrow=c(2,2))
par(cex=1.5)
par(mar = c(0, 0, 0.9,0), oma = c(6, 6, 3, 3),mgp=c(1,1,0))
#dev.set(dev.prev())
for (j in cluster_select)
{
  gene_index<-j
  #cvfit <- cv.grpreg(X_big_int_exp, data_clustered[,gene_index], group=group,penalty="grLasso",nfolds=nrow(data_clustered),seed=23556,alpha=alpha)
  #beta_ini<-grpreg(X_big_int_exp, data_clustered[,gene_index], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha)$beta
  #group_id<-c(1,2,rep(3:(2+cluster),each=num_cov))
  #beta_each_g<-unlist(tapply(beta_ini,group_id,function(x) x))
  #weight_ini<-rep()
  #for(i in 1:cluster)
  #{
  #  m<-mean((X_big_int_exp_intcep[,(3+(i-1)*num_cov):(3+(i-1)*num_cov+num_cov-1)]%*%beta_each_g[(3+(i-1)*num_cov):(3+(i-1)*num_cov+num_cov-1)])^2)
  #  weight_ini<-c(weight_ini,m)
  #}
  #weight<-(1/(sqrt(weight_ini)+0.0001))^0.3
  cvfit <- cv.grpreg(X_big_int_exp, data_clustered[,gene_index], group=group,penalty="grLasso",nfolds=10,lambda.min=lambda.min,seed=23556,alpha=alpha)
  fit<-grpreg(X_big_int_exp, data_clustered[,gene_index], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha)
  beta_select<-fit$beta
  index_nonzero<-which(beta_select!=0)
  beta_select[index_nonzero]
  
  fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
  fitted_self<-X_big_int_exp_intcep[,1:2]%*%beta_select[1:2]
  fitted_self<-t(fitted_self)
  self_size<-rbind(self_size,fitted_self)
  
  #par(mfrow = c(1, 2))
  ##########j=1
  if (j %in% (cluster_select[1]))
  {
    plot(x=exp_index[high],y=data_clustered[high,gene_index],axes=FALSE,ylim=c(-1.4,max(data_clustered[,gene_index])+0.2),pch=19,col="navyblue",xlim=c(min(exp_index),max(exp_index)+100))
    points(x=exp_index[low],y=data_clustered[low,gene_index],pch=19,cex=1.2,col="mediumorchid1")
    lines(x=exp_index,y=fitted_self,col="blue",lwd=2.5)
    lines(x=exp_index,y=fitted,col="darkorange",lwd=2.5)
    text(x=rep(max(exp_index)+3, 1), y=fitted[ncol(data_original)], pos=4, labels=j)
    
    find_index<-index_nonzero[-(1:2)]-2
    influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
    num_gene<-length(find_index)/num_cov
    gene_cluster<-rep()
    fitted_gene_matrix<-rep()
    if (num_gene>0)
    {
      for (i in 1:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
      }
      influ_order<-influ[order(fitted_gene_matrix[ncol(data_original),])]
      text(x=rep(max(exp_index)+3, length(influ_order)), y=seq(-1,2,length=length(influ_order)), pos=4, labels=influ_order)
      axis(2,col = "grey40", col.axis = "grey20", at = seq(-1,7, 1),cex.axis=1.2)
      box(col = "grey60",lwd=2)
      gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
      gene_whole<-rbind(gene_whole,gene_cluster)
    }
    
    par(new = TRUE)
    par("plt" = c(0.49, 0.89, 0.28, 0.57),mgp=c(1,0.2,0))
    
    if (num_gene>0)
    {
      gene_relate<-find_index[(1-1)*num_cov+1]%/%num_cov+1
      gene_cluster<-c(gene_cluster,gene_relate)
      fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]%*%beta_select[index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]
      fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
      plot(x=exp_index,y=fitted_gene_matrix[,1],col="green",ylab='',xlab='',axes=FALSE,lwd=1.5,type="l",ylim=c(-0.1,0.3))
      
      for (i in 2:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
      }
      axis(2,col = "grey40", col.axis = "grey20", at = c(0,0.3),cex.axis=0.7)
      #axis(1,col = "grey40", col.axis = "grey20", at = exp_index, labels=round(exp_index/10000,digits=1),cex.axis=0.7)
      box(col = "grey60",lwd=2)
    }
  }
  
  ######j=2
  par(mar = c(0, 0, 0.9,0),mgp=c(1,1,0))
  if (j %in% (cluster_select[2]))
  {
    plot(x=exp_index[high],y=data_clustered[high,gene_index],axes=FALSE,ylim=c(-1,max(data_clustered[,gene_index])+0.2),pch=19,col="navyblue",xlim=c(min(exp_index),max(exp_index)+100))
    points(x=exp_index[low],y=data_clustered[low,gene_index],pch=19,cex=1.2,col="mediumorchid1")
    lines(x=exp_index,y=fitted_self,col="blue",lwd=2.5)
    lines(x=exp_index,y=fitted,col="darkorange",lwd=2.5)
    text(x=rep(max(exp_index)+3, 1), y=fitted[ncol(data_original)], pos=4, labels=j)
    
    
    find_index<-index_nonzero[-(1:2)]-2
    influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
    num_gene<-length(find_index)/num_cov
    gene_cluster<-rep()
    fitted_gene_matrix<-rep()
    if (num_gene>0)
    {
      for (i in 1:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
      }
      influ_order<-influ[order(fitted_gene_matrix[ncol(data_original),])]
      text(x=rep(max(exp_index)+3, length(influ_order)), y=seq(-1,2,length=length(influ_order)), pos=4, labels=influ_order)
      box(col = "grey60",lwd=2)
      gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
      gene_whole<-rbind(gene_whole,gene_cluster)
    }
    
    par(new = TRUE)
    par("plt" = c(0.49, 0.89, 0.28, 0.57),mgp=c(1,0.2,0))
    if (num_gene>0)
    {
      gene_relate<-find_index[(1-1)*num_cov+1]%/%num_cov+1
      gene_cluster<-c(gene_cluster,gene_relate)
      fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]%*%beta_select[index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]
      fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
      plot(x=exp_index,y=fitted_gene_matrix[,1],ylab='',xlab='',axes=FALSE,col="green",lwd=1.5,ylim=c(-0.1,0.05),type="l")
      for (i in 2:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
      }
    }
    axis(2,col = "grey40", col.axis = "grey20", at = c(-0.1,0),cex.axis=0.7)
    #axis(1,col = "grey40", col.axis = "grey20", at = exp_index, labels=round(exp_index/10000,digits=1),cex.axis=0.7)
    box(col = "grey60",lwd=2)
  }
  
  
  #######j=3
  par(mar = c(0, 0, 0.9,0),mgp=c(1,1,0))
  
  if (j %in% (cluster_select[3]))
  {
    plot(x=exp_index[high],y=data_clustered[high,gene_index],axes=FALSE,ylim=c(-1,max(data_clustered[,gene_index])+0.6),pch=19,col="navyblue",xlim=c(min(exp_index),max(exp_index)+100))
    points(x=exp_index[low],y=data_clustered[low,gene_index],pch=19,cex=1.2,col="mediumorchid1")
    lines(x=exp_index,y=fitted_self,col="blue",lwd=2.5)
    lines(x=exp_index,y=fitted,col="darkorange",lwd=2.5)
    text(x=rep(max(exp_index)+3, 1), y=fitted[ncol(data_original)], pos=4, labels=j)
    
    
    find_index<-index_nonzero[-(1:2)]-2
    influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
    num_gene<-length(find_index)/num_cov
    gene_cluster<-rep()
    fitted_gene_matrix<-rep()
    if (num_gene>0)
    {
      for (i in 1:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
        influ_order<-influ[order(fitted_gene_matrix[ncol(data_original),])]
      }
      influ_order<-influ[order(fitted_gene_matrix[ncol(data_original),])]
      text(x=rep(max(exp_index)+3, length(influ_order)), y=seq(-0.5,0.5,length=length(influ_order)), pos=4, labels=influ_order)
      axis(2,col = "grey40", col.axis = "grey20", at = seq(-1,5, 1),cex.axis=1.2)
      axis(1, col = "grey40", col.axis = "grey20", at = exp_index, labels=round(exp_index/10000,digits=2),cex.axis=1.2)
      box(col = "grey60",lwd=2)
      gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
      gene_whole<-rbind(gene_whole,gene_cluster)
    }
    
    par(new = TRUE)
    par("plt" = c(0.49, 0.89, 0.28, 0.57),mgp=c(1,0.2,0))
    if (num_gene>0)
    {
      gene_relate<-find_index[(1-1)*num_cov+1]%/%num_cov+1
      gene_cluster<-c(gene_cluster,gene_relate)
      fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]%*%beta_select[index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]
      fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
      plot(x=exp_index,y=fitted_gene_matrix[,1],ylab='',xlab='',col="green",lwd=1.5,axes=FALSE,ylim=c(-0.05,0.3),type="l")
      for (i in 2:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
        influ_order<-influ[order(fitted_gene_matrix[ncol(data_original),])]
      }
      axis(2,col = "grey40", col.axis = "grey20", at = c(0,0.3),cex.axis=0.7)
      #axis(1,col = "grey40", col.axis = "grey20", at = exp_index, labels=round(exp_index/10000,digits=1),cex.axis=0.7)
      box(col = "grey60",lwd=2)
    }
  } 
  
  ##########j=4
  par(mar = c(0, 0, 0.9,0),mgp=c(1,1,0))
  
  if (j %in% (cluster_select[4]))
  {
    plot(x=exp_index[high],y=data_clustered[high,gene_index],axes=FALSE,ylim=c(-1,max(data_clustered[,gene_index])+0.2),pch=19,col="navyblue",xlim=c(min(exp_index),max(exp_index)+100))
    points(x=exp_index[low],y=data_clustered[low,gene_index],pch=19,cex=1.2,col="mediumorchid1")
    lines(x=exp_index,y=fitted_self,col="blue",lwd=2.5)
    lines(x=exp_index,y=fitted,col="darkorange",lwd=2.5)
    text(x=rep(max(exp_index)+3, 1), y=fitted[ncol(data_original)], pos=4, labels=j)
    
    find_index<-index_nonzero[-(1:2)]-2
    influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
    num_gene<-length(find_index)/num_cov
    gene_cluster<-rep()
    fitted_gene_matrix<-rep()
    if (num_gene>0)
    {
      for (i in 1:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
      }
      influ_order<-influ[order(fitted_gene_matrix[ncol(data_original),])]
      text(x=rep(max(exp_index)+3, length(influ_order)), y=seq(-1,1,length=length(influ_order)), pos=4, labels=influ_order)
      axis(1, col = "grey40", col.axis = "grey20", at = exp_index,labels=round(exp_index/10000,digits=2),cex.axis=1.2)
      box(col = "grey60",lwd=2)
      gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
      gene_whole<-rbind(gene_whole,gene_cluster)
    }
    
    par(new = TRUE)
    par("plt" = c(0.49, 0.89, 0.28, 0.57),mgp=c(1,0.2,0))
    if (num_gene>0)
    {
      gene_relate<-find_index[(1-1)*num_cov+1]%/%num_cov+1
      gene_cluster<-c(gene_cluster,gene_relate)
      fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]%*%beta_select[index_nonzero[(3+(1-1)*num_cov):(3+1*num_cov-1)]]
      fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
      plot(x=exp_index,y=fitted_gene_matrix[,1],ylab='',xlab='', axes=FALSE, col="green",lwd=1.5,ylim=c(-0.05,0.2),type="l")
      
      for (i in 2:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster<-c(gene_cluster,gene_relate)
        fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
        fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
        lines(x=exp_index,y=fitted_gene_matrix[,i],ylab="y",col="green",lwd=1.5)
      }
      axis(2,col = "grey40", col.axis = "grey20", at = c(0,0.2),cex.axis=0.7)
      #axis(1,col = "grey40", col.axis = "grey20", at = exp_index, labels=round(exp_index/10000,digits=1),cex.axis=0.7)
      box(col = "grey60",lwd=2)
    }
  }
  print(j)
}
mtext("Expression Index (1e+04)", side = 1, outer = TRUE, cex = 2, line = 2.2,
      col = "grey20",font=1)
mtext("Expression of Individual Modules", side = 2, outer = TRUE, cex = 2, line = 2.2,
      col = "grey20",font=1)
