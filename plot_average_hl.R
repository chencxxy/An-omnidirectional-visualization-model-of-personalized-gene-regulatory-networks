#######################################################
#####continue of network construction
#######################################################

#################################
##average effect of h & l
#################################

###order phenotype
dim(data_phenotype)
total<-apply(data_gene_exp,2,sum)
data_phenotype_order<-data_phenotype[,order(total)]
head(data_phenotype_order)
subject_id_h<-data_phenotype_order[2,]==1
subject_id_l<-data_phenotype_order[2,]==0
value_h<-rep()
value_l<-rep()
for (time_status in c(0.08,1,30,180))
{
    i_subject_time_status<-which(subject_id_h & data_phenotype_order[1,]==time_status)
    affect_value<-gene_whole[,i_subject_time_status+2]
  value_h<-cbind(value_h,apply(affect_value,1,mean))
  
    i_subject_time_status<-which(subject_id_l & data_phenotype_order[1,]==time_status)
    affect_value<-gene_whole[,i_subject_time_status+2]
  value_l<-cbind(value_l,apply(affect_value,1,mean))
}
dim(gene_whole)
dim(value_h)
dim(value_l)
value_diff<-value_h-value_l
dim(value_diff)

dim(self_size)
ave_self_0.08<-apply(self_size[,data_phenotype_order[1,]==0.08],1,mean)
ave_self_1<-apply(self_size[,data_phenotype_order[1,]==1],1,mean)
ave_self_30<-apply(self_size[,data_phenotype_order[1,]==30],1,mean)
ave_self_180<-apply(self_size[,data_phenotype_order[1,]==180],1,mean)
averages_self_size<-self_size

####################################################
####construct the data we need to plot for  l->h at time=0.08  
####################################################

###try not circle
#lines(x=exp_index_scale,y=fitted,ylab="y")
par(mfrow = c(2, 2))
par(cex=0.8)
#par(mar = c(0, 0, 0, 0), oma = c(6, 6, 3, 3))
par(mar = c(0.15, 0.15, 0.15, 0.15),oma = c(1, 1, 1, 1))
#plot time=0
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,1]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.001,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_0.08
for (i in 1)
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
  nodes_gene<-data.frame(averages_self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(47,49,44),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  V(net)$size <- (V(net)$self_weight)^0.4*4
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.6*60
  
  #change coordinate to separate the nodes
  set.seed(4343)
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
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.15,mark.groups=list(c(1:50)),mark.col=c("palegoldenrod"),mark.border=NA)
}

###plot l->h at time=1
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,2]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.001,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_1
for (i in 1)
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
  nodes_gene<-data.frame(averages_self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(47,46,38,34,49),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  V(net)$size <- (V(net)$self_weight)^0.4*4
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.6*60
  
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
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.15,mark.groups=list(c(1:50)),mark.col=c("lightcyan"),mark.border=NA)
}

#plot l->h at time=30
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,3]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.001,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_30
for (i in 1)
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
  nodes_gene<-data.frame(averages_self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(1,38,46,44,34,49,13),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  V(net)$size <- (V(net)$self_weight)^0.4*4
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.6*60
  
  #change coordinate to separate the nodes
  set.seed(4343)
  net.bg<-net
  l <- layout.fruchterman.reingold(net.bg)
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
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  l[15,2]<-0.8
  l[11,2]<-0.83
  l[41,2]<-0.85
  l[40,2]<-0.6
  l[14,2]<-0.5
  l[30,2]<-0.75
  l[2,1]<-0.5
  l[35,1]<-0.7
  l[21,1]<-0.72
  l[19,1]<-0.7
  l[27,1]<-0.5
  l[31,1]<-0.74
  l[3,2]<-0.45
  l[1,]<-c(0.15,-0.2)
  l[33,]<-c(-0.7,-0.15)

  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #V(net) #total 34 clusters included in the network
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.28,mark.groups=list(c(1:50)),mark.col=c("moccasin"),mark.border=NA)
}

#plot l->h at time=180
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,4]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.001,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_30
for (i in 1)
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
  nodes_gene<-data.frame(averages_self_size[,c(1,(i+1))])
  colnames(nodes_gene)<-c("id","self_weight")
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(44,46,49,7),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  V(net)$size <- (V(net)$self_weight)^0.4*4
  #initial_weight<-abs(gene_whole[,c(i_40_1+2)]-E(net)$tissue_name)
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.6*60
  
  #change coordinate to separate the nodes
  set.seed(4343)
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
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.15,mark.groups=list(c(1:50)),mark.col=c("darkolivegreen1"),mark.border=NA)
}

