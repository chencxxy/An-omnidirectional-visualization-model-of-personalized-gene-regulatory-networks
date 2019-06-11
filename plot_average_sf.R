#######################################################
#####continue of network construction
#######################################################

#################################
##average effect of s & f
#################################

###order phenotype
dim(data_phenotype)
total<-apply(data_gene_exp,2,sum)
data_phenotype_order<-data_phenotype[,order(total)]
head(data_phenotype_order)
subject_id_s<-data_phenotype_order[1,data_phenotype_order[3,]==1]
subject_id_f<-data_phenotype_order[1,data_phenotype_order[3,]==0]
value_s<-rep()
value_f<-rep()
for (time_status in c(0,1,7,28))
{
  affect_value<-rep()
  for (subject_id in subject_id_s)
  {
    i_subject_time_status<-which(data_phenotype_order[1,]==subject_id & data_phenotype_order[2,]==time_status)
    affect_value<-cbind(affect_value,gene_whole[,i_subject_time_status+2])
  }
  value_s<-cbind(value_s,apply(affect_value,1,mean))
  
  affect_value<-rep()
  for (subject_id in subject_id_f)
  {
    i_subject_time_status<-which(data_phenotype_order[1,]==subject_id & data_phenotype_order[2,]==time_status)
    affect_value<-cbind(affect_value,gene_whole[,i_subject_time_status+2])
  }
  value_f<-cbind(value_f,apply(affect_value,1,mean))
}
dim(gene_whole)
dim(value_s)
dim(value_f)
value_diff<-value_s-value_f
dim(value_diff)

dim(self_size)
ave_self_0<-apply(self_size[,data_phenotype_order[2,]==0],1,mean)
ave_self_1<-apply(self_size[,data_phenotype_order[2,]==1],1,mean)
ave_self_7<-apply(self_size[,data_phenotype_order[2,]==7],1,mean)
ave_self_28<-apply(self_size[,data_phenotype_order[2,]==28],1,mean)
averages_self_size<-self_size

####################################################
####construct the data we need to plot for  f->s at time=0  
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
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.003,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_0
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
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(124,53,59),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.7*25
  
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
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("palegoldenrod"),mark.border=NA)
}

###plot time=1
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,2]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.003,]
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
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(109,59,5,53,124),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.7*25
  
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
  l[5,2]<--0.3
  l[4,1]<--0.6
  l[109,]<-c(-0.47,-0.58)
  l[2,]<-c(-0,-0.62)
  l[9,1]<-0.3
  l[60,2]<--0.9
  l[120,1]<-0.34
  l[52,1]<-0.4
  l[143,1]<-0.55
  l[78,2]<- l[91,2]<-0.65
  l[83,2]<-0.45
  l[16,2]<-0.88
  l[117,1]<-0.3
  l[123,1]<-0.4
  l[77,1]<-0.55
  l[124,1]<-0.25
  l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  #V(net) #total 34 clusters included in the network
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("lightcyan"),mark.border=NA)
}

#plot time=7
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,3]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.003,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_7
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
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(59,98),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.7*25
  
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
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("moccasin"),mark.border=NA)
}

#plot time=28
gene_whole_threhold<-gene_whole
gene_whole_threhold[,1+2]<-value_diff[,4]
gene_whole_threhold<-gene_whole_threhold[abs(gene_whole_threhold[,1+2])>0.003,]
dim(gene_whole_threhold)
averages_self_size[,2]<-ave_self_28
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
  nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(53,109,98,5),1,ifelse(nodes_gene$id %in% (setdiff(averages_self_size[,1],self_index)),2,3))
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
  E(net)$width <- abs(gene_whole_threhold[,i+2])^0.7*25
  
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
  plot(net.bg, rescale=F, layout=l,edge.arrow.size=.07,mark.groups=list(c(1:145)),mark.col=c("darkolivegreen1"),mark.border=NA)
}

