# Ali Tafti
# Optional Advanced Lab 2: SAP Community Network
#Read in the hs0 data over the internet using the read.table() function.
getwd()
# Save the data file to a location on your hard drive and specify the path here (Windows systems use forward slashes)
#dir_path <-"~/YourWorkingDirectoryFilePath"
#setwd(dir_path)

# clear everything out of memory
rm(list=ls()) 

# This is a 10% random sample for class exercises
infile_sub<-"C:\\R\\SAPFull_SubGraph_EdgeList.csv"

## Load package
library(igraph)
library(formattable)
el=read.csv(infile_sub, header = TRUE, sep = ",")
class(el)
# ---
# [1] "data.frame"
# ---
# Describe the data frame
str(el)

# Create the directed graph object
g_SAPSub=graph.data.frame(el, directed = TRUE, vertices= NULL)


# Edges
ecount(g_SAPSub)
## Vertices
vcount(g_SAPSub)


## Check whether Self_loops exist, as do multiple edges
is.simple(g_SAPSub)
#Is it a simple graph? No!
# ---
#[1] FALSE
# ---

# Create edge weights
E(g_SAPSub)$weight <-1
E(g_SAPSub)$weight 
g_SAPSub_simpl<-simplify(g_SAPSub, edge.attr.comb="sum")
is.simple(g_SAPSub_simpl)

# Edges
ecount(g_SAPSub_simpl)
## Vertices
vcount(g_SAPSub_simpl)

# Use the inverse of log weight for some of the network measure calculations
inv_weight<-1/log(E(g_SAPSub_simpl)$weight  + 1)
num_weight<-E(g_SAPSub_simpl)$weight 
length(inv_weight)
E(g_SAPSub_simpl)$weight <-inv_weight

# Plotting the networks-------------------------->>>>>>>>>>>>>>>>>>>>>

par(mfrow=c(1, 1))
layout=layout_with_fr(g_SAPSub_simpl,weights =inv_weight,dim=3,niter=150 )

plot(g_SAPSub_simpl,weight=inv_weight,edge.size=5,edge.width=0.7,edge.arrow.width=0.5,
     layout=layout,vertex.size=5,vertex.color='lightblue',edge.color='black',vertex.label=NA)
title("Simplified SAP Network using Fruchterman-Reingold")

#----------------------------------------------------------------------------------
# You can see the neighbors of some selected nodes
neighbors(g_SAPSub_simpl, v=c('900'))
neighbors(g_SAPSub_simpl, v=c('592540'))

## For example, there is only a path to 814 from 511; not from it. So there are outpaths from 511 to 814, but no in paths. (or the other vectors)
E(g_SAPSub_simpl)$weight <- inv_weight

reciprocity(g_SAPSub_simpl)
is.connected(g_SAPSub_simpl)
is.connected(g_SAPSub_simpl, mode="strong")
is.connected(g_SAPSub_simpl, mode="weak")

##------------------------------------------------------------------->>>>>>>>>>>>>>>>>>>>>
# Diameter with both kinds of weights

# Avg. path length and diameter
average.path.length(g_SAPSub_simpl, directed=TRUE)
diameter(g_SAPSub_simpl)
diameter(g_SAPSub_simpl, weights= num_weight)
diameter(g_SAPSub_simpl, weights= inv_weight)


# Summarize the graph structure
summary(g_SAPSub_simpl)

# Clique structure: 5 cliques of size 5, 39 cliques of size 4, 335 triangles
table(sapply(maximal.cliques(g_SAPSub_simpl), length))
#A <- get.adjacency(g_SAPSub_simpl, sparse=FALSE)

# Can try either of these weighting schemes for various measures; they change the interpretation of the measures
# Inverse weight
E(g_SAPSub_simpl)$weight <- inv_weight # Using inverse weighing schems
# Regular weight
#E(g_SAPSub_simpl)$weight <- num_weight

# Embeddedness/ inverse of structural hole access (see Burt 2004)
constraints_SAP <- round(constraint(g_SAPSub_simpl, nodes=V(g_SAPSub_simpl)), digits=4)
# Degree centrality
  degree_sap <- degree(g_SAPSub_simpl)
# Node betweenness
betweens_SAP <- round(betweenness(g_SAPSub_simpl, v=V(g_SAPSub_simpl), directed = TRUE, nobigint =TRUE, normalized = FALSE))
# Edge betwenness
edgebetweens_SAP<-edge.betweenness(g_SAPSub_simpl, e=E(g_SAPSub_simpl), directed = TRUE)
# Local clustering coefficients
clustering_SAP <- transitivity(g_SAPSub_simpl, type="local", vids=V(g_SAPSub_simpl)) 

# Plots 1 and 2: Can run them together
par(mfrow=c(1, 2))
edge_frame<-data.frame(edgebetweens_SAP, num_weight, inv_weight)
a_edge<-aggregate(edgebetweens_SAP ~ inv_weight, data=edge_frame, mean)
plot(a_edge, col="blue", log="xy", xlab="Weight (inverse) of edge", ylab="Average Betweenness of edges")
node_frame<-data.frame(betweens_SAP, constraints_SAP, clustering_SAP, degree_sap)
a_node<-aggregate(betweens_SAP ~ clustering_SAP, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Clustering", ylab="Average Betweenness of nodes")


# Plot set 2: Four plots 
par(mfrow=c(2, 2))
a_node<-aggregate(betweens_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Betweenness")
a_edge<-aggregate(edgebetweens_SAP ~ num_weight, data=edge_frame, mean)
plot(a_edge, col="blue", log="xy", xlab="Weight of edge", ylab="Average Betweenness of edges")
a_node<-aggregate(clustering_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Clustering")
a_node<-aggregate(constraints_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Constraint (Embeddedness)")

# Log-log degree distributino
par(mfrow=c(1, 2))
d.net <-degree(g_SAPSub_simpl)
dd.net <- degree.distribution(g_SAPSub_simpl)
d <- 1:max(d.net)-1
ind <- (dd.net != 0)
plot(d[ind], dd.net[ind], log="xy", col="blue",
     xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")

# CHUNK 8# Average neighbor degree versus vertex degree
a.nn.deg <- graph.knn(g_SAPSub_simpl,V(g_SAPSub_simpl))$knn
plot(d.net, a.nn.deg, log="xy", 
     col="goldenrod", xlab=c("Log Vertex Degree"),
     ylab=c("Log Average Neighbor Degree"))

#--------------Additional Analysis------------------------------------


Metric_Val<-c(round(reciprocity(g_SAPSub_simpl),digits=6),round(transitivity(g_SAPSub_simpl, weights = inv_weight),digits=6),is.connected(g_SAPSub_simpl),is.connected(g_SAPSub_simpl, mode="strong"),is.connected(g_SAPSub_simpl, mode="weak"))
Metric_Name<-c('Reciprocity','Transitivity','IsConnected','IsConnected-Strong','IsConnected-Weak')

df<-data.frame(Metric_Name,Metric_Val)

formattable(df,align = c("l", rep("r", NCOL(df) - 1)),
            list(`Metric_Name` = formatter("span", style = ~ style(color = "black")),
                 `Metric_Val`= formatter("span", 
                                         x ~ icontext(ifelse(x > 0, x, "remove"), ifelse(x > 0, x, "No")), 
                                         style = x ~ style(color = ifelse(x == 0, "red","green" )))))


Path_Metric<-c('Avg. Path Length','Diameter','Diam.- numeric weights','Diam.- Inv Weights')
Path_Metric_Vals<-c(average.path.length(g_SAPSub_simpl, directed=TRUE),
                    diameter(g_SAPSub_simpl),
                    diameter(g_SAPSub_simpl, weights= num_weight),
                    diameter(g_SAPSub_simpl, weights= inv_weight))

df1<-data.frame(Path_Metric,Path_Metric_Vals)
formattable(df1,align = c("l",rep("r", NCOL(df1) - 1)),
            list(`Path_Metric` = formatter("span", style = ~ style(color = "black")), 
                 `Path_Metric_Vals` = color_bar("#add8e6")))



t1<-table((sapply(maximal.cliques(g_SAPSub_simpl), length)))
clique_number<-c(3320,335,39,5)
clique_size<-c(2,3,4,5)
df3<-data.frame(clique_size,clique_number)
formattable(df3,align = c("l",rep("r", NCOL(df3) - 1)),
            list(`Path_Metric` = formatter("span", style = ~ style(color = "black")), 
                 `Path_Metric_Vals` = color_bar("#add8e6")))


par(mfrow=c(1, 1))
fourcliques <- max_cliques(g_SAPSub_simpl, min = 4, max = 4, subset = NULL,file = NULL)
fourcliquesall<- c()

for (i in c(1:39)) {
  fourcliquesall <- c(fourcliquesall,fourcliques[[i]])
}


fourcliquesall

g4 <- induced.subgraph(graph=g_SAPSub_simpl,vids=(fourcliquesall))
plot(g4,vertex.label=V(g4)$name,vertex.size=15,vertex.color='lightblue',edge.color='black',
     edge.arrow.width=0.3,vertex.label.cex=0.7, main="Cliques of size 4")




# Embeddedness/ inverse of structural hole access (see Burt 2004)
constraints_SAP <- round(constraint(g_SAPSub_simpl, nodes=V(g_SAPSub_simpl)), digits=4)
constraint_df<-data.frame(constraints_SAP)
constraint_df<-data.frame(rownames(constraint_df),constraints_SAP)
colnames(constraint_df)<-c('Vertex','constraints_SAP')

top25_str_holes<-constraint_df[order(constraint_df$constraints_SAP),][1:25,]
v1<-as.character((top25_str_holes$Vertex))
top25_str_holes_subgraph<-induced_subgraph(g_SAPSub_simpl,vids=v1)

plot(top25_str_holes_subgraph,vertex.size=15,vertex.color='lightblue',edge.color='black',
     edge.arrow.width=0.3,vertex.label.cex=0.7) + title("Top-25 Nodes with Highest Structural Hole Access ")

# Degree centrality
degree_sap <- degree(g_SAPSub_simpl)
degree_df<-data.frame(degree_sap)
degree_df<-data.frame(rownames(degree_df),degree_sap)
colnames(degree_df)<-c('Vertex','Degree')

top20_degree<-degree_df[order(degree_df$Degree,decreasing = T),][1:20,]
d1<-as.character(top20_degree$Vertex)
top20_degree_subgraph<-induced_subgraph(g_SAPSub_simpl,vids=d1)

plot(top20_degree_subgraph,vertex.size=degree(top20_degree_subgraph)*3,vertex.color='lightblue',edge.color='black',
     edge.arrow.width=0.3,vertex.label.cex=0.7) + title("Top-20 Nodes with Highest Degree")


# Node betweenness
betweens_SAP <- round(betweenness(g_SAPSub_simpl, v=V(g_SAPSub_simpl), directed = TRUE, nobigint =TRUE, normalized = FALSE))
between_sap_df<-data.frame(betweens_SAP)
between_sap_df<-data.frame(rownames(between_sap_df),betweens_SAP)
colnames(between_sap_df)<-c('Vertex','Betweenness')

top10_betweenness<-between_sap_df[order(between_sap_df$Betweenness,decreasing = T),][1:10,]
b1<-as.character(top10_betweenness$Vertex)
top10_betweenness_subgraph<-induced_subgraph(g_SAPSub_simpl,vids=b1)

plot(top10_betweenness_subgraph,vertex.color='lightblue',edge.color='black',
     edge.arrow.width=0.3,vertex.label.cex=0.7) + title("Top-10 Nodes with Highest Betweenness")

# Edge betwenness
edgebetweens_SAP<-edge.betweenness(g_SAPSub_simpl, e=E(g_SAPSub_simpl), directed = TRUE)
order(edgebetweens_SAP,decreasing = T)

# Local clustering coefficients
clustering_SAP <- transitivity(g_SAPSub_simpl, type="local", vids=V(g_SAPSub_simpl)) 



agirvan.comm <- edge.betweenness.community(g_SAPSub_simpl, weights=inv_weight)
Girvan_size <- sizes(girvan.comm)
c.m.edge <- membership(girvan.comm)

plot(girvan.comm, g_SAPSub_simpl, vertex.edge= NA, vertex.size=10,vertex.label=NA)+ title(main = "Girvan Newman Algorithm")

comm.det.walktrap=walktrap.community(g_SAPSub_simpl,weights=E(g_SAPSub_simpl)$weight)
c.m.walktrap <- membership(comm.det.walktrap)
plot(comm.det.walktrap,g_SAPSub_simpl,vertex.label=NA, vertex.size=10)+ title(main = "Walk-Trap Algorithm")






