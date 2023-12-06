# Libraries ---------------------------------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(DESeq2)
library(ggplot2)
library(psych)
library(NetworkToolbox)
library(ggnet)
library(GGally)
library(sna)
library(network)
library(dplyr)
library(gridExtra)


# 1. Data download -----------------------------------------------------------

proj <- "TCGA-KIRP"
dir.create(file.path(proj))

# Cancer
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")
rna.expr.data.C <- assay(rna.data.C)
genes.info.C <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

# Normal
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata")
rna.expr.data.N <- assay(rna.data.N)
genes.info.N <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))

# Clinical
clinical.query <- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)

# Check
all(na.omit(genes.info.N) == na.omit(genes.info.C))
dim(rna.expr.data.C)
dim(rna.expr.data.N)
length(unique(clinical.query$submitter_id))


save(rna.query.C, rna.data.C, rna.expr.data.C, genes.info.C, 
     rna.query.N, rna.data.N, rna.expr.data.N, genes.info.N,
     clinical.query, file="MyData.RData")

# 1.1 Data cleaning -----------------------------------------------------------

load("MyData.RData")

# Initial checks
ncol(rna.expr.data.C)
head(colnames(rna.expr.data.C))
head(substr(colnames(rna.expr.data.N), 1, 12))

dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1, 12))) # no duplicates 
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1, 12))) # no duplicates 

# Keeping patients with one sample
expr.C <- as.data.frame(rna.expr.data.C)
expr.N <- as.data.frame(rna.expr.data.N)

# Rename patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1, 12)
unique(colnames(expr.C))

colnames(expr.N) <- substr(colnames(expr.N), 1, 12)
length(unique(colnames(expr.N)))

intersect(colnames(expr.N), colnames(expr.C))
setdiff(colnames(expr.N), colnames(expr.C))

# Check the actual counts
typeof(expr.C[1,1])            # ok
any(is.na(expr.C))             # ok
any(is.nan(as.matrix(expr.C))) # ok

typeof(expr.N[1,1])            # ok
any(is.na(expr.N))             # ok
any(is.nan(as.matrix(expr.N))) # ok

# Consider only patients for which we have both normal and cancer samples
expr.C <- expr.C[, colnames(expr.N)]


# 1.2 Normalizing data with Deseq2 --------------------------------------------

all(rownames(expr.C) == rownames(expr.N))
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)
n_col     <- ncol(full.data)
dim(full.data)

metad <- c(rep("normal", n_col/2), rep("cancer", n_col/2))
metad
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])
full.data <- cbind(rownames(full.data), full.data)
ncol(full.data)
ncol(metad)

dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)

dim(counts(dds))

# Filtering: at least ten counts on 90% of patients?
threshold <- round((n_col*90) / 100)
keep <- rowSums(counts(dds) >= 10) >= threshold
dds  <- dds[keep,]
dim(counts(dds))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == n_col) # no null rows

filtr.expr.n <- as.data.frame(normalized_counts[, 1:(n_col/2)])
filtr.expr.c <- as.data.frame(normalized_counts[, 20:n_col])
# Cancerous sample names were added a ".1" in full.data because they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1, 12)




# 1.3 Gene selection ----------------------------------------------------------

genes <- read.csv2("targets.csv", row.names = 1)
genes <- genes[,1]
head(genes) 

genes.c <- intersect(rownames(filtr.expr.c), 
                     genes.info.C[ genes.info.C$gene_name %in% genes , "gene_id"]   ) 
genes.n <- intersect(rownames(filtr.expr.n),  
                     genes.info.N[ genes.info.N$gene_name %in% genes , "gene_id"]   )  

setdiff(genes.c, genes.n)

length(genes)
length(genes.c)
length(genes.n)

filtr.expr.n <- filtr.expr.n[genes.n, ]
filtr.expr.c <- filtr.expr.c[genes.c, ]

rownames(filtr.expr.c) <- genes.info.C[genes.c, "gene_name"]
rownames(filtr.expr.n) <- genes.info.N[genes.n, "gene_name"]



# 2. Differentially expressed genes (DEGs) -----------------------------------

fc <- log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n)) 
names(fc) <- rownames(filtr.expr.c)
head(fc)

dim(filtr.expr.c)
dim(filtr.expr.n)

pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i, ], filtr.expr.n[i, ] ))$p.value)
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")

expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[,1] <- round(expr.table[, 1], 2)

fc.threshold   <- 1
pval.threshold <- 0.05
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= fc.threshold & expr.table$pval.fc.fdr <= pval.threshold,]) 
length(deg.genes)

head(expr.table[deg.genes, ], 10)
#write.table(expr.table[deg.genes,], file = "DEG.csv", sep = ";")

#------------------#
#-- Volcano plot --#
#------------------#

expr.table$diffexpressed <- "NO";
expr.table$diffexpressed[expr.table$fc >= fc.threshold & expr.table$pval.fc.fdr <= pval.threshold] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -fc.threshold & expr.table$pval.fc.fdr <= pval.threshold] <- "DOWN"
head(expr.table)

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
summary(expr.table$diffexpressed)

custom_colors <- c("NO" = "#CFCFCF", "UP" = "#FF9999", "DOWN" = "#99FF99")
pastel_orange <- "#FFD1A1"

ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed)) +  
  geom_point() +
  scale_color_manual(values=custom_colors, name="Differential expressed") +
  xlab("Fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(pval.threshold), col=pastel_orange) +
  geom_vline(xintercept=fc.threshold, col=pastel_orange) +
  geom_vline(xintercept=-fc.threshold, col=pastel_orange) +
  ggtitle("Volcano Plot") +
  theme_bw() +
  theme(legend.position=c(0.95, 0.95), 
        legend.justification=c(1, 1), 
        legend.background = element_rect(color="black", linewidth=0.2, linetype="solid"), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

# print and enrichment 
cat(deg.genes, sep = "\n")


# 3. Co-expression networks ----------------------------

#-----------------#
#-- Computation --#
#-----------------#

rho.threshold <- 0.65

# Cancer network
filtr.expr.c <- filtr.expr.c[deg.genes, ] # taking only the significant genes
log2.filtr.expr.c <- log2(filtr.expr.c + 1)
cor.mat.c <- cor(t(log2.filtr.expr.c), method="pearson")
diag(cor.mat.c) <- 0
adj.mat.c  <- ifelse(abs(cor.mat.c) < rho.threshold, 0, 1)

# Normal network
filtr.expr.n <- filtr.expr.n[deg.genes, ] # taking only the significant genes
log2.filtr.expr.n <- log2(filtr.expr.n + 1)
cor.mat.n <- cor(t(log2.filtr.expr.n), method="pearson")
diag(cor.mat.n) <- 0
adj.mat.n <- ifelse(abs(cor.mat.n) < rho.threshold, 0, 1)

#--------------#
#-- Analysis --#
#--------------#

#-- Compute the degree index and check if the network is a scale free network --#

# Cancer network 
net.c <- network(adj.mat.c, matrix.type="adjacency", 
                 ignore.eval=FALSE, names.eval="weights")

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)

# Normal network
net.n <- network(adj.mat.n, matrix.type="adjacency",
                 ignore.eval=FALSE, names.eval = "weights")
degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing=T)

#------------------------------#
#-- Degree distribution plot --#
#------------------------------#

# Preparing the data for ggplot
cancer_data <- data.frame(degree = degree.c, type = "Cancer")
normal_data <- data.frame(degree = degree.n, type = "Normal")
combined_data <- rbind(cancer_data, normal_data)
percentiles <- combined_data %>% 
  group_by(type) %>% 
  summarize(p95 = quantile(degree[degree > 0], 0.95))

# Plotting in a grid
plot_cancer <- ggplot(subset(combined_data, type == "Cancer"), aes(x=degree)) +
  geom_histogram(binwidth = 9, fill="#FF9999", color="black") +
  geom_vline(data=subset(percentiles, type == "Cancer"), aes(xintercept=p95),
             color="darkred", linetype="dashed", linewidth=0.7) +
  ggtitle("Degree Distribution: Cancer Network") +
  xlab("Degree") + ylab("Frequency") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

plot_normal <- ggplot(subset(combined_data, type == "Normal"), aes(x=degree)) +
  geom_histogram(binwidth = 5, fill="#99CCFF", color="black") +
  geom_vline(data=subset(percentiles, type == "Normal"), aes(xintercept=p95),
             color="darkblue", linetype="dashed", linewidth=0.5) +
  ggtitle("Degree Distribution: Normal Network") +
  xlab("Degree") + ylab("Frequency") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.7))

grid.arrange(plot_cancer, plot_normal, ncol=2)

#-- If so, find the hubs (5% of the nodes with highest degree values), compare hubs sets related to the two condition (cancer, normal) and identify the hubs selectively characterizing each network --#

# Custom colors
custom_colors <- c("hub" = "#FF9999", "non-hub" = "#CFCFCF")
pastel_blue <- "#ADD8E6"
pastel_green <- "#99FF99"

# Cancer network

# Identify hubs
x <- quantile(degree.c[degree.c > 0], 0.95)
hubs.c <- degree.c[degree.c >= x]
names(hubs.c)

# Set vertex and edge attributes
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c), "hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net.c %v% "size" = (rowSums(adj.mat.c != 0) / max(rowSums(adj.mat.c != 0)) + 0.1) * 10
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, pastel_blue, pastel_green))

# Compute the layout
coord.c <- gplot.layout.fruchtermanreingold(net.c, NULL)
net.c %v% "x" = coord.c[, 1]
net.c %v% "y" = coord.c[, 2]

# Plotting the network
ggnet2(net.c, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Cancer Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))

library(network)
library(GGally)
library(sna)
library(ggplot2)
hub_names <- names(hubs.c)

# Set vertex and edge attributes
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% hub_names, "hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net.c %v% "size" = (rowSums(adj.mat.c != 0) / max(rowSums(adj.mat.c != 0)) + 0.1) * 10
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, pastel_blue, pastel_green))

# Compute the layout
coord.c <- gplot.layout.fruchtermanreingold(net.c, NULL)
net.c %v% "x" = coord.c[, 1]
net.c %v% "y" = coord.c[, 2]

net_df <- fortify(net.c)

# Filter for hub nodes
hub_df <- net_df[net_df$label %in% hub_names, ]

# Plotting the network with labels for hubs
gg <- ggnet2(net.c, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
             edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5) +
  geom_text(data = hub_df, aes(label = label, x = x, y = y), vjust = -1, size = 3) +
  scale_color_manual(values = custom_colors, name = "Node Type") +
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("Cancer Network") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(gg)


# Normal network 

# Identify hubs
y <- quantile(degree.n[degree.n > 0], 0.95)
hubs.n <- degree.n[degree.n >= y]
names(hubs.n)

# Set vertex and edge attributes
net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n), "hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net.n %v% "size" = (rowSums(adj.mat.n != 0) / max(rowSums(adj.mat.n != 0)) + 0.1) * 10
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, pastel_blue, pastel_green))

# Compute the layout
coord.n <- gplot.layout.fruchtermanreingold(net.n, NULL)
net.n %v% "x" = coord.n[, 1]
net.n %v% "y" = coord.n[, 2]

# Plotting the network
ggnet2(net.n, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Differential Co-Expressed Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))

intersect(names(hubs.c), names(hubs.n))



# 4. Differential Co-expressed Network ---------------------------------------

# Fisher z-transformation
z.mat.c <- atanh(cor.mat.c)
z.mat.n <- atanh(cor.mat.n)

# Sample sizes for each condition
n.c <- ncol(filtr.expr.c)
n.n <- ncol(filtr.expr.n)

# Compute z-scores for differential correlation
z.scores <- (z.mat.c-z.mat.n) / sqrt((1/(n.c-3)) + (1/(n.n-3)))

# Apply the threshold |Z| < 3
z.threshold <- 3
diff.cor.mat <- ifelse(abs(z.scores) < z.threshold, 0, 1)

# Custom colors
custom_colors <- c("hub" = "#FF9999", "non-hub" = "#CFCFCF")
pastel_orange <- "#FFD1A1"
pastel_blue <- "#ADD8E6"
pastel_green <- "#99FF99"

# Creating the network
net <- network(diff.cor.mat, matrix.type="adjacency", 
               ignore.eval=FALSE, names.eval="weights")
degree <- rowSums(adj.mat.c != 0)
names(degree) <- rownames(diff.cor.mat)
degree <- sort(degree, decreasing = T)

# Identify hubs
z <- quantile(degree[degree > 0], 0.95)
hubs <- degree[degree >= z]
names(hubs)

# Set vertex and edge attributes
net %v% "type" = ifelse(network.vertex.names(net) %in% names(hubs), "hub", "non-hub")
net %v% "color" = ifelse(net %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net %v% "size" = (rowSums(adj.mat.c != 0) / max(rowSums(adj.mat.c != 0)) + 0.1) * 10
set.edge.attribute(net, "edgecolor", ifelse(net %e% "weights" > 0, pastel_blue, pastel_green))

coord <- gplot.layout.fruchtermanreingold(net, NULL)
net %v% "x" = coord[, 1]
net %v% "y" = coord[, 2]

# Plotting the network
ggnet2(net, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Differential Co-Expressed Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))

# Identify hubs: top 5% nodes with highest degree
num.hubs <- ceiling(length(degree) * 0.05)
hubs <- names(degree[order(-degree)[1:num.hubs]])
hubs


# 5. Patient Similarity Network (PSN) -------------------------------------

library(igraph)

# Compute the Euclidean distance matrix for all significant genes
distance_matrix.c <- as.matrix(dist(t(filtr.expr.c)))

# Convert distances to similarities (normalized and inversed)
similarity_matrix.c <- 1 - (distance_matrix.c / max(distance_matrix.c))
diag(similarity_matrix.c) <- 0

psn.c <- graph_from_adjacency_matrix(similarity_matrix.c, mode =  "undirected", weighted = TRUE)

clusterlouvain.c <- cluster_louvain(psn.c, resolution=1)
num_communities.c <- length(unique(clusterlouvain.c$membership))

coords.c <- layout_with_fr(psn.c)
plot(psn.c, vertex.color=rainbow(num_communities.c, alpha=0.6)[clusterlouvain.c$membership],
     layout=coords.c, vertex.label=NA, edge.width=(E(psn.c)$weight), edge.color='black')
title('Community detection of patients based\non the similarity network\n(all genes with euclidian)')

f.c <- subgraph.edges(graph=psn.c, eids = E(psn.c)[.from(which(clusterlouvain$membership==1)[5])])
la.c <- layout_with_fr(f.c)
plot(f.c, vertex.color=c('red','cyan'), vertex.label=NA, 
     edge.width=(E(f.c)$weight*5), edge.color='black', layout=la.c)
title('Visualization of the difference in link weights
      from one selected node to all the others')

# Bonus points ------------------------------------------------------------


# Bonus 1- Compute a different centrality index (CI) ----------------------------

install.packages("igraph")  # Install igraph if you haven't already
library(igraph)             # Load igraph package

top5pct_degree.c <- head(degree.c, length(degree.c) * 0.05)

# betweenness
betweenness.c <- betweenness(net.c)#, directed = FALSE)
names(betweenness.c) <- rownames(adj.mat.c)
betweenness.c <- sort(betweenness.c, decreasing = TRUE)
top5pct_betweenness.c <- head(betweenness.c, length(betweenness.c) * 0.05)

# closeness
closeness.c <- closeness(net.c)#, mode = "all")
names(closeness.c) <- rownames(adj.mat.c)
closeness.c <- sort(closeness.c, decreasing = TRUE)
top5pct_closeness.c <- head(closeness.c, length(closeness.c) * 0.05)

# pagerank
# Assuming adj.mat.c is your adjacency matrix
net.c.graph <- graph_from_adjacency_matrix(adj.mat.c, mode = "undirected")
pagerank.c <- page.rank(net.c.graph)$vector
names(pagerank.c) <- rownames(adj.mat.c)
pagerank.c <- sort(pagerank.c, decreasing = TRUE)
top5pct_pagerank.c <- head(pagerank.c, length(pagerank.c) * 0.05)

# overlap with degree based hubs
overlap_betweenness.c <- intersect(names(top5pct_betweenness.c), names(top5pct_degree.c))
overlap_betweenness.c
overlap_closeness.c <- intersect(names(top5pct_closeness.c), names(top5pct_degree.c))
overlap_closeness.c
overlap_pagerank.c <- intersect(names(top5pct_pagerank.c), names(top5pct_degree.c))
overlap_pagerank.c

# overlap between all the metrics
intersect(intersect(overlap_betweenness.c, overlap_closeness.c), overlap_pagerank.c)


# Bonus 2- Perform the study using a different similarity measure ---------------

#--------------------------------#
#-- Using Spearman Correlation --#
#--------------------------------#

# Bonus 2.3- Co-expression networks ----------------------------

#-----------------#
#-- Computation --#
#-----------------#

rho.threshold <- 0.65

# Cancer network
filtr.expr.c <- filtr.expr.c[deg.genes, ] # taking only the significant genes
log2.filtr.expr.c <- log2(filtr.expr.c + 1)
cor.mat.c <- cor(t(log2.filtr.expr.c), method="spearman")
diag(cor.mat.c) <- 0
adj.mat.c  <- ifelse(abs(cor.mat.c) < rho.threshold, 0, 1)

# Normal network
filtr.expr.n <- filtr.expr.n[deg.genes, ] # taking only the significant genes
log2.filtr.expr.n <- log2(filtr.expr.n + 1)
cor.mat.n <- cor(t(log2.filtr.expr.n), method="spearman")
diag(cor.mat.n) <- 0
adj.mat.n <- ifelse(abs(cor.mat.n) < rho.threshold, 0, 1)

#--------------#
#-- Analysis --#
#--------------#

#-- Compute the degree index and check if the network is a scale free network --#

# Cancer network 
net.c <- network(adj.mat.c, matrix.type="adjacency", 
                 ignore.eval=FALSE, names.eval="weights")

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)

# Normal network
net.n <- network(adj.mat.n, matrix.type="adjacency",
                 ignore.eval=FALSE, names.eval = "weights")
degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing=T)

#------------------------------#
#-- Degree distribution plot --#
#------------------------------#

# Preparing the data for ggplot
cancer_data <- data.frame(degree = degree.c, type = "Cancer")
normal_data <- data.frame(degree = degree.n, type = "Normal")
combined_data <- rbind(cancer_data, normal_data)
percentiles <- combined_data %>% 
  group_by(type) %>% 
  summarize(p95 = quantile(degree[degree > 0], 0.95))

# Plotting in a grid
plot_cancer <- ggplot(subset(combined_data, type == "Cancer"), aes(x=degree)) +
  geom_histogram(binwidth = 7, fill="#FF9999", color="black") +
  geom_vline(data=subset(percentiles, type == "Cancer"), aes(xintercept=p95),
             color="darkred", linetype="dashed", linewidth=0.7) +
  ggtitle("Degree Distribution: Cancer Network") +
  xlab("Degree") + ylab("Frequency") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

plot_normal <- ggplot(subset(combined_data, type == "Normal"), aes(x=degree)) +
  geom_histogram(binwidth = 5, fill="#99CCFF", color="black") +
  geom_vline(data=subset(percentiles, type == "Normal"), aes(xintercept=p95),
             color="darkblue", linetype="dashed", linewidth=0.5) +
  ggtitle("Degree Distribution: Normal Network") +
  xlab("Degree") + ylab("Frequency") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.7))

grid.arrange(plot_cancer, plot_normal, ncol=2)

#-- If so, find the hubs (5% of the nodes with highest degree values), compare hubs sets related to the two condition (cancer, normal) and identify the hubs selectively characterizing each network --#

# Custom colors
custom_colors <- c("hub" = "#FF9999", "non-hub" = "#CFCFCF")
pastel_blue <- "#ADD8E6"
pastel_green <- "#99FF99"

# Cancer network

# Identify hubs
x <- quantile(degree.c[degree.c > 0], 0.95)
hubs.c <- degree.c[degree.c >= x]
names(hubs.c)

# Set vertex and edge attributes
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c), "hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net.c %v% "size" = (rowSums(adj.mat.c != 0) / max(rowSums(adj.mat.c != 0)) + 0.1) * 10
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, pastel_blue, pastel_green))

# Compute the layout
coord.c <- gplot.layout.fruchtermanreingold(net.c, NULL)
net.c %v% "x" = coord.c[, 1]
net.c %v% "y" = coord.c[, 2]

# Plotting the network
ggnet2(net.c, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Cancer Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))


# Normal network 

# Identify hubs
y <- quantile(degree.n[degree.n > 0], 0.95)
hubs.n <- degree.n[degree.n >= y]
names(hubs.n)

# Set vertex and edge attributes
net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n), "hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net.n %v% "size" = (rowSums(adj.mat.n != 0) / max(rowSums(adj.mat.n != 0)) + 0.1) * 10
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, pastel_blue, pastel_green))

# Compute the layout
coord.n <- gplot.layout.fruchtermanreingold(net.n, NULL)
net.n %v% "x" = coord.n[, 1]
net.n %v% "y" = coord.n[, 2]

# Plotting the network
ggnet2(net.n, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Differential Co-Expressed Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))

intersect(names(hubs.c), names(hubs.n))


# Bonus 3- Perform gene set enrichment analysis ---------------------------------



# Bonus 4- Task 5 using gene expression profiles related to normal --------------

# Compute the Euclidean distance matrix for all significant genes
distance_matrix.n <- as.matrix(dist(t(filtr.expr.n)))

# Convert distances to similarities (normalized and inversed)
similarity_matrix.n <- 1 - (distance_matrix.n / max(distance_matrix.n))
diag(similarity_matrix.n) <- 0

psn.n <- graph_from_adjacency_matrix(similarity_matrix.n, mode="undirected", weighted = TRUE)

clusterlouvain.n <- cluster_louvain(psn, resolution=1)
num_communities.n <- length(unique(clusterlouvain.n$membership))

coords.n <- layout_with_fr(psn.n)
plot(psn.n, vertex.color=rainbow(num_communities.n, alpha=0.6)[clusterlouvain.n$membership],
     layout=coords.n, vertex.label=NA, edge.width=(E(psn.n)$weight), edge.color='black')
title('Community detection of patients based\non the similarity network\n(all genes with euclidian)')

#plotting a subset to highlight the edges width based on their weights (similarity measure)
f.n  <- subgraph.edges(graph=psn.n, eids = E(psn.n)[.from(which(clusterlouvain.n$membership==1)[5])])
la.n <- layout_with_fr(f.n)
plot(f.n, vertex.color=c('red', 'cyan'), vertex.label=NA, edge.width=(E(f)$weight*5),
     edge.color='black', layout=la.n)
title('Visualization of the difference in link weights
      from one selected node to all the others')

# Bonus 5- Perform PSN communities characterization -----------------------------

