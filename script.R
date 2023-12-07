# Libraries ---------------------------------------------------------------

if ("package:igraph" %in% search()) detach("package:igraph", unload = TRUE, character.only = TRUE)
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
library(RColorBrewer)

set.seed(123)

# Define a palette
pal <- c("olivedrab3", "olivedrab4", "dodgerblue3", "dodgerblue4",
         "darkorange1", "darkorange3", "firebrick1", "firebrick3",
         "gray75", "gray65")


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
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1, 12))) # no duplicates 

dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1, 12))) # no duplicates 

# Keeping patients with one sample
expr.C <- as.data.frame(rna.expr.data.C)
expr.N <- as.data.frame(rna.expr.data.N)

# Rename patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1, 12)
length(unique(colnames(expr.C))) == length(unique(substr(colnames(rna.expr.data.C), 1, 12)))

colnames(expr.N) <- substr(colnames(expr.N), 1, 12)
length(unique(colnames(expr.N))) == length(unique(substr(colnames(rna.expr.data.N), 1, 12)))

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
expr.C <- expr.C[, intersect(colnames(expr.N), colnames(expr.C))]
expr.N <- expr.N[, intersect(colnames(expr.N), colnames(expr.C))]
dim(expr.C)
dim(expr.N)

# 1.2 Normalizing data with Deseq2 --------------------------------------------

all(rownames(expr.C) == rownames(expr.N))
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)
n_col     <- ncol(full.data)
dim(full.data)

metad <- c(rep("normal", n_col/2), rep("cancer", n_col/2))
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

# Filtering: at least ten counts on 90% of patients
threshold <- round((n_col*90) / 100)
keep <- rowSums(counts(dds) >= 10) >= threshold
dds  <- dds[keep,]
dim(counts(dds))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == n_col) # no null rows

filtr.expr.n <- as.data.frame(normalized_counts[, 1:(n_col/2)])
filtr.expr.c <- as.data.frame(normalized_counts[, (n_col/2+1):n_col])
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

fc.threshold   <- 1.2
pval.threshold <- 0.05
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= fc.threshold & expr.table$pval.fc.fdr <= pval.threshold,]) 
length(deg.genes) # number of significant genes

head(expr.table[deg.genes, ], 10)

#------------------#
#-- Volcano plot --#
#------------------#

expr.table$diffexpressed <- "NO";
expr.table$diffexpressed[expr.table$fc >= fc.threshold & expr.table$pval.fc.fdr <= pval.threshold] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -fc.threshold & expr.table$pval.fc.fdr <= pval.threshold] <- "DOWN"
head(expr.table)

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
summary(expr.table$diffexpressed)


custom_colors1 <- c("NO" = pal[9], "UP" = pal[3], "DOWN" = pal[7])
custom_colors2 <- c("NO" = pal[10], "UP" = pal[4], "DOWN" = pal[8])

ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed)) +  
  geom_point(size = 1.85,
             alpha = 0.6,
             stroke = 0.7,
             aes(color = diffexpressed, fill = diffexpressed)) +  # Fill color for points
  scale_color_manual(values = custom_colors1, name = "Differential expressed") +
  scale_fill_manual(values = custom_colors2) +  # Darker colors for outline
  xlab("Fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept = -log10(pval.threshold), col = pal[5]) +
  geom_vline(xintercept = fc.threshold, col = pal[5]) +
  geom_vline(xintercept = -fc.threshold, col = pal[5]) +
  ggtitle("Volcano Plot") +
  xlim(-3, 5) +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95), 
        legend.justification = c(1, 1), 
        legend.background = element_rect(color = "black", linewidth = 0.2, linetype = "solid"), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")) +
  guides(fill = FALSE)



# 3. Co-expression networks ----------------------------

#-----------------#
#-- Computation --#
#-----------------#

rho.threshold <- 0.675

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

plot_cancer <- ggplot(subset(combined_data, type == "Cancer"), aes(x=degree)) +
  geom_histogram(binwidth = 4, fill=pal[7], color="darkred") +
  geom_vline(data=subset(percentiles, type == "Cancer"), aes(xintercept=p95),
             color=pal[6], linetype="dashed", linewidth=0.6) +
  ggtitle("Degree Distribution: Cancer Network") +
  xlab("Degree") + ylab("Frequency") +
  geom_text(data=subset(percentiles, type=="Cancer"), aes(x=p95, y=60, label="Quantile\nat 95%"),
            color="black", vjust=-0.5, hjust=1.2) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

plot_normal <- ggplot(subset(combined_data, type == "Normal"), aes(x=degree)) +
  geom_histogram(binwidth = 5, fill=pal[3], color="darkblue") +
  geom_vline(data=subset(percentiles, type == "Normal"), aes(xintercept=p95),
             color=pal[6], linetype="dashed", linewidth=0.6) +
  ggtitle("Degree Distribution: Normal Network") +
  xlab("Degree") + ylab("Frequency") +
  geom_text(data=subset(percentiles, type == "Normal"), aes(x=p95, y=30.1, label="Quantile\nat 95%"),
            color="black", vjust=-0.5, hjust=1.2) +
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

label.nodes.c <- ifelse(network.vertex.names(net.c) %in% names(hubs.c), names(hubs.c), "")

# Plotting the network
ggnet2(net.c, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5,
       label = label.nodes.c, label.size = 2.5) +
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

label.nodes.n <- ifelse(network.vertex.names(net.n) %in% names(hubs.n), names(hubs.n), "")

# Plotting the network
ggnet2(net.n, color = "color", alpha = 0.7, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.8, edge.size = 0.5,
       label = label.nodes.n, label.size = 2.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Normal Network") +
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
custom_colors <- c("hub"=pal[5], "non-hub"=pal[9])

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


plot_data <- data.frame(degree = degree)
percentiles <- plot_data %>% 
  summarize(p95 = quantile(degree[degree > 0], 0.95))

ggplot(plot_data, aes(x=degree)) +
        geom_histogram(binwidth = 4, fill=pal[7], color="darkred") +
        geom_vline(data=subset(percentiles), aes(xintercept=p95),
                   color=pal[6], linetype="dashed", linewidth=0.6) +
        ggtitle("Degree Distribution of the Co-Expressed Network") +
        xlab("Degree") + ylab("Frequency") +
        geom_text(data=subset(percentiles), aes(x=p95, y=60, label="Quantile\nat 95%"),
                  color="black", vjust=-0.5, hjust=1.2) +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5))


# Set vertex and edge attributes
net %v% "type" = ifelse(network.vertex.names(net) %in% names(hubs), "hub", "non-hub")
net %v% "color" = ifelse(net %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net %v% "size" = (rowSums(adj.mat.c != 0) / max(rowSums(adj.mat.c != 0)) + 0.1) * 10
set.edge.attribute(net, "edgecolor", ifelse(net %e% "weights" > 0, pal[3], pal[7]))

coord <- gplot.layout.fruchtermanreingold(net, NULL)
net %v% "x" = coord[, 1]
net %v% "y" = coord[, 2]

label.nodes <- ifelse(network.vertex.names(net) %in% names(hubs), network.vertex.names(net), "")

# Plotting the network
ggnet2(net, color = "color", alpha = 0.8, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.2, edge.size = 0.4,
       label = label.nodes, label.size = 2.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Differential Co-Expressed Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))


# 5. Patient Similarity Network (PSN) -------------------------------------

#-- 5.1 Compute the Patient Similarity Network using cancer gene expression profile --#

# Compute the Euclidean distance matrix for all significant genes
distance_matrix.c <- as.matrix(dist(t(filtr.expr.c)))

# Convert distances to similarities (normalized and inversed)
similarity_matrix.c <- 1 - (distance_matrix.c / max(distance_matrix.c))
diag(similarity_matrix.c) <- 0

psn.c <- network(similarity_matrix.c, matrix.type="adjacency", 
                ignore.eval=FALSE, names.eval="weights")


#-- 5.2 Perform the community detection (e.g. apply Louvain algorithm to the PSN) --#

l.comp.c <- component.largest(psn.c, result = "graph")
l.comp.c <- similarity_matrix.c[rownames(l.comp.c), rownames(l.comp.c)]

write.csv2(l.comp.c, "input-matrix-c.csv")

# Let's open the terminal 
# pip install bctpy
# python3 btc-community-c.py input-matrix-c.csv

comm.res.c <- read.csv2("output-c.txt", header = FALSE)
rownames(comm.res.c) <- rownames(l.comp.c)

net.final.c <- network(l.comp.c, matrix.type="adjacency", ignore.eval=FALSE, names.eval="weights", directed=F)
set.edge.attribute(net.final.c, "edgecolor", ifelse(net.final.c %e% "weights" > 0, pal[9], pal[10]))

all(net.final.c %v% "vertex.names" == rownames(comm.res.c)) # ok
net.final.c %v% "community" <-  as.character(comm.res.c[,1])

#let's reverse the membership information
comm.res.c <- cbind(comm.res.c, rownames(comm.res.c))
colnames(comm.res.c) <- c("com", "gene") 

n_comm.c <- length(unique(net.final.c %v% "community"))
communities.c <- vector("list", n_comm.c)

for (i in 1:n_comm.c) communities.c[[i]] <- comm.res.c[comm.res.c$com == i, 2]

# Force positioning 
copy.comp.c <- l.comp.c*0

# Add strong weight between members of the same community 
for (cindex in 1:n_comm.c){
  copy.comp.c[communities.c[[cindex]], 
              communities.c[[cindex]]] <- 1}

copy.net.c <- network(copy.comp.c,  matrix.type="adjacency", ignore.eval=FALSE,
                      names.eval="weights", directed=F)
coord.comm.c <- gplot.layout.fruchtermanreingold(copy.net.c, NULL)

net.final.c %v% "x" = coord.comm.c[, 1]
net.final.c %v% "y" = coord.comm.c[, 2]

names(pal) <- c(4, 6, 5, 2, 1, 7, 8, 3, 9, 10)

ggnet2(net.final.c, color="community", alpha = 0.85, palette=pal,
       mode=c("x", "y"), edge.color="edgecolor", edge.alpha=0.6, edge.size=0.25) +
  guides(size="none") +
  theme_minimal() +
  ggtitle("Patient Similarity Network\nusing Cancer's Gene") +
  labs(color = "Community") +
  theme(plot.title = element_text(hjust = 0.5))


# Bonus points ------------------------------------------------------------


# Bonus 1- Compute a different centrality index (CI) ----------------------------

top5pct_degree.c <- head(degree.c, length(degree.c) * 0.05)

# Betweenness
betweenness.c <- betweenness(net.c)#, directed = FALSE)
names(betweenness.c) <- rownames(adj.mat.c)
betweenness.c <- sort(betweenness.c, decreasing = TRUE)
top5pct_betweenness.c <- head(betweenness.c, length(betweenness.c) * 0.05)
names(top5pct_betweenness.c)

# Closeness
closeness.c <- closeness(net.c)#, mode = "all")
names(closeness.c) <- rownames(adj.mat.c)
closeness.c <- sort(closeness.c, decreasing = TRUE)
top5pct_closeness.c <- head(closeness.c, length(closeness.c) * 0.05)
names(top5pct_closeness.c)

# Overlap with degree based hubs
overlap_betweenness.c <- intersect(names(top5pct_betweenness.c), names(top5pct_degree.c))
overlap_betweenness.c
overlap_closeness.c <- intersect(names(top5pct_closeness.c), names(top5pct_degree.c))
overlap_closeness.c

# overlap between all the metrics
intersect(overlap_betweenness.c, overlap_closeness.c)


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
  geom_histogram(binwidth = 4, fill=pal[7], color="darkred") +
  geom_vline(data=subset(percentiles, type == "Cancer"), aes(xintercept=p95),
             color=pal[6], linetype="dashed", linewidth=0.6) +
  ggtitle("Degree Distribution: Cancer Network") +
  geom_text(data=subset(percentiles, type=="Cancer"), aes(x=p95, y=60, label="Quantile\nat 95%"),
            color="black", vjust=-0.5, hjust=1.2) +
  xlab("Degree") + ylab("Frequency") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

plot_normal <- ggplot(subset(combined_data, type == "Normal"), aes(x=degree)) +
  geom_histogram(binwidth = 5, fill=pal[3], color="darkblue") +
  geom_vline(data=subset(percentiles, type == "Normal"), aes(xintercept=p95),
             color=pal[6], linetype="dashed", linewidth=0.6) +
  ggtitle("Degree Distribution: Normal Network") +
  geom_text(data=subset(percentiles, type == "Normal"), aes(x=p95, y=25, label="Quantile\nat 95%"),
            color="black", vjust=-0.5, hjust=1.2) +
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
hubs.c1 <- degree.c[degree.c >= x]
names(hubs.c1)

# Set vertex and edge attributes
net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c1), "hub", "non-hub")
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

intersect(names(hubs.c1), names(hubs.n))

# Bonus 2.4- Differential Co-expressed Network ---------------------------------------

pal <- c("olivedrab3", "olivedrab4", "dodgerblue3", "dodgerblue4",
         "darkorange1", "darkorange3", "firebrick1", "firebrick3",
         "gray75", "gray65")

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
custom_colors <- c("hub"=pal[5], "non-hub"=pal[9])

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

plot_data <- data.frame(degree = degree)
percentiles <- plot_data %>% 
  summarize(p95 = quantile(degree[degree > 0], 0.95))

ggplot(plot_data, aes(x=degree)) +
  geom_histogram(binwidth = 4, fill=pal[3], color="darkblue") +
  geom_vline(data=subset(percentiles), aes(xintercept=p95),
             color=pal[6], linetype="dashed", linewidth=0.6) +
  ggtitle("Degree Distribution of the Co-Expressed Network") +
  xlab("Degree") + ylab("Frequency") +
  geom_text(data=subset(percentiles), aes(x=p95, y=60, label="Quantile\nat 95%"),
            color="black", vjust=-0.5, hjust=1.2) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

# Set vertex and edge attributes
net %v% "type" = ifelse(network.vertex.names(net) %in% names(hubs), "hub", "non-hub")
net %v% "color" = ifelse(net %v% "type" == "hub", custom_colors["hub"], custom_colors["non-hub"])
net %v% "size" = (rowSums(adj.mat.c != 0) / max(rowSums(adj.mat.c != 0)) + 0.1) * 10
set.edge.attribute(net, "edgecolor", ifelse(net %e% "weights" > 0, pal[3], pal[7]))

coord <- gplot.layout.fruchtermanreingold(net, NULL)
net %v% "x" = coord[, 1]
net %v% "y" = coord[, 2]

label.nodes <- ifelse(network.vertex.names(net) %in% names(hubs), network.vertex.names(net), "")

# Plotting the network
ggnet2(net, color = "color", alpha = 0.8, size = "size", mode = c("x", "y"),
       edge.color = "edgecolor", edge.alpha = 0.2, edge.size = 0.4,
       label = label.nodes, label.size = 2.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Differential Co-Expressed Network") +
  labs(color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5))

# Bonus 3- Perform gene set enrichment analysis ---------------------------------

selected.genes <- names(hubs.c)

# enrichR on R 
library(enrichR)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)
if (websiteLive) dbs$libraryName

dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
if (websiteLive) {enriched <- enrichr(selected.genes, dbs)}
if (websiteLive) View(enriched[["GO_Biological_Process_2023"]])

# Plotting the results
# GO
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")}
# KEGG
if (websiteLive) {
  plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", 
             orderBy = "P.value") }


# Bonus 4- Task 5 using gene expression profiles related to normal --------------

#-- Bonus 4.1 Compute the Patient Similarity Network using cancer gene expression profile --#

# Compute the Euclidean distance matrix for all significant genes
distance_matrix.n <- as.matrix(dist(t(filtr.expr.n)))

# Convert distances to similarities (normalized and inversed)
similarity_matrix.n <- 1 - (distance_matrix.n / max(distance_matrix.n))
diag(similarity_matrix.n) <- 0

psn.n <- network(similarity_matrix.n, matrix.type="adjacency", 
                 ignore.eval=FALSE, names.eval="weights")

#-- Bonus 4.2 Perform the community detection (e.g. apply Louvain algorithm to the PSN) --#

l.comp.n <- component.largest(psn.n, result = "graph")
l.comp.n <- similarity_matrix.n[rownames(l.comp.n), rownames(l.comp.n)]

write.csv2(l.comp.n, "input-matrix-n.csv")

# Let's open the terminal 
# pip install bctpy
# python3 btc-community-n.py input-matrix-n.csv

comm.res.n <- read.csv2("output-n.txt", header = FALSE)
rownames(comm.res.n) <- rownames(l.comp.n)

net.final.n <- network(l.comp.n, matrix.type="adjacency", ignore.eval=FALSE, names.eval="weights", directed=F)
set.edge.attribute(net.final.n, "edgecolor", ifelse(net.final.n %e% "weights" > 0, pal[9], pal[10]))

all(net.final.n %v% "vertex.names" == rownames(comm.res.n)) # ok
net.final.n %v% "community" <-  as.character(comm.res.n[,1])

#let's reverse the membership information
comm.res.n <- cbind(comm.res.n, rownames(comm.res.n))
colnames(comm.res.n) <- c("com", "gene") 

n_comm.n <- length(unique(net.final.n %v% "community"))
communities.n <- vector("list", n_comm.n)

for (i in 1:n_comm.n) communities.n[[i]] <- comm.res.n[comm.res.n$com == i, 2]

# Force positioning 
copy.comp.n <- l.comp.n*0

# Add strong weight between members of the same community 
for (cindex in 1:n_comm.n){
  copy.comp.n[communities.n[[cindex]], 
              communities.n[[cindex]]] <- 1}

copy.net.n <- network(copy.comp.n,  matrix.type="adjacency", ignore.eval=FALSE,
                      names.eval="weights", directed=F)
coord.comm.n <- gplot.layout.fruchtermanreingold(copy.net.n, NULL)

net.final.n %v% "x" = coord.comm.n[, 1]
net.final.n %v% "y" = coord.comm.n[, 2]

pal <- c("olivedrab3", "olivedrab4", "dodgerblue3", "dodgerblue4",
         "darkorange1", "darkorange3", "firebrick1", "firebrick3",
         "gray75", "gray65")
names(pal) <- c(4, 6, 5, 2, 3, 7, 8, 1, 9, 10)

ggnet2(net.final.n, color="community", alpha = 0.85, palette=pal,
       mode=c("x", "y"), edge.color="edgecolor", edge.alpha=0.6, edge.size=0.25) +
  guides(size="none") +
  theme_minimal() +
  ggtitle("Patient Similarity Network\nusing Normal's Gene") +
  labs(color = "Community") +
  theme(plot.title = element_text(hjust = 0.5))

# Bonus 5- Perform PSN communities characterization -----------------------------

rownames(clinical.query) <- clinical.query$submitter_id

# Status
table(clinical.query$vital_status)

# Reformat it to 0=censored, 1=dead
clinical.query$status <- ifelse(clinical.query$vital_status == "Alive", 0, 1)

# What is overall survival? 
clinical.query$os <- ifelse(clinical.query$vital_status == "Alive", 
                            clinical.query$days_to_last_follow_up,
                            clinical.query$days_to_death)

library(survival)
library(ggsurvfit)
library(dplyr)

survfit2(Surv(os, status) ~ 1, data=clinical.query) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability")+ 
  add_confidence_interval()

#risk table
survfit2(Surv(os, status) ~ 1, data = clinical.query) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability")+ 
  add_confidence_interval()+
  add_risktable()

library(survminer)
ggsurvplot(survfit(Surv(os, status) ~ 1, data = clinical.query),
           conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() )


### Groups comparison ###

# Sex
ggsurvplot(survfit(Surv(os, status) ~ gender, data = clinical.query),
           pval = TRUE, conf.int = T,
           risk.table.col = "strata",  break.time.by = 500, 
           size = 1, legend.labs =  c("Females", "Males"),
           risk.table.height = 0.25, 
           ggtheme = theme_bw(), xlab = "Time in days", risk.table = "abs_pct",
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE,
           #ncensor.plot = TRUE,
           palette = c("darkorange1", "dodgerblue3"))

table(clinical.query$prior_malignancy)
table(clinical.query$race)

# Race
ggsurvplot(survfit(Surv(os, status) ~ race, data = clinical.query),
           pval = TRUE, conf.int = F,
           risk.table.col = "strata",  break.time.by = 500, 
           size = 1, risk.table.height = 0.25, 
           ggtheme = theme_bw(), xlab = "Time in days", risk.table = "abs_pct",
           risk.table.y.text.col = T, risk.table.y.text = FALSE,
           palette = c("#E7B800", "#2E9FDF", "#f75752", "#75a716", "#ba62fe"))

# Age (over 65)
clinical.query$over65 <- ifelse(clinical.query$age_at_index >= 65, T, F)

ggsurvplot(survfit(Surv(os, status) ~ over65, data = clinical.query),
           pval = TRUE, conf.int = T,
           risk.table.col = "strata",  break.time.by = 500, 
           size = 1, legend.labs =  c("Under 65", "Over 65"),
           risk.table.height = 0.25, 
           ggtheme = theme_bw(), xlab = "Time in days", risk.table = "abs_pct",
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE,
           #ncensor.plot = TRUE,
           palette = c("darkorange1", "dodgerblue3"))
