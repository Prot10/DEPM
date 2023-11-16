# Data download -----------------------------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)

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



# Data cleaning -----------------------------------------------------------

ncol(rna.expr.data.C)
head(colnames(rna.expr.data.C))
head(substr(colnames(rna.expr.data.N), 1,12))

dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #no duplicates 
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #duplicates!

patients.C <- substr(colnames(rna.expr.data.C), 1,12)
sort(table(patients.C)) 

unique.patients.C <- names(which(table(patients.C) == 1))
# Let's get their index in the list of patients
idx.unique.pats <- match(unique.patients.C, substr(colnames(rna.expr.data.C), 1,12) )

# Keeping patients with one sample
expr.C <- as.data.frame(rna.expr.data.C[,idx.unique.pats])
expr.N <- as.data.frame(rna.expr.data.N)

#let's rename patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1,12)
unique(colnames(expr.C))

colnames(expr.N) <- substr(colnames(expr.N), 1,12)
unique(colnames(expr.N))

intersect(colnames(expr.N), colnames(expr.C))
setdiff(colnames(expr.N), colnames(expr.C))

# 3 normal sample do not have a cancerous sample to compare to. Let's remove them
match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N)) #idx to remove
expr.N <- expr.N[,-c(10,15,20)]

length(intersect(colnames(expr.N), colnames(expr.C)))

#let's check the actual counts
typeof(expr.C[1,1]) #ok
any(is.na(expr.C)) #ok
any(is.nan(as.matrix(expr.C))) #ok

typeof(expr.N[1,1]) #ok
any(is.na(expr.N)) #ok
any(is.nan(as.matrix(expr.N))) #ok

#let's consider only patients for which we have both normal and cancer samples
expr.C <- expr.C[, colnames(expr.N)]

