# Load the dataset
file_name <- "repository-files-table.2023-11-14.tsv"
data <- read.table(file = file_name, header = TRUE, sep = "\t")

# View the first rows
head(data)
