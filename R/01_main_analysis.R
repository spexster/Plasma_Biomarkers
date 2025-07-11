r
# Load necessary libraries
library(Seurat)
library(tidyverse)

r
setwd("C:/Users/farza/Documents/Plasma_Biomarkers")

r
#set seed
set.seed(333)  # Classic choice, but any integer works

library(SeuratDisk)

r
install.packages("SeuratDisk")

r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")

# Metadata Inspection Without SeuratDisk
r
install.packages("hdf5r")
library(hdf5r)

# Inspect H5AD structure
h5file <- H5File$new("your_file.h5ad", mode = "r")
print(h5file$ls(recursive = TRUE))  # Show all groups/datasets

# Extract metadata (assuming standard AnnData format)
obs <- h5file[["obs"]][]
rownames(obs) <- h5file[["obs/_index"]][]
head(obs)  # View first 5 rows of metadata

# Close file
h5file$close()install.packages("hdf5r")
library(hdf5r)

# Inspect H5AD structure
h5file <- H5File$new("your_file.h5ad", mode = "r")
print(h5file$ls(recursive = TRUE))  # Show all groups/datasets

# Extract metadata (assuming standard AnnData format)
obs <- h5file[["obs"]][]
rownames(obs) <- h5file[["obs/_index"]][]
head(obs)  # View first 5 rows of metadata

# Close file
h5file$close()

save(file = "my_analysis.R", list = c("important_function", "other_object"))save(file = "my_analysis.R", list = c("important_function", "other_object"))

r
system("git add .")
system('git commit -m "Add analysis and setup files"')
system("git push origin main")

system("git init")

system("git add .")  # Add all files to the staging area

system('git commit -m "Initial commit with project setup and scripts"')

system('git remote add origin https://github.com/spexster/Plasma_Biomarkers.git')

system('git push -u origin main')

https://github.com/spexster/Plasma_Biomarkers

git init