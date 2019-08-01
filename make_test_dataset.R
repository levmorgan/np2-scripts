# Create a small test dataset

# Optionally clear out existing variables
rm(list=ls())


N_REGULATORS <- 100
N_GENES <- 300

options(stringsAsFactors = FALSE)

library(seqinr)
setwd("C:/Users/morga/Projects/brent_lab/np2/YEAST/RESOURCES")

# Read in data
conditions <- read.table("conditions")
data.expr <- read.table("data.expr", sep="\t")
data.fc.tsv <- read.table("data.fc.tsv", sep="\t", header=TRUE)
genes <- read.table("genes")
promoter.fasta <- read.fasta("promoter.fasta")
regulators <- read.table("regulators")
signed.de.adj <- read.table("signed.de.adj", sep="\t")

# Cut down data
test_regs_i <- sample(seq_along(regulators[,1]), N_REGULATORS)
test_regs_i <- sort(test_regs_i)
test_regs <- regulators[test_regs_i,]

# Our sample genes will be the regulators and N_GENES other genes
test_genes <- c(test_regs, sample(setdiff(genes[,1], test_regs), N_GENES))
test_genes_i <- sort(match(test_genes, genes[,1]))
test_genes <- genes[test_genes_i,]

test.data.expr <- data.expr[test_genes_i,]
test.data.fc.tsv <- data.fc.tsv[,test_genes_i]
test.promoter.fasta <- promoter.fasta[test_genes]
test.signed.de.adj <- signed.de.adj[test_regs_i,test_genes_i]

# Write data back out
dir.create("TEST_RESOURCES")
write.table(data.frame(conditions), file="TEST_RESOURCES/conditions", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(test.data.expr, file="TEST_RESOURCES/data.expr", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(test.data.fc.tsv, file="TEST_RESOURCES/data.fc.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=TRUE)
write.table(data.frame(test_genes), file="TEST_RESOURCES/genes", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(data.frame(test_regs), file="TEST_RESOURCES/regulators", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(test.signed.de.adj, file="TEST_RESOURCES/signed.de.adj", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.fasta(test.promoter.fasta, names(test.promoter.fasta), "TEST_RESOURCES/promoter.fasta")
