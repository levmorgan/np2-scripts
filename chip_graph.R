setwd("C:/Users/morga/Projects/brent_lab/np2/YEAST")
# Optionally clear out existing variables
rm(list=ls())

# Read in data
chip_benchmark <- read.table("BENCHMARKS/ChIP_netowrk.adjlst", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pwm_benchmark <- read.table("BENCHMARKS/PWM_netowrk.adjlst", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genes <- read.delim("OUTPUT/genes", header = FALSE, stringsAsFactors=FALSE)
regulators <- read.delim("OUTPUT/regulators", header = FALSE, stringsAsFactors=FALSE)
np2_output_base <- read.delim("OUTPUT/netprophet2_network.adjmtr", header = FALSE, sep = "\t", 
			 skipNul = TRUE, stringsAsFactors = FALSE)

benchmark <- chip_benchmark
np2_output <- np2_output_base

# Remove NAs
np2_output <- np2_output[,colSums(is.na(np2_output)) < nrow(np2_output)]

# Add gene names
colnames(np2_output) <- genes[,1]
# Remove genes missing from benchmark data
np2_output <- np2_output[,colnames(np2_output) %in% benchmark[,2]]

# Add regulator names
num_cols <- dim(np2_output)[2] + 1
np2_output[,num_cols] <- regulators[,1]
# Remove regulators missing from benchmark data
np2_output <- np2_output[np2_output[,num_cols] %in% benchmark[,1],]

# Eval whether edges have benchmark data support
np2_long <- melt(np2_output, id.vars = c(num_cols), measure.vars = 1:(num_cols - 1))
colnames(np2_long) <- c("regulators", "genes", "weight")
np2_long <- np2_long[np2_long$weight > 0,]
np2_long <- np2_long[order(np2_long$regulators, np2_long$weight, decreasing=TRUE),]


eval_benchmark_support <- function(benchmark_data, weights) {
	weights_key <- paste(as.character(weights[,1]), as.character(weights[,2]))
	benchmark_key <- paste(as.character(benchmark_data[,1]), as.character(benchmark_data[,2]))
	return(weights_key %in% benchmark_key)
}

np2_long$has_support <- eval_benchmark_support(benchmark, np2_long)

# Find % benchmark data support for networks with 10-100 edges per TF
np2_split <- split(np2_long, np2_long$regulators)

plot_data = c()
for (i in 1:100) {
	support_sums <- Map(function (genes) {
				     sum(genes$has_support[1:i])}, 
				     np2_split)
	support_sums <- as.vector(support_sums, mode="numeric")
	support_sum <- sum(support_sums)
	plot_data[i] <- support_sum/(i*length(np2_split))
}

plot(plot_data)

#create_evaluation_plot(np2_output, genes, regulators, chip_benchmark)
#create_evaluation_plot(np2_output, genes, regulators, pwm_benchmark)
