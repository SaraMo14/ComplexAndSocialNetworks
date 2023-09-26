library(igraph)

write("## Running script",stderr())

write("## Creating output directory",stderr())
mainDir <- dirname(rstudioapi::getSourceEditorContext()$path) # where this script is stored
dir.create(file.path(mainDir, "output"),F)
# setwd(file.path(mainDir, "output"))

#################################################################################################################################
# (A) Plot of the clustering coefficient and the average shortest-path as a function of the parameter p of the WS model ----
#################################################################################################################################


write("## Plotting of (A) statred",stderr())


#p <- seq(0.0001, 0.001, 0.01, 0.1, 1 )
p <- 10^(seq(-4,0,0.2))

compute_metrics_A <- function(p){
  transitivity_values <- numeric()
  path_values <- numeric()
  #for each p, iterate 20 times
  for (j in 1:20) {
    ws_graph <- sample_smallworld(1, 1000, 4, p) 
    transitivity_values <-c(transitivity_values, transitivity(ws_graph, type = "undirected", vids = NULL,
                                    weights = NULL)) #computes transitivity
    path_values <- c(path_values,average.path.length(ws_graph)) #computes path
  }
  # Compute the average of transitivity and path values for this p
  avg_trans <- mean(transitivity_values)
  avg_path <- mean(path_values)
  
  return(list(p = p, avg_trans = avg_trans, avg_path = avg_path))
  
}
# Use sapply to iterate over p_values and compute metrics
results <- lapply(p, compute_metrics_A)

# Extract the results into separate vectors
p_results <- sapply(results, function(x) x$p)
avg_trans_results <- sapply(results, function(x) x$avg_trans)
avg_path_results <- sapply(results, function(x) x$avg_path)

#normalize
avg_trans_results <-avg_trans_results/avg_trans_results[1]
avg_path_results <-avg_path_results/avg_path_results[1]

# Plot the results
pdf(file = file.path(getwd(), "output/A.pdf"), width = 6, height = 4)
options(scipen=100)
plot(p,avg_trans_results, ylim = c(0,1), log='x', pch=15, col=2, ylab = "Values")
points(p,avg_path_results, ylim = c(0,1),pch=15, col = 'blue')
#legend
legend( x="bottomleft",  legend=c("C(p)/C(0)", "L(p)/L(0)"), col=c("red", "blue"), pch=15) 
dev.off()



#################################################################################################################################
# (B) Plot of the average shortest-path length as a function of the network size of the ER model ----
#################################################################################################################################

#NB: for the graph to be almost surely connected, we need: p > (1+eps) ln/n


network_sizes <- seq(1, 16001,5000)
eps =  0.05


# Preallocate a numeric vector to store results
avg_shortest_paths <- numeric(length(network_sizes))

# Iterate over network_sizes
for (i in seq_along(network_sizes)) {
  n <- network_sizes[i]
  print(n)
  threshold <- (1 + eps) * log(n) / n
  prob <- threshold + eps 
  avg_shortest_paths[i] <- mean(replicate(20, average.path.length(sample_gnp(n, prob))))
}

# Create a data frame with the results (if needed)
result_b <- data.frame(n = network_sizes, avg = avg_shortest_paths)


#converted to data frame so i can extract as follows:
n_results <- result_b$n
avg_path_results <- result_b$avg
#replace NaN values with 0
avg_path_results[is.nan(avg_path_results)] <- 0


pdf(file = file.path(getwd(), "output/B.pdf"), width = 6,  height = 4)
plot(n_results,avg_path_results, xlab = "n", ylab='Average Shortest Path',  type = "b")

dev.off()