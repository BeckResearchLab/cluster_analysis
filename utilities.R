# to be sourced

# return the directory path for a k / cluster pair
dir.k.cluster <- function(dir.output, k, cluster, make.dir = F) {
	dir <- sprintf("%s/k_%d.dir/cluster_%d.dir", dir.output, k, cluster)
	# equivalent to to mkdir -p
	if (make.dir) {
		dir.create(dir, recursive = T, showWarnings = F)
	}
	return(dir)
}

# the sum of member to centroid distances
get.distsum <- function() {
	x <- env$cluster.ensemble
	X <- x@k
	Y <- sapply(x@models, function(z) info(z, "distsum"))
	Z <- c(NA, diff(Y))
	df = data.frame(k=X, distsum=Y, distsum.delta=Z)
	return(df);
}
