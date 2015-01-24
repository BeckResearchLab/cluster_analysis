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
