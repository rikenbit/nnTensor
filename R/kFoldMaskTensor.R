kFoldMaskTensor <- function(X, k=5, avoid.zero=TRUE, seeds=123){
	# Setting
	set.seed(seeds)
	if(length(dim(X)) == 2){
		X <- as.tensor(X)
	}
	if(avoid.zero){
		target <- which(X@data != 0)
	}else{
		target <- seq(length(X@data))
	}
	# Shuffle
	target <- sample(target, length(target))
	index <- suppressWarnings(split(target, seq(k)))
	# Output
	out <- lapply(seq(k), function(x){
		X@data[] <- 1
		X@data[index[[x]]] <- 0
		X
	})
	set.seed(NULL)
	out
}
