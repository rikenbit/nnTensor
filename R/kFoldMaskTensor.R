kFoldMaskTensor <- function(X, k=3, seeds=123){
	# Setting
	set.seed(seeds)
	# Input Type
	checkMatrix <- length(dim(X)) == 2
	if(checkMatrix){
		X <- as.tensor(X)
	}
	# Non-zero
	position_nz <- which(X@data != 0)
	# NA
	position_na <- which(is.na(X@data))
	# Shuffle
	position_nz <- sample(position_nz, length(position_nz))
	index <- suppressWarnings(split(position_nz, seq(k)))
	# Output
	out <- lapply(seq(k), function(x){
		X@data[] <- 1
		X@data[index[[x]]] <- 0
		X@data[position_na] <- 0
		if(checkMatrix){
			X <- X@data
		}
		X
	})
	set.seed(NULL)
	out
}
