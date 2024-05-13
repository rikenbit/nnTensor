kFoldMaskTensor <- function(X, k=3, seeds=123, sym=FALSE){
	# Setting
	set.seed(seeds)
	# Input Type
	checkMatrix <- length(dim(X)) == 2
	if(checkMatrix){
		X <- as.tensor(X)
	}
	# Symmetric-Mode
	if(sym){
		if(!checkMatrix || !isSymmetric(X@data)){
			stop("sym option is available only when matrix is specified")
		}
		# Non-zero in Upper Triangle Parts
		position_nz_ut <- which(X@data[which(upper.tri(X@data))] != 0)
		# NA
		position_na_ut <- which(is.na(X@data[which(upper.tri(X@data))]))
		# Shuffle
		position_nz_ut <- sample(position_nz_ut, length(position_nz_ut))
		index <- suppressWarnings(split(position_nz_ut, seq(k)))
		# Output
		out <- lapply(seq(k), function(x){
			X@data[] <- 1
			X@data[upper.tri(X@data)][index[[x]]] <- 0
			X@data[upper.tri(X@data)][position_na_ut] <- 0
			tX <- t(X@data)
			X@data[lower.tri(X@data)] <- tX[lower.tri(tX)]
			X@data
		})
	# Asymmetric-Mode
	}else{
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
	}
	set.seed(NULL)
	out
}
