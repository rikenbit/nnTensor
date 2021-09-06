plotTensor2D <- function(X = NULL, method=c("sd", "mad"),
    sign=c("positive", "negative", "both"), thr=2){
    # Argument check
    method <- match.arg(method)
    sign <- match.arg(sign)
    # Array -> Tensor
    if (is.array(X)){
        X <- as.tensor(X)
    }
    # Setting
    coordinate <- expand.grid(1:dim(X)[1], 1:dim(X)[2])
    allpoints <- vec(X)
    # Color value
    # SD
    if(method == "sd"){
        sdval <- (allpoints - mean(allpoints)) / sd(allpoints)
        if(sign == "positive"){
            usingpoints <- which(sdval >= thr)
        }
        if(sign == "negative"){
            usingpoints <- which(sdval <= - thr)
        }
        if(sign == "both"){
            usingpoints <- which(abs(sdval) >= thr)
        }
    }
    # MAD
    if(method == "mad"){
        madval <- median(abs(allpoints - median(allpoints)))
        madnorm <- (allpoints - median(allpoints)) / madval
        if(sign == "positive"){
            usingpoints <- which(madnorm >= thr)
        }
        if(sign == "negative"){
            usingpoints <- which(madnorm <= - thr)
        }
        if(sign == "both"){
            usingpoints <- which(abs(madnorm) >= thr)
        }
    }
    # Color schema
    col <- rep(0, length=length(allpoints))
    col[usingpoints] <- 6
    # Plot
    scatter2D(coordinate[,1], coordinate[,2],
        col=col, bty = "g", colkey = FALSE, pch=15)
}