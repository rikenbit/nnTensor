plotTensor3D <-
function (X = NULL) 
{
    if (is(X)[1] != "Tensor") {
        stop("Please specify the Tensor\n")
    }
    coordinate <- expand.grid(1:dim(X)[1], 1:dim(X)[2], 1:dim(X)[3])
    colnames(coordinate) <- c("A1", "A2", "A3")
    allpoints <- as.vector(cs_unfold(X, m = 3)@data)
    usingpoints <- which(abs((allpoints - mean(allpoints))/sd(allpoints)) > 
        2)
    value <- smoothPalette(-allpoints, pal = "RdBu")
    if (length(usingpoints) != 0) {
        scatter3D(coordinate[usingpoints, 1], coordinate[usingpoints, 
            2], coordinate[usingpoints, 3], bty = "g", col = "blue", 
            colkey = FALSE)
    }
    else {
        scatter3D(coordinate[, 1], coordinate[, 2], coordinate[, 
            3], bty = "g", col = "blue", colkey = FALSE)
    }
}
