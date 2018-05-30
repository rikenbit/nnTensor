recTensor <-
function (S = NULL, A = NULL, idx = 1:3, reverse = FALSE) 
{
    if (is(S)[1] == "Tensor") {
        X_bar <- S
    }
    else if ("numeric" %in% is(S) && length(S) == 1) {
        X_bar <- rand_tensor(c(1, 1, 1))
        X_bar@data[] <- S
    }
    else if ("numeric" %in% is(S) && length(S) != 1) {
        X_bar <- rand_tensor(sapply(A, nrow))
        X_bar@data[] <- 0
        for (x in 1:length(A)) {
            X_bar@data[x, x, x] <- S[x]
        }
    }
    if (!reverse) {
        for (n in idx) {
            X_bar <- ttm(X_bar, t(A[[n]]), m = n)
        }
    }
    else {
        for (n in idx) {
            X_bar <- ttm(X_bar, A[[n]], m = n)
        }
    }
    return(X_bar)
}
