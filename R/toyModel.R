toyModel <- function(model = "CP", seeds=123){
    if (model == "CP") {
        set.seed(seeds)
        X <- rand_tensor(c(30, 30, 30))
        X@data[,,] <- rpois(30*30*30, lambda=1)

        X[1:5, 1:5, 1:5] <- rpois(5*5*5, lambda=100)
        X[10:14, 10:14, 10:14] <- rpois(5*5*5, lambda=100)
        X[18:22, 18:22, 18:22] <- rpois(5*5*5, lambda=100)
        X[26:30, 26:30, 26:30] <- rpois(5*5*5, lambda=100)
        set.seed(NULL)
    }
    else if (model == "Tucker") {
        set.seed(seeds)
        X <- rand_tensor(c(50, 50, 50))
        X@data[,,] <- rpois(50*50*50, lambda=1)

        X[1:5, 1:5, 1:50] <- rpois(5*5*50, lambda=100)
        X[46:50, 46:50, 1:50] <- rpois(5*5*50, lambda=100)
        X[1:50, 16:20, 1:5] <- rpois(50*5*5, lambda=100)
        X[1:50, 30:35, 46:50] <- rpois(50*6*5, lambda=100)
        X[10:15, 1:50, 16:20] <- rpois(6*50*5, lambda=100)
        X[30:35, 1:50, 30:35] <- rpois(6*50*6, lambda=100)
        set.seed(NULL)
    }
    else if (model == "NMF") {
        set.seed(seeds)
        X <- matrix(rpois(100 * 300, lambda=1),
            nrow = 100, ncol = 300)

        X[1:50, 1:50] <- rpois(50*50, lambda=100)
        X[50:70, 51:100] <- rpois(21*50, lambda=200)
        X[60:65, 151:170] <- rpois(6*20, lambda=150)
        X[25:35, 200:220] <- rpois(11*21, lambda=150)
        X[51:100, 220:300] <- rpois(50*81, lambda=150)
        set.seed(NULL)
    }
    else if (model == "siNMF_Easy") {
        set.seed(seeds)
        X1 <- matrix(rpois(100 * 300, lambda=1),
            nrow = 100, ncol = 300)
        X2 <- matrix(rpois(100 * 200, lambda=1),
            nrow = 100, ncol = 200)
        X3 <- matrix(rpois(100 * 150, lambda=1),
            nrow = 100, ncol = 150)

        X1[1:30, 1:90] <- rpois(30*90, lambda=100)
        X1[31:60, 91:180] <- rpois(30*90, lambda=100)
        X1[61:90, 181:270] <- rpois(30*90, lambda=100)

        X2[1:30, 61:120] <- rpois(30*60, lambda=100)
        X2[31:60, 121:180] <- rpois(30*60, lambda=100)
        X2[61:90, 1:60] <- rpois(30*60, lambda=100)

        X3[1:30, 81:120] <- rpois(30*40, lambda=100)
        X3[31:60, 1:40] <- rpois(30*40, lambda=100)
        X3[61:90, 41:80] <- rpois(30*40, lambda=100)

        X <- list(X1=X1, X2=X2, X3=X3)
        set.seed(NULL)
    }
    else if (model == "siNMF_Hard") {
        set.seed(seeds)
        X1 <- matrix(rpois(100 * 300, lambda=1),
            nrow = 100, ncol = 300)
        X2 <- matrix(rpois(100 * 200, lambda=1),
            nrow = 100, ncol = 200)
        X3 <- matrix(rpois(100 * 150, lambda=1),
            nrow = 100, ncol = 150)

        X1[1:30, 1:90] <- rpois(30*90, lambda=100)
        X1[31:60, 91:180] <- rpois(30*90, lambda=100)
        X1[61:90, 181:270] <- rpois(30*90, lambda=100)

        X1[61:75, 1:90] <- rpois(15*90, lambda=50)
        X1[1:15, 181:270] <- rpois(15*90, lambda=50)
        X1[31:45, 181:270] <- rpois(15*90, lambda=50)

        X2[1:30, 61:120] <- rpois(30*60, lambda=100)
        X2[31:60, 121:180] <- rpois(30*60, lambda=100)
        X2[61:90, 1:60] <- rpois(30*60, lambda=100)

        X2[76:90, 61:120] <- rpois(15*60, lambda=50)
        X2[16:30, 121:180] <- rpois(15*60, lambda=50)

        X3[1:30, 81:120] <- rpois(30*40, lambda=100)
        X3[31:60, 1:40] <- rpois(30*40, lambda=100)
        X3[61:90, 41:80] <- rpois(30*40, lambda=100)

        X3[1:15, 1:40] <- rpois(15*40, lambda=50)
        X3[76:90, 1:40] <- rpois(15*40, lambda=50)
        X3[16:30, 41:80] <- rpois(15*40, lambda=50)
        X3[46:60, 41:80] <- rpois(15*40, lambda=50)
        X3[61:75, 81:120] <- rpois(15*40, lambda=50)

        X <- list(X1=X1, X2=X2, X3=X3)
        set.seed(NULL)
    }
    else {
        stop("Please specify model parameter as CP, Tucker, or NMF")
    }
    X
}