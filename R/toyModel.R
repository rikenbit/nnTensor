toyModel <-
function (model = "CP", seeds=123)
{
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
        X <- matrix(rpois(100 * 300, lambda=1), nrow = 100, ncol = 300)
        X[1:50, 1:50] <- rpois(50*50, lambda=100)
        X[50:70, 51:100] <- rpois(21*50, lambda=200)
        X[60:65, 151:170] <- rpois(6*20, lambda=150)
        X[25:35, 200:220] <- rpois(11*21, lambda=150)
        X[51:100, 220:300] <- rpois(50*81, lambda=150)
        set.seed(NULL)
    }
    else {
        stop("Please specify model parameter as CP, Tucker, or NMF")
    }
    X
}
