S <- rand_tensor(c(10,20,30))
A <- list(
	A1=matrix(runif(4*10), nrow=10, ncol=4),
	A2=matrix(runif(5*20), nrow=20, ncol=5),
	A3=matrix(runif(6*30), nrow=30, ncol=6))

t_A <- list(
	A1=matrix(runif(10*4), nrow=4, ncol=10),
	A2=matrix(runif(20*5), nrow=5, ncol=20),
	A3=matrix(runif(30*6), nrow=6, ncol=30))

out1 <- recTensor(S, A)
out2 <- recTensor(S, t_A, reverse=TRUE)

expect_equivalent(dim(out1), c(4,5,6))
expect_equivalent(dim(out2), c(4,5,6))
