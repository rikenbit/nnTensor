Xs <- list(
	toyModel(model="CP"),
	toyModel(model="Tucker"),
	toyModel(model="NMF")
)

for(i in seq_along(Xs)){
	X <- Xs[[i]]
	Ms <- kFoldMaskTensor(X=X)
	lapply(seq_along(Ms), function(x){
		expect_true(identical(dim(Ms[[x]]), dim(X)))
	})
}
