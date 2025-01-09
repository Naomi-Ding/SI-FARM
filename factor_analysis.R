Library(R.matlab)
require(nFactors)

### settings
args <- commandArgs(TRUE)
# print(args)
n <- eval(parse(text = args[[1]])) # 100 or 200
p <- eval(parse(text = args[[2]])) # p = 3
nv <- eval(parse(text = args[[3]])) # 100 or 200
nsimu <- eval(parse(text = args[[4]]))
ul <- eval(parse(text = args[[5]])) # choose from (0,1,2,3) for different level of heterogeneity
index.g <- eval(parse(text = args[[6]]))
print(paste0("n=", n, ", p=", p, ", nv=", nv, ", nsimu=", nsimu, ", ul=", ul, ", index_g=", index.g))


### load the residual matrix
data <- readMat(sprintf("residuals_n%d_s%d_nsimu%d_g%d_ul%d.mat", n, nv, nsimu, index_g, ul))
all.PR <- data$all.PR
q <- data$q
rm(data)

### consider 4 methods for factor analysis,
### including OC & AF & Parallel & Kaiser respectively
# q.est <- array(0, dim=c(4, nsimu))
q.fa.est <- array(0, dim = c(4, nsimu))
start <- Sys.time()
for (i in 1:nsimu) {
  print(i)
  data <- data.frame(all.PR[, , i])
  eigenvalues <- eigenComputes(x = data)
  aparallel <- eigenBootParallel(x = data, quantile = 0.99, nboot = 100, option = "permutation", model = "factors")$quantile
  aparallel.fa <- parallel(var = length(eigenvalues), rep = 100, quantile = 0.99, model = "factors")$eigen$qevpea
  # results <- nScree(x = eigenvalues, aparallel = aparallel, model = "factors")
  results.fa <- nScree(x = eigenvalues, aparallel = aparallel.fa, model = "factors")
  # q.est[,i] <- as.numeric(results$Components)
  q.fa.est[, i] <- as.numeric(results.fa$Components)
  end <- Sys.time()
  print(end - start)
}
# rownames(q.est) <- colnames(results$Components)
# rowMeans(q.est==q)
rownames(q.fa.est) <- colnames(results.fa$Components)
accuracy <- rowMeans(q.fa.est == q)
print(accuracy)

# ### AF turns out to be the best method
# q.est <- q.fa.est['naf', ]

### save the results
resname <- sprintf("factor_analysis_n%d_s%d_nsimu%d_g%d_ul%d.mat", n, nv, nsimu, index.g, ul)
writeMat(resname, qest = q.fa.est, accuracy = accuracy)

print(paste0("estimated number of latent factors: ", median(q.fa.est)))