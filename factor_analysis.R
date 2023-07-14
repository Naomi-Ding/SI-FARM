# install.packages("R.matlab")
library(R.matlab)
##### Real data #####
# data <- readMat('realdata_p6_residual.mat')
# data <- readMat('D:/OneDrive - Florida State Students/PhD/paper_2/real_data/COVID19_data/est_result_residual.mat')
data <- readMat('realdata_p6_sivc_ini_residual.mat')
R = data$R
PR = data$PR
R.smooth = data$R.smooth
PR.smooth = data$PR.smooth
library(lattice)
library(nFactors)
data = data.frame(R.smooth)
# data = data.frame(R)
# data = data.frame(PR)
# data = data.frame(PR.smooth)
eigenvalues <- eigenComputes(x = data)
aparallel <- eigenBootParallel(x = data, quantile=0.99, nboot = 100, option = "permutation", model = "factors")$quantile
aparallel.fa <- parallel(var = length(eigenvalues), rep = 100, quantile = 0.99, model = "factors")$eigen$qevpea
results <- nScree(x = eigenvalues, aparallel = aparallel, model = "factors")
results.fa <- nScree(x = eigenvalues, aparallel = aparallel.fa, model = "factors")
p.oc <- results$Components$noc
p.af <- results$Components$naf
p.parallel <- results$Components$nparallel
p.kaiser <- results$Components$nkaiser
p.oc.fa <- results.fa$Components$noc
p.af.fa <- results.fa$Components$naf
p.parallel.fa <- results.fa$Components$nparallel
p.kaiser.fa <- results.fa$Components$nkaiser
#####

# data <- readMat('residual_matrix_g5_u0.mat')
# # str(data)
# # head(data)
# all.R <- data$all.R
# all.U <- data$all.U
# all.PR <- data$all.PR


# install.packages('cate')
library(cate)

p.ed <- 0
p.bcv <- 0
for (i in 1:100)
{
  p.ed[i] <- est.factor.num(all.PR[,,i], method = "ed")
  #p.bcv[i] <- est.factor.num(all.R[,,i], method = "bcv")
}
print(p.ed)
print(p.bcv)


# install.packages('jackstraw')
# install.packages('ClusterR')
library(ClusterR)
library(jackstraw)
# install.packages("devtools")
library("devtools")
# install_github("ncchung/jackstraw")


# install.packages("paran")
library(MASS)
library(paran)
p.pc <- 0
for (i in 1:100)
{
  print(i)
  tmp.R <- paran(all.R[,,i], centile = 95, cfa = TRUE)
  p.pc[i] <- tmp.R$Retained
}
print(p.pc)


# install.packages("nFactors")
library(lattice)
library(nFactors)
p.oc <- 0
p.af <- 0
p.parallel <- 0
p.kaiser <- 0
p.oc.fa <- 0
p.af.fa <- 0
p.parallel.fa <- 0
p.kaiser.fa <- 0
for (i in 1:100)
{
  print(i)
  # data = data.frame(all.R[,,i])
  data = data.frame(all.PR[,,i])
  eigenvalues <- eigenComputes(x = data)
  aparallel <- eigenBootParallel(x = data, quantile=0.99, nboot = 100, option = "permutation", model = "factors")$quantile
  aparallel.fa <- parallel(var = length(eigenvalues), rep = 100, quantile = 0.99, model = "factors")$eigen$qevpea
  results <- nScree(x = eigenvalues, aparallel = aparallel, model = "factors")
  results.fa <- nScree(x = eigenvalues, aparallel = aparallel.fa, model = "factors")
  p.oc[i] <- results$Components$noc
  p.af[i] <- results$Components$naf
  p.parallel[i] <- results$Components$nparallel
  p.kaiser[i] <- results$Components$nkaiser
  p.oc.fa[i] <- results.fa$Components$noc
  p.af.fa[i] <- results.fa$Components$naf
  p.parallel.fa[i] <- results.fa$Components$nparallel
  p.kaiser.fa[i] <- results.fa$Components$nkaiser
}

print(p.oc)
print(p.oc.fa)
mean(p.oc == 1)
mean(p.oc.fa == 1)
print(p.af)
print(p.af.fa)
mean(p.af == 1)
mean(p.af.fa == 1)
print(p.parallel)
print(p.parallel.fa)
mean(p.parallel == 1)
mean(p.parallel.fa == 1)
print(p.kaiser)
print(p.kaiser.fa)
mean(p.kaiser == 1)
mean(p.kaiser.fa == 1)