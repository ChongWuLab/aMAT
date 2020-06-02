# aMAT

aMAT is the software for our proposed method, termed adaptive Multi-trait Association Test (aMAT). It is computationally efficient and only takes a few hours to get the whole genome results in a single core. In other words, it takes a few minutes if we use 100 or more cores.

We only need to prepare the following two inputs:

* GWAS summary results (SNPs by traits), i.e., each row is the Z scores for a particular SNP of p traits of interest. The order of traits should be the same as that in the trait correlation matrix.
* Trait correlation matrix: We and others recommend using [LD score regression](https://github.com/bulik/ldsc) to  calculate the correlation between traits. LD score regression takes many biases into account and yield a more accurate estimation oon the trait correlation matrix.



## Examples

We use a toy example for illustration.

* Generate the data

```R
ar1_cor <- function(n, rho) {
exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
    (1:n - 1))
rho^exponent
}

# generate trait correlation matrix
n.trait = 50
trait.cor = ar1_cor(n.trait,0.7)

# generate Z scores
tmp <- svd(trait.cor)
eigenvalue <- tmp$d
D <- diag(eigenvalue)
CovSsqrt <- tmp$u %*% sqrt(D) %*% t(tmp$v)

n.snp = 10000

Z = matrix(NA,n.snp,n.trait)
for(i in 1:n.snp) {
  Z[i,] <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
}

```

* Run aMAT for multi-trait association analysis

```R
source("multitrait_test_support.R")
start.time = proc.time()[3]
res = aMAT(Z,trait.cor)
head(res)

#running time
proc.time()[3] - start.time #10 seconds in a single core

# empricial type 1 error rate is well controlled
colSums(res<0.05)/dim(res)[1]
```

