load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/realdata/HCP_afterQC_chr21.RData")
load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/realdata/brain_volumes.rda")

X <- as.matrix(genotmp); X[is.na(X)] <- 0; X <- X[rownames(X) %in% d[,1],] 
#X <- cbind(as.numeric(rownames(X)), X);colnames(X)[1] <- "ID"; rownames(X) <- NULL
dat <- as.matrix(d); dat <- dat[dat[,1] %in% rownames(X),]; dat <- dat[,-1]

n <- dim(dat)[1]
nom <- sum(dat^2)/n
Y <- as.matrix(t(apply(dat, 1, function(x) scale(x, center=FALSE, scale = TRUE)))) 
#f <- sqrt(nom/sum(Y[1,]^2))
#Y <- Y*f
Y <- Y*100
Y <- scale(Y, scale = FALSE)


save(Y, file = "/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/realdata/phenotype.RData")

