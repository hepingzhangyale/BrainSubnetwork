#!/usr/bin/env Rscript
#SBATCH --mem-per-cpu=4g  --time=2:00:00 --mail-type=ALL --mail-user=yisha.yao@yale.edu


library(plyr)
library(Matrix)
load("/gpfs/ysm/project/zhang_heping/wd278/Yisha/Data/Genotype/chr21/HCP_afterQC_chr21.RData")
load("/gpfs/ysm/project/zhang_heping/wd278/Yisha/Data/Phenotype/brain_volumes.rda")
X <- as.matrix(genotmp); X[is.na(X)] <- 0; X <- X[rownames(X) %in% d[,1],] #X <- cbind(as.numeric(rownames(X)), X);colnames(X)[1] <- "ID"; rownames(X) <- NULL
#Y <- as.matrix(d); Y <- Y[Y[,1] %in% rownames(X),]; Y <- Y[,-1]
X <- scale(X, scale = FALSE)
#Y <- scale(Y, scale = FALSE)
n <- nrow(X)
d <- ncol(X)

load("/gpfs/ysm/project/zhang_heping/wd278/Yisha/Data/Genotype/chr21/HCP_snp_group_chr21.RData")
all(colnames(genotmp) == rownames(BIM_snp_group))  #for chr2, colnames(genotmp)[93039] != rownames(BIM_snp_group); these are also FALSE: chr3, chr4, chr14, chr17
SNP_group <- BIM_snp_group
freq <- as.data.frame(table(SNP_group$group))
colnames(freq) <- c("group", "Freq")
groups <- as.data.frame(unique(SNP_group$group))#the order of the groups are the same as in the original BIM_snp_group; correct order
colnames(groups) <- c("group")
group_freq <- join(groups, freq) #the order of the groups are the same as in the original BIM_snp_group; correct order
q <- dim(group_freq)[1]  #number of groups in X
SNP_group_freq <- join(SNP_group, freq) #contains SNP names, groups belong to, group sizes; correct order of groups
rm(freq, groups)
all(colnames(genotmp) == SNP_group_freq$SNP)
rm(chr_groups, BIM_snp_group, BIMtmp, genotmp)

X_group_sizes <- group_freq$Freq
low <- rep(0, q)
upp <- rep(0, q)
count <- 0
for (i in 1:q){
  low[i] <- count + 1
  upp[i] <- count+ X_group_sizes[i]
  count <- count + X_group_sizes[i]
}
X_grouping <- rep(1:q, X_group_sizes)
rm(i)
#to obtain the PCs of the groups
PC_design <- matrix(0, n, 1)
PC_group_sizes <- rep(0, q)
for (i in 1:q){
  pca <- prcomp(X[, low[i]:upp[i]], scale. = FALSE)
  props <- (pca$sdev)^2/sum((pca$sdev)^2)
  cumprops <- cumsum(props)
  cutoff <- which(cumprops > 0.8)[1] #set a low 0.8 to make sure that top PCs of different groups 
  PC <- pca$x[, 1:cutoff]
  PC_group_sizes[i] <- cutoff
  PC_design <- cbind(PC_design, PC)
}
rm(i, pca, props, cumprops, cutoff, PC)
PC_design <- as.matrix(PC_design[, -1])
PC_grouping <- rep(1:q, PC_group_sizes)
pc_d <- length(PC_grouping)
pc_low <- rep(0, q)
pc_upp <- rep(0, q)
count <- 0
for (i in 1:q){
  pc_low[i] <- count + 1
  pc_upp[i] <- count+ PC_group_sizes[i]
  count <- count + PC_group_sizes[i]
}
rm(i, count)
s_PC_design <- scale(PC_design, scale = TRUE) #should used normalized design since the coefficients are largely affected by the scale of the corresponding covariates 
PC_cor <- t(s_PC_design) %*% s_PC_design/(n-1) - diag(pc_d)

# to combine groups with strong correlations
maxcors <- matrix(0, q, q)
for (i in 1:(q-1)){
  for (l in (i+1):q){
    cor <- PC_cor[pc_low[i]:pc_upp[i], pc_low[l]:pc_upp[l]]
    maxcors[i,l] <- max(max(cor), abs(min(cor)))
  }
}
maxcor_temp <- maxcors
maxc <- max(maxcor_temp)
edges <- matrix(0, 1,3)
while (maxc >0.3 | maxc ==0.3){
  maxind <- which(maxcor_temp==maxc, arr.ind=TRUE)
  edges <- rbind(edges, c(maxind[1,],maxc))
  maxcor_temp[maxind[1,1], maxind[1,2]] <- 0
  maxc <- max(maxcor_temp)
}
edges <- edges[-1,]
colnames(edges) <- c("row", "col", "correlation")
rm(i,l, maxcor_temp, maxc, maxind, cor)
cliques <- list()
cliques[[1]] <- edges[1, 1:2]
m <- 1 #the number of cliques
U <- cliques[[1]]

cut8 <- max(which(edges[,3] >0.8))
for (i in 2:cut8){
  mem <- 0; mem1 <- 0; mem2 <- 0
  if (length(intersect(edges[i,1:2], U)) ==2){
    log1 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1])) !=0)
    mem1 <- which(log1 == TRUE)
    log2 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,2])) !=0)
    mem2 <- which(log2 == TRUE)
    if (mem1 != mem2){
      cliques[[min(mem1, mem2)]]<- union(cliques[[mem1]], cliques[[mem2]])
      cliques[[max(mem1, mem2)]]<- NULL
      m <- m-1
    }
  } else if (length(intersect(edges[i,1:2], U)) ==1){
    logg <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1:2])) !=0)
    mem <- which(logg == TRUE)
    cliques[[mem]] <- union(cliques[[mem]], edges[i,1:2])
    U <- union(U, edges[i,1:2])
  } else {
    m <- m+1
    cliques[[m]] <- edges[i,1:2]
    U <- union(U, edges[i,1:2])
  }
}

cut7 <- max(which(edges[,3] >0.65)) #the 0.6 can be adjusted in each chromosome
for (i in (cut8+1):cut7){
  mem <- 0; mem1 <- 0; mem2 <- 0
  if (length(intersect(edges[i,1:2], U)) ==2){
    log1 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1])) !=0)
    mem1 <- which(log1 == TRUE)
    log2 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,2])) !=0)
    mem2 <- which(log2 == TRUE)
    if (mem1 != mem2){
      meancor <- sum(maxcors[cliques[[mem1]], cliques[[mem2]]])/(length(cliques[[mem1]])*length(cliques[[mem2]]))
      if (meancor > 0.25){#meancor > 0.25 can be adjusted
        cliques[[min(mem1, mem2)]]<- union(cliques[[mem1]], cliques[[mem2]])
        cliques[[max(mem1, mem2)]]<- NULL
        m <- m-1
      }
    }
  } else if (length(intersect(edges[i,1:2], U)) ==1){
    logg <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1:2])) !=0)
    mem <- which(logg == TRUE)
    if (length(intersect(cliques[[mem]], edges[i,1])) == 0){loge <- edges[i,1]}else{loge <- edges[i,2]}
    meancor <- sum(maxcors[loge, cliques[[mem]]])/length(cliques[[mem]])
    if (meancor > 0.25){
      cliques[[mem]] <- union(cliques[[mem]], edges[i,1:2])
      U <- union(U, edges[i,1:2])
    }
  } else {
    m <- m+1
    cliques[[m]] <- edges[i,1:2]
    U <- union(U, edges[i,1:2])
  }
}

cut5 <- max(which(edges[,3] >0.5))
for (i in (cut7+1):cut5){
  mem <- 0; mem1 <- 0; mem2 <- 0
  if (length(intersect(edges[i,1:2], U)) ==2){
    log1 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1])) !=0)
    mem1 <- which(log1 == TRUE)
    log2 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,2])) !=0)
    mem2 <- which(log2 == TRUE)
    if (mem1 != mem2){
      meancor <- sum(maxcors[cliques[[mem1]], cliques[[mem2]]])/(length(cliques[[mem1]])*length(cliques[[mem2]]))
      if (meancor > 0.3){#meancor > 0.3 can be adjusted
        cliques[[min(mem1, mem2)]]<- union(cliques[[mem1]], cliques[[mem2]])
        cliques[[max(mem1, mem2)]]<- NULL
        m <- m-1
      }
    }
  } else if (length(intersect(edges[i,1:2], U)) ==1){
    logg <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1:2])) !=0)
    mem <- which(logg == TRUE)
    if (length(intersect(cliques[[mem]], edges[i,1])) == 0){loge <- edges[i,1]}else{loge <- edges[i,2]}
    meancor <- sum(maxcors[loge, cliques[[mem]]])/length(cliques[[mem]])
    if (meancor > 0.3){
      cliques[[mem]] <- union(cliques[[mem]], edges[i,1:2])
      U <- union(U, edges[i,1:2])
    }
  } else {#m <- m+1;cliques[[m]] <- edges[i,1:2];U <- union(U, edges[i,1:2])
    }
}
rm(loge, logg, mem, log1, mem1, log2, mem2, meancor)
#to record the old-new grouping correspondence
o_cliques <- sapply(seq_along(cliques), function(k) unname(sort(cliques[[k]])))
first_ele <- rep(0, length(o_cliques))
first_ele <- sapply(seq_along(o_cliques), function(k) first_ele[k] <- o_cliques[[k]][1])
o_first_ele <- sort(first_ele, index.return=T)
oo_cliques <- unname(o_cliques[o_first_ele$ix])
rm(o_cliques, first_ele, o_first_ele)
reorganized <- unlist(oo_cliques)
remain <- c(1:q)[is.na(pmatch(c(1:q),reorganized))]
n_q <- length(remain) + length(oo_cliques) #this is the new group size
n_X_group_sizes <- rep(0, n_q)
old_new_groups <- matrix(0,length(remain),2)
reordering_ind <- rep(0, d)
n_low <- rep(0, n_q)
n_upp <- rep(0, n_q)
count <- 0
for (l in 1:length(remain)){
  n_X_group_sizes[l] <- X_group_sizes[remain[l]]
  old_new_groups[l,] <- c(remain[l], l)
  n_low[l] <- count + 1
  n_upp[l] <- count + X_group_sizes[remain[l]]
  count <- n_upp[l]
  reordering_ind[n_low[l]:n_upp[l]] <- c(low[remain[l]]:upp[remain[l]]) 
}
reordering_ind <- reordering_ind[-which(reordering_ind==0)]
for (i in 1:length(oo_cliques)){
  n_X_group_sizes[i+length(remain)] <- sum(X_group_sizes[oo_cliques[[i]]])
  block <- cbind(oo_cliques[[i]], rep((i+length(remain)), length(oo_cliques[[i]])))
  old_new_groups <- rbind(old_new_groups, block)
  n_low[i+length(remain)] <- count + 1
  n_upp[i+length(remain)] <- count + n_X_group_sizes[i+length(remain)]
  count <- count + n_X_group_sizes[i+length(remain)]
  v <- oo_cliques[[i]]
  for (k in 1:length(v)){
    reordering_ind <- c(reordering_ind, low[v[k]]:upp[v[k]])
  }
}
rm(i, l, k, block, v, remain, reorganized, count, U)
colnames(old_new_groups) <- c("old_group_label", "new_group_label")

#the group-combined new X with groups switched
n_X <- X[, reordering_ind]
n_X_grouping <- rep(1:n_q, n_X_group_sizes)
#to obtain the new PCs
n_PC_design <- matrix(0, n, 1)
n_PC_group_sizes <- rep(0, n_q)
for (i in 1:n_q){
  pca <- prcomp(n_X[, n_low[i]:n_upp[i]], scale. = FALSE)
  props <- (pca$sdev)^2/sum((pca$sdev)^2)
  cumprops <- cumsum(props)
  cutoff <- which(cumprops > 0.8)[1] #set a low 0.8 to make sure that top PCs of different groups 
  PC <- pca$x[, 1:cutoff]
  n_PC_group_sizes[i] <- cutoff
  n_PC_design <- cbind(n_PC_design, PC)
}
rm(i, pca, props, cumprops, cutoff, PC)
n_PC_design <- as.matrix(n_PC_design[, -1])
n_PC_grouping <- rep(1:n_q, n_PC_group_sizes)
n_pc_d <- length(n_PC_grouping)
n_pc_low <- rep(0, n_q)
n_pc_upp <- rep(0, n_q)
count <- 0
for (i in 1:n_q){
  n_pc_low[i] <- count + 1
  n_pc_upp[i] <- count+ n_PC_group_sizes[i]
  count <- count + n_PC_group_sizes[i]
}
rm(i, count, reordering_ind)
n_s_PC_design <- scale(n_PC_design, scale = TRUE) 
n_PC_cor <- t(n_s_PC_design) %*% n_s_PC_design/(n-1) - diag(n_pc_d)
length(which(n_PC_cor > 0.6, arr.ind = TRUE)); length(which(n_PC_cor > 0.8, arr.ind = TRUE))

save(n_X, n_X_group_sizes, old_new_groups, n_PC_group_sizes, 
     n_PC_grouping, n_pc_d, n_pc_low, n_pc_upp, n_s_PC_design, n, q, n_q,
     file = "/gpfs/ysm/project/zhang_heping/yy634/cluster_gene_association/real_data/chr21/chr21_design.RData")

