#!/usr/bin/env Rscript
#SBATCH --mem-per-cpu=4g  --time=24:00:00 --mail-type=ALL --mail-user=yisha.yao@yale.edu


#the input of sort_grouping is a dataframe with two variables, "index" and "label".
sort_grouping <- function(df){
  attach(df)
  temp <- df[order(label, index),]
  labels <- sort(unique(temp$label))
  num <- length(labels)
  L <- list()
  for (i in 1:num){
    ind <- which(temp$label == labels[i])
    L[[i]] <- temp$index[ind]
  }
  L <- L[order(sapply(L, function(x) x[1], simplify=TRUE))]
  group_sizes <- rep(0, num)
  for(l in 1:num){
    group_sizes[l] <- length(L[[l]])
  }
  ordered_index <- unlist(L)
  new_labels <- rep(1:num, group_sizes)
  new_df <- as.data.frame(cbind(ordered_index, new_labels))
  colnames(new_df) <- c("index", "label")
  return(new_df)
}


TPR <- function(tr_sparsity, hat_sparsity){
  real_pos <- length(which(tr_sparsity != 0))
  detected_pos <- length(which(tr_sparsity !=0 & hat_sparsity != 0))
  rate <- detected_pos/real_pos
  return(rate)
}
TNR <- function(tr_sparsity, hat_sparsity){
  real_neg <- length(which(tr_sparsity == 0))
  detected_neg <- length(which(tr_sparsity ==0 & hat_sparsity == 0))
  rate <- detected_neg/real_neg
  return(rate)
}
ACC <- function(tr_sparsity, hat_sparsity){
  detected_pos <- length(which(tr_sparsity !=0 & hat_sparsity != 0))
  detected_neg <- length(which(tr_sparsity ==0 & hat_sparsity == 0))
  total <- dim(as.matrix(tr_sparsity))[1] * dim(as.matrix(tr_sparsity))[2]
  rate <- (detected_pos+detected_neg)/total
  return(rate)
}

comparetuning <- matrix(0, 9, 5)

load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha8_lamthred2800.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
data1 <- matrix(0, T, 5)
for (t in 1:T){
  data1[t,1] <- sum((Beta-Betas[[t]])^2)/(d*p)
  data1[t,2] <- ACC(Beta_labels, Blabels[[t]])
  data1[t,3] <- length(which(true_labels != cluster_labels[[t]]))/(q*p) 
  data1[t,4] <- length(which(Beta_cluster_num != cluster_nums[[t]]))/q
  data1[t,5] <- sum((true_means - means[[t]])^2)/(d*p)
}
comparetuning[1,] <- data1[T,]
#save(data2, file = "/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation2/setting2_results.RData")


load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha8_lamthred3000.RData")
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
data1 <- matrix(0, T, 5)
for (t in 1:T){
  data1[t,1] <- sum((Beta-Betas[[t]])^2)/(d*p)
  data1[t,2] <- ACC(Beta_labels, Blabels[[t]])
  data1[t,3] <- length(which(true_labels != cluster_labels[[t]]))/(q*p) 
  data1[t,4] <- length(which(Beta_cluster_num != cluster_nums[[t]]))/q
  data1[t,5] <- sum((true_means - means[[t]])^2)/(d*p)
}
comparetuning[2,] <- data2[nr,]


load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha8_lamthred3200.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data3 <- matrix(0, nr, 5)
for (t in 1:nr){
  data3[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data3[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data3[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data3[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data3[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[3,] <- data3[nr,]


load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha12_lamthred2800.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data4 <- matrix(0, nr, 5)
for (t in 1:nr){
  data4[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data4[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data4[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data4[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data4[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[4,] <- data4[nr,]



load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha12_lamthred3000.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data5 <- matrix(0, nr, 5)
for (t in 1:nr){
  data5[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data5[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data5[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data5[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data5[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[5,] <- data5[nr,]



load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha12_lamthred3200.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data6 <- matrix(0, nr, 5)
for (t in 1:nr){
  data6[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data6[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data6[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data6[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data6[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[6,] <- data6[nr,]



load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha16_lamthred2800.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data7 <- matrix(0, nr, 5)
for (t in 1:nr){
  data7[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data7[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data7[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data7[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data7[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[7,] <- data7[nr,]


load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha16_lamthred3000.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data8 <- matrix(0, nr, 5)
for (t in 1:nr){
  data8[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data8[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data8[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data8[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data8[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[8,] <- data8[nr,]


load("/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/setting6_alpha16_lamthred3200.RData")
d <- n_pc_d;  q <- n_q
#the true labels are stored in Beta_labels; the estimated labels are in classifications.
true_labels <- matrix(0, q, p)
for (i in 1:q){
  index4 <- which(Beta_labels[i,] !=0)
  if (length(index4)==0){}else{
    index3 <- which(Beta_labels[i,] ==0)
    part3 <- as.data.frame(cbind(index3, 0))
    colnames(part3) <- c("index", "label")
    part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
    colnames(part4_temp) <- c("index", "label")
    part4 <- sort_grouping(part4_temp)
    colnames(part4) <- c("index", "label")
    tr <- rbind(part3, part4)
    tr <- tr[order(tr$index),]
    true_labels[i,] <- as.vector(tr$label)
  }
}
rm(i, index3, index4, part3, part4, part4_temp, tr)
#to standardize/reorder classifications (the estimated labels of Betas)
cluster_labels <- list()
for (l in 1:length(classifications)){
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[l]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  cluster_labels[[l]] <- labels_temp
  rm(i, index1, index2, part1, part2, part2_temp, est)
}
#group_sizes, Beta, Beta_labels, Beta_cluster_num, t_means are the truths. Beta_labels records the true classifications
T <- length(Betas)
nr <- min(T, 10)
data9 <- matrix(0, nr, 5)
for (t in 1:nr){
  data9[t,1] <- sum((Beta-Betas[[T-nr+t]])^2)/(d*p)
  data9[t,2] <- ACC(Beta_labels, Blabels[[T-nr+t]])
  data9[t,3] <- length(which(true_labels != cluster_labels[[T-nr+t]]))/(q*p) 
  data9[t,4] <- length(which(Beta_cluster_num != cluster_nums[[T-nr+t]]))/q
  data9[t,5] <- sum((true_means - means[[T-nr+t]])^2)/(d*p)
}
comparetuning[9,] <- data9[nr,]

save(comparetuning, file = "/gpfs/gibbs/project/zhang_heping/yy634/cluster_gene_association/submit2/simulation/tuning_setting6.RData")


#sum((Beta-Betas[[6]])^2)/(d*p)   #MSE
#ACC(Beta_labels, Blabels[[6]])   #accuracy
#length(which(true_labels != cluster_labels[[6]]))/(q*p)   #classification error rate
#length(which(Beta_cluster_num != cluster_nums[[6]]))/q    #cluster number error rate
#sum((true_means - means[[6]])^2)/(d*p)   #gaussian center MSE   
