# put all the score_chr.xlsx and clustering_chr.xlsx files and Yeo_DK_key.csv,  yeo7_NetworkNames.csv in the same folder.
# replace "/path of results/" by the path to this folder 
rm(list = ls())
library(readxl)
library(dplyr)

path = "/path of results/"
dk_to_yeo7_key = read.csv(paste0(path, "Yeo_DK_key.csv"), row.names = 1)
yeo7 = read.csv(paste0(path, "yeo7_NetworkNames.csv"))
colnames(yeo7) = c("yeo_krienen", "NetworkNames")
dk_to_yeo7_key = merge(dk_to_yeo7_key, yeo7, by = "yeo_krienen")
rm(yeo7)


setwd(path)

yeo7_overlap <- list()
for(chr in 1:22){
  score <- read_xlsx(paste0("score_chr", chr, ".xlsx"))
  brain_names = colnames(score)[c(-1, -2, -3)]
  brain_names = paste0(sapply(sapply(brain_names, function(x) strsplit(x, "[_]")[[1]][1]), function(xx) ifelse(xx == "L", "lh", "rh")), "_",
         sapply(brain_names, function(x) strsplit(x, "[_]")[[1]][2]))
  
  clustering <- read_xlsx(paste0("clustering_chr", chr, ".xlsx"))
  clustering <- clustering[, c("new_group_label", "cluster0")]
  clustering$gene = sapply(clustering$new_group_label, function(x){
    tmp = score$group[score$new_group_label == x]
    tmp = tmp[grep("ENSG", tmp)]
    paste(tmp, collapse = ",")
  })
  clustering <- clustering %>% filter(!is.na(cluster0) & (gene != "") & !grepl("huge", cluster0))
  yeo7_overlap_chr <- list()
  for(j in 1:nrow(clustering)){
    cc = strsplit(clustering$cluster0[j], "[,]")[[1]]
    cc = lapply(cc, function(x){
      as.numeric(strsplit(x, "[-]")[[1]])
    })
    tmp = dk_to_yeo7_key
    for(c in 1:length(cc)){
      v = rep(0, nrow(tmp))
      v[sapply(brain_names[cc[[c]]], function(x) which(tmp$roi == x))] = 1
      tmp = cbind(tmp, v)
      colnames(tmp)[ncol(tmp)] = paste0("cluster", c)
    }
    d <- tmp %>% group_by(NetworkNames) %>%
      summarise(across(starts_with("cluster"), mean))
    d2 <- tmp %>% group_by(NetworkNames) %>%
      summarise(across(starts_with("cluster"), sum))
    s = apply(d2[, -1], 2, sum)
    d <- cbind(d, t(apply(d2[, -1], 1, function(x) x/s)))
    colnames(d)[(ncol(d)-length(cc)+1) : ncol(d)] = paste0("across_cluster", 1:length(cc))
    d$gene = clustering$gene[j]
    d$new_group_label = clustering$new_group_label[j]
    d$chr = chr
    n_all = sapply(cc, function(x) paste(brain_names[x], collapse = ","))
    for(c in 1:length(cc)){
      d = cbind(d, n_all[c])
      colnames(d)[ncol(d)] = paste0("rois_cluster", c)
    }
    yeo7_overlap_chr[[j]] = d
  }
  yeo7_overlap[[chr]] = yeo7_overlap_chr
}

yeo7_overlap_df = NULL
for(chr in 1:22){
  yeo7_overlap_chr = yeo7_overlap[[chr]]
  tmp = lapply(yeo7_overlap_chr, function(x){
    xtmp <- x %>% group_by(NetworkNames, gene, new_group_label, chr) %>%
      mutate(max_prop = max(c_across(starts_with("cluster"))),
             across_max_prop = max(c_across(starts_with("across"))),
             cluster_id = which.max(c_across(starts_with("cluster"))),
             across_cluster_id = which.max(c_across(starts_with("across"))))
    xx = xtmp[, c("NetworkNames", "gene", "new_group_label", "chr", "max_prop", "cluster_id", "across_max_prop", "across_cluster_id")]
    xx$rois = unlist(xtmp[1, paste0("rois_cluster", xtmp$cluster_id)])
    return(xx)
  })
  yeo7_overlap_chr = do.call("rbind", tmp)
  yeo7_overlap_df = rbind(yeo7_overlap_df, yeo7_overlap_chr)
}
yeo7_overlap_df <- yeo7_overlap_df %>% arrange(desc(max_prop), desc(across_max_prop))


region<- read.csv("Yeo_DK_key.csv")

region$network <- NA
region$network[region$yeo_krienen==1]= "Visual"
region$network[region$yeo_krienen==2]= "Somatomotor"
region$network[region$yeo_krienen==3]= "Dorsal Attention"
region$network[region$yeo_krienen==4]= "Salience / Ventral Attention"
region$network[region$yeo_krienen==5]= "Limbic"
region$network[region$yeo_krienen==6]= "Control"
region$network[region$yeo_krienen==7]= "Default"

num_region <- c(4, 20,3, 13, 8, 8, 12 )

ii=0
cluster_summary_all <- NULL
for(chr in 1:22)
{
  chr_list <- yeo7_overlap[[chr]]
  print(length(chr_list))
  for(j in 1:length(chr_list))
  {
    cluster_list <- chr_list[[j]]
    cluster_num <- length(grep("across_cluster", colnames(cluster_list)))
    
    for( k in 1:cluster_num)
    {
      prop_cluster<-cluster_list[, paste0("cluster",k)]
      num_region_cluster <- prop_cluster*num_region
      cluster_list[, paste0("across_cluster",k)] <- num_region_cluster/sum(num_region_cluster)
      
      df <- cluster_list[,c("gene", "new_group_label", "chr",
                            "NetworkNames",paste0("cluster",k),
                            paste0("across_cluster",k),
                            paste0("rois_cluster",k))]
      
      df$cluster_id <- k 
      df$cluster_length <- length(unlist( strsplit(as.character(df[1,paste0("rois_cluster",k)]),"[,]") ))
      
      chi_matrix <- matrix(NA, nrow=2, ncol=7)
      chi_matrix[1,] <- num_region_cluster
      chi_matrix[2,] <- num_region-num_region_cluster
      test_result <- fisher.test(chi_matrix)
      df$chi_pvalue <- test_result$p.value
      
      colnames(df)[c(5,6,7)]<- c("cluster","across_cluster","rois_cluster")
      
      cluster_summary_all <- rbind(cluster_summary_all, df)
    }
    
  }
}

sum(cluster_summary_all$chi_pvalue< 0.05)
select_cluster_pvalue <- cluster_summary_all[cluster_summary_all$chi_pvalue< 0.05,]
head(cluster_summary_all)


cluster_gene <- cluster_summary_all$gene
cluster_gene <- strsplit(cluster_gene, "[,]")
cluster_gene<- unlist(cluster_gene)
cluster_gene <- unique(cluster_gene)
length(cluster_gene)

gene_label <- read.table("magma.genes.out", header=T)

length(intersect(cluster_gene,gene_label$GENE ))


cluster_summary_all$gene_label=NA
for(i in 1:nrow(cluster_summary_all))
{
  gene_list <- cluster_summary_all$gene[i]
  gene_list <- unlist(strsplit(gene_list, "[,]"))
  
  
  if(length(gene_label$SYMBOL[gene_label$GENE==gene_list[1]])==0)
  {
    cluster_summary_all$gene_label[i] = gene_list[1]
  }else{
  cluster_summary_all$gene_label[i]=  gene_label$SYMBOL[gene_label$GENE==gene_list[1]]
  }
  if(length(gene_list)>1){
  for(j in 2:length(gene_list)){
    
    if(length(gene_label$SYMBOL[gene_label$GENE==gene_list[j]])==0)
    {
      cluster_summary_all$gene_label[i] = paste0(cluster_summary_all$gene_label[i], ",", gene_list[j])
      
      
    }else{
      cluster_summary_all$gene_label[i] <- paste0(cluster_summary_all$gene_label[i], ",", gene_label$SYMBOL[gene_label$GENE==gene_list[j]])
      
    }
    
  }
  }
  
}


head(cluster_summary_all)

# make table 1 
table1_long <- cluster_summary_all[,c("gene_label","cluster_length","NetworkNames","across_cluster","cluster_id")]

table1_wide <- reshape(table1_long, direction = "wide", timevar="NetworkNames", idvar =c("gene_label","cluster_length","cluster_id"))

colnames(table1_wide) <- c("Gene", "# of regions","cluster_id","Frontoparietal", "Default",  "Dorsal Attention", "Limbic",  "Salience", "Somatomotor",  "Visual")


table1 <- table1_wide[, -3]

  
save(cluster_summary_all,yeo7_overlap,yeo7_overlap_df,table1_wide, file="yeo7_overlap_summary_result.RData")
