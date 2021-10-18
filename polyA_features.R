### Script for analyzing poly(A) features in
args = commandArgs(trailingOnly=TRUE)
tpm = read.table(args[1], header=T, stringsAsFactors = F) # TPM.gene.tsv
features = read.table(args[2], header=F, stringsAsFactors = F) # polyA_features.bed
gff = read.table(args[3], stringsAsFactors = F, sep='\t', quote='"') # Arabidopsis_thaliana.TAIR10.46.gff3.gz

# snoRNAs = gsub('ID=transcript:([^\\.]+).*$','\\1',gff[which(gff[,3] == 'snoRNA'),9])

gff[gff[,1] == 'Mt',1] = 'm'
gff[gff[,1] == 'Pt',1] = 'c'
gff[,1] = paste('Ath_chr', gff[,1], sep='')



utrs = gff[gff[,3] == 'three_prime_UTR',]
utrs[utrs[,7] == '+', 5] = utrs[utrs[,7] == '+', 5] + 200
utrs[utrs[,7] == '-', 4] = utrs[utrs[,7] == '-', 4] - 200
utrs[,9] = gsub('Parent=transcript:','',utrs[,9])

exons = gff[gff[,3] == 'exon',]
exons[,9] = gsub('Parent=transcript:([^;]+).*$','\\1',exons[,9])
no_utr_exons = exons[! exons[,9] %in% utrs[,9],]
last_exons = unlist(sapply(unique(no_utr_exons[,9]),function(x){
  these_exons = no_utr_exons[no_utr_exons[,9]==x,]
  if(these_exons[1,7] == '+'){
    line = these_exons[nrow(these_exons),]
    line[,5] = line[,5] + 200
  }else{
    line = these_exons[1,]
    line[,4] = line[,4] - 200
  }
  return(line)
}))
last_exons = matrix(last_exons, ncol = 9, byrow = T)

overlapping_utrs = apply(features, 1, function(element){
  FL = as.numeric(element[2])
  FR = as.numeric(element[3])
  chromutrs = utrs[utrs[,1] == as.character(element[1]) & utrs[,7] == as.character(element[6]),]
  overlaps = chromutrs[as.numeric(chromutrs[,4]) <= FR & as.numeric(chromutrs[,5]) >= FL,9]
  overlapping_genes = paste(unique(gsub('\\.[0-9]+$','',overlaps)), collapse = ',')
  if(overlapping_genes != ''){
    return(overlapping_genes)
  }else{
    chromutrs = last_exons[last_exons[,1] == as.character(element[1]) & last_exons[,7] == as.character(element[6]),]
    overlaps = chromutrs[as.numeric(chromutrs[,4]) <= FR & as.numeric(chromutrs[,5]) >= FL,9]
    return(paste(unique(gsub('\\.[0-9]+$','',overlaps)), collapse = ','))
  }
})
names(overlapping_utrs) = features[,4]

noncanonical = features[features[,4] %in% names(which(overlapping_utrs == '')),]
genes = gff[gff[,3] == "gene",]
genes[,9] = gsub('ID=gene:([^;]+).*$','\\1',genes[,9])
overlapping_genes = apply(noncanonical, 1, function(element){
  FL = as.numeric(element[2])
  FR = as.numeric(element[3])
  chromgenes = genes[genes[,1] == as.character(element[1]) & genes[,7] == as.character(element[6]),]
  overlaps = chromgenes[as.numeric(chromgenes[,4]) <= FR & as.numeric(chromgenes[,5]) >= FL,9]
  overlapping_genes = paste(unique(gsub('\\.[0-9]+$','',overlaps)), collapse = ',')
  return(overlapping_genes)
})
names(overlapping_genes) = noncanonical[,4]


excluded_regions = c(
  features[features[,1] == 'Ath_chr2' & features[,3] < 15000,4], # Chr2 rDNA
  features[features[,1] == 'Ath_chr3' & features[,3] < 14205000 & features[,4] > 14190000, 4], # Chr3 rDNA
  features[features[,1] == 'Ath_chrc',4],
  features[features[,1] == 'Ath_chrm',4]
)

filtered_features = features[!features[,4] %in% excluded_regions,]


cov_list = list()
for(i in c('wt_1','wt_2','wt_3','ect234_1','ect234_2','ect234_3')){
  cov_list[[i]] = list()
  cat(i,'\n')
  ct = read.table(paste('feature_coverage/',i,'.coverage.tsv',sep=''), sep='\t', stringsAsFactors = F)
  for(n in 1:nrow(ct)){
    cov_list[[i]][[ct[n,1]]] = as.numeric(strsplit(ct[n,2],',')[[1]])
  }
}

cov_false_list = list()
for(i in c('wt_1','wt_2','wt_3','ect234_1','ect234_2','ect234_3')){
  cov_false_list[[i]] = list()
  cat(i,'\n')
  ct = read.table(paste('feature_coverage/',i,'.false.coverage.tsv',sep=''), sep='\t', stringsAsFactors = F)
  for(n in 1:nrow(ct)){
    cov_false_list[[i]][[ct[n,1]]] = as.numeric(strsplit(ct[n,2],',')[[1]])
  }
}


cov_wt = sapply(features[,4], function(x){
  return(
    cov_list[['wt_1']][[x]] + 
    cov_list[['wt_2']][[x]] + 
    cov_list[['wt_3']][[x]] 
  )
})
cov_ect = sapply(features[,4], function(x){
  return(
    cov_list[['ect234_1']][[x]] + 
    cov_list[['ect234_2']][[x]] + 
    cov_list[['ect234_3']][[x]] 
  )
})
cov_total = sapply(features[,4], function(x){
  return(cov_wt[[x]] + cov_ect[[x]])
})

cov_false_total = sapply(features[,4], function(x){
  return(
    cov_false_list[['wt_1']][[x]] + 
    cov_false_list[['wt_2']][[x]] + 
    cov_false_list[['wt_3']][[x]] + 
    cov_false_list[['ect234_1']][[x]] + 
    cov_false_list[['ect234_2']][[x]] + 
    cov_false_list[['ect234_3']][[x]] 
  )
})

false_clusters = names(which(unlist(lapply(cov_false_total, sum)) > unlist(lapply(cov_total, sum))))
filtered_features = filtered_features[!filtered_features[,4] %in% false_clusters,]


canonical_sites = intersect(filtered_features[,4], names(which(overlapping_utrs != "")))
noncanonical_sites = intersect(filtered_features[,4], names(which(overlapping_utrs == "")))
APA = names(which(overlapping_genes[noncanonical_sites] != ""))

output_features = filtered_features
output_features[,4] = overlapping_utrs[filtered_features[,4]]
output_features[overlapping_utrs[filtered_features[,4]]=="",4] = overlapping_genes[names(which(overlapping_utrs[filtered_features[,4]] == ""))]
output_features[,4] = make.names(output_features[,4], unique = T)
output_features[,4] = gsub("^X", "polyA", output_features[,4])
write.table(output_features, 'dec14/polyA_clusters.bed', quote=F, sep='\t', row.names = F, col.names = F)

lapply(genesets, function(x)c('num'=sum(x%in%multisite), 'pct'=mean(x%in%multisite)))

multisite = names(which(table(gsub('\\.[0-9]+$', '', output_features[,4])) > 1))
multisite = intersect(multisite, c(genesets[["A"]], genesets[["C"]]))

fishertable = rbind(
  'unchanged' = c( sum(dominant_sites_agree[intersect(multisite, genesets[['B']])]), sum(dominant_sites_agree[intersect(multisite, genesets[['C']])]) ),
  'changed' = c( sum(!dominant_sites_agree[intersect(multisite, genesets[['B']])]), sum(!dominant_sites_agree[intersect(multisite, genesets[['C']])]) )
)
colnames(fishertable) = c('targets', 'non-targets')
fisher.test(fishertable)

top2200 = names(sort(mean_TPM[genesets$C],decreasing = T)[1:2200])
sum(top2200 %in% multisite)

fishertable = rbind(
  'unchanged' = c( sum(dominant_sites_agree[intersect(multisite, top2200)]), sum(dominant_sites_agree[intersect(multisite, genesets[['C']])]) ),
  'changed' = c( sum(!dominant_sites_agree[intersect(multisite, top2200)]), sum(!dominant_sites_agree[intersect(multisite, genesets[['C']])]) )
)

gene = multisite[1]
dominant_sites_agree = sapply(multisite, function(gene){
  clusternames = filtered_features[grep(gene, output_features[,4]),4]
  wt_abund = sapply(clusternames, function(cluster)ifelse(cluster%in%names(cov_wt), sum(cov_wt[[cluster]]), 0))
  ect_abund = sapply(clusternames, function(cluster)ifelse(cluster%in%names(cov_ect), sum(cov_ect[[cluster]]), 0))
  return(names(wt_abund)[order(wt_abund, decreasing = T)[1]] == names(ect_abund)[order(ect_abund, decreasing = T)[1]])
})

cov_diff = sapply(features[,4], function(x){
  return(cov_wt[[x]] - cov_ect[[x]])
})
peak_position = lapply(cov_total, function(x){
  order(x, decreasing = T)[1]
})

feature_count_list = lapply(cov_list, function(x)unlist(lapply(x, sum)))
feature_counts = cbind(
  feature_count_list[['wt_1']],
  feature_count_list[['wt_2']],
  feature_count_list[['wt_3']],
  feature_count_list[['ect234_1']],
  feature_count_list[['ect234_2']],
  feature_count_list[['ect234_3']]
)
colnames(feature_counts) = c('wt_1','wt_2','wt_3','ect234_1','ect234_2','ect234_3')


count_data = feature_counts[filtered_features[,4],]
rownames(count_data) = output_features[,4]

APA_proportion = t(sapply(multisite, function(x){
  subtable = count_data[grep(x, output_features[,4], value = T),,drop=FALSE]
  dominant_site = rownames(subtable)[order(rowSums(subtable[,1:3]), decreasing = T)[1]]
  nondom = colSums(subtable[rownames(subtable) != dominant_site,,drop=FALSE])
  tot = colSums(subtable)
  return(c(wt=sum(nondom[1:3])/sum(tot[1:3]), ect234=sum(nondom[4:6])/sum(tot[4:6])))
}))

library(beeswarm)
beeswarm(as.data.frame(APA_proportion[intersect(targets.stringent,multisite),]), pch='.', spacing = .4)


samples = data.frame(
  sample=c('wt_1','wt_2','wt_3','ect234_1','ect234_2','ect234_3'),
  condition=c('wt','wt','wt','mut','mut','mut')
)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(feature_counts[c(canonical_sites, noncanonical_sites),]), colData = samples, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "wt")
# keep = rowSums(counts(dds)) >= 3
# dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, independentFiltering = F)

downregulated = rownames(res)[res$log2FoldChange < (-1) & res$padj < 0.05]
downregulated = downregulated[!is.na(downregulated)]

upregulated = rownames(res)[res$log2FoldChange > (1) & res$padj < 0.05]
upregulated = upregulated[!is.na(upregulated)]


DEcolors = rep('gray70', nrow(res))
names(DEcolors) = rownames(res)
DEcolors[downregulated] = 'cadetblue'
DEcolors[upregulated] = 'firebrick2'

lfc = res$log2FoldChange
lfc[lfc > 5] = 5
lfc[lfc < (-5)] = -5
pv = -log10(res$padj)
pv[pv > 40] = 40


plot(x = lfc, y=pv, pch='.', cex=2, col=DEcolors, xlab='Log2 Fold Change', ylab='-log10(adj. p-value)')


pdf('polyA_cluster_details.pdf', 6, 4, useDingbats = F, pointsize = 8)
par(mfrow=c(1,2), lwd=1, lend='butt', ljoin='mitre')
  feature_table = table(filtered_features[,3] - filtered_features[,2])
  # feature_table = table(ceiling((filtered_features[,3] - filtered_features[,2])/10)*10)
  plot(NA, xlim=c(0,500), ylim=c(0,max(feature_table)), main = "Width of poly(A) clusters", xlab="poly(A) cluster length (nt)", ylab="Frequency", las=1)
  polygon(x=c(as.numeric(names(feature_table)),0), y=c(feature_table,0), border = NA, col='firebrick3', lty=1)
  lines(x=c(0,as.numeric(names(feature_table))), y=c(0,feature_table), lwd=1, col='black', lty=1)
  abline(h=0)
  
  number_of_sites = unlist(lapply(cov_total, function(x)sum(x>0)))
  positions_table = table(unlist(number_of_sites))
  # positions_table = table(ceiling(unlist(number_of_sites)/10)*10)
  plot(NA, xlim=c(0,150), ylim=c(0,max(positions_table)), xlab="Number of poly(A) positions", ylab="Frequency", main="Poly(A) sites per cluster")
  polygon(x=c(as.numeric(names(positions_table)),0), y=c(positions_table,0), border = NA, col='firebrick3', lty=1)
  lines(x=c(0,as.numeric(names(positions_table))), y=c(0,positions_table), lwd=1, col='black', lty=1)
  abline(h=0)
dev.off()


library(pheatmap)
heatmat = matrix(data=unlist(lapply(cov_total, function(x)x[1:400])), ncol = 400, byrow = T)
rownames(heatmat) = names(cov_total)

pheatmap(
  heatmat[names(overlapping_utrs)[overlapping_utrs %in% targets.stringent],],
  cluster_cols = F, cluster_rows = F, 
  scale = 'row'
)
