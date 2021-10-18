#setwd("C:/Users/michael.schon/Documents/HOME_WORKSPACE/brodersen/")

samplenames = c('wt', 'ect234')

polyA_consensus = list()
# polyA_median = list()
for(name in samplenames){
  # polyA_consensus[[name]] = read.table(paste(wd,'/smartseq2_data/polyA_bedgraphs/',name,'.consensus.bed',sep = ''), header=F, stringsAsFactors = F, col.names = c('chrom', 'start', 'end', 'gene', 'score', 'strand'))
  polyA_consensus[[name]] = read.table(paste(name,'.consensus.bed',sep = ''), header=F, stringsAsFactors = F, col.names = c('chrom', 'start', 'end', 'gene', 'score', 'strand'))
  rownames(polyA_consensus[[name]]) = polyA_consensus[[name]][,'gene']
}

shared = intersect(polyA_consensus[["wt"]][,"gene"], polyA_consensus[["ect234"]][,"gene"])

oldgenesets = list()
oldgenesets[['A']] = scan('targets_permissive.txt', 'character')
oldgenesets[['B']] = scan('targets_stringent.txt', 'character')
oldgenesets[['C']] = scan('nontargets.txt', 'character') 

load('genes.list_proto.Rdat')
genesets = list()
genesets[['A']] = genes.list.proto$genes.roots
genesets[['B']] = genes.list.proto$genes.stringent
genesets[['C']] = genes.list.proto$genes.nontargets.roots




geneTPM = read.table('TPM.gene.tsv', sep='\t', header=T, stringsAsFactors = F)
test_set = intersect(as.character(unique(unlist(genesets))), rownames(geneTPM))
class_by_gene = rep(NA,length(test_set))
names(class_by_gene) = test_set
class_by_gene[genesets[['C']]] = 'C'
class_by_gene[genesets[['A']]] = 'A'
class_by_gene[genesets[['B']]] = 'B'


wt = c('wt_1','wt_2','wt_3')
ect = c('ect234_1','ect234_2','ect234_3')

TPM_wt = rowMeans(geneTPM[test_set,wt])
TPM_ect = rowMeans(geneTPM[test_set,ect])

mean_TPM = (TPM_wt + TPM_ect)/2
sorted_TPM = sort(mean_TPM)
cutoff = 0
shared = intersect(shared,names(mean_TPM[mean_TPM > cutoff]))

geneset_colors = c('#40B35C', '#006F00', '#595959')
names(geneset_colors) = c('A','B','C')

subsets = lapply(genesets, function(x)intersect(x,shared))

polyA_diffs = list()
diff_tables = list()
diff_bins = list()
# bins = rbind(
#   c(-10000, -1001),
#   c(-1000, -101),
#   c(-100, -51),
#   c(-50, -11),
#   c(-10, -2),
#   c(-1, 1),
#   c(2, 10),
#   c(11, 50),
#   c(51, 100),
#   c(101, 1000),
#   c(1001, 10000)
# )
bins = rbind(
  c(-10000, -1001),
  c(-1000, -101),
  c(-100, -11),
  c(-10, 10),
  c(11, 100),
  c(101, 1000),
  c(1001, 10000)
)

for(geneset in c('A','B','C')){
  control_set = polyA_consensus[["wt"]][subsets[[geneset]],]
  control_plus = control_set[control_set[,"strand"]=='+',]
  control_minus = control_set[control_set[,"strand"]=='-',]
  
  test_set = polyA_consensus[["ect234"]][subsets[[geneset]],]
  test_plus = test_set[test_set[,"strand"]=='+',]
  test_minus = test_set[test_set[,"strand"]=='-',]
  
  polyA_diffs[[geneset]] = c(
    test_plus[,"start"] - control_plus[,"start"],
    -(test_minus[,"start"] - control_minus[,"start"])
  )
  names(polyA_diffs[[geneset]]) = c(rownames(test_plus), rownames(test_minus))
  diff_tables[[geneset]] = table(polyA_diffs[[geneset]])/length(polyA_diffs[[geneset]])
  diff_bins[[geneset]] = apply(bins, 1, function(x){sum(polyA_diffs[[geneset]] >= x[1] & polyA_diffs[[geneset]] <= x[2])})
}

pdf('polyA_site_distribution_consensus.highest.pdf', useDingbats = F, width = 5, height = 4, pointsize = 8)
  layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow = TRUE))
  par(mar=c(0,4,1,1), oma=c(4,0,0,0), lend='butt', ljoin='mitre')
  # bin_names = c('<-1000', '-1000:-101', '-100:-51', '-50:-11', '-10:-2', '-1:1', '2:10', '11:50', '51:100', '101:1000', '>1000')
  bin_names = c('<-1000', '-1000:-101', '-100:-11', '-10:10', '11:100', '101:1000', '>1000')
  binned_matrix = rbind(
    diff_bins[['C']]/length(subsets[["C"]]),
    diff_bins[['A']]/length(subsets[["A"]]),
    diff_bins[['B']]/length(subsets[["B"]])
    
  )
  barplot(
    binned_matrix, 
    col = geneset_colors[c('C','A','B')], 
    border=NA, names.arg = bin_names, 
    las=2, ylim=c(0,0.9), 
    ylab='Frequency', beside = T, 
    xlab='Distance from wild-type poly(A) site'
  )
  
  # plot(ecdf(polyA_diffs[["C"]]),col=geneset_colors['C'],xlim=c(-200,200),ylim=c(0,1), verticals = T, pch = '.', lwd=2, ylab='Cumulative frequency', xlab='Distance from wild-type poly(A) site', main=NA)
  # lines(ecdf(polyA_diffs[['A']]),col=geneset_colors['A'],xlim=c(-200,200),ylim=c(0,1), verticals = T, pch = '.', lwd=2)
  # lines(ecdf(polyA_diffs[['B']]),col=geneset_colors['B'],xlim=c(-200,200),ylim=c(0,1), verticals = T, pch = '.', lwd=2)
  ymax = .75
  xvals = as.numeric(names(diff_tables[['C']])) 
  yvals = as.numeric(diff_tables[['C']])
  plot(NA,xlim=c(-50,50), type='l', ylim=c(0,ymax), col=geneset_colors['C'], ylab=NA, las=1)
  abline(v=0,lty=2)
  polygon(x=c(xvals,xvals[length(xvals)],xvals[1]), y=c(yvals,0,0), col=geneset_colors['C'], border=NA)
  abline(h=0)
  
  xvals = as.numeric(names(diff_tables[['A']])) 
  yvals = as.numeric(diff_tables[['A']])
  plot(NA,xlim=c(-50,50), type='l', ylim=c(0,ymax), col=geneset_colors['A'], ylab=NA, las=1)
  abline(v=0,lty=2)
  polygon(x=c(xvals,xvals[length(xvals)],xvals[1]), y=c(yvals,0,0), col=geneset_colors['A'], border=NA)
  abline(h=0)
  
  xvals = as.numeric(names(diff_tables[['B']])) 
  yvals = as.numeric(diff_tables[['B']])
  plot(NA,xlim=c(-50,50), type='l', ylim=c(0,ymax), col=geneset_colors['B'], ylab='Density', las=1)
  abline(v=0,lty=2)
  polygon(x=c(xvals,xvals[length(xvals)],xvals[1]), y=c(yvals,0,0), col=geneset_colors['B'], border=NA)
  abline(h=0)

dev.off()




