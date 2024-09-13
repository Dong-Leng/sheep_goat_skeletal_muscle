# Analysis in R platform
# R v4.3.2
# Tximport import quantitative files of kallisto output
library(tximport)
library(tximeta)
base_dir = "./kallisto_out"
sample_id = dir(file.path(base_dir))
files = file.path(base_dir, sample_id, "abundance.h5")
names(files) = paste0(sample_id)
tx2gene = read.table(file = "tx2gene", header = TRUE)
# "tx2gene" : Documents corresponding to transcripts id and genes id
all_sample_quant.tsv = tximport(files, type="kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene=tx2gene, ignoreTxVersion = T)


# Quantitative tables for sheep and goats were merged
library(dplyr)
goat_all_TPM <- read.csv("./goat_all_quant_TPM.csv")
sheep_all_TPM <- read.csv("./sheep_all_quant_TPM.csv")
merge_TPM <- merge(goat_all_TPM,sheep_all_TPM,by = "Gene")
unique_TPM <- distinct(merge_TPM, Gene, .keep_all = TRUE)
goat_all_counts <- read.csv("./goat_all_counts.csv")
sheep_all_counts <- read.csv("./sheep_all_counts.csv")
merge_counts <- merge(goat_all_counts,sheep_all_counts,by = "Gene")
unique_counts <- distinct(merge_counts, Gene, .keep_all = TRUE)
write.csv(unique_TPM,"all_TPM.csv",row.names = F)
write.csv(unique_counts,"all_counts.csv",row.names = F)


# T-SNE analysis and plot
library(ggpubr)
library(ggthemes)
library(Rtsne)
all_TPM <- read.csv("all_TPM.csv",row.names = 1)
group <- colnames(all_TPM) %>% as.data.frame()
colnames(group) <- "Sample"
group <- separate(group1, Sample, into = c("Individual", "Muscle"), sep = "_")
group$Species <- gsub("[0-9]", "", group$Individual)
group <- group[order(group1$Species),]
set.seed(150)
tsne <- Rtsne(t(all_TPM), perplexity = 3)
colnames(tsne$Y) <- c("TSNE1","TSNE2")
tsne_data <- data.frame(sample=colnames(all_TPM),
                         Type = group$Species,
                        tsne$Y)
p1 <- ggplot(tsne_data, aes(TSNE1, TSNE2, color = Type, shape = Type))+ 
  geom_point(size = 3.2)+ 
  geom_hline(yintercept = 0,linetype="dashed") + 
  geom_vline(xintercept = 0,linetype="dashed") + 
  theme_bw() + 
  geom_mark_ellipse(data = tsne_data, aes(x = TSNE1, y = TSNE2, color = Type), expand = unit(0, "mm"))+ 
  scale_color_manual(values=c("#ECA8A9","#6CA9D1"))+
  scale_shape_manual(values = c(1:2))
ggsave("t-SNE plot.pdf", plot=p1, width = 6.7, height =5.5)


# Spearman correlation heatmap (no clustering)
library(tidyr)
all_TPM_1 <- all_TPM[, match(group1$Sample, colnames(all_TPM))]
sample_cor <- round(cor(all_TPM_1,method = "spearman"),digits = 5)
anno <- group1[,c(2,4)]
musclecolor <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF","#197EC0FF","#F05C3BFF","#46732EFF","#71D0F5FF")
names(musclecolor) <- unique(group1$muscle)
breedcolor <- c("#ECA8A9","#6CA9D1")
names(breedcolor) <- unique(group1$breed)
ann_colors <- list(muscle=musclecolor, breed = breedcolor)
rownames(anno) <- rownames(sample_cor)
library(pheatmap)
p2 <- pheatmap(sample_cor,
         cluster_rows = F,
         cluster_cols = F,
         border = F,  
         annotation_col = anno,
         annotation_colors = ann_colors,
         show_rownames =F,
         show_colnames =F,
         scale = "none", 
         color=colorRampPalette(c("#4575B4","white","#D73027"))(200)
   ) 
ggsave("TPM_Spearman_cor_heatmap.pdf", plot=p2, width = 7, height =6)


# Analysis of DEG between muscles of different parts of sheep and goats
library(DESeq2)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(dplyr)
all_counts <- read.csv("sheep_all_counts.csv",header = T,check.names = F)
unique_data <- distinct(all_counts, Gene, .keep_all = TRUE)
rownames(unique_data) <- unique_data$Gene
unique_data <- unique_data[,-1]
unique_data <- round(unique_data)
info = data.frame(sample = colnames(unique_data), condition = as.character(lapply(strsplit(colnames(unique_data), "_"), tail, 1)))
condition = as.factor(info[["condition"]])
infoData = data.frame(row.names = info[["sample"]], condition)
dds = DESeqDataSetFromMatrix(unique_data, infoData, design= ~ condition)
dds = DESeq(dds)
sample_combine = combn(unique(info$condition), 2, simplify=F)
diff_num_all <- as.data.frame(matrix(nrow=10,ncol=10))
rownames(diff_num_all) = colnames(diff_num_all) = unique(info$condition)
for (i in sample_combine) { 
  result01 = results(dds, contrast = c("condition", i[1], i[2]))
  result02 = result01[order(result01$pvalue),]
  filename1 = paste("diff_", i[1], "vs", i[2], ".csv", sep = "")
  write.csv(result02, file = filename1)
  diff_gene_deseq2 = subset(result02, padj < 0.01 & (log2FoldChange > 2 | log2FoldChange < -2))
  filename2 = paste("signification_diff_", i[1], "vs", i[2], ".csv", sep = "")
  write.csv(diff_gene_deseq2, file = filename2)
  diff_num_all[i[2],i[1]] = length(rownames(diff_gene_deseq2))
}
all_counts <- read.csv("goat_all_counts.csv",header = T,check.names = F)
unique_data <- distinct(all_counts, Gene, .keep_all = TRUE)
rownames(unique_data) <- unique_data$Gene
unique_data <- unique_data[,-1]
unique_data <- round(unique_data)
info = data.frame(sample = colnames(unique_data), condition = as.character(lapply(strsplit(colnames(unique_data), "_"), tail, 1)))
condition = as.factor(info[["condition"]])
infoData = data.frame(row.names = info[["sample"]], condition)
dds = DESeqDataSetFromMatrix(unique_data, infoData, design= ~ condition)
dds = DESeq(dds)
sample_combine = combn(unique(info$condition), 2, simplify=F)
diff_num_all <- as.data.frame(matrix(nrow=10,ncol=10))
rownames(diff_num_all) = colnames(diff_num_all) = unique(info$condition)
for (i in sample_combine) { 
  result01 = results(dds, contrast = c("condition", i[1], i[2]))
  result02 = result01[order(result01$pvalue),]
  filename1 = paste("diff_", i[1], "vs", i[2], ".csv", sep = "")
  write.csv(result02, file = filename1)
  diff_gene_deseq2 = subset(result02, padj < 0.01 & (log2FoldChange > 2 | log2FoldChange < -2))
  filename2 = paste("signification_diff_", i[1], "vs", i[2], ".csv", sep = "")
  write.csv(diff_gene_deseq2, file = filename2)
  diff_num_all[i[1],i[2]] = length(rownames(diff_gene_deseq2))
}
library(pheatmap)
diff_num_all[is.na(diff_num_all)] <- 0
p3 <- pheatmap(diff_num_all, 
               cluster_rows = F, 
               cluster_cols = F, 
               border = F,
               show_rownames = T,
               legend = T,
               display_numbers=T,
               number_format = "%d", 
               fontsize=13,
               number_color = "black")
ggsave("each_muscle_DEG_all_heatmap.pdf", plot=p3, width = 6.7, height =5.8)


# DEG analysis of the same muscle sites in sheep and goats
goat_all_TPM <- read.csv("goat_all_counts.csv",header = T,check.names = F)
unique_data <- distinct(goat_all_TPM, Gene, .keep_all = TRUE)
rownames(unique_data) <- unique_data$Gene
goat_all_TPM1 <- unique_data
sheep_all_TPM <- read.csv("sheep_all_counts.csv",header = T,check.names = F)
unique_data1 <- distinct(sheep_all_TPM, Gene, .keep_all = TRUE)
rownames(unique_data1) <- unique_data1$Gene
sheep_all_TPM1 <- unique_data1
goat_sheep_merge <- merge(goat_all_TPM1, sheep_all_TPM1, by = "Gene")
rownames(goat_sheep_merge) <- goat_sheep_merge$Gene
goat_sheep_merge <- goat_sheep_merge[,-1]
goat_sheep_merge <- round(goat_sheep_merge)
info = data.frame(sample = colnames(goat_sheep_merge), condition = as.character(gsub("\\d+", "", colnames(goat_sheep_merge))))
condition = as.factor(info[["condition"]])
infoData = data.frame(row.names = info[["sample"]], condition)
dds = DESeqDataSetFromMatrix(goat_sheep_merge, infoData, design= ~ condition)
dds = DESeq(dds)
sample_combine = combn(unique(info$condition), 2, simplify=F)
sample_combine1 = sample_combine[c(10,29,47,64,80,95,109,122,134,145)]
goat_sheep_alltime_diff <- data.frame()
for (i in sample_combine1) { 
  result01 = results(dds, contrast = c("condition", i[1], i[2]))
  result02 = result01[order(result01$pvalue),]
  filename1 = paste("diff_", i[1], "vs", i[2], ".csv", sep = "")
  write.csv(result02, file = filename1)
  diff_gene_deseq2 = subset(result02, padj < 0.01 & (log2FoldChange > 2 | log2FoldChange < -2))
  filename2 = paste("signification_diff_", i[1], "vs", i[2], ".csv", sep = "")
  write.csv(diff_gene_deseq2, file = filename2)
  diff_gene_deseq3 <- as.data.frame(diff_gene_deseq2)
  diff_gene_deseq3$cluster = i[1]
  diff_gene_deseq3$gene = rownames(diff_gene_deseq3)
  goat_sheep_alltime_diff = rbind(goat_sheep_alltime_diff, diff_gene_deseq3)
}
goat_sheep_alltime_diff2 <- goat_sheep_alltime_diff[,c(8,2,5,6,7)]
colnames(goat_sheep_alltime_diff2) <- c("gene","avg_log2FC","p_val","p_val_adj","cluster")
my_cols <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF","#197EC0FF","#F05C3BFF","#46732EFF","#71D0F5FF")
library(scRNAtoolVis)
p4 <- jjVolcano(diffData = goat_sheep_alltime_diff2,tile.col = my_cols,topGeneN = 0,size = 3.5,fontface = 'italic')
ggsave("goat_vs_sheep.pdf", plot=p4, width = 9, height = 7.5)


# Genes specifically expressed in different muscle sites of sheep and goats and their functions
library(pheatmap)
library(ggplot2)
library(dplyr)
library(gprofiler2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)
tauData = read.table("goat_all_quant_TPM.TAU.txt", header = TRUE)
rownames(tauData) = tauData[['ID']]
tauData = tauData[tauData[["Specific...0.75."]] == TRUE, ]
GO_info <- data.frame()
for (part in unique(tauData[['Tissue']])) { 
  one_part = tauData[tauData[['Tissue']]==part, ]
  orth_info = gorth(
    query = one_part[['ID']],
    source_organism = "chircus", 
    target_organism = "hsapiens",
    numeric_ns = "",
    mthreshold = Inf,
    filter_na = FALSE
  )
  go_ifo <- data.frame(Gene = orth_info$ortholog_name,
                       Muscle = part)
  GO_info = rbind(GO_info, go_ifo)
}
GO_info1 <- na.omit(GO_info)
GO_info1 <- GO_info1[GO_info1$Gene != "N/A",]
fenshu2xiaoshu <- function(ratio){
  sapply(ratio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
}
group <- unique(GO_info1$Muscle)
for(i in group){ 
  each_markers <- filter(GO_info1, Muscle %in% i)
  ego_ALL <- enrichGO(gene = each_markers$Gene,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = 'SYMBOL',
                      ont           = "ALL",  
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1)    
  ego_all <- data.frame(ego_ALL)
  enrichment_fold = fenshu2xiaoshu(ego_all$GeneRatio)/fenshu2xiaoshu(ego_all$BgRatio)
  enrichment_fold = as.numeric(enrichment_fold)
  ego_all$enrichment_fold = enrichment_fold
  log10P = -log10(ego_all$pvalue)
  log10P = as.numeric(log10P)
  ego_all$log10P = log10P
  write.csv(ego_all, paste0(i,"_GO.csv"))
}
goat_go$Description = factor(goat_go$Description, levels = rev(unique(goat_go$Description)), ordered = TRUE)
goat_go$Group = factor(goat_go$Group, levels = unique(goat_go$Group), ordered = TRUE)
p5 = ggplot(goat_go, aes(x = Group, y = Description, size = Count, color = log10P))+ 
  geom_point(stat ="identity")+
  scale_colour_gradient(low="#d90424",high="#374a89")+ 
  labs(
    color=expression(-Log10P),
    size="Count Number",
    x=""
  )+
  theme_bw()+
  scale_size_continuous(range = c(3.5, 7))+  
  theme(
    axis.text.y = element_text(size = rel(1.2)),
    axis.title.x = element_text(size=rel(1.2)),
    axis.title.y = element_blank()
  )
ggsave(p5, filename = "goat_go.pdf", width = 13, height = 7.5)
gene_table = read.csv('goat_all_quant_TPM.csv', header = TRUE)
rownames(gene_table) = gene_table[['Gene']]
merge_table = merge(x = tauData,
                    y = gene_table,
                    by.x = 'ID',
                    by.y = 'Gene')
rownames(merge_table) = merge_table[['ID']]
merge_table = merge_table %>% arrange(Tissue)
plot_data = merge_table[ , (ncol(tauData)+1) : ncol(merge_table)]
plot_data = as.data.frame(t(plot_data))
plot_data[['sampleID']] = rownames(plot_data)
plot_data$group_1 <- as.character(lapply(strsplit(plot_data$sampleID, "_"), tail, 1))
plot_data = plot_data %>% arrange(group_1)
rownames(plot_data) = plot_data[['sampleID']]
plot_data_2 = plot_data[ , 1:(nrow(merge_table))]
for (gene_id in colnames(plot_data_2)) {
  plot_data_2[[gene_id]] = scale(as.numeric(plot_data_2[[gene_id]]), center=T, scale=T)  
}
data = as.data.frame(t(plot_data_2))
annotation_col = data.frame(Group = factor(plot_data[["group_1"]]))
annotation_row = data.frame(Group = factor(merge_table[["Tissue"]]))
rownames(annotation_col) = rownames(plot_data_2)
rownames(annotation_row) = colnames(plot_data_2)
ann_colors = list(Group = c("BJJ" = '#FED439FF',
                            "BMJ" = '#709AE1FF',
                            "BSTJ" = '#8A9197FF',
                            "BZCJ" = '#D2AF81FF',
                            "FCJ" = '#FD7446FF',
                            "FWXJ" = '#D5E4A2FF',
                            "GETJ" = '#197EC0FF',
                            "JXFJ" = '#F05C3BFF',
                            "SJJ" = '#46732EFF',
                            "STJ" = '#71D0F5FF'))
pdf("z-score_heatmap_of_tissue-specific_genes.pdf",height = 7, width = 4.5)
pheatmap(as.matrix(data),
         breaks=seq(-1,1,0.02),
         color = colorRampPalette(c("#295f32", "#a0c769", "#fbeeb0", "#d86e43", "#942320"))(100),
         cluster_cols = F,
         cluster_rows = F,
         border_color= NA,
         show_rownames = F,
         show_colnames = F,
         main = "",
         display_numbers=F,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors[1],
         angle_col = 90,
         fontsize_row = 12,
         fontsize_col = 12
)+theme_test()
dev.off()


# Heatmap of meat quality characteristics of different muscle sites in sheep and goats
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
all_TPM <- read.csv("amino_acid_content.csv",row.names = 1)
row_group <- data.frame(Name = rownames(all_TPM),
                        Class = all_TPM$Class)
row_group <- row_group[order(row_group$Class),]
annotation_row <- data.frame(Class = row_group[["Class"]])
classcolor <- c("#8E7FB8","#A2C9AE")
names(classcolor) <- unique(row_group$Class)
all_TPM_1 <- all_TPM[,-1]
group1 <- colnames(all_TPM_1) %>%
  as.data.frame()
colnames(group1) <- "Sample"
group1 <- separate(group1, Sample, into = c("Individual", "Muscle"), sep = "_")
group1$Sample <- paste(group1$Individual, group1$Muscle, sep = "_")
group1$Species <- gsub("[0-9]", "", group1$Individual)
group1 <- group1[order(group1$Species, group1$Muscle),]
all_TPM_2 <- all_TPM_1[, match(group1$Sample, colnames(all_TPM_1))]
all_TPM_2 <- all_TPM_2[row_group$Name, ]
all_TPM_3 <- as.data.frame(apply(all_TPM_2, 2, function(x) log2(x + 1)))
annotation_col <- group1[,c(2,4)]
Musclecolor <- c("#709AE1FF","#FD7446FF","#D5E4A2FF","#197EC0FF","#F05C3BFF","#46732EFF","#71D0F5FF")
names(Musclecolor) <- unique(group1$Muscle)
Speciescolor <- c("#ECA8A9","#6CA9D1")
names(Speciescolor) <- unique(group1$Species)
ann_colors <- list(Muscle=Musclecolor, Species = Speciescolor, Class = classcolor)
rownames(annotation_col) <- colnames(all_TPM_3)
rownames(annotation_row) <- rownames(all_TPM_3)
p6 <- pheatmap(all_TPM_3,
         cluster_rows = F,
         cluster_cols = F,
         border = F,  
         border_color= NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames =T,
         show_colnames =F,
         scale = "none",
         color=colorRampPalette(c("#4575B4","white","#D73027"))(200)
   ) 
ggsave("heatmap.pdf", plot=p6, width = 7, height =4)


# Weighted gene co-expression network analysis (WGCNA)
library(WGCNA)
library(reshape2)
library(stringr)
library(dplyr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
all_TPM <- read.csv("goat_all_quant_TPM.csv")
unique_data <- distinct(all_TPM, Gene, .keep_all = TRUE)
rownames(unique_data) <- unique_data$Gene
unique_data <- unique_data[,-1]
patterns <- c("BJJ", "BSTJ", "BZCJ")
pattern <- paste(patterns, collapse = "|")
unique_data_filtered <- unique_data[, !grepl(pattern, colnames(unique_data))]
RNAseq_voom <- unique_data_filtered 
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr0 <- WGCNA_matrix
datExpr <- datExpr0
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
datTraits <- read.csv("goat_phenotype.csv")
rownames(datTraits) <- datTraits$Sample
datTraits$Sample <- factor(datTraits$Sample,levels = unique(datTraits$Sample), ordered = TRUE)
sampleNames = rownames(datExpr)
traitRows = match(sampleNames, datTraits$Sample)
rownames(datTraits) = datTraits[traitRows, 1]
datTraits <- datTraits[,-1]
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("goat_soft_threshold.pdf",width = 8.3,height = 4.2)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
mergedColors = labels2colors(net$colors)
pdf("goat_tree+module.pdf",width = 7,height = 4.5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Dynamic tree cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  moduleColors <- labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  sizeGrWindow(10,6)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  pdf("goat_Module-meat_quality-relationships.pdf",width = 6.5,height = 8)
  par(mar = c(6, 8.5, 3, 3));
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(moduleTraitCor),
                 yLabels = rownames(moduleTraitCor),
                 xSymbols = colnames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-meat quality relationships"))
  dev.off()
  write.csv(rownames(moduleTraitCor),"goat_module_color.csv",row.names = F) 
}
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, method = "spearman"));
write.csv(geneModuleMembership,"goat_each_module_MM_value.csv")