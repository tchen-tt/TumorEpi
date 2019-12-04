library(DESeq2)
library(magrittr)

#--------------------- combind gene expression file------------------------
files <- list.files("../", pattern = "*genes.results")
expressions <- lapply(files, function(x) {
  data <- read.table(paste0("../", x), header = TRUE, stringsAsFactors = FALSE,
                     row.names = 1)
  data <- data[, 4, drop = FALSE]
  names(data) <- gsub(pattern = ".genes.results", replacement = "", x)
  data
})
expressions <- as.data.frame(expressions)

#--------------------- start gene expression analysis ---------------------
coldata <- data.frame(row.names = colnames(expressions),
                      condition = rep(c("LLC", "MLE"), each = 4))
expressions <- round(expressions, digits = 0)
expressions <- expressions[rowMeans(expressions) > 1,]
dds <- DESeqDataSetFromMatrix(countData = expressions,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, format = "DataFrame",
               contrast = c("condition", "LLC", "MLE"))
t <- counts(dds, normalized = TRUE)

meads.llc <- rowMeans(t[,1:4])
meads.mle <- rowMeans(t[,5:8])

diff_t <- res %>% as.data.frame %>%
  rownames_to_column("gene_id") %>%
  filter(abs(log2FoldChange)>1 & padj < 0.001)

#------------------- gene annotations using clusterProfilerr---------------
library(clusterProfiler)
library(org.Mm.eg.db)
diffgene <- diff_t$gene_id
gene.df <- bitr(diffgene, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
# outputs
write.csv(gene.df, file = "../diff_gene_ensembl_symnol.csv", row.names = FALSE, quote = FALSE)
ego <- enrichGO(gene = gene.df$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)
#-------------llc cell line and lung norrmal tissue----------------------------
llc <- read.csv("../LLC_cell_line_counts.csv", header = TRUE, stringsAsFactors = FALSE,
                row.names = 1)
lung <- read.csv("../lung_normal_tissue.csv", header = TRUE, stringsAsFactors = FALSE,
                 row.names = 1)
mouse_tissue <- read.table("../GSE138103_annotated_combined.counts.txt",
                         row.names = 1, header = TRUE, stringsAsFactors = FALSE)
llc_lung_gene <- intersect(rownames(llc), rownames(lung))
llc <- llc[llc_lung_gene,]
colnames(llc) <- paste(rep("LLC", 3), 1:3, sep = "_")
colnames(expressions) <- paste(rep(c("LLC", "MLE"), each = 4), c(4:7, 1:4), sep = "_")
lung <- lung[llc_lung_gene,]
cell_line <- cbind.data.frame(llc, expressions[llc_lung_gene, ])
mouse <- cbind.data.frame(cell_line, 
                          mouse_tissue[llc_lung_gene, 1:(ncol(mouse_tissue)-1)])
mouse <- round(mouse, digits = 0)
#----------load GFP LLC cell lines ----------------------------------------------
GFP.file <- list.files("../LLC_GFP/", pattern = "genes.result")
GFP.llc <- lapply(GFP.file, function(x) {
  da <- read.table(paste0("../LLC_GFP/", x), header = TRUE, row.names = 1)
  da <- da[, 4, drop = FALSE]
  colnames(da) <- gsub(pattern = ".genes.results", replacement = "", x)
  da
})
GFP.llc <- as.data.frame(GFP.llc)
GFP.llc <- GFP.llc[llc_lung_gene,]
GFP.llc <- round(GFP.llc, digits = 0)
colnames(GFP.llc) <- paste(rep(c("LLCVitro", "LLCVivo"), c(3, 3)), c(1:3, 1:3), sep = "_")
mouse <- cbind.data.frame(GFP.llc, mouse)
conditions = data.frame(row.names = colnames(mouse),
                        conditions = 1:ncol(mouse))
dds_mouse <- DESeqDataSetFromMatrix(countData = mouse, 
                                    colData = conditions,
                                    design = ~conditions)
# load exome
exomes <- readRDS("../mouse_exons_by_genes.Rds")
exomes <- exomes[rownames(mouse)]
rowRanges(dds_mouse) <- exomes
mouse_fpkm <- fpkm(dds_mouse, robust = FALSE)
mouse_fpkm <- round(mouse_fpkm, digits = 2)


#---------------------loaad gene variant and mhc bing data ------------------
binding <- read.csv("../gene_mhc_I_II_expre.csv", header = TRUE, stringsAsFactors = FALSE)
testing_expr <- mouse_fpkm[unique(binding$gene_id), 1:3]

tt_testing <- log2(mouse_fpkm[unique(binding$gene_id),]+1)
tt_testing[tt_testing > 8] = 8

names <- str_split(colnames(tt_testing), "_", simplify = TRUE)[,c(1:2)] %>% apply(., 1, function(x) paste(x, collapse = "_"))
colnames(tt_testing) <- names
pheatmap::pheatmap(tt_testing, color = rev(brewer.pal(11,"RdBu")[1:6]),
                   show_rownames = FALSE, cluster_cols = TRUE)
express_gene <- unique(binding$gene_id)[rowMeans(tt_testing[unique(binding$gene_id),7:13])>1]
not_expression <- unique(binding$gene_id)[rowMeans(tt_testing[unique(binding$gene_id),7:13])<1]

pheatmap::pheatmap(tt_testing[c(express_gene),], color = rev(brewer.pal(11,"RdBu")[1:6]),
                   show_rownames = FALSE, cluster_cols = TRUE, cluster_rows = TRUE)
h_express <- hclust(dist(tt_testing[express_gene,]))
h_notexpress <- hclust(dist(tt_testing[not_expression,]))

pheatmap::pheatmap(tt_testing[h_notexpress$labels[h_notexpress$order],], color = colorRampPalette(rev(brewer.pal(11,"RdBu")[1:6]))(9),
                   show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
colorRampPalette(rev(brewer.pal(11,"RdBu")[1:6]))(9)

t_testing <- tt_testing %>% 
  as.data.frame %>%
  rownames_to_column("id") %>% 
  gather(., -id, key = "sample", value = "expression") %>%
  mutate(id = factor(id, levels = c(h_express$labels[h_express$order],
                                    h_notexpress$labels[h_notexpress$order])),
         sample = factor(sample, levels = colnames(tt_testing)))
t_testing$class <- str_split(t_testing$sample, pattern = "_", simplify = TRUE)[,1]
t_testing$class <- ifelse(t_testing$class == "LLCVitro", "LLC", t_testing$class)


P.plot <- t_testing %>% 
  group_by(id, class) %>% 
  summarise(mean.expression = mean(expression)) %>% 
  filter(class != "LLCVivo") %>%
  mutate(mean.expression = ifelse(mean.expression>8, 8, mean.expression)) %>%
  mutate(mean.expression = floor(mean.expression)) %>%
  ungroup %>%
  mutate(id = factor(id, levels = c(h_express$labels[h_express$order],
                                    h_notexpress$labels[h_notexpress$order])),
         class = factor(class, levels = c("LLC", "MLE", "Lung", "Brain",
                                          "Colon", "FRT", "Heart", "iLN", 
                                          "Kidney", "Liver", "mesLN",
                                          "PBMC", "SI", "Skin", "Spleen",
                                          "Thymus","BM"))) %>% 
  mutate(mean.expression = as.character(mean.expression)) %>%
  ggplot(., aes(x = class, y = id, fill = mean.expression)) +
  geom_tile() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  scale_fill_manual("log2(FPKM+1)",
                    labels = c("[0, 1)", "[1, 2)", "[2, 3)",
                               "[3, 4)", "[4, 5)", "[5, 6)",
                               "[6, 7)", "[7, 8)",  ">=8"),
                    values = colorRampPalette(rev(brewer.pal(11,"RdBu")[1:6]))(9))
#---------------------not use it -------------------------------------------------
P.plot <- t_testing %>%
  mutate(expression = ifelse(expression>8, 8, expression))%>%
  mutate(expression = floor(expression)) %>%
  mutate(expression = as.character(expression)) %>%
  ggplot(., aes(x = sample, y = id, fill = expression)) +
  geom_tile() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0),
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  scale_fill_manual("log2(FPKM+1)",
                    labels = c("[0, 1)", "[1, 2)", "[2, 3)",
                               "[3, 4)", "[4, 5)", "[5, 6)",
                               "[6, 7)", "[7, 8)",  ">=8"),
    values = colorRampPalette(rev(brewer.pal(11,"RdBu")[1:6]))(9))
#-----------------------------------not use it-----------------------------------
P.build <- ggplot_build(plot = P.plot)
P.build$layout$panel_params[[1]]$x.range
x.pos <- min(P.build$layout$panel_params[[1]]$x.range)
x.min <- x.pos - diff(P.build$layout$panel_params[[1]]$x.range) * 0.025
y.range <- P.build$layout$panel_params[[1]]$x.range
P.plot + annotation_raster(rep(brewer.pal(12, "Set3")[c(2, 4)], c(length(not_expression), length(express_gene))), 
                           xmin = x.min, xmax = x.pos, 
                           ymin = -Inf, ymax = Inf) +
  coord_cartesian(xlim = c(x.min, max(P.build$layout$panel_params[[1]]$x.range)),
                  expand = TRUE)


 #----------------MHC I and II gene express levels ----------------
# H-2D H-2K H-2L(mouse C56BTL/6 none) H-2Ia
mhc.gene <- c("ENSMUSG00000073411", "ENSMUSG00000061232",
              "ENSMUSG00000036594")
mhc.gene.expression <- mouse_fpkm[mhc.gene,]

mhc.expression <- mhc.gene.expression %>%
  as.data.frame %>%
  rownames_to_column("id") %>% 
  gather(-id, key="sample", value = "expression")
mhc.expression$group1 <- str_split(mhc.expression$sample, pattern = "_", simplify = TRUE)[,1]

mhc.expression %>% 
  ggplot(., aes(x = group1, y = expression)) +
  geom_boxplot() + theme_classic() + facet_grid(id ~.)



mhc.expression %>% 
  filter(group1 == "LLC") %>%
  ggplot(., aes(x = group1, y = expression, colour = id)) +
  geom_boxplot() + theme_classic() + geom_point()

mhc.expression %>% group_by(., id, group1) %>%
  summarise(mean.expr = mean(expression),
            len = length(expression),
            sem.expre = sd(expression)/sqrt(length(expression)))
mhc.expression$tissue <- "Mouse"
mhc.expression$tissue <- ifelse(mhc.expression$group1 == "LLC", "LLC", mhc.expression$tissue)
mhc.expression$tissue <- ifelse(mhc.expression$group1 == "MLE", "MLE", mhc.expression$tissue)
mhc.expression$tissue <- ifelse(mhc.expression$group1 == "LLCVitro", "LLC", mhc.expression$tissue)
mhc.expression$tissue <- ifelse(mhc.expression$group1 == "LLCVivo", "LLCVivo", mhc.expression$tissue)
mhc.expression$tissue <- ifelse(mhc.expression$group1 %in% c("iLN", "mesLN"), "lymph node", mhc.expression$tissue)
mhc.expression$tissue <- ifelse(mhc.expression$group1 == "Spleen", "Spleen", mhc.expression$tissue)
mhc.expression$tissue <- ifelse(mhc.expression$group1 == "PBMC", "PBMC", mhc.expression$tissue)


P7 <- mhc.expression %>%
  filter(id == "ENSMUSG00000073411",
         tissue != "LLCVivo") %>%
  mutate(tissue = factor(tissue, levels = c("LLC", "MLE", "PBMC", "lymph node", "Spleen", "Mouse"))) %>%
  ggplot(., aes(x = tissue, y = expression)) +
  geom_boxplot() + 
  theme(axis.title.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, size = 1)) +
  ylab("Gene expression FPKM") + 
  scale_x_discrete(labels = c("LLC", "MLE", "PBMC", "LymphNode", "Spleen", "Tissue"))

P8<-mhc.expression %>%
  filter(id == "ENSMUSG00000061232",
         tissue != "LLCVivo") %>%
  mutate(tissue = factor(tissue, levels = c("LLC", "MLE", "PBMC", "lymph node", "Spleen", "Mouse"))) %>%
  ggplot(., aes(x = tissue, y = expression)) +
  geom_boxplot() + 
  theme(axis.title.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, size = 1)) +
  ylab("Gene expression FPKM") + 
  scale_x_discrete(labels = c("LLC", "MLE", "PBMC", "LymphNode", "Spleen", "Tissue"))

P9<-mhc.expression %>%
  filter(id == "ENSMUSG00000036594",
         tissue != "LLCVivo") %>%
  mutate(tissue = factor(tissue, levels = c("LLC", "MLE", "PBMC", "lymph node", "Spleen", "Mouse"))) %>%
  ggplot(., aes(x = tissue, y = expression)) +
  geom_boxplot() + 
  theme(axis.title.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, size = 1)) +
  ylab("Gene expression FPKM") + 
  scale_x_discrete(labels = c("LLC", "MLE", "PBMC", "LymphNode", "Spleen", "Tissue"))

P10 <- mhc.expression %>%
  filter(group1 %in% c("LLCVitro", "LLCVivo")) %>% 
  group_by(id, group1) %>%
  summarise(mean.expression = mean(expression),
            sem.expression = sd(expression)/sqrt(n()))%>%
  ggplot(., aes(x = id, y = mean.expression, fill = group1)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) + theme_classic() +
  scale_y_continuous(expand = c(0, 2)) + 
  geom_errorbar(aes(ymin=mean.expression-sem.expression, 
                    ymax=mean.expression+sem.expression),
                position = position_dodge(0.9), width = 0.2) +
  scale_x_discrete(labels = c("H2-Aa", "H2-K1", "H2-D1")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic")) +
  ylab("Gene Expression FPKM") +
  scale_fill_discrete("") + theme(legend.position = c(0.2, 0.8)) + 
  scale_fill_manual("",values = brewer.pal(12, "Paired")[c(8, 10)])
files <- list.files(FilePath, pattern = "sorted.sorted.txt")
counts <- lapply(files, function(x) {
  read.table(x, row.names = 1,
             col.names = c("id", gsub("//.sorted//.sorted//.txt", replacement = "", x)))
})
expression_counts <- as.data.frame(counts)
