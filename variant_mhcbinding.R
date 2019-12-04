vcfFile <- read.table("/Users/chentao/Desktop/mouse_LLC_MLE/immue_paper/extract.vcf",
                      header = FALSE, stringsAsFactors = FALSE)
binding <- read.table("/Users/chentao/Desktop/mouse_LLC_MLE/immue_paper/trans_exom_ep.csv",
                      header = TRUE, stringsAsFactors = FALSE,
                      sep = ",")
vcfFile[1, 10] %>% str_split(":") %>% unlist %>% .[3] %>% as.numeric
vcfFile$mutation_type <- lapply(vcfFile$V8, function(x) {
  x %>% str_split("\\|") %>% unlist %>% .[2]
}) %>% unlist
vcfFile$gene_type <- lapply(vcfFile$V8, function(x) {
  x %>% str_split("\\|") %>% unlist %>% .[8]
}) %>% unlist
vcfFile$AF1 <- lapply(vcfFile$V10, function(x) {
  x %>% str_split(":") %>% unlist %>% .[3] %>% as.numeric
}) %>% unlist
vcfFile$AF2 <- lapply(vcfFile$V11, function(x) {
  x %>% str_split(":") %>% unlist %>% .[3] %>% as.numeric
}) %>% unlist

#vcfFile$AF <- lapply(1:nrow(vcfFile), function(x) {
mean(vcfFile[x, c(18, 19), drop = TRUE], na.rm = T)
})

vcfFile$SNV_indel <- ifelse(nchar(vcfFile$V4) == nchar(vcfFile$V5), "SNP", "InDel")

result <- vcfFile[, c(1, 2, 4, 5, 8, 16, 17, 18, 20)]
immu_result <- list()
immu_result$mutation_summary <- table(result$mutation_type, result$gene_type) %>% as.matrix

# single nucleotide mutations
t <- data.frame(type = c("UTRs", "synonymous", "missense", "premature stop"),
                Number = c(651, 1608, 1883, 71))
t$type <- factor(t$type, levels = t$type)
library(ggplot2)
windowsFonts(Arial = windowsFont("Arial"))
P1 <- ggplot(t, aes(x = type, y = Number)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  ylab("Number of variant") + scale_y_continuous(expand = c(0, 20))
P2 <- vcfFile %>% group_by(SNV_indel) %>% dplyr::count() %>% 
  ggplot(., aes( x = SNV_indel, y = n)) +
  geom_bar(stat = "identity") + 
  ylab("Number of variant") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_y_continuous(expand = c(0, 60)) +
  xlab("")
plot_grid(P2, P1, labels = c("A", "B"))

#---------------------extract missvariant and protein coding ---------------------------------
result_mis <- result[result$mutation_type == "missense_variant" & result$gene_type == "protein_coding",]
mis_index <- paste(result_mis$V1, result_mis$V2, sep = "_")
bin_index <- paste(binding$CHROM, binding$POS, sep = "_")
binding_extract <- binding[bin_index %in% mis_index, ]

#---------------------build fasta file papre for IEDB data base------------------------------
#MHC I
for(i in 1:nrow(binding_extract)) {
  paste0(">", binding_extract[i,]$MHCI) %>%
    paste(., str_to_upper(binding_extract[i,]$MHCI),
          sep = "\n") %>% 
    cat(., "\n", file = "/Users/taochen/Desktop/immue_paper/mouse_mhci.fa",
        append = TRUE)
}
#MHC II
for(i in 1:nrow(binding_extract)) {
  paste0(">", binding_extract[i,]$MHCII) %>%
    paste(., str_to_upper(binding_extract[i,]$MHCII),
          sep = "\n") %>% 
    cat(., "\n", file = "/Users/taochen/Desktop/immue_paper/mouse_mhcii.fa",
        append = TRUE)
}

class1 <- read.csv("/Users/chentao/Desktop/mouse_LLC_MLE/immue_paper/mhci.csv",
                   header = TRUE, stringsAsFactors = FALSE)
class2 <- read.csv("/Users/chentao/Desktop/mouse_LLC_MLE/immue_paper/predict_result.csv",
                   header = TRUE, stringsAsFactors = FALSE)
get.min.rank <- function(alleles, binding.data) {
  if(!(alleles %in% c("H-2-Db", "H-2-Kb", "H2-IAb"))) {
    sotp("check your inputs alleles.")
  }
  data <- binding.data %>% 
    filter(., allele %in% alleles) %>%
    group_by(seq_num) %>% 
    top_n(n = 1, wt = - Percentile.Rank) %>%
    ungroup
  tt <- lapply(unique(data$seq_num), function(x) {
    dd <- data[data$seq_num == x,]
    dd[1, "peptide"] = paste(dd$peptide, collapse = ";")
    dd[1,]
  })
  do.call(rbind, tt)
}f

class1.H2Db <- get.min.rank(alleles = "H-2-Db", binding.data = class1)
class1.H2Kb <- get.min.rank(alleles = "H-2-Kb", binding.data = class1)
class1.bind <- left_join(class1.H2Db, class1.H2Kb, by = "seq_num",
                         suffix = c(".H2Db", ".H2Kb"))


get.min.rank2 <- function(alleles, binding.data) {
  if(!(alleles %in% c("H-2-Db", "H-2-Kb", "H2-IAb"))) {
    sotp("check your inputs alleles.")
  }
  data <- binding.data %>% 
    filter(., allele %in% alleles) %>%
    group_by(seq_num) %>% 
    top_n(n = 1, wt = - percentile_rank) %>%
    ungroup
  tt <- lapply(unique(data$seq_num), function(x) {
    dd <- data[data$seq_num == x,]
    dd[1, "peptide"] = paste(dd$peptide, collapse = ";")
    dd[1,]
  })
  do.call(rbind, tt)
}
class2.H2IAb <- get.min.rank2(alleles = "H2-IAb", binding.data = class2)
mhc <- left_join(class1.bind, class2.H2IAb, by = "seq_num")

library(ggtern)
ggtern(mhc, aes(x = Percentile.Rank.H2Db, y = Percentile.Rank.H2Kb,
                z = percentile_rank, colour = mark)) + geom_point(size = 0.7) +
  theme_arrowdefault() + xlab("H2-Db") +
  ylab("H2-Kb") + zlab("H2-IAb") + 
  scale_colour_discrete(labels = c("A", "B", "C", "D")) + 
  theme_classic() + 
  scale_colour_manual(values = c(brewer.pal(11, "RdBu")[6], brewer.pal(12, "Paired")[c(8, 10, 12)]))
mhc$mark <- 0
mhc$mark <- ifelse(mhc$percentile_rank < 3, mhc$mark + 1, mhc$mark)
mhc$mark <- ifelse(mhc$Percentile.Rank.H2Db < 3, mhc$mark + 2, mhc$mark)
mhc$mark <- ifelse(mhc$Percentile.Rank.H2Kb < 3, mhc$mark + 4, mhc$mark)
mhc$mark <- ifelse(mhc$mark < 3, 0, mhc$mark)
mhc$mark <- ifelse(mhc$mark == 6, 0, mhc$mark)
mhc$mark <- ifelse(mhc$mark == 4, 0, mhc$mark)
mhc$mark <- as.character(mhc$mark)


#------------------------------------- annotation mutation genes--------------------------
library(clusterProfiler)
library(org.Mm.eg.db)

gene <- binding_extract$gene_id
gene.df <- bitr(gene, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
mutation.gene <- binding_extract$gene_id
ggo <- groupGO(gene = gene.df$ENTREZID,
               OrgDb = org.Mm.eg.db,
               ont = "cc",
               level = 3,
               readable = TRUE)
tumor_path_way <- c("EGFR tyrosine kinase inhibitor resistance",
                    "cGMP-PKG signaling pathway",
                    "VEGF signaling pathway",
                    "MAPK signaling pathway",
                    "PI3K-Akt signaling pathway",
                    "cAMP signaling pathway",
                    "NF-kappa B signaling pathway",
                    "p53 signaling pathway",
                    "Wnt signaling pathway",
                    "DNA replication",
                    "Cell cycle",
                    "HIF-1 signaling pathway",
                    "TGF-beta signaling pathway",
                    "Notch signaling pathway")




ggos <- enrichGO(gene = gene.df$ENTREZID,
                 keyType = "ENTREZID",
                 ont = "CC",
                 readable = TRUE,
                 OrgDb = org.Mm.eg.db,
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 1)
ggosMF <- enrichGO(gene = gene.df$ENTREZID,
                   keyType = "ENTREZID",
                   ont = "MF",
                   readable = TRUE,
                   OrgDb = org.Mm.eg.db,
                   pvalueCutoff = 0.5,
                   qvalueCutoff = 1)

kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = "mmu",
                 pvalueCutoff = 0.5
)
barplot(ggo)
kegg_enrich <- kk@result

P3 <- kegg_enrich[kegg_enrich$Description %in% tumor_path_way,] %>% 
  ggplot(., aes(x = reorder(Description, -Count), y = Count)) + 
  geom_bar(stat = "identity") + coord_flip() +
  theme_classic() + scale_y_continuous(expand = c(0,0.2)) +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  ylab("Number of gene")

gocc <- ggos@result

gocc$type <- "CC"
goMF <- ggosMF@result
goMF$type <- "MF"
go <- rbind(gocc, goMF)
go$type <- factor(go$type, levels = c("MF", "CC"))


gos <- go %>% dplyr::group_by(type) %>%
  dplyr::top_n(n = 7, wt = Count) %>% 
  arrange(., type, Count)
gos$Description <- factor(gos$Description, levels = gos$Description)
P4 <- gos%>%
  ggplot(., aes(x = Description, y = Count)) +
  geom_bar(stat = "identity") + 
  coord_flip() + theme_classic() + 
  scale_y_continuous(expand = c(0, 0.4)) + 
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  ylab("Number of gene")
Pf1<- plot_grid(P2, P3, P1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)), P4, 
          labels = c("A", "C", "B", "D"), nrow = 2,
          rel_widths = c(1, 2))
plot_grid(P2, 
          P1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)), nrow = 2, labels = c("A", "B"))


result <- mhc %>% arrange(seq_num)
result <- cbind(binding_extract, result)
write.csv(result, file = "/Users/taochen/Desktop/immue_paper/gene_mhc_I_II_expre.csv",
          row.names = F, quote = F)
tiff(filename = "ggtern_figure2_ternF.tiff", width = 700*5, height = 600*5, res = 72*5)
ggtern(mhc, aes(x = Percentile.Rank.H2Db, y = Percentile.Rank.H2Kb,
                z = percentile_rank, colour = mark)) + geom_point(size = 0.7) +
  theme_arrowdefault() + xlab("H2-Db") + geom_point(size = 0.5) +
  ylab("H2-Kb") + zlab("H2-IAb") + 
  theme_classic() + 
  theme(legend.title  = element_blank()) +
  scale_colour_manual(labels =c("all greater than 3", "H2-IAb and H2-Db are less than 3",
                        "H2-IAb and H2-Kb are less than 3",
                        "all less than 3"),values = c(brewer.pal(11, "RdBu")[6], brewer.pal(12, "Paired")[c(8, 10, 6)]))
dev.off()

ggtern() 
  geom_point(data = mhc[mhc$mark == "0"], aes(x = Percentile.Rank.H2Db, y = Percentile.Rank.H2Kb,
                                              z = percentile_rank))


