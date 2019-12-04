# Created on: 2018/08/27
# Author:chentao
rm(list=ls())
setwd("/home/taotao/work/project_kagnti/exom/result/")
library(tidyverse)
library(magrittr)
library(doParallel)
library(foreach)
library(biomaRt)

ensembl=useMart(biomart = "ensembl")
#listDatasets(ensembl)
mart <- useMart(biomart = "ensembl",dataset="mmusculus_gene_ensembl")
getseq <- function(x){
  acids <- intersect_test$vars[x] %>% str_extract_all('[a-zA-Z]{3}') %>% .[[1]]
  locat <- intersect_test$vars[x] %>% str_extract('[0-9]+') %>% .[[1]] %>% as.numeric
  POS <- intersect_test$POS[x]
  x=intersect_test$trans_id[x]
  seq=getSequence(id=x,
                  mart=mart,
                  type = "ensembl_transcript_id",
                  seqType = "peptide")
  snv <- seq$peptide %>% substring(locat, locat)
  if(length(snv)!=0){
    if(acids[1] %in% acid$from & acids[2] %in% acid$from){
      if(acid[acid$from==acids[1],]$to==snv){
        peptite_ref = seq$peptide %>% substring(locat-13, locat+13)
        peptite_var = str_c(substring(seq$peptide,locat-13,locat-1), acid[acid$from==acids[2],]$to,
                            substring(seq$peptide,locat+1,locat+13))
        result <- data.frame(POS=POS,trans_id=x, 
                            peptite_ref=peptite_ref, peptite_var=peptite_var)
        return(result)
      }
    }
  }
}

vcf1 <- read.table("sample_rs1.vcf", header=FALSE, stringsAsFactors=FALSE)
vcf1 <- vcf1 %>% as.tibble %>% dplyr::select(V1:V8) %>% filter(V1!="MT")

vcf2 <- read.table("sample_rs2.vcf", header=FALSE, stringsAsFactors=FALSE)
vcf2 <- vcf2 %>% as.tibble %>% dplyr::select(V1:V8) %>% filter(V1!="MT")

vcf3 <- read.table("sample_rs3.vcf", header=FALSE, stringsAsFactors=FALSE)
vcf3 <- vcf3 %>% as.tibble %>% dplyr::select(V1:V8) %>% filter(V1!="MT")

acid <- read.csv("acid.csv", header=TRUE, stringsAsFactors=FALSE)


cl <- makeCluster(20)
registerDoParallel(cl)

func <- function(x,vcf){
  pre <- vcf[x,] %>% dplyr::select(1:5)
  t <- vcf[x,] %>% dplyr::select(8) %>% as.vector %>% str_split(";") %>% .[[1]]
  ann <- t %>% str_detect("protein_coding")
  if(sum(ann)>0){
    t <- t[ann] %>% str_split(',') %>% .[[1]]
    pro_cod <- str_detect(t, "protein_coding")
    t <- t[pro_cod]
    gene_id <- str_extract(t, "ENSMUSG[0-9]{11}")
    trans_id <- str_extract(t, "ENSMUST[0-9]{11}")
    vars <- str_extract(t, "[a-z]{1}\\.[a-zA-Z]{3}[0-9]+[a-zA-Z]{3}")
    rr <- data.frame(gene_id,trans_id,vars)
    n <- sum(pro_cod)
    if(n > 0){
      pre <- pre %>% as.vector %>% list %>% rep(each=n) %>% unlist %>% matrix(nrow=n, byrow=TRUE)
      pre <- pre %>% data.frame
      colnames(pre) <- c("CHROM","POS","ID","REF","ALT")
      rr <- cbind.data.frame(pre, rr)
      return(rr)
    }
  }
}








test1 <- foreach(x=1:nrow(vcf1), .combine="rbind", .packages=c("magrittr","dplyr","stringr")) %dopar% func(x,vcf=vcf1)
test2 <- foreach(x=1:nrow(vcf2), .combine="rbind", .packages=c("magrittr","dplyr","stringr")) %dopar% func(x,vcf=vcf2)
test3 <- foreach(x=1:nrow(vcf3), .combine="rbind", .packages=c("magrittr","dplyr","stringr")) %dopar% func(x,vcf=vcf3)



test1_1 <- test1 %>% as.tibble %>% unite(index, POS, trans_id, sep='_')
test2_2 <- test2 %>% as.tibble %>% unite(index, POS, trans_id, sep='_')
test3_3 <- test3 %>% as.tibble %>% unite(index, POS, trans_id, sep='_')


all_index <- intersect(intersect(test1_1$index, test2_2$index), test3_3$index)
intersect_test <- test1_1 %>% filter(index %in% all_index)

intersect_test <- intersect_test %>% separate(index, into=c("POS","trans_id")) %>% na.omit
#peptite <- lapply(intersect_test$trans_id, getseq)
peptite <- foreach(x=1:nrow(intersect_test), .combine="rbind", .packages=c("magrittr","dplyr","stringr","biomaRt")) %dopar% getseq(x)

peptite1 <- peptite %>% unite(index, POS, trans_id, sep='_')
intersect_test1 <- intersect_test %>% unite(index, POS, trans_id, sep="_")
result_test <- left_join(intersect_test1, peptite1, by="index")


result_test1 <- na.omit(result_test)
result_test1 <- apply(result_test1, 2, as.character) %>% as.tibble

write.csv(result_test1, "result_trans_varpeptite.csv", row.names=FALSE, quote=FALSE)

result_test11 <- result_test1 %>%
  separate(index, into=c("POS","trans_id"), sep="_") %>%
  unite(gene_id_pep_var,gene_id, peptite_var, sep="_")


result_test2 <- result_test11 %>% 
  dplyr::select(CHROM,POS,REF,ALT,peptite_ref,gene_id_pep_var) %>%
  unique

pep_gene_da <- result_test11 %>% filter(gene_id_pep_var %in% result_test2$gene_id_pep_var) %>%
  dplyr::select(gene_id_pep_var, ID) %>%
  unique


tt <-left_join(result_test2,pep_gene_da) %>%
  separate(gene_id_pep_var, into=c("gene_id","peptite_var"), sep = "_") %>%
  dplyr::select("CHROM", "POS", "ID", "gene_id", "REF", "ALT", "peptite_ref", "peptite_var")


write.csv(tt, "result_var_peptite.csv", row.names=FALSE, quote=F)

transcrips <- read.csv("/home/taotao/work/project_kagnti/transcript/cleandata/ReadCount_Fpkm_value.csv",header=T,
                       stringsAsFactors=FALSE)
trans_exom_combine <- transcrips %>% 
  transmute(gene_id, gene_name, gene_biotype,mean.fpkm=(LLC4.fpkm+LLC5.fpkm+LLC4.fpkm)/3) %>% 
  filter(gene_biotype=="protein_coding") %>% 
  merge(tt) %>% 
  dplyr::select(gene_id,gene_name,CHROM,POS,ID,REF,ALT,peptite_ref,peptite_var,mean.fpkm)
write.csv(trans_exom_combine,"trans_exom_combine.csv",
          row.names=FALSE)
stopCluster(cl)

write_fasta <- function(x){
  line <- paste0(">",x,"\n", x,"\n")
  cat(line,file="peptite.fasta", append=TRUE)
}
