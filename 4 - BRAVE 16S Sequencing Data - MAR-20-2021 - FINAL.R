# Analysis of BRAVE Kids 16S rRNA sequencing data 

remove(list=ls())
setwd("C:/Users/msk37/Google Drive/Research/BRAVE_Kids/16S_Sequencing/") 
library(phyloseq)
library(dplyr)
library(plyr)
library(data.table)
library(gridExtra)
library(ggplot2)
library(metagenomeSeq)
library(readr)
library(tidyverse)
library(vegan)
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("ANCOM_v2.1/scripts/ancom_v2.1.R")
library(DESeq2)

set.seed(1234)

phy.brave <- readRDS("phy.brave.alex.rds")
nsamples(phy.brave)
ntaxa(phy.brave)

# Number ASVs and store the exact sequences in the reference sequence slot - can access with refseq(ps)
dna <- Biostrings::DNAStringSet(taxa_names(phy.brave))
names(dna) <- taxa_names(phy.brave)
phy.brave <- merge_phyloseq(phy.brave, dna)
taxa_names(phy.brave) <- paste0(seq(ntaxa(phy.brave)))
phy.brave
remove(dna)

# Prune NP samples that don't have at least 1000 reads
phy.brave.1000 <- prune_samples(sample_sums(phy.brave)>=1000, phy.brave)
nsamples(phy.brave.1000)

# Remove samples from children not exposed to or infected with SARS-CoV-2 or with indeterminate testing
phy.brave.1000 <- subset_samples(phy.brave.1000, corona %in% c("Positive", "Negative"))
nsamples(phy.brave.1000)
# Remove nasal samples
phy.pruned <- subset_samples(phy.brave.1000, method_np!="Nasal")
remove(phy.brave, phy.brave.1000)
phy.pruned <- prune_taxa(taxa_sums(phy.pruned) > 0, phy.pruned)
nsamples(phy.pruned)
ntaxa(phy.pruned)

# **********************************************************************************************************************
# CLEANUP THE TAXONOMY TABLE
# **********************************************************************************************************************

# Create dataframe for filtered taxtable
taxtable <- data.frame(as(tax_table(phy.pruned),"matrix"),stringsAsFactors=FALSE)

# Replaces NA with data from higher taxonomic rank (identifies taxonomic rank with letter)
rplc<-which(taxtable$phylum=="" | is.na(taxtable$class))
taxtable$phylum[rplc]<-paste(taxtable$domain[rplc],";p",sep="")
rplc<-which(taxtable$class=="" | is.na(taxtable$class))
taxtable$class[rplc]<-paste(taxtable$phylum[rplc],";c",sep="")
rplc<-which(taxtable$order=="" | is.na(taxtable$order))
taxtable$order[rplc]<-paste(taxtable$class[rplc],";o",sep="")
rplc<-which(taxtable$family=="" | is.na(taxtable$family))
taxtable$family[rplc]<-paste(taxtable$order[rplc],";f",sep="")
rplc<-which(taxtable$genus=="" | is.na(taxtable$genus) | taxtable$genus=="")
taxtable$genus[rplc]<-paste(taxtable$family[rplc],";g",sep="")
rplc<-which(taxtable$genus=="NA;c;o;f;g")
taxtable$genus[rplc]<-"Unassigned"
rplc<-which(taxtable$species=="" | is.na(taxtable$species) | taxtable$species=="")
taxtable$species[rplc]<-paste(taxtable$genus[rplc],";s",sep="")
remove(rplc)
# Use this as replacement for phyloseq object taxonomy table
tax_table(phy.pruned) <- as(taxtable,"matrix")
taxtable <- data.frame(taxtable)
head(tax_table(phy.pruned))

# **********************************************************************************************************************
# CREATE PHYLOSEQ OBJECTS
# **********************************************************************************************************************

# Create phyloseq object and metadata file for analysis of SARS-CoV-2 exposure
phy.brave.cvd1 <- phy.pruned
metadata_cvd1 <- data.frame(sample_data(phy.brave.cvd1))
metadata_cvd1$brave_id <- as.character(metadata_cvd1$brave_id)
metadata_cvd1$age_cat <- as.factor(metadata_cvd1$age_cat)
metadata_cvd1$group2[metadata_cvd1$corona=="Negative"] <- "NEG"
metadata_cvd1$group2[metadata_cvd1$corona=="Positive" & metadata_cvd1$resp_sx_any=="1"] <- "POS_RESP_SX"
metadata_cvd1$group2[metadata_cvd1$corona=="Positive" & metadata_cvd1$resp_sx_any=="0"] <- "POS_ASX"
metadata_cvd1$group2 <- as.factor(metadata_cvd1$group2)
metadata_cvd1$age_cat2[metadata_cvd1$age<3] <- "0-2"
metadata_cvd1$age_cat2[metadata_cvd1$age>=3 & metadata_cvd1$age<6] <- "3-5"
metadata_cvd1$age_cat2[metadata_cvd1$age>=6 & metadata_cvd1$age<9] <- "6-8"
metadata_cvd1$age_cat2[metadata_cvd1$age>=9 & metadata_cvd1$age<12] <- "9-11"
metadata_cvd1$age_cat2[metadata_cvd1$age>=12 & metadata_cvd1$age<15] <- "12-14"
metadata_cvd1$age_cat2[metadata_cvd1$age>=15 & metadata_cvd1$age<18] <- "15-17"
metadata_cvd1$age_cat2[metadata_cvd1$age>=18] <- "18-20"
metadata_cvd1$age_cat2 <- as.factor(metadata_cvd1$age_cat2)
sample_data(phy.brave.cvd1) <- metadata_cvd1
summary(sample_sums(phy.brave.cvd1))
length(unique(metadata_cvd1$brave_id))
ntaxa(phy.brave.cvd1)

# Create file with reference sequences for ASVs
taxa_sums <- data.frame(sort(taxa_sums(phy.brave.cvd1), TRUE))
taxa <- merge(taxa_sums, taxtable, by="row.names", all.x=TRUE)
remove(taxtable)
colnames(taxa)[2] <- "Abundance"
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-1]
taxa_seqs <- data.frame(refseq(phy.brave.cvd1))
taxtable <- merge(taxa, taxa_seqs, by="row.names")
names(taxtable)[names(taxtable)=="Row.names"] <- "ASV"
taxtable <- taxtable[,-2]
write.csv(taxtable, "blast_asvs.csv")
remove(taxa_seqs, taxa_sums, taxa)

# Create phyloseq object and metadata file for analysis of respiratory symptoms in SARS-CoV-2
phy.brave.cvd2 <- subset_samples(phy.brave.cvd1, corona!="Negative")
nsamples(phy.brave.cvd2)
metadata_cvd2 <- data.frame(sample_data(phy.brave.cvd2))
metadata_cvd2$brave_id <- as.character(metadata_cvd2$brave_id)
metadata_cvd2$resp_sx_any <- as.factor(metadata_cvd2$resp_sx_any)
metadata_cvd2$age_cat <- as.factor(metadata_cvd2$age_cat)
sample_data(phy.brave.cvd2) <- metadata_cvd2
summary(sample_sums(phy.brave.cvd2))
length(unique(metadata_cvd2$brave_id))
ntaxa(phy.brave.cvd2)

# Create phyloseq object and metadata file for analysis of respiratory symptoms in SARS-CoV-2 compared to healthy controls
phy.brave.cvd3 <- subset_samples(phy.brave.cvd1, group2!="POS_ASX")
nsamples(phy.brave.cvd3)
metadata_cvd3 <- data.frame(sample_data(phy.brave.cvd3))
metadata_cvd3$brave_id <- as.character(metadata_cvd3$brave_id)
metadata_cvd3$resp_sx_any <- as.factor(metadata_cvd3$resp_sx_any)
metadata_cvd3$age_cat <- as.factor(metadata_cvd3$age_cat)
sample_data(phy.brave.cvd3) <- metadata_cvd3
summary(sample_sums(phy.brave.cvd3))
length(unique(metadata_cvd3$brave_id))
ntaxa(phy.brave.cvd3)

# **********************************************************************************************************************
# FIRST COMPARE DEMOGRAPHICS BY SARS-COV-2 INFECTION STATUS
# **********************************************************************************************************************

table(metadata_cvd1$group2)
tapply(metadata_cvd1$age, metadata_cvd1$group2, summary)
summary(aov(age ~ group2, data=metadata_cvd1))

table(metadata_cvd1$group2, metadata_cvd1$sex)
prop.table(table(metadata_cvd1$group2, metadata_cvd1$sex), 1)
chisq.test(metadata_cvd1$group2, metadata_cvd1$sex)

table(metadata_cvd1$group2, metadata_cvd1$race, useNA="always")
prop.table(table(metadata_cvd1$group2, metadata_cvd1$race), 1)
fisher.test(metadata_cvd1$group2, metadata_cvd1$race)

table(metadata_cvd1$group2, metadata_cvd1$asthma)
prop.table(table(metadata_cvd1$group2, metadata_cvd1$asthma), 1)
chisq.test(metadata_cvd1$group2, metadata_cvd1$asthma)

table(metadata_cvd1$group2, metadata_cvd1$obesity)
prop.table(table(metadata_cvd1$group2, metadata_cvd1$obesity), 1)
chisq.test(metadata_cvd1$group2, metadata_cvd1$obesity)

table(metadata_cvd1$group2, metadata_cvd1$abx_30d)
prop.table(table(metadata_cvd1$group2, metadata_cvd1$abx_30d), 1)
fisher.test(metadata_cvd1$group2, metadata_cvd1$abx_30d)

table(metadata_cvd1$group2, metadata_cvd1$probx_30d)
prop.table(table(metadata_cvd1$group2, metadata_cvd1$probx_30d), 1)
fisher.test(metadata_cvd1$group2, metadata_cvd1$probx_30d)

summary(metadata_cvd2$timing_np_dx)

# **********************************************************************************************************************
# CREATE DATAFRAMES FOR MEASURING ALPHA DIVERSITY PRIOR TO FILTERING
# **********************************************************************************************************************

# Create dataframe with alpha diversity measures for each specimen
diversity_tmp <- estimate_richness(phy.brave.cvd1, measures = c("Shannon", "Chao1"))
setDT(diversity_tmp, keep.rownames = TRUE)[]
diversity_tmp$brave_id <- diversity_tmp$rn
diversity_tmp$rn <- NULL
diversity_tmp$brave_id <- gsub("PCOV1", "PCOV1-", diversity_tmp$brave_id)
diversity_cvd1 <- merge(metadata_cvd1, diversity_tmp, by="brave_id")
remove(diversity_tmp, metadata_cvd1)

diversity_tmp <- estimate_richness(phy.brave.cvd2, measures = c("Shannon", "Chao1"))
setDT(diversity_tmp, keep.rownames = TRUE)[]
diversity_tmp$brave_id <- diversity_tmp$rn
diversity_tmp$rn <- NULL
diversity_tmp$brave_id <- gsub("PCOV1", "PCOV1-", diversity_tmp$brave_id)
diversity_cvd2 <- merge(metadata_cvd2, diversity_tmp, by="brave_id")
remove(diversity_tmp, metadata_cvd2)

diversity_tmp <- estimate_richness(phy.brave.cvd3, measures = c("Shannon", "Chao1"))
setDT(diversity_tmp, keep.rownames = TRUE)[]
diversity_tmp$brave_id <- diversity_tmp$rn
diversity_tmp$rn <- NULL
diversity_tmp$brave_id <- gsub("PCOV1", "PCOV1-", diversity_tmp$brave_id)
diversity_cvd3 <- merge(metadata_cvd3, diversity_tmp, by="brave_id")
remove(diversity_tmp, metadata_cvd3)

# **********************************************************************************************************************
# ANALYSES BY AGE
# **********************************************************************************************************************

shapiro.test(diversity_cvd1$Shannon)
qqnorm(diversity_cvd1$Shannon); qqline(diversity_cvd1$Shannon)
# Shannon diversity is approximately normally distributed
summary(lm(Shannon ~ age + group2, data=diversity_cvd1))
# Increasing Shannon diversity with increasing age (P<0.0001)

shapiro.test(diversity_cvd1$Chao1)
qqnorm(diversity_cvd1$Chao1); qqline(diversity_cvd1$Chao1)
diversity_cvd1$Chao1_log <- log(diversity_cvd1$Chao1)
summary(lm(Chao1_log ~ age + group2, data=diversity_cvd1))
# Decreasing Chao1 richness with increasing age (p=0.02)

adonis(phyloseq::distance(phy.brave.cvd1, method="bray") ~ age + group2,
       data = diversity_cvd1)
# Overall microbiome composition differs by age (P<0.001)

# *****************
# ANCOM-II ANALYSES
# *****************

meta_data <- get_variable(phy.brave.cvd1)
otu_data <- data.frame(otu_table(phy.brave.cvd1))
otu_id <- row.names(otu_data)
otu_id <- gsub("PCOV1", "PCOV1-", otu_id)
otu_data <- data.frame(t(otu_table(phy.brave.cvd1)))
colnames(otu_data) <- otu_id

feature_table = otu_data; sample_var = "brave_id"; group_var = "age_cat"
out_cut = 0.05; zero_cut = 0.95; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "age"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "group2"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
output <- data.frame(res$out)
output <- output[order(-output$W),]
write_csv(output, "Results_ancom2/age_ancom2.csv")
remove(adj_formula, alpha, group_var, lib_cut, main_var, neg_lb, otu_id, out_cut, p_adj_method, rand_formula, sample_var,
       struc_zero, zero_cut, res, feature_table, meta_data, otu_data, prepro, output)

age_ancom2 <- read.csv("Results_ancom2/age_ancom2.csv")
age_ancom2 <- subset(age_ancom2, detected_0.8=="TRUE")
age_ancom2 <- subset(age_ancom2, W!="Inf")  # these represent structural zeros
names(age_ancom2)[names(age_ancom2)=="taxa_id"] <- "ASV"
age_ancom2 <- merge(taxtable, age_ancom2, by="ASV")
age_ancom2 <- age_ancom2[,-c(9,13,14)]

phy.deseq.cvd1 <- phy.brave.cvd1
filter <- phyloseq::genefilter_sample(phy.deseq.cvd1, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd1))
phy.deseq.cvd1 <- prune_taxa(filter, phy.deseq.cvd1) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd1, ~ group2 + age)
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd1)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
age_ancom2 <- merge(age_ancom2, sigtab, by="ASV")
names(age_ancom2)[names(age_ancom2)=="log2FoldChange"] <- "DESeq2_log2FoldChange"
remove(filter, phy.deseq.cvd1, ds, dds, res, sigtab)
write.csv(age_ancom2, "Results_ancom2/age_ancom2_final.csv")

# **********************************************************************************************************************
# ANALYSES BY CORONA
# **********************************************************************************************************************

tapply(diversity_cvd1$Shannon, diversity_cvd1$corona, summary)
wilcox.test(diversity_cvd1$Shannon~diversity_cvd1$corona)
summary(lm(Shannon ~ corona + age, data=diversity_cvd1))
# No difference in Shannon diversity by SARS-CoV-2 status (p=0.36)

tapply(diversity_cvd1$Chao1, diversity_cvd1$corona, summary)
diversity_cvd1$log_Chao1 <- log(diversity_cvd1$Chao1)
wilcox.test(diversity_cvd1$log_Chao1~diversity_cvd1$corona)
summary(lm(log_Chao1 ~ corona + age, data=diversity_cvd1))
# Higher Chao1 richness in children with SARS-CoV-2 infection adjusting for age (p=0.008)

adonis(phyloseq::distance(phy.brave.cvd1, method="bray") ~ corona + age,
       data = diversity_cvd1)
# Overall microbiome composition differs by SARS-CoV-2 status (p=0.007)

# *****************
# ANCOM-II ANALYSES
# *****************

meta_data <- get_variable(phy.brave.cvd1)
otu_data <- data.frame(otu_table(phy.brave.cvd1))
otu_id <- row.names(otu_data)
otu_id <- gsub("PCOV1", "PCOV1-", otu_id)
otu_data <- data.frame(t(otu_table(phy.brave.cvd1)))
colnames(otu_data) <- otu_id

feature_table = otu_data; sample_var = "brave_id"; group_var = "corona"
out_cut = 0.05; zero_cut = 0.95; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "corona"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "age"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
output <- data.frame(res$out)
output <- output[order(-output$W),]
write_csv(output, "Results_ancom2/corona_ancom2.csv")
remove(adj_formula, alpha, group_var, lib_cut, main_var, neg_lb, otu_id, out_cut, p_adj_method, rand_formula, sample_var,
       struc_zero, zero_cut, res, feature_table, meta_data, otu_data, prepro, output)

corona_ancom2 <- read.csv("Results_ancom2/corona_ancom2.csv")
corona_ancom2 <- subset(corona_ancom2, detected_0.8=="TRUE")
corona_ancom2 <- subset(corona_ancom2, W!="Inf")  # these represent structural zeros
names(corona_ancom2)[names(corona_ancom2)=="taxa_id"] <- "ASV"
corona_ancom2 <- merge(taxtable, corona_ancom2, by="ASV")
corona_ancom2 <- corona_ancom2[,-c(9,13,14)]

ASVs_all <- c("1030")
phy.deseq.cvd1 <- phy.brave.cvd1
filter <- phyloseq::genefilter_sample(phy.deseq.cvd1, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd1))
phy.deseq.cvd1 <- prune_taxa(filter, phy.deseq.cvd1) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd1, ~ age_cat2 + corona)
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd1)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab1 <- subset(sigtab, ASV %in% ASVs_all)

ASVs_9yo <- c("1017")
phy.deseq.cvd2 <- subset_samples(phy.brave.cvd1, age>=9)
filter <- phyloseq::genefilter_sample(phy.deseq.cvd1, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd1))
phy.deseq.cvd1 <- prune_taxa(filter, phy.deseq.cvd1) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd1, ~ age + corona)
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd1)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab2 <- subset(sigtab, ASV %in% ASVs_9yo)
remove(filter, phy.deseq.cvd2, ds, dds, res, sigtab)

sigtab_final <- rbind(sigtab1, sigtab2)
corona_ancom2 <- merge(corona_ancom2, sigtab_final, by="ASV")
names(corona_ancom2)[names(corona_ancom2)=="log2FoldChange"] <- "DESeq2_log2FoldChange"
write.csv(corona_ancom2, "C:/Users/msk37/Google Drive/Research/BRAVE_Kids/16S_Sequencing/Results_ancom2/corona_ancom2_final.csv")

# **********************************************************************************************************************
# ANALYSES BY RESP_SX_ANY
# **********************************************************************************************************************

tapply(diversity_cvd2$Shannon, diversity_cvd2$resp_sx_any, summary)
wilcox.test(diversity_cvd2$Shannon~diversity_cvd2$resp_sx_any)
summary(lm(Shannon ~ resp_sx_any + age, data=diversity_cvd2))
# No difference in Shannon diversity by respiratory symptoms (p=0.42)

tapply(diversity_cvd2$Chao1, diversity_cvd2$resp_sx_any, summary)
diversity_cvd2$log_Chao1 <- log(diversity_cvd2$Chao1)
wilcox.test(diversity_cvd2$log_Chao1~diversity_cvd2$resp_sx_any)
summary(lm(log_Chao1 ~ resp_sx_any + age, data=diversity_cvd2))
# No difference in Chao1 richness by respiratory symptoms (p=0.14)

adonis(phyloseq::distance(phy.brave.cvd2, method="bray") ~ resp_sx_any + age,
       data = diversity_cvd2)
# Overall microbiome composition differs by respiratory symptoms (p=0.008)

# *****************
# ANCOM-II ANALYSES
# *****************

meta_data <- get_variable(phy.brave.cvd2)
otu_data <- data.frame(otu_table(phy.brave.cvd2))
otu_id <- row.names(otu_data)
otu_id <- gsub("PCOV1", "PCOV1-", otu_id)
otu_data <- data.frame(t(otu_table(phy.brave.cvd2)))
colnames(otu_data) <- otu_id

feature_table = otu_data; sample_var = "brave_id"; group_var = "resp_sx_any"
out_cut = 0.05; zero_cut = 0.95; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "resp_sx_any"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "age"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
output <- data.frame(res$out)
output <- output[order(-output$W),]
write_csv(output, "Results_ancom2/resp_ancom2.csv")
remove(adj_formula, alpha, group_var, lib_cut, main_var, neg_lb, otu_id, out_cut, p_adj_method, rand_formula, sample_var,
       struc_zero, zero_cut, res, feature_table, meta_data, otu_data, prepro, output)

resp_ancom2 <- read.csv("Results_ancom2/resp_ancom2.csv")
resp_ancom2 <- subset(resp_ancom2, detected_0.8=="TRUE")
resp_ancom2 <- subset(resp_ancom2, W!="Inf")  # these represent structural zeros
names(resp_ancom2)[names(resp_ancom2)=="taxa_id"] <- "ASV"
resp_ancom2 <- merge(taxtable, resp_ancom2, by="ASV")
resp_ancom2 <- resp_ancom2[,-c(9,13,14)]

ASVs_all <- c("1014", "1030", "1491")
phy.deseq.cvd2 <- phy.brave.cvd2
filter <- phyloseq::genefilter_sample(phy.deseq.cvd2, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd2))
phy.deseq.cvd2 <- prune_taxa(filter, phy.deseq.cvd2) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd2, ~ age_cat2 + resp_sx_any)
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd2)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab1 <- subset(sigtab, ASV %in% ASVs_all)
remove(filter, phy.deseq.cvd2, ds, dds, res, sigtab)

ASVs_9yo <- c("1017", "1032")
phy.deseq.cvd2 <- subset_samples(phy.brave.cvd2, age>=9)
filter <- phyloseq::genefilter_sample(phy.deseq.cvd2, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd2))
phy.deseq.cvd2 <- prune_taxa(filter, phy.deseq.cvd2) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd2, ~ age_cat2 + resp_sx_any)
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd2)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab2 <- subset(sigtab, ASV %in% ASVs_9yo)
remove(filter, phy.deseq.cvd2, ds, dds, res, sigtab)

ASVs_12yo <- c("1060", "1535", "1593", "1618")
phy.deseq.cvd2 <- subset_samples(phy.brave.cvd2, age>=12)
nsamples(phy.deseq.cvd2)
filter <- phyloseq::genefilter_sample(phy.deseq.cvd2, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd2))
phy.deseq.cvd2 <- transform_sample_counts(phy.deseq.cvd2, function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd2, ~ age + resp_sx_any)
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd2)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab3 <- subset(sigtab, ASV %in% ASVs_12yo)

sigtab_final <- rbind(sigtab1, sigtab2, sigtab3)
resp_ancom2 <- merge(resp_ancom2, sigtab_final, by="ASV")
names(resp_ancom2)[names(resp_ancom2)=="log2FoldChange"] <- "DESeq2_log2FoldChange"
remove(filter, phy.deseq.cvd2, ds, dds, res, sigtab, sigtab1, sigtab2, sigtab3, sigtab_final, ASVs_all, ASVs_12yo)
write.csv(resp_ancom2, "Results_ancom2/resp_ancom2_final.csv")

# **********************************************************************************************************************
# ANALYSES COMPARING SARS-COV-2+ PARTICIPANTS WITH RESPIRATORY SYMPTOMS TO UNINFECTED PARTICIPANTS (CONTROLS)
# **********************************************************************************************************************

tapply(diversity_cvd3$Shannon, diversity_cvd3$corona, summary)
wilcox.test(diversity_cvd3$Shannon~diversity_cvd3$corona)
summary(lm(Shannon ~ corona + age, data=diversity_cvd3))
# No difference in Shannon diversity by SARS-CoV-2 status (p=0.55)

tapply(diversity_cvd3$Chao1, diversity_cvd3$corona, summary)
diversity_cvd3$log_Chao1 <- log(diversity_cvd3$Chao1)
wilcox.test(diversity_cvd3$log_Chao1~diversity_cvd3$corona)
summary(lm(log_Chao1 ~ corona + age, data=diversity_cvd3))
# Higher Chao1 richness in children with SARS-CoV-2 infection with respiratory symptoms adjusting for age (p=0.008)

adonis(phyloseq::distance(phy.brave.cvd3, method="bray") ~ corona + age,
       data = diversity_cvd3)
# Overall microbiome composition differs by respiratory symptoms (p=0.002)

# *****************
# ANCOM-II ANALYSES
# *****************

meta_data <- get_variable(phy.brave.cvd3)
otu_data <- data.frame(otu_table(phy.brave.cvd3))
otu_id <- row.names(otu_data)
otu_id <- gsub("PCOV1", "PCOV1-", otu_id)
otu_data <- data.frame(t(otu_table(phy.brave.cvd3)))
colnames(otu_data) <- otu_id

feature_table = otu_data; sample_var = "brave_id"; group_var = "corona"
out_cut = 0.05; zero_cut = 0.95; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "corona"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "age"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
output <- data.frame(res$out)
output <- output[order(-output$W),]
write_csv(output, "Results_ancom2/controls_ancom2.csv")
remove(adj_formula, alpha, group_var, lib_cut, main_var, neg_lb, otu_id, out_cut, p_adj_method, rand_formula, sample_var,
       struc_zero, zero_cut, res, feature_table, meta_data, otu_data, prepro, output)

controls_ancom2 <- read.csv("Results_ancom2/controls_ancom2.csv")
controls_ancom2 <- subset(controls_ancom2, detected_0.8=="TRUE")
controls_ancom2 <- subset(controls_ancom2, W!="Inf")  # these represent structural zeros
names(controls_ancom2)[names(controls_ancom2)=="taxa_id"] <- "ASV"
controls_ancom2 <- merge(taxtable, controls_ancom2, by="ASV")
controls_ancom2 <- controls_ancom2[,-c(9,13,14)]

ASVs_all <- c("1030")
phy.deseq.cvd3 <- phy.brave.cvd3
filter <- phyloseq::genefilter_sample(phy.deseq.cvd3, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd3))
phy.deseq.cvd3 <- prune_taxa(filter, phy.deseq.cvd3) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd3, ~ age_cat2 + corona)
colData(ds)$corona <- factor(colData(ds)$corona,
                                  levels=c("Negative","Positive"))
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd3)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab1 <- subset(sigtab, ASV %in% ASVs_all)
remove(filter, phy.deseq.cvd2, ds, dds, res, sigtab)

ASVs_9yo <- c("1017", "1032")
phy.deseq.cvd3 <- subset_samples(phy.brave.cvd3, age>=9)
nsamples(phy.deseq.cvd3)
filter <- phyloseq::genefilter_sample(phy.deseq.cvd3, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd3))
phy.deseq.cvd3 <- prune_taxa(filter, phy.deseq.cvd3) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd3, ~ age + corona)
colData(ds)$corona <- factor(colData(ds)$corona,
                             levels=c("Negative","Positive"))
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd3)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab2 <- subset(sigtab, ASV %in% ASVs_9yo)

ASVs_12yo <- c("1060", "1582", "1618")
phy.deseq.cvd3 <- subset_samples(phy.brave.cvd3, age>=12)
nsamples(phy.deseq.cvd3)
filter <- phyloseq::genefilter_sample(phy.deseq.cvd3, filterfun_sample(function(x) x >= 1), A = 0.05*nsamples(phy.deseq.cvd3))
phy.deseq.cvd3 <- prune_taxa(filter, phy.deseq.cvd3) %>%
    transform_sample_counts(., function(x) x+1)
ds <- phyloseq_to_deseq2(phy.deseq.cvd3, ~ age + corona)
colData(ds)$corona <- factor(colData(ds)$corona,
                             levels=c("Negative","Positive"))
dds <- DESeq(ds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- cbind(as(res, "data.frame"), 
                as(tax_table(phy.deseq.cvd3)[rownames(res), ], "matrix")) 
sigtab$ASV <- rownames(sigtab)
sigtab <- sigtab[,c(14,2)]
sigtab3 <- subset(sigtab, ASV %in% ASVs_12yo)
sigtab_final <- rbind(sigtab1, sigtab2, sigtab3)
controls_ancom2 <- merge(controls_ancom2, sigtab_final, by="ASV")
names(controls_ancom2)[names(controls_ancom2)=="log2FoldChange"] <- "DESeq2_log2FoldChange"
remove(filter, phy.deseq.cvd3, ds, dds, res, sigtab, sigtab1, sigtab2, sigtab3, sigtab_final, ASVs_all, ASVs_9yo, ASVs_12yo)
write.csv(controls_ancom2, "Results_ancom2/controls_ancom2_final.csv")