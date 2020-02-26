library(dplyr)
library(magrittr)
library(pheatmap)
library(sleuth)
library(topGO)
library(GO.db)
library(reutils)
library(GSEABase)
library(KEGGREST)
library(pathview)
library(ggplot2)
library(lsa)
library("RColorBrewer")
library("PoiClaClu")
library(rlang)
library(kohonen)
library(oposSOM)
library(mefa)
library(VennDiagram)
library(UpSetR)
library(OmicCircos)
library(stats)

options(scipen = 1000)

options(reutils.email = '271296251017a@gmail.com', reutils.show.headlines = 10, reutils.verbose.queries = T, reutils.test.remote = F)

# abundance_full <- read.table("~/yeast_transcriptomics/technical/rh_all_abundance.tsv", header=TRUE, row.names = 1)
# poisd <- PoissonDistance(t(abundance_full))
# samplePoisDistMatrix <- as.matrix(poisd$dd)
# rownames(samplePoisDistMatrix) <- colnames(abundance_full)
# colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(samplePoisDistMatrix,  clustering_distance_rows=poisd$dd,clustering_distance_cols=poisd$dd, col=colors)
write.table(annot_full, file.path('~', 'Pea_transcriptomics', 'TableS1.tsv'), col.names = T, row.names = F, sep = '\t')


sample_id <- dir(file.path('~', 'Pea_transcriptomics', '25052019_kallisto+trinotate', 'kallisto_dir'))
metadata <- read.table('/home/reverend_casy/Pea_transcriptomics/25052019_kallisto+trinotate/metadata.txt', header = T) %>%
  mutate(path = file.path('~', 'Pea_transcriptomics', '25052019_kallisto+trinotate', 'kallisto_dir', sample_id)) %>% mutate(condition = as.factor(condition))

annot_full %>% filter(!is.na(GOs)) %>% nrow()

so_pea <- sleuth_prep(metadata, extra_bootstrap_summary = T)
so_pea <- sleuth_fit(so_pea, ~1, 'reduced')
so_pea <- sleuth_fit(so_pea, ~condition, 'full')
so_pea <- sleuth_wt(so_pea, 'condition30_days')
so_pea <- sleuth_lrt(so_pea, 'reduced', 'full')

lrt <- sleuth_results(so_pea, 'reduced:full', 'lrt', show_all = F)
wt <- sleuth_results(so_pea, 'condition30_days', 'wt', show_all = F)

lrt_filtered <- lrt %>% filter(qval < 0.001)
wt_filtered <- wt %>% filter(qval < 0.001)

joint_results <- wt_filtered[wt_filtered$target_id %in% lrt_filtered$target_id,]

#### Functional annotation

BLAST_decoder <- function(data, column) {
  column <- deparse(substitute(column))
  return(sapply(data[[column]][which(!is.na(data[[column]]))], function(x) x %>% strsplit('\\^') %>% unlist() %>% extract2(1)))
}

annot_full <-  read.table('final_report_processed_3.tsv', header = TRUE, sep = '\t', fill = TRUE, quote = "") %>% mutate(GOs = as.character(GOs),
                                                                                                                           gene_ontology_pfam = as.character(gene_ontology_pfam),
                                                                                                                           plants_BLASTP = as.character(plants_BLASTP),
                                                                                                                           plants_BLASTX = as.character(plants_BLASTX))



annot_full[annot_full == ""] <- NA
annot_full %>% View()
annot_full$GOs <- ifelse(is.na(annot_full$GOs), annot_full$gene_ontology_pfam, annot_full$GOs)

full_blastp_queries <- sapply(annot_full$plants_BLASTP[which(!is.na(annot_full$plants_BLASTP))], function(x)
  x %>% strsplit('\\^') %>% unlist() %>% extract2(1))
full_blastx_queries <- sapply(annot_full$plants_BLASTX[which(!is.na(annot_full$plants_BLASTX))], function(x)
  x %>% strsplit('\\^') %>% unlist() %>% extract2(1))

full_blastp_targets <- c()
full_blastp_targets_1 <- efetch(full_blastp_queries[27001:27005], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
full_blastp_targets_1 <- full_blastp_targets_1[grep('DEFINITION', full_blastp_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
full_blastp_targets <- c(full_blastp_targets, full_blastp_targets_1)
annot_full %>% View()
names(full_blastp_targets) <- names(full_blastp_queries)
annot_full$BLASTP_names <- NA
annot_full$BLASTP_names[which(!is.na(annot_full$plants_BLASTP))] <- full_blastp_targets


full_blastx_targets <- c()
full_blastx_targets_1 <- efetch(full_blastx_queries[27001:27237], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
full_blastx_targets_1 <- full_blastx_targets_1[grep('DEFINITION', full_blastx_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
full_blastx_targets <- c(full_blastx_targets, full_blastx_targets_1)
names(full_blastx_targets) <- names(full_blastx_queries)
annot_full$BLASTX_names <- NA
annot_full$BLASTX_names[which(!is.na(annot_full$plants_BLASTX))] <- full_blastx_targets

#### MapMan
mercator_res %>% View()
mercator_res <- read.table(file.path('~', 'Pea_transcriptomics','mercator.results.GO.tsv'), header = T, sep = '\t', fill = T) %>% mutate_all(as.character())
mercator_res <- mercator_res[match(annot_full$transcript_id, mercator_res$Contig),]
annot_full <- cbind(annot_full, mercator_res %>% dplyr::select(GO_MapMan, MapMan_terms))
annot_full$GOs <- ifelse(is.na(annot_full$GOs), annot_full$GO_MapMan, annot_full$GOs)

annot_full$GO_MapMan %<>% as.character()

for (i in 1:nrow(annot_full)) {
  if (!is.na(annot_full$gene_ontology_pfam[i]) & !is.na(annot_full$GOs[i])) {
    if (length(grep(annot_full$gene_ontology_pfam[i], annot_full$GOs[i])) == 0) {
      # print(grep(annot_full$gene_ontology_pfam[i], annot_full$GOs[i]))
      annot_full$GOs[i] <- paste(annot_full$GOs[i], annot_full$gene_ontology_pfam[i], sep = ', ')
    }
  }
}
annot_full %>% View()
length(which(!is.na(annot_full$GOs)))
#### KEGG
annot_trans <-  read.table(file.path('~', 'Pea_transcriptomics', 'final_report_processed_3.tsv'), header = TRUE, sep = '\t', fill = TRUE, quote = "") %>% mutate_all(as.character())
annot_trans <- annot_trans[match(annot_full$prot_id, annot_trans$prot_id),]
kegg_mapping <- read.table(file.path('~', 'Pea_transcriptomics', 'BlastKOALA.merged.filtered.tsv'), sep = '\t', header = T) %>% mutate_all(as.character())
annot_full$KEGG_Pathway <- annot_trans$KEGG_Pathway
annot_full$KEGG_Pathway[annot_full$KEGG_Pathway == ''] <- NA
annot_full$KEGG_Pathway %<>% as.character()
kegg_mapping <- kegg_mapping[match(annot_full$prot_id, kegg_mapping$target_id),]
annot_full$KEGG_Pathway_Blast_KOALA <- kegg_mapping$KEGG_Pathway %>% as.character()
# annot_full$KEGG_Pathway <- ifelse(is.na(annot_full$KEGG_Pathway), annot_full$KEGG_Pathway_Blast_KOALA, annot_full$KEGG_Pathway)
annot_full$KEGG_Pathway[which(is.na(annot_full$KEGG_Pathway))] <- annot_full$KEGG_Pathway_Blast_KOALA[which(is.na(annot_full$KEGG_Pathway))]

View(annot_full$KEGG_Pathway)


#### KEGGREST
View(keggList('brite'))




annot_full %>% dplyr::select(KEGG_Pathway) %>% View()

annot_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, seed_eggNOG_ortholog, GOs,
                                               BLASTX_names, BLASTP_names, GO_MapMan, MapMan_terms) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>%
  filter(transcript_id %in% joint_results$target_id)
annot_to_merge <- annot_to_merge[match(joint_results$target_id, annot_to_merge$transcript_id),]


# mercator_mapping <- read.table('pce12231-sup-0002-si.csv', header = T, sep = '\t', fill = T) %>% mutate_all(as.character())



annotated_joint <- cbind(joint_results, annot_to_merge) %>% dplyr::select(-c(transcript_id)) %>% mutate_all(as.character())
annotated_joint %>% View()
lrt_unique <- lrt_filtered[!(lrt_filtered$tartarget_id %in% wt_filtered$target_id),]
wt_unique <- wt_filtered[!(wt_filtered$target_id %in% lrt_filtered$target_id),]
annotated_joint %>% nrow()
upregulated_30_days <- annotated_joint %>% filter(b > 0)
downregulated_30_days <- annotated_joint %>% filter(b < 0)

write.table(joint_results, file.path('~', 'Pea_transcriptomics', 'outputs', 'full_table.tsv'), sep = '\t', row.names = F, col.names = T)
write.table(downregulated_30_days, file.path('~', 'Pea_transcriptomics', 'outputs', 'downreguated_30_days.tsv'), sep = '\t', row.names = F, col.names = T)
write.table(upregulated_30_days, file.path('~', 'Pea_transcriptomics', 'outputs', 'upregulated_30_days.tsv'), sep = '\t', row.names = F, col.names = T)

### GO Annotation
full_mapping <- annot_full %>% dplyr::select(transcript_id, GOs) %>% mutate(transcript_id = as.character(transcript_id),
                                                                            GOs = as.character(GOs))
full_mapping <- full_mapping %>% filter(is.na(GOs) == F)
full_mapping$GOs <- sapply(full_mapping$GOs, function(x) strsplit(x, ', '))
names(full_mapping$GOs) <- full_mapping$transcript_id
full_mapping <- full_mapping$GOs
# full_mapping %>% filter(is.na(GOs) == F) %>% nrow()
# full_mapping %>% head(1)

up30_GO_dummy <- upregulated_30_days %>% dplyr::select(target_id) %>% mutate(target_id = as.character(target_id)) %>% unlist()
down30_GO_dummy <- downregulated_30_days %>% dplyr::select(target_id) %>% mutate(target_id = as.character(target_id)) %>% unlist()

newFuncadelicGo <- function(x, database, mode, intGenes, cutoff = 0.01) {
  myGOdata <- new("topGOdata", 
                  description = '', ontology = database,
                  allGenes=x,
                  geneSel=selector,
                  annot=annFUN.gene2GO, gene2GO=full_mapping, nodeSize=10)
  resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")

  allRes <- GenTable(myGOdata, weight01Fisher = resultFisher,
                     orderBy = "weight01Fisher", ranksOf = "weight01Fisher", topNodes = 35)
  allRes$weight01Fisher <- as.double(allRes$weight01Fisher)
  allRes[is.na(allRes)] <- 0
  if (mode=='table') {
    allRes$Percentage_of_Differentially_Expressed_Genes_from_number_of_DEGs <- ((allRes$Significant / nrow(intGenes))*100) %>% round(2)
    allRes$Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes <- ((allRes$Significant / allRes$Annotated)*100) %>% round(2)
    allRes <- allRes %>% filter(weight01Fisher < cutoff)
    return(allRes)
  } else if (mode=='pvals'){
    retlist <- (allRes %>% filter(weight01Fisher < cutoff))$weight01Fisher
    names(retlist) <- (allRes %>% filter(weight01Fisher < cutoff))$GO.ID
    return(retlist)
  }
  else if(mode=='graph') {
    GoGraph <- showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = allRes %>% filter(weight01Fisher < cutoff) %>% nrow(), useInfo = "all")$dag
    return(GoGraph)
  }
}
full_mapping %>% head()

full_mapping_2 <- annot_full %>% dplyr::select(transcript_id, GOs) %>% mutate(transcript_id = as.character(transcript_id),
                                                                            GOs = as.character(GOs))
full_mapping_2 <- full_mapping_2 %>% filter(is.na(GOs) == F)
full_mapping_2$GOs <- sapply(full_mapping_2$GOs, function(x) strsplit(x, ', '))
names(full_mapping_2$GOs) <- full_mapping_2$transcript_id
full_mapping_2 <- full_mapping_2$GOs
full_mapping_2 %>% length()

geneUniverse <- full_mapping %>% names()
geneUniverse_2 <- full_mapping %>% names()

geneList_up <- factor(as.integer(geneUniverse %in% up30_GO_dummy))
names(geneList_up) <- geneUniverse
geneList_down <- factor(as.integer(geneUniverse %in% down30_GO_dummy))
names(geneList_down) <- geneUniverse

iterating_function_GO <- function(subset, ont, cutoff, make_table = T){
  dummy <- subset %>% dplyr::select(target_id) %>% mutate(target_id = as.character(target_id)) %>% unlist()
  geneList <- factor(as.integer(geneUniverse %in% dummy))
  names(geneList) <- geneUniverse_2
  name_val <- paste(ont, substitute(subset) %>% as.character(), sep = '_')
  if (make_table) assign(name_val, newFuncadelicGo(geneList, ont, mode = 'table', intGenes=subset, cutoff = cutoff), envir = .GlobalEnv)
  }

subsets <- list(upregulated_30_days, downregulated_30_days)
ontologies <- c('BP', 'CC', 'MF')

sapply(ontologies, function(x) iterating_function_GO(upregulated_30_days, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(downregulated_30_days, x, 0.001))



sapply(ls()[grep('(BP|CC|MF)_(up)|(down)regulated', ls())], function(x) get(x) %>% write.table(., file = file.path('~', 'Pea_transcriptomics', 'comparisons',
                                                                              paste(x, 'tsv', sep = '.')),
                                                                        sep = '\t', col.names = T, row.names = F))

# gsub("DEFINITION  ", "", unlist(strsplit(content(efetch(unlist(strsplit(upregulated_30_days$plants_BLASTX, '\\^'))[1], db = 'protein', rettype = 'gp', retmode = 'text')), '\n'))[2],)
# unlist(strsplit(content(efetch(unlist(strsplit(upregulated_30_days$plants_BLASTX, '\\^'))[1], db = 'protein', rettype = 'gp', retmode = 'text')), '\n'))[2]


# strsplit(upregulated_30_days$plants_BLASTX[377], '\\^') %>% unlist() %>% extract2(1) %>% efetch(db = 'protein', rettype = 'gp', retmode = 'text') %>% 
#   content() %>% strsplit('\n') %>% unlist() %>% extract2(2) %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.)
# 
# # annotated_joint %>% View()
# # efetch(unlist(strsplit(upregulated_30_days$plants_BLASTX[1], '\\^'))[1])
# 
# blastx_queries <- sapply(annotated_joint$plants_BLASTX[which(!is.na(annotated_joint$plants_BLASTX))], function(x)
#   x %>% strsplit('\\^') %>% unlist() %>% extract2(1))
# 
# blastp_queries <- sapply(annotated_joint$plants_BLASTP[which(!is.na(annotated_joint$plants_BLASTP))], function(x)
#   x %>% strsplit('\\^') %>% unlist() %>% extract2(1))


### Some bruteforce magic
# blastp_targets <- с()
# blastp_targets_1 <- efetch(blastp_queries[10001:10370], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
# blastp_targets_1 <- blastp_targets_1[grep('DEFINITION', blastp_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
# blastp_targets <- c(blastp_targets, blastp_targets_1)
# length(blastp_targets_1)
# length(blastp_targets)
# 
# names(blastp_targets) <- names(blastp_queries)
# blastp_queries %>% length()
# annotated_joint %>% View()
# 
# # annotated_joint$plants_BLASTP %>% is.na() %>% which() %>% length()
# annotated_joint$BLASTP_names <- NA
# annotated_joint$BLASTP_names[which(!is.na(annotated_joint$plants_BLASTP))] <- blastp_targets
# # blastx_queries %>% length()
# # 
# # blastx_targets <- c()
# # 
# # blastx_targets_1 <- efetch(blastx_queries[10001:10400], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
# # blastx_targets_1 <- blastx_targets_1[grep('DEFINITION', blastx_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
# # blastx_targets <- c(blastx_targets, blastx_targets_1)
# # 
# # names(blastx_targets) <- names(blastx_queries)
# # 
# annotated_joint$BLASTX_names <- NA
# annotated_joint$BLASTX_names[which(!is.na(annotated_joint$plants_BLASTX))] <- blastx_targets
# 
# 
# upregulated_30_days <- annotated_joint %>% filter(b > 0)
# downregulated_30_days <- annotated_joint %>% filter(b < 0)

write.table(upregulated_30_days, file.path('~', 'Pea_transcriptomics', 'sleuth_outputs', 'upregulated_30_days.tsv'), sep = '\t', col.names = T,
            row.names = F)
write.table(downregulated_30_days, file.path('~', 'Pea_transcriptomics', 'sleuth_outputs', 'downregulated_30_days.tsv'), sep = '\t', col.names = T,
            row.names = F)

####### WELCOME TO THE RICE FIELDS
chsample_id <- dir(file.path('~', 'Pea_transcriptomics', 'kallisto_subsamples_chinese'))
chmetadata <- read.table(file.path('~', 'Pea_transcriptomics', 'chinese_metadata.txt'), sep = '\t', header = T) %>% 
  mutate(path = file.path('~', 'Pea_transcriptomics', 'kallisto_subsamples_chinese', chsample_id))
chmetadata_zhewan <- chmetadata %>% filter(line == 'Zhewan_1') %>% mutate(condition = as.factor(condition))
chmetadata_zhongwan <- chmetadata %>% filter(line == 'Zhongwan_6')%>% mutate(condition = as.factor(condition))

### ZHEWAN
zhe_de <- sleuth_prep(chmetadata_zhewan, extra_bootstrap_summary = T)
zhe_de <- sleuth_fit(zhe_de, ~1, 'reduced')
zhe_de <- sleuth_fit(zhe_de, ~condition, 'full')
zhe_de <- sleuth_wt(zhe_de, 'condition25')
zhe_de <- sleuth_lrt(zhe_de, 'reduced', 'full')
zhe_lrt <- sleuth_results(zhe_de, 'reduced:full', 'lrt', show_all = F) %>% filter(qval < 0.001)
zhe_wt <- sleuth_results(zhe_de, 'condition25', 'wt', show_all = F) %>% filter(qval < 0.001)
zhe_merged <- zhe_wt %>% filter(target_id %in% zhe_lrt$target_id)
annot_full %>% View()
zhe_merged %>% filter(target_id %in% annotated_joint$target_id) %>% nrow()
zhe_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, seed_eggNOG_ortholog, GOs, BLASTX_names,
                                             BLASTP_names, KEGG_Pathway, MapMan_terms) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% zhe_merged$target_id)
zhe_to_merge <- zhe_to_merge[match(zhe_merged$target_id, zhe_to_merge$transcript_id),]
zhe_merged <- cbind(zhe_merged, zhe_to_merge) %>% dplyr::select(-c(transcript_id))
zhe_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))
zhe_merged %>% View()
# zhe_merged$BLASTX_names <- NA
# zhe_merged$BLASTP_names <- NA
# 
# zhe_blastx_queries <- BLAST_decoder(zhe_merged, plants_BLASTX)
# zhe_blastx_queries %>% length()
# zhe_blastx_targets <- c()
# zhe_blastx_targets_1 <- efetch(zhe_blastx_queries[6501:6609], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
# zhe_blastx_targets_1 <- zhe_blastx_targets_1[grep('DEFINITION', zhe_blastx_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
# zhe_blastx_targets <- c(zhe_blastx_targets, zhe_blastx_targets_1)
# names(zhe_blastx_targets) <- names(zhe_blastx_queries)
# 
# zhe_merged$BLASTX_names <- ifelse(!is.na(zhe_merged$plants_BLASTX), zhe_blastx_targets, zhe_merged$BLASTX_names)
# 
# zhe_blastp_queries <- BLAST_decoder(zhe_merged, plants_BLASTP)
# zhe_blastp_queries %>% length()
# zhe_blastp_targets <- c()
# 
# zhe_blastp_targets_1 <- efetch(zhe_blastp_queries[6501:6616], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
# zhe_blastp_targets_1 <- zhe_blastp_targets_1[grep('DEFINITION', zhe_blastp_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
# zhe_blastp_targets <- c(zhe_blastp_targets, zhe_blastp_targets_1)
# zhe_blastp_targets_1 %>% length()
# zhe_blastp_targets %>% length()
# names(zhe_blastp_targets) <- names(zhe_blastp_queries)
# 
# zhe_merged$BLASTP_names[which(!is.na(zhe_merged$plants_BLASTP))] <- zhe_blastp_targets
# zhe_merged$BLASTX_names[which(!is.na(zhe_merged$plants_BLASTX))] <- zhe_blastx_targets

zhe_up <- zhe_merged %>% filter(b > 0)
zhe_down <- zhe_merged %>% filter(b < 0)

write.table(zhe_up, file.path('~', 'Pea_transcriptomics', 'chinese_sleuth', 'zhewan_up.tsv'), row.names = F, sep = '\t')
write.table(zhe_down, file.path('~', 'Pea_transcriptomics', 'chinese_sleuth', 'zhewan_down.tsv'), row.names = F, sep = '\t')

### ZHONGWAN
zhong_de <- sleuth_prep(chmetadata_zhongwan, extra_bootstrap_summary = T)
zhong_de <- sleuth_fit(zhong_de, ~1, 'reduced')
zhong_de <- sleuth_fit(zhong_de, ~condition, 'full')
zhong_de <- sleuth_wt(zhong_de, 'condition25')
zhong_de <- sleuth_lrt(zhong_de, 'reduced', 'full')
zhong_lrt <- sleuth_results(zhong_de, 'reduced:full', 'lrt', show_all = F) %>% filter(qval < 0.001)
zhong_wt <- sleuth_results(zhong_de, 'condition25', 'wt', show_all = F) %>% filter(qval < 0.001)
zhong_merged <- zhong_wt %>% filter(target_id %in% zhong_lrt$target_id)
sleuth_results(zhong_de, 'reduced:full', 'lrt', show_all = F) %>% filter(target_id == 'TRINITY_10_DN52145_c0_g1_i1')
zhong_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                               seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names, KEGG_Pathway, MapMan_terms) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% zhong_merged$target_id)
zhong_to_merge <- zhong_to_merge[match(zhong_merged$target_id, zhong_to_merge$transcript_id),]
zhong_merged <- cbind(zhong_merged, zhong_to_merge) %>% dplyr::select(-c(transcript_id))
zhong_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

# zhong_merged$BLASTX_names <- NA
# zhong_merged$BLASTP_names <- NA
# 
# zhong_blastx_queries <- BLAST_decoder(zhong_merged, plants_BLASTX)
# zhong_blastx_queries %>% length()
# zhong_blastx_targets <- c()
# zhong_blastx_targets_1 <- efetch(zhong_blastx_queries[3001:3416], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
# zhong_blastx_targets_1 <- zhong_blastx_targets_1[grep('DEFINITION', zhong_blastx_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
# zhong_blastx_targets <- c(zhong_blastx_targets, zhong_blastx_targets_1)
# zhong_blastx_targets_1 %>% length()
# zhong_blastx_targets %>% length()
# 
# names(zhong_blastx_targets) <- names(zhong_blastx_queries)
# 
# zhong_blastp_queries <- BLAST_decoder(zhong_merged, plants_BLASTP)
# zhong_blastp_queries %>% length()
# zhong_blastp_targets <- c()
# zhong_blastp_targets_1 <- efetch(zhong_blastp_queries[3001:3424], db = 'protein', rettype = 'gp', retmode = 'text') %>% content() %>% strsplit('\n') %>% unlist()
# zhong_blastp_targets_1 <- zhong_blastp_targets_1[grep('DEFINITION', zhong_blastp_targets_1)] %>% sapply(function(x) x %>%  gsub("DEFINITION  ", "", x=.) %>% gsub("\\.", "", x=.))
# zhong_blastp_targets <- c(zhong_blastp_targets, zhong_blastp_targets_1)
# zhong_blastp_targets_1 %>% length()
# zhong_blastp_targets %>% length()
# 
# names(zhong_blastp_targets) <- names(zhong_blastp_queries)
# 
# zhong_merged$BLASTP_names[which(!is.na(zhong_merged$plants_BLASTP))] <- zhong_blastp_targets
# zhong_merged$BLASTX_names[which(!is.na(zhong_merged$plants_BLASTX))] <- zhong_blastx_targets

zhong_up <- zhong_merged %>% filter(b > 0)
zhong_down <- zhong_merged %>% filter(b < 0)

write.table(zhong_up, file.path('~', 'Pea_transcriptomics', 'chinese_sleuth', 'zhongwan_up.tsv'), row.names = F, sep = '\t')
write.table(zhong_down, file.path('~', 'Pea_transcriptomics', 'chinese_sleuth', 'zhongwan_down.tsv'), row.names = F, sep = '\t')

### GO Annotation
sapply(ontologies, function(x) iterating_function_GO(zhe_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(zhe_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(zhong_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(zhong_down, x, 0.001))

sapply(ls()[grep('(BP)|(CC)|(MF)_(zhe|zhong)', ls())], function(x) get(x) %>% write.table(., file = file.path('~', 'Pea_transcriptomics', 'chinese_sleuth', 'GO_tables',
                                                                                            paste(x, 'tsv', sep = '.')),
                                                                        sep = '\t', col.names = T, row.names = F))


#### RICE FIELDS VS SPRING-2
chmetadata$condition %<>% as.factor()
metadata$sample %<>% as.factor()

######## 1 ################

#### Sprint2-10 vs Zhewan1-10
meta_1 <- rbind(chmetadata %>% filter(line == 'Zhewan_1', condition == '10') %>% dplyr::select(-c(line)), metadata %>% filter(condition == '10_days'))

sleuth_cross_1 <- sleuth_prep(meta_1, extra_bootrstrap_summary = T)
sleuth_cross_1 <- sleuth_fit(sleuth_cross_1, ~1, 'reduced')
sleuth_cross_1 <- sleuth_fit(sleuth_cross_1, ~condition, 'full')
models(sleuth_cross_1)
sleuth_cross_1 <- sleuth_lrt(sleuth_cross_1, 'reduced', 'full')
sleuth_cross_1 <- sleuth_wt(sleuth_cross_1, 'condition10_days')

cross_1_lrt <- sleuth_results(sleuth_cross_1, 'reduced:full', 'lrt', show_all = F)
cross_1_wt <- sleuth_results(sleuth_cross_1, 'condition10_days', 'wt', show_all = F)

cross_1_merged <- cross_1_wt %>% filter(target_id %in% cross_1_lrt$target_id) %>% filter(qval < 0.001)
cross_1_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                               seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names, KEGG_Pathway, MapMan_terms) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_1_merged$target_id)
cross_1_to_merge <- cross_1_to_merge[match(cross_1_merged$target_id, cross_1_to_merge$transcript_id),]
cross_1_merged <- cbind(cross_1_merged, cross_1_to_merge) %>% dplyr::select(-c(transcript_id))
cross_1_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_1_up <- cross_1_merged %>% filter(b > 0)
cross_1_down <- cross_1_merged %>% filter(b < 0)

cross_1_up %>% nrow()
#### Sprint2-10 vs Zhongwan6-10

meta_2 <- rbind(chmetadata %>% filter(line == 'Zhongwan_6', condition == '10') %>% dplyr::select(-c(line)), metadata %>% filter(condition == '10_days'))

sleuth_cross_2 <- sleuth_prep(meta_2, extra_bootrstrap_summary = T)
sleuth_cross_2 <- sleuth_fit(sleuth_cross_2, ~1, 'reduced')
sleuth_cross_2 <- sleuth_fit(sleuth_cross_2, ~condition, 'full')
models(sleuth_cross_2)
sleuth_cross_2 <- sleuth_lrt(sleuth_cross_2, 'reduced', 'full')
sleuth_cross_2 <- sleuth_wt(sleuth_cross_2, 'condition10_days')

cross_2_lrt <- sleuth_results(sleuth_cross_2, 'reduced:full', 'lrt', show_all = F)
cross_2_wt <- sleuth_results(sleuth_cross_2, 'condition10_days', 'wt', show_all = F)

cross_2_merged <- cross_2_wt %>% filter(target_id %in% cross_2_lrt$target_id) %>% filter(qval < 0.001)
cross_2_merged %>% nrow()

cross_2_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_2_merged$target_id)
cross_2_to_merge <- cross_2_to_merge[match(cross_2_merged$target_id, cross_2_to_merge$transcript_id),]
cross_2_merged <- cbind(cross_2_merged, cross_2_to_merge) %>% dplyr::select(-c(transcript_id))
cross_2_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_2_up <- cross_2_merged %>% filter(b > 0)
cross_2_down <- cross_2_merged %>% filter(b < 0)

cross_2_down %>% filter(target_id %in% cross_1_down$target_id) %>% nrow()
cross_2_up %>% filter(target_id %in% cross_1_up$target_id) %>% nrow()

cross_2_up %>% nrow()
#### Shared genes

shared_point1_down <- cross_2_down$target_id[cross_2_down$target_id %in% cross_1_down$target_id]
shared_point1_up <- cross_2_up$target_id[cross_2_up$target_id %in% cross_1_up$target_id]

length(shared_point1_down) / length(cross_2_down$target_id)
length(shared_point1_down) / length(cross_1_down$target_id)


######### 3 ###################
### Zhewan1-10 vs Zhongwan6-10
meta_3 <- chmetadata %>% filter(condition == '10')
sleuth_cross_3 <- sleuth_prep(meta_3, extra_bootrstrap_summary = T)
sleuth_cross_3 <- sleuth_fit(sleuth_cross_3, ~1, 'reduced')
sleuth_cross_3 <- sleuth_fit(sleuth_cross_3, ~line, 'full')
models(sleuth_cross_3)
sleuth_cross_3 <- sleuth_lrt(sleuth_cross_3, 'reduced', 'full')
sleuth_cross_3 <- sleuth_wt(sleuth_cross_3, 'lineZhongwan_6')

cross_3_lrt <- sleuth_results(sleuth_cross_3, 'reduced:full', 'lrt', show_all = F)
cross_3_wt <- sleuth_results(sleuth_cross_3, 'lineZhongwan_6', 'wt', show_all = F)

cross_3_merged <- cross_3_wt %>% filter(target_id %in% cross_3_lrt$target_id) %>% filter(qval < 0.001)

cross_3_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_3_merged$target_id)
cross_3_to_merge <- cross_3_to_merge[match(cross_3_merged$target_id, cross_3_to_merge$transcript_id),]
cross_3_merged <- cbind(cross_3_merged, cross_3_to_merge) %>% dplyr::select(-c(transcript_id))
cross_3_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_3_up <- cross_3_merged %>% filter(b > 0)
cross_3_down <- cross_3_merged %>% filter(b < 0)

cross_3_down %>% filter(target_id %in% cross_2_down$target_id) %>% nrow()
cross_3_up %>% filter(target_id %in% cross_2_up$target_id) %>% nrow()

cross_3_down %>% filter(target_id %in% cross_2_up$target_id) %>% nrow()
cross_3_up %>% filter(target_id %in% cross_2_down$target_id) %>% nrow()

cross_3_down %>% filter(target_id %in% cross_1_down$target_id) %>% nrow()
cross_3_up %>% filter(target_id %in% cross_1_up$target_id) %>% nrow()

cross_3_down %>% filter(target_id %in% cross_1_up$target_id) %>% nrow()
cross_3_up %>% filter(target_id %in% cross_1_down$target_id) %>% nrow()


########## 5 ###################
### Sprint2-30 vs Zhewan1-25

meta_4 <- rbind(chmetadata %>% filter(line == 'Zhewan_1', condition == '25') %>% dplyr::select(-c(line)), metadata %>% filter(condition == '30_days'))
sleuth_cross_4 <- sleuth_prep(meta_4, extra_bootrstrap_summary = T)
sleuth_cross_4 <- sleuth_fit(sleuth_cross_4, ~1, 'reduced')
sleuth_cross_4 <- sleuth_fit(sleuth_cross_4, ~condition, 'full')
models(sleuth_cross_4)
sleuth_cross_4 <- sleuth_lrt(sleuth_cross_4, 'reduced', 'full')
sleuth_cross_4 <- sleuth_wt(sleuth_cross_4, 'condition30_days')

cross_4_lrt <- sleuth_results(sleuth_cross_4, 'reduced:full', 'lrt', show_all = F)
cross_4_wt <- sleuth_results(sleuth_cross_4, 'condition30_days', 'wt', show_all = F)

cross_4_merged <- cross_4_wt %>% filter(target_id %in% cross_4_lrt$target_id) %>% filter(qval < 0.001)

cross_4_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_4_merged$target_id)
cross_4_to_merge <- cross_4_to_merge[match(cross_4_merged$target_id, cross_4_to_merge$transcript_id),]
cross_4_merged <- cbind(cross_4_merged, cross_4_to_merge) %>% dplyr::select(-c(transcript_id))
cross_4_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_4_up <- cross_4_merged %>% filter(b > 0)
cross_4_down <- cross_4_merged %>% filter(b < 0)
sessionInfo()

cross_4_up %>% nrow()
cross_4_down %>% nrow()

### Sprint2-30 vs Zhongwan6-25
meta_5 <- rbind(chmetadata %>% filter(line == 'Zhongwan_6', condition == '25') %>% dplyr::select(-c(line)), metadata %>% filter(condition == '30_days'))
sleuth_cross_5 <- sleuth_prep(meta_5, extra_bootrstrap_summary = T)
sleuth_cross_5 <- sleuth_fit(sleuth_cross_5, ~1, 'reduced')
sleuth_cross_5 <- sleuth_fit(sleuth_cross_5, ~condition, 'full')
models(sleuth_cross_5)
sleuth_cross_5 <- sleuth_lrt(sleuth_cross_5, 'reduced', 'full')
sleuth_cross_5 <- sleuth_wt(sleuth_cross_5, 'condition30_days')

cross_5_lrt <- sleuth_results(sleuth_cross_5, 'reduced:full', 'lrt', show_all = F)
cross_5_wt <- sleuth_results(sleuth_cross_5, 'condition30_days', 'wt', show_all = F)

cross_5_merged <- cross_5_wt %>% filter(target_id %in% cross_5_lrt$target_id) %>% filter(qval < 0.001)

cross_5_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_5_merged$target_id)
cross_5_to_merge <- cross_5_to_merge[match(cross_5_merged$target_id, cross_5_to_merge$transcript_id),]
cross_5_merged <- cbind(cross_5_merged, cross_5_to_merge) %>% dplyr::select(-c(transcript_id))
cross_5_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_5_up <- cross_5_merged %>% filter(b > 0)
cross_5_down <- cross_5_merged %>% filter(b < 0)

cross_5_up %>% nrow()
cross_5_down %>% nrow()


##### 7 ###################
### Sprint2-10 vs Zhewan1-25
meta_6 <- rbind(metadata %>% filter(condition == '10_days'), chmetadata %>% filter(line == 'Zhewan_1', condition == '25') %>% dplyr::select(-c(line)))
sleuth_cross_6 <- sleuth_prep(meta_6, extra_bootrstrap_summary = T)
sleuth_cross_6 <- sleuth_fit(sleuth_cross_6, ~1, 'reduced')
sleuth_cross_6 <- sleuth_fit(sleuth_cross_6, ~condition, 'full')
models(sleuth_cross_6)

sleuth_cross_6 <- sleuth_lrt(sleuth_cross_6, 'reduced', 'full')
sleuth_cross_6 <- sleuth_wt(sleuth_cross_6, 'condition25')

cross_6_lrt <- sleuth_results(sleuth_cross_6, 'reduced:full', 'lrt', show_all = F)
cross_6_wt <- sleuth_results(sleuth_cross_6, 'condition10_days', 'wt', show_all = F)

cross_6_merged <- cross_6_wt %>% filter(target_id %in% cross_6_lrt$target_id) %>% filter(qval < 0.001)
cross_6_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_6_merged$target_id)
cross_6_to_merge <- cross_6_to_merge[match(cross_6_merged$target_id, cross_6_to_merge$transcript_id),]
cross_6_merged <- cbind(cross_6_merged, cross_6_to_merge) %>% dplyr::select(-c(transcript_id))
cross_6_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_6_up <- cross_6_merged %>% filter(b > 0)
cross_6_down <- cross_6_merged %>% filter(b < 0)

# ann <- ann %>% mutate(ann[,1] = as.character(ann[,1])

### Sprint2-10 vs Zhongwan6-25
meta_7 <- rbind(metadata %>% filter(condition == '10_days'), chmetadata %>% filter(line == 'Zhongwan_6', condition == '25') %>% dplyr::select(-c(line)))
sleuth_cross_7 <- sleuth_prep(meta_7, extra_bootrstrap_summary = T)
sleuth_cross_7 <- sleuth_fit(sleuth_cross_7, ~1, 'reduced')
sleuth_cross_7 <- sleuth_fit(sleuth_cross_7, ~condition, 'full')
models(sleuth_cross_7)

sleuth_cross_7 <- sleuth_lrt(sleuth_cross_7, 'reduced', 'full')
sleuth_cross_7 <- sleuth_wt(sleuth_cross_7, 'condition25')
models(sleuth_cross_7)
cross_7_lrt <- sleuth_results(sleuth_cross_7, 'reduced:full', 'lrt', show_all = F)
cross_7_wt <- sleuth_results(sleuth_cross_7, 'condition25', 'wt', show_all = F)

cross_7_merged <- cross_7_wt %>% filter(target_id %in% cross_7_lrt$target_id) %>% filter(qval < 0.001)
cross_7_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_7_merged$target_id)
cross_7_to_merge <- cross_7_to_merge[match(cross_7_merged$target_id, cross_7_to_merge$transcript_id),]
cross_7_merged <- cbind(cross_7_merged, cross_7_to_merge) %>% dplyr::select(-c(transcript_id))
cross_7_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_7_up <- cross_7_merged %>% filter(b > 0)
cross_7_down <- cross_7_merged %>% filter(b < 0)

cross_7_up %>% nrow()
cross_7_down %>% nrow()
##################### Переделать GO

###### 8 #########
#### Zhewan1-25 vs Zhongwan6-25
meta_8 <- chmetadata %>% filter(condition == '25') %>% dplyr::select(-c(condition))
sleuth_cross_8 <- sleuth_prep(meta_8, extra_bootrstrap_summary = T)
sleuth_cross_8 <- sleuth_fit(sleuth_cross_8, ~1, 'reduced')
sleuth_cross_8 <- sleuth_fit(sleuth_cross_8, ~line, 'full')
models(sleuth_cross_8)

sleuth_cross_8 <- sleuth_lrt(sleuth_cross_8, 'reduced', 'full')
sleuth_cross_8 <- sleuth_wt(sleuth_cross_8, 'lineZhongwan_6')

cross_8_lrt <- sleuth_results(sleuth_cross_8, 'reduced:full', 'lrt', show_all = F)
cross_8_wt <- sleuth_results(sleuth_cross_8, 'lineZhongwan_6', 'wt', show_all = F)

cross_8_merged <- cross_8_wt %>% filter(target_id %in% cross_8_lrt$target_id) %>% filter(qval < 0.001)
cross_8_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_8_merged$target_id)
cross_8_to_merge <- cross_8_to_merge[match(cross_8_merged$target_id, cross_8_to_merge$transcript_id),]
cross_8_merged <- cbind(cross_8_merged, cross_8_to_merge) %>% dplyr::select(-c(transcript_id))
cross_8_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_8_up <- cross_8_merged %>% filter(b > 0)
cross_8_down <- cross_8_merged %>% filter(b < 0)

cross_8_up %>% nrow()
cross_8_down %>% nrow()

#### 9 ######
### Zhewan1-10 vs Zhongwan6-25
meta_9 <- rbind(chmetadata %>% filter(condition == '10', line == 'Zhewan_1') %>% dplyr::select(-c(line)),
                chmetadata %>% filter(condition == '25', line == 'Zhongwan_6') %>% dplyr::select(-c(line)))
sleuth_cross_9 <- sleuth_prep(meta_9, extra_bootrstrap_summary = T)
sleuth_cross_9 <- sleuth_fit(sleuth_cross_9, ~1, 'reduced')
sleuth_cross_9 <- sleuth_fit(sleuth_cross_9, ~condition, 'full')

sleuth_cross_9 <- sleuth_lrt(sleuth_cross_9, 'reduced', 'full')
sleuth_cross_9 <- sleuth_wt(sleuth_cross_9, 'condition25')

cross_9_lrt <- sleuth_results(sleuth_cross_9, 'reduced:full', 'lrt', show_all = F)
cross_9_wt <- sleuth_results(sleuth_cross_9, 'condition25', 'wt', show_all = F)

cross_9_merged <- cross_9_wt %>% filter(target_id %in% cross_9_lrt$target_id) %>% filter(qval < 0.001)
cross_9_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_9_merged$target_id)
cross_9_to_merge <- cross_9_to_merge[match(cross_9_merged$target_id, cross_9_to_merge$transcript_id),]
cross_9_merged <- cbind(cross_9_merged, cross_9_to_merge) %>% dplyr::select(-c(transcript_id))
cross_9_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_9_up <- cross_9_merged %>% filter(b > 0)
cross_9_down <- cross_9_merged %>% filter(b < 0)

##### 10
### Zhongwan6-10 vs Zhewan1-25
meta_10 <- rbind(chmetadata %>% filter(condition == '10', line == 'Zhongwan_6') %>% dplyr::select(-c(line)),
                 chmetadata %>% filter(condition == '25', line == 'Zhewan_1') %>% dplyr::select(-c(line)))

sleuth_cross_10 <- sleuth_prep(meta_10, extra_bootrstrap_summary = T)
sleuth_cross_10 <- sleuth_fit(sleuth_cross_10, ~1, 'reduced')
sleuth_cross_10 <- sleuth_fit(sleuth_cross_10, ~condition, 'full')

sleuth_cross_10 <- sleuth_lrt(sleuth_cross_10, 'reduced', 'full')
sleuth_cross_10 <- sleuth_wt(sleuth_cross_10, 'condition25')

cross_10_lrt <- sleuth_results(sleuth_cross_10, 'reduced:full', 'lrt', show_all = F)
cross_10_wt <- sleuth_results(sleuth_cross_10, 'condition25', 'wt', show_all = F)

cross_10_merged <- cross_10_wt %>% filter(target_id %in% cross_10_lrt$target_id) %>% filter(qval < 0.001)
cross_10_to_merge <- annot_full %>% dplyr::select(transcript_id, plants_BLASTX, plants_BLASTP, Pfam, gene_ontology_pfam, 
                                                 seed_eggNOG_ortholog, GOs, BLASTX_names, BLASTP_names) %>% 
  as.data.frame() %>% mutate_all(as.character()) %>% filter(transcript_id %in% cross_10_merged$target_id)
cross_10_to_merge <- cross_10_to_merge[match(cross_10_merged$target_id, cross_10_to_merge$transcript_id),]
cross_10_merged <- cbind(cross_10_merged, cross_10_to_merge) %>% dplyr::select(-c(transcript_id))
cross_10_merged %<>% mutate(plants_BLASTX = as.character(plants_BLASTX), plants_BLASTP = as.character(plants_BLASTP))

cross_10_up <- cross_10_merged %>% filter(b > 0)
cross_10_down <- cross_10_merged %>% filter(b < 0)



#### GO annotation
ls()[grep('cross_[0-9]_(up|down)', ls())]
# sapply(ontologies, function(x) sapply(ls()[grep('cross_[0-9]_(up|down)', ls())], function(y) iterating_function_GO(y, x, 0.001)))
sapply(ontologies, function(x) iterating_function_GO(cross_1_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_1_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_2_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_2_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_3_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_3_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_4_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_4_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_5_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_5_down, x, 0.001))

sapply(ontologies, function(x) iterating_function_GO(cross_6_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_6_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_7_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_7_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_8_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_8_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_9_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_9_down, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_10_up, x, 0.001))
sapply(ontologies, function(x) iterating_function_GO(cross_10_down, x, 0.001))


# sapply(ls()[grep('(BP|CC|MF)_cross_([0-9])', ls())], function(x) get(x) %>% write.table(., file = file.path('~', 'Pea_transcriptomics', 'comparisons',
#                                                                                             paste(x, 'tsv', sep = '.')), sep = '\t'))
# sapply(ls()[grep('^cross_([0-9]|10)_(up|down)', ls())], function(x) get(x) %>% 
#   write.table(file.path('~', 'Pea_transcriptomics', 'comparisons', paste(x, 'tsv', sep = '.')), sep = '\t', row.names = F, col.names = T))

######PCA
### Plot PCA by tpm metrics
alab_sprint <- read.table(file.path('~', 'Pea_transcriptomics', '25052019_kallisto+trinotate', 'all_trinity_abundance_sprint.tsv'),
                          sep = '\t', header = T)
alab_chinese <- read.table(file.path('~', 'Pea_transcriptomics', 'kallisto_subsamples_chinese', 'all_trinity_abundance_chinese.tsv'),
                           sep = '\t', header = T)
alab_chinese <- alab_chinese[match(alab_sprint$target_id, alab_chinese$target_id),]
alab_full <- cbind(alab_sprint, alab_chinese %>% dplyr::select(-c(target_id)))


library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")

colnames(alab_full)[5] <- 'Sprint-2,10 (1)'
poisd <- PoissonDistance(t(alab_full[,2:25]))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- colnames(alab_full[2:25])
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(samplePoisDistMatrix,  clustering_distance_rows=poisd$dd,clustering_distance_cols=poisd$dd, col=colors,
         labels_col = rownames(samplePoisDistMatrix))

alab_for_pca = t(alab_full[,2:25])
rownames(alab_for_pca) <- colnames(alab_full)[2:25]
colnames(alab_for_pca) <- alab_full$target_id
alab_pca <- prcomp(alab_for_pca[,-c(which(apply(alab_for_pca, 2, var) == 0))], center = T, scale. = T)
summary(alab_pca)
str(alab_pca)
abs(alab_pca$rotation[,c(1:5)]) %>% head()
alab_pca
annot_full %>% filter(transcript_id == 'TRINITY_10_DN47076_c0_g1_i1') %>% dplyr::select(BLASTP_names, BLASTX_names)
pca_abs_top5 <- data.frame('PC1' = names(abs(alab_pca$rotation[,1]) %>% sort(decreasing = T)) %>% as.character(),
                           'PC1_values' = abs(alab_pca$rotation[,1]) %>% sort(decreasing = T),
                            'PC2' = names(abs(alab_pca$rotation[,2]) %>% sort(decreasing = T)) %>% as.character(),
                           'PC2_values' = abs(alab_pca$rotation[,2]) %>% sort(decreasing = T),
                           'PC3' = names(abs(alab_pca$rotation[,3]) %>% sort(decreasing = T)) %>% as.character(),
                           'PC3_values' = abs(alab_pca$rotation[,3]) %>% sort(decreasing = T),
                           'PC4' = names(abs(alab_pca$rotation[,4]) %>% sort(decreasing = T)) %>% as.character(),
                           'PC4_values' = abs(alab_pca$rotation[,4]) %>% sort(decreasing = T),
                           'PC5' = names(abs(alab_pca$rotation[,5]) %>% sort(decreasing = T)) %>% as.character(),
                           'PC5_values' = abs(alab_pca$rotation[,5]) %>% sort(decreasing = T))
annot_full %>% filter(transcript_id == 'TRINITY_10_DN1570_c0_g1_i1') %>% dplyr::select(BLASTP_names)
annot_full[annot_full$transcript_id %in% pca_abs_top5$PC1,][match(pca_abs_top5$PC1, annot_full[annot_full$transcript_id %in% pca_abs_top5$PC1,]$transcript_id),] %>% dplyr::select(transcript_id, BLASTP_names) %>% head()
pca_abs_top5 %<>% mutate('PC1_names' = 
                           annot_full[annot_full$transcript_id %in% pca_abs_top5$PC1,][match(pca_abs_top5$PC1, annot_full[annot_full$transcript_id %in% pca_abs_top5$PC1,]$transcript_id),] %>% dplyr::select(BLASTP_names) %>% unlist(),
                         'PC2_names' = 
                           annot_full[annot_full$transcript_id %in% pca_abs_top5$PC2,][match(pca_abs_top5$PC2, annot_full[annot_full$transcript_id %in% pca_abs_top5$PC2,]$transcript_id),] %>% dplyr::select(BLASTP_names) %>% unlist(),
                         'PC3_names' = 
                           annot_full[annot_full$transcript_id %in% pca_abs_top5$PC3,][match(pca_abs_top5$PC3, annot_full[annot_full$transcript_id %in% pca_abs_top5$PC3,]$transcript_id),] %>% dplyr::select(BLASTP_names) %>% unlist(),
                         'PC4_names' = 
                           annot_full[annot_full$transcript_id %in% pca_abs_top5$PC4,][match(pca_abs_top5$PC4, annot_full[annot_full$transcript_id %in% pca_abs_top5$PC4,]$transcript_id),] %>% dplyr::select(BLASTP_names) %>% unlist(),
                         'PC5_names'= 
                           annot_full[annot_full$transcript_id %in% pca_abs_top5$PC5,][match(pca_abs_top5$PC5, annot_full[annot_full$transcript_id %in% pca_abs_top5$PC5,]$transcript_id),] %>% dplyr::select(BLASTP_names) %>% unlist())
pca_abs_top5 <- pca_abs_top5[,c(1,11,2,3,12,4,5,13,6,7,14,8,9,15,10)]
write.table(pca_abs_top5, file.path('~', 'Pea_transcriptomics', 'PC_prominent_values_cont_names.tsv'), col.names = T, row.names = F, sep = '\t')

pca_abs_top5_alt <- data.frame('PC1_top' = names(alab_pca$rotation[,1] %>% sort(decreasing = T) %>% head(20)),
                           'PC1_bottom' = names(alab_pca$rotation[,1] %>% sort(decreasing = T) %>% tail(20)),
                           'PC2_top' = names(alab_pca$rotation[,2] %>% sort(decreasing = T) %>% head(20)),
                           'PC2_bottom' = names(alab_pca$rotation[,2] %>% sort(decreasing = T) %>% tail(20)),
                           'PC3_top' = names(alab_pca$rotation[,3] %>% sort(decreasing = T) %>% head(20)),
                           'PC3_bottom' = names(alab_pca$rotation[,3] %>% sort(decreasing = T) %>% tail(20)),
                           'PC4_top' = names(alab_pca$rotation[,4] %>% sort(decreasing = T) %>% head(20)),
                           'PC4_bottom' = names(alab_pca$rotation[,4] %>% sort(decreasing = T) %>% tail(20)),
                           'PC5_top' = names(alab_pca$rotation[,5] %>% sort(decreasing = T) %>% head(20)),
                           'PC5_bottom' = names(alab_pca$rotation[,5] %>% sort(decreasing = T) %>% tail(20)))
pca_abs_top5_alt %<>% transmute('PC1_top' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC1_top) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC1_top, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                              'PC1_bottom' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC1_bottom) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC1_bottom, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC2_top' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC2_top) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC2_top, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC2_bottom' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC2_bottom) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC2_bottom, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC3_top' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC3_top) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC3_top, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC3_bottom' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC3_bottom) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC3_bottom, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC4_top' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC4_top) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC4_top, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC4_bottom' = 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC4_bottom) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC4_bottom, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC5_top'= 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC5_top) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC5_top, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist(),
                            'PC5_bottom'= 
                              annot_full %>% filter(transcript_id %in% pca_abs_top5_alt$PC5_bottom) %>% dplyr::select(transcript_id, BLASTP_names) %>% arrange(match(pca_abs_top5_alt$PC5_bottom, transcript_id)) %>% dplyr::select(BLASTP_names) %>% unlist())

write.table(pca_abs_top5_alt, file.path('~', 'Pea_transcriptomics', 'PC_prominent_values_original.tsv'), col.names = T, row.names = F, sep = '\t')

library("ggbiplot")

pca_data <- data.frame('Sample' = c('Sprint-2,10_1', 'Sprint-2,10_2', 'Sprint-2,10_3', 'Sprint-2,10_4',
                                    'Sprint-2,30_5', 'Sprint-2,30_6', 'Sprint-2,30_7', 'Sprint-2,30_8',
                                    as.character(chmetadata$sample)),
                       'Line' = c(rep('Sprint-2', 8), rep('Zhewan-1', 8), rep('Zhongwan-6', 8)),
                       'Condition' = c(rep(10, 4), rep(20, 4), rep(10, 4), rep(25, 4), rep(10, 4), rep(25, 4)))
pca_data$Color <- paste(pca_data$Line, pca_data$Condition, sep = ', ')



for (i in c(1:10)) {

g <- ggbiplot(alab_pca, choices = gtools::combinations(5, 2, c(1,2,3,4,5))[i,], varname.size = 0, labels.size = 0, alpha = 0, var.axes = F, groups = pca_data$Color) +
  expand_limits(x=c(-3,1)) + theme_bw() +
  theme(text = element_text(size=13, family = 'arial', color='black'), axis.text.x = element_text(color = 'black', family = 'arial'),
        axis.text.y = element_text(color = 'black', family = 'arial',),
        axis.line.x = element_line(color = 'black'),
        legend.text = element_text(size = 13),
        legend.title =
          element_text( size = 15)) + geom_point(aes(colour = pca_data$Color), size = 4, alpha = 0.3,
          color = c('#068013',  '#f1cb20', '')) + labs(colour = 'Cultivar and seed age')

g
ggsave(filename = file.path('~/Pea_transcriptomics/publication_pics/intext_figures/Figure6-7', paste(i, 'png', sep = '.')),
      plot =  last_plot(), device = png(), width = 20, height = 15, unit = 'cm')
}
g

gtools::combinations(5, 2, c(1,2,3,4,5))[1,] %>% as.vector()
#### Hacking plot_pca

library(tcR)
cosine.similarity(V1, V2V4, .do.norm = T)
cosine.similarity(V1, V3V5, .do.norm = T)



cosine(V1, V2V4)
cosine(V1, V3V5)
sleuth_results(sleuth_cross_4, 'reduced:full', 'lrt', show_all = F) %>% nrow()
cross_4_up %>% nrow()

######## Graph summary
alab_chinese %>% head()
alab_chinese$Zhewan1_10_mean <- apply(alab_chinese[,c(2:5)], 1, mean)
alab_chinese$Zhongwan6_10_mean <- apply(alab_chinese[,c(10:13)], 1, mean)
alab_chinese$Zhewan1_25_mean <- apply(alab_chinese[,c(6:9)], 1, mean)
alab_chinese$Zhongwan1_25_mean <- apply(alab_chinese[,c(14:17)], 1, mean)
alab_chinese$BLASTP_names <- annot_full[match(alab_chinese$target_id, annot_full$transcript_id),]$BLASTP_names
chinese_splot <- alab_chinese %>% select(target_id, Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean, Zhongwan1_25_mean)

alab_sprint$BLASTP_names <- annot_full[match(alab_sprint$target_id, annot_full$transcript_id),]$BLASTP_names
alab_sprint$Sprint2_10_mean <- apply(alab_sprint[,c(2:5)], 1, mean)
alab_sprint$Sprint2_30_mean <- apply(alab_sprint[,c(6:9)], 1, mean)
chinese_splot %<>% cbind(alab_sprint[,c(10,11,12)])
chinese_splot %>% head()

alab_chinese[match(alab_sprint$target_id, alab_chinese$target_id),]

##### getting rid of pesky legalbumins
chinese_splot %<>% filter(target_id != 'TRINITY_10_DN3139_c0_g1_i3')
chinese_splot %>% head()
chinese_splot_5k <- chinese_splot %>% filter_at(vars(-c(target_id, BLASTP_names, is_storage, MapMan_Category)), all_vars(.<=5000))
chinese_splot %>% head()
ggplot(chinese_splot_5k, aes(x = value, y = value, alpha  = 0.1, size = I(1.5))) + 
  geom_text(aes(label = ifelse(Zhewan1_10_mean >= 75000 & Zhewan1_25_mean >= 75000, as.character(BLASTP_names), '' )))
  geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan1_25_mean, col = 'Zhongwan6', alpha = 0.05)) +
  geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, col = 'Sprint2', alpha = 0.05)) +
  geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean, col = 'Zhewan1', alpha = 0.05)) +
  xlab(label = 'Earlier condition') + ylab('Latter condition') +
  guides(alpha = 'none') +
  scale_y_continuous(labels = scales::comma)

  
storage_mapman <- c('TRINITY_10_DN12311_c0_g2_i1', 'TRINITY_10_DN1695_c0_g1_i1', 'TRINITY_10_DN20343_c0_g1_i1', 'TRINITY_10_DN1695_c0_g1_i2',
                      'TRINITY_10_DN22874_c0_g1_i1', 'TRINITY_10_DN22921_c0_g1_i1', 'TRINITY_10_DN26967_c0_g1_i1', 'TRINITY_10_DN3482_c0_g1_i3',
                      'TRINITY_10_DN35816_c0_g1_i1', 'TRINITY_10_DN3926_c0_g1_i1', 'TRINITY_10_DN3926_c0_g1_i3', 'TRINITY_10_DN3926_c1_g1_i1',
                      'TRINITY_10_DN3985_c0_g1_i1', 'TRINITY_10_DN4215_c0_g1_i1', 'TRINITY_10_DN4215_c0_g1_i2', 'TRINITY_10_DN50438_c0_g1_i1',
                      'TRINITY_10_DN563_c0_g2_i2', 'TRINITY_10_DN563_c0_g2_i4 ', 'TRINITY_10_DN6344_c0_g1_i1', 'TRINITY_10_DN822_c1_g1_i5',
                      'TRINITY_10_DN907_c0_g1_i2', 'TRINITY_30_DN1328_c0_g2_i1', 'TRINITY_30_DN1478_c0_g1_i11', 'TRINITY_30_DN1478_c0_g1_i6',
                      'TRINITY_30_DN1776_c0_g1_i2', 'TRINITY_30_DN244_c0_g1_i4', 'TRINITY_30_DN244_c0_g1_i7', 'TRINITY_30_DN552_c0_g1_i18',
                      'TRINITY_30_DN552_c0_g1_i3', 'TRINITY_30_DN664_c0_g1_i13', 'TRINITY_30_DN8039_c0_g1_i1')
  
storage_mapman <- storage_mapman[!(storage_mapman %in% storage_proteins$target_id)]  

storage_proteins <- c(annot_full$BLASTP_names[grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTP_names)],
                      annot_full$BLASTX_names[grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTX_names)[
                        !(grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTX_names) %in% grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTP_names))
                        ]]
                      )
storage_proteins <- c(storage_proteins,
                      (chinese_splot %>% filter(Sprint2_10_mean < 1000 & Sprint2_30_mean > 10000) %>% select(BLASTP_names) %>% unlist() %>% unname())[
                        !(chinese_splot %>% filter(Sprint2_10_mean < 1000 & Sprint2_30_mean > 10000) %>% select(BLASTP_names) %>% unlist() %>% unname() %in% storage_proteins)
                      ],
                      annot_full$BLASTP_names[annot_full$transcript_id %in% storage_mapman])

names(storage_proteins) <- c(annot_full$transcript_id[grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTP_names)] %>% as.character(),
                             annot_full$transcript_id[grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTX_names)[
                               !(grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTX_names) %in% grep('([Vv]icilin)|([Ss]torage)|([Ll]eg[(umin)A])', annot_full$BLASTP_names))
                               ]] %>% as.character(),   
  chinese_splot %>% filter(Sprint2_10_mean < 1000 & Sprint2_30_mean > 10000) %>% filter(!(BLASTP_names %in% storage_proteins[c(1:31)])) %>% 
    select(target_id) %>% unlist() %>% unname() %>% as.character(),
  annot_full$transcript_id[annot_full$transcript_id %in% storage_mapman] %>% as.character())



storage_proteins <- data.frame('target_id' = names(storage_proteins), 'BLASTP_name' = storage_proteins)

storage_proteins %<>% arrange(match(target_id, chinese_splot$target_id))
chinese_splot$is_storage <- ifelse(chinese_splot$target_id %in% storage_proteins$target_id, 'Yes', 'No') %>% as.factor()

mapman_for_alab <- annot_full %>% select(transcript_id) %>% mutate(MapMan_Category = 
                                                                     sapply(annot_full$MapMan_terms %>% as.character, function(x) strsplit(x, '[.]') %>% unlist() %>%
                                                                              unname() %>% first() %>% strsplit('[,]') %>% unlist() %>% unname() %>% first())) %>%
  mutate(MapMan_Category_2 = sapply(annot_full$MapMan_terms %>% as.character, function(x) strsplit(x, '[.]') %>% unlist() %>%
                                      unname() %>% nth(2) %>% strsplit('[,]') %>% unlist() %>% unname() %>% first())) %>%
  mutate(MapMan_Category_3 = sapply(annot_full$MapMan_terms %>% as.character, function(x) strsplit(x, '[.]') %>% unlist() %>%
                                      unname() %>% nth(3) %>% strsplit('[,]') %>% unlist() %>% unname() %>% first())) %>%
  arrange(match(chinese_splot$target_id, mapman_for_alab$transcript_id))



chinese_splot %<>% mutate(MapMan_Category = mapman_for_alab$MapMan_Category, MapMan_Category2 = mapman_for_alab$MapMan_Category_2,
                          MapMan_Category3 = mapman_for_alab$MapMan_Category_3, relation_to_folding = sapply(annot_full$MapMan_terms, function(x) grepl('folding', x)),
                          HSF = sapply(annot_full$MapMan_terms, function(x) grepl('Heat-shock', x)))
chinese_splot$Heat_shock[chinese_splot$relation_to_folding == T] <- 'folding-related'
chinese_splot$Heat_shock[chinese_splot$HSF == T] <- 'Heat-shock_related_TF'
chinese_splot$Heat_shock %<>% as.factor()

chinese_splot$MapMan_Category[chinese_splot$is_storage == 'Yes'] <- 'storage protein'
chinese_splot$MapMan_Category[chinese_splot$MapMan_Category == 'ell wall'] <- 'cell wall'
chinese_splot$MapMan_Category[chinese_splot$MapMan_Category == 'mino acid metabolism'] <- 'amino acid metabolism'
chinese_splot$MapMan_Category[chinese_splot$MapMan_Category == 'major CHO metabolism'] <- 'CHO metabolism'
chinese_splot$MapMan_Category[chinese_splot$MapMan_Category == 'minor CHO metabolism'] <- 'CHO metabolism'
chinese_splot$MapMan_Category[chinese_splot$MapMan_Category == 'cell'] <- 'misc'
chinese_splot$MapMan_Category %<>% as.factor()
chinese_splot$MapMan_Category2 %<>% as.factor()
chinese_splot$MapMan_Category3 %<>% as.factor()

length(chinese_splot$MapMan_Category)

chinese_splot %>% filter(MapMan_Category == 'protein') %>% filter(Zhewan1_10_mean > 2000)

ggplot(chinese_splot %>% filter(MapMan_Category == 'protein'), aes(x = value, y = value, alpha  = 0.1, size = I(3), color = MapMan_Category2)) + 
  # geom_text(aes(label = ifelse(Zhewan1_10_mean >= 75000 & Zhewan1_25_mean >= 75000, as.character(BLASTP_names), '' )))
  # geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan1_25_mean, shape = 'Zhongwan6', alpha = 0.05)) +
  geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, shape = 'Sprint2'), alpha = 0.4) +
  # geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean, shape = 'Zhewan1'), alpha = 0.4) +
  xlab(label = 'Earlier condition') + ylab('Latter condition') +
  guides(alpha = 'none') +
  scale_y_continuous(labels = scales::comma)

####### Gotta plot some Δtpm
delta_splot <- chinese_splot %>% select(target_id, Sprint2_10_mean, Sprint2_30_mean, Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean, Zhongwan1_25_mean,
                                        BLASTP_names, MapMan_Category, MapMan_Category2, MapMan_Category3) %>% mutate(deltaSZhe_early = Zhewan1_10_mean - Sprint2_10_mean,
                                        deltaSZhe_late = Zhewan1_25_mean - Sprint2_30_mean, deltaSZho_early = Zhongwan6_10_mean - Sprint2_10_mean,
                                        deltaSZho_late = Zhongwan1_25_mean - Sprint2_30_mean)

delta_splot %<>% mutate(mdeltaSZhe_early = deltaSZhe_early / (Zhewan1_10_mean + Sprint2_10_mean),
                        mdeltaSZhe_late = deltaSZhe_late / (Zhewan1_25_mean + Sprint2_30_mean),
                        mdeltaSZho_early =  deltaSZho_early / (Zhongwan6_10_mean + Sprint2_10_mean),
                        mdeltaSZho_late = deltaSZho_late / (Zhongwan1_25_mean + Sprint2_30_mean)
                        )
delta_splot %>% head()

ggplot(delta_splot, aes(x = deltaSZhe_early, y = deltaSZhe_late, col = MapMan_Category)) + 
  # geom_text(aes(label = ifelse(Zhewan1_10_mean >= 75000 & Zhewan1_25_mean >= 75000, as.character(BLASTP_names), '' )))
  geom_point(alpha = 0.6, size = I(2.5)) +
  labs(title = 'Distribution of Δtpm regarding the seed maturation stage', x = 'Zhongwan6_10_mean - Sprint2_10_mean / Zhongwan6_10_mean + Sprint2_10_mean',
       y ='Zhongwan6_25_mean - Sprint2_30_mean / Zhongwan6_25_mean + Sprint2_10_mean') +
  guides(alpha = 'none') +
  scale_y_continuous(labels = scales::comma)
storage_proteins$BLASTX_names <- annot_full %>% filter(transcript_id %in% storage_proteins$target_id) %>% arrange(match(transcript_id, storage_proteins$target_id)) %>% select(BLASTX_names)
write.table(storage_proteins, file = file.path('~', 'Pea_transcriptomics', '130719_pics', 'storage_proteins.tsv'), col.names = T, row.names = F, sep = '\t')
storage_proteins
# / Zhewan1_10_mean + Sprint2_10_mean
# / Zhewan1_25_mean + Sprint2_30_mean
delta_splot %>% filter(is.nan(mdeltaSZhe_late)) %>% head()

PoisHeatMap <- function(df, pal){
  poisd <- PoissonDistance(t(df[,2:ncol(df)]))
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  rownames(samplePoisDistMatrix) <- colnames(df[2:ncol(df)])
  colors <- colorRampPalette(rev(brewer.pal(9, pal)) )(255)
  pheatmap(samplePoisDistMatrix,  clustering_distance_rows=poisd$dd,clustering_distance_cols=poisd$dd, col=colors,
           labels_col = rownames(samplePoisDistMatrix))
}
alab_no_storage %>% nrow()
alab_full %>% colnames()
colnames(alab_full)[2:25] <- c("Sprint-2,10 (1)", "Sprint-2,10 (2)", "Sprint-2,10 (3)", "Sprint-2, 10(1)",
                               "Sprint-2,20 (1)", "Sprint-2,20 (2)", "Sprint-2,20 (3)", "Sprint-2,20 (4)",
                               "Zhewan-1,10 (1)", "Zhewan-1,10 (2)", "Zhewan-1,10 (3)", "Zhewan-1,10 (4)",
                               "Zhewan-1,25 (1)", "Zhewan-1,25 (2)", "Zhewan-1,25 (3)", "Zhewan-1,25 (4)",
                               "Zhongwan-6,10 (1)", "Zhongwan-6,10 (2)", "Zhongwan-6,10 (3)", "Zhongwan-6,10 (4)",
                               "Zhongwan-6,25 (1)", "Zhongwan-6,25 (2)", "Zhongwan-6,25 (3)", "Zhongwan-6,25 (4)")

PoisHeatMap(alab_full[1:25] %>% filter(!(target_id %in% annot_full[grepl('stress', annot_full$MapMan_terms),]$transcript_id)), 'Greens')
PoisHeatMap(alab_full[1:25], 'Greens')

alab_no_storage <- alab_full %>% filter(!(target_id %in% storage_proteins$target_id)) %>%  filter(target_id != 'TRINITY_10_DN3139_c0_g1_i3')
alab_full$MapMan_Category <- chinese_splot$MapMan_Category
alab_full %>% select(target_id, MapMan_Category) %>% tail(10)
chinese_splot %>% select(target_id, MapMan_Category) %>% tail(10)

alab_full %>% filter(MapMan_Category == 'storage protein') %>% nrow()
alab_full %>% filter(MapMan_Category == 'protein') %>% nrow()
PoisHeatMap(alab_no_storage, "Greens")
PoisHeatMap(alab_full %>% filter(MapMan_Category == 'storage protein' | MapMan_Category == 'protein') %>% select(-MapMan_Category), "Greens")
alab_full %>% nrow()
chinese_splot %>% nrow()

#### Group-wise plots
group_1 <- c('amino acid metabolism', 'N-metabolism', 'nucleotide metabolism', 'polyamine metabolism', 'protein')
group_2 <- c('Biodegradation of Xenobiotics', 'metal handling', 'stress')
group_3 <- c('C1-metabolism', 'CHO metabolism', 'fermentation', 'glycolysis', 'gluconeogenesis / glyoxylate cycle', 'OPP')
group_4 <- c('mitochondrial electron transport / ATP synthesis', 'PS', 'redox', 'TCA / org')
group_5 <- c('cell wall', 'development', 'signalling', 'hormone metabolism', 'transport')
group_6 <- c('DNA', 'RNA', 'micro RNA')
group_7 <- c('tetrapyrrole synthesis', 'secondary metabolism', 'S-assimilation', 'lipid metabolism', 'Co-factor and vitamine metabolism')

blues_1507 <- c('mitochondrial electron transport / ATP synthesis', 'N-metabolism', 'not assigned', '')

delta_splot %>% head()

ggplot(chinese_splot[chinese_splot$MapMan_Category == 'storage protein',], aes(x = value, y = value, alpha  = I(0.005), size = I(3))) + 
  # geom_text(aes(label = ifelse(Zhewan1_10_mean >= 75000 & Zhewan1_25_mean >= 75000, as.character(BLASTP_names), '' )))
  # geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan1_25_mean, shape = 'Zhongwan6'), alpha = 0.3) +
  geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, shape = 'Sprint2'), alpha = 0.4, col = '#C853F1') +
  # geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean, shape = 'Zhewan1'), alpha = 0.4) +
  xlab(label = 'Earlier condition') + ylab('Latter condition') +
  guides(alpha = 'none') +
  scale_y_continuous(labels = scales::comma)

# %>% filter_at(vars(-c(target_id, MapMan_Category, BLASTP_names)), all_vars(.>=-1000 & .<=1000))
delta_splot %>% head()


ggplot(delta_splot %>% filter(MapMan_Category == 'protein', MapMan_Category2 == 'targeting'),
       aes(x = deltaSZhe_early, y = deltaSZhe_late)) + 
  # geom_text(aes(label = ifelse(Zhewan1_10_mean >= 75000 & Zhewan1_25_mean >= 75000, as.character(BLASTP_names), '' )))
  geom_point(alpha = 0.4, size = I(2.5), col = '#C853F1') +
  labs(title = 'Δtpm of protein targeting-related transcripts between Zhewan1 and Sprint2', x = 'Zhewan1_10_mean - Sprint2_10_mean',
       y ='Zhewan1_25_mean - Sprint2_30_mean') +
  guides(alpha = 'none') +
  scale_y_continuous(labels = scales::comma) + theme(plot.title = element_text(hjust = 0.5))

delta_splot %>% filter(MapMan_Category == 'amino acid metabolism') %>% filter(deltaSZho_early < -3000)

alab_full$MapMan_Category2 <- chinese_splot$MapMan_Category2
chinese_splot %>% select(target_id, MapMan_Category2) %>% head()
alab_full %>% select(target_id, MapMan_Category2) %>% head()
PoisHeatMap(alab_full %>% filter(MapMan_Category == 'protein', MapMan_Category2 == 'targeting') %>% select(-c(MapMan_Category,
                                                                                                                  MapMan_Category2)), "Greens")


#### Adding regression lines
Sprint2_lm <- lm(Sprint2_30_mean ~ Sprint2_10_mean, chinese_splot %>% filter(MapMan_Category == 'storage protein', Sprint2_10_mean < 100,
                                                                             Sprint2_30_mean < 20000))
Sprint2_line <- Sprint2_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
Sprint2_lower <- Sprint2_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
Sprint2_upper <- Sprint2_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()

Zhewan_lm <- lm(Zhewan1_25_mean ~ Zhewan1_10_mean, chinese_splot %>% filter(MapMan_Category == 'protein', MapMan_Category2 == 'folding',
                Zhewan1_25_mean < 200, Zhewan1_10_mean < 200))
Zhewan_line <- Zhewan_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
Zhewan_lower <- Zhewan_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
Zhewan_upper <- Zhewan_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()

Zhongwan_lm <- lm(Zhongwan1_25_mean ~ Zhongwan6_10_mean, chinese_splot %>% filter(MapMan_Category == 'stress', MapMan_Category3 == 'heat'))
Zhongwan_line <- Zhongwan_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
Zhongwan_lower <- Zhongwan_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
Zhongwan_upper <- Zhongwan_lm %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()

ggplot(chinese_splot %>% filter(MapMan_Category == 'storage protein', Sprint2_10_mean < 100, Sprint2_30_mean < 20000),
       aes(size = I(3))) + 
  # geom_text(aes(label = ifelse(Zhewan1_10_mean >= 75000 & Zhewan1_25_mean >= 75000, as.character(BLASTP_names), '' )))
  # geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean), col = '#C853F1', alpha = 0.4) +
  geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, col = 'Sprint2'), col = '#C853F1', alpha = 0.4) +
  # geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean), col = '#C853F1', alpha = 0.4) +
  guides(alpha = 'none') +
  scale_y_continuous(labels = scales::comma) +
  geom_line(aes(Sprint2_10_mean, Sprint2_line), size = 1, col = '#228AD0') + geom_ribbon(aes(
    Sprint2_10_mean, ymin = Sprint2_lower, ymax = Sprint2_upper), fill = '#61CDFF', alpha = 0.2) + labs(title='Folding-related transcript distribution in Zhewan1',
  xlab = 'Earlier condition', ylab = 'Latter condition') + theme(plot.title = element_text(hjust = 0.5))

alab_full$MapMan_Category3 <- chinese_splot$MapMan_Category3
PoisHeatMap(alab_full %>% filter(MapMan_Category == 'stress', MapMan_Category3 == 'heat' | MapMan_Category2 == 'heat') %>% select(-c(MapMan_Category,
                                                                                                              MapMan_Category2,MapMan_Category3)), "Greens")

cosine(V2V4, V3V5)
annot_full[grepl(' transcription factor', annot_full$Pfam),] %>% dplyr::select(BLASTP_names) %>% nrow()

##### Building a new MapMan mapping
MapMan_split <- function(string) {
  strsplit(string, damn_this_pattern, perl = T)
}

strsplit('DNA.synthesis/chromatin structure,protein.postranslational modification', damn_this_pattern, perl = T)
?strsplit
damn_this_pattern <- '(?<=,)(?=amino|Biodegradation|C1|cell|Co-factor|development|DNA|ell|fermentation|gluconeogenesis|glycolysis|hormone|lipid|major|metal handling|micro RNA|mino acid|minor|misc|mitochondrial|N-metabolism|not assigned|nucleotide|OPP|polyamine|protein|PS|redox|RNA|S-assimilation|secondary|signalling|stress|TCA|tetrapyrrole|transport)'

# damn_this_pattern <- ',(?=amino)(?=Biodegradation)(?=C1)(?=cell)(?=Co-factor)(?=development)(?=DNA)(?=ell)(?=fermentation)(?=gluconeogenesis)(?=glycolysis)(?=hormone)
#            (?=lipid)(?=major)(?=metal)(?=micro)(?=mino)(?=minor)(?=misc)(?=mitochondrial)(?=N-metabolism)(?=not assigned)(?=nucleotide)
#            (?=OPP)(?=polyamine)(?=protein)(?=PS)(?=redox)(?=RNA)(?=S-assimilation)(?=secondary)(?=signalling)(?=stress)(?=TCA)(?=tetrapyrrole)(?=transport)'


MapMan_df <- data.frame('transcript_id' = sapply(c(1:nrow(annot_full)), function(x) rep(as.character(annot_full$transcript_id[x]), lengths(regmatches(as.character(annot_full$MapMan_terms[x]),
                        gregexpr(damn_this_pattern, as.character(annot_full$MapMan_terms[x]), perl = T))) + 1)) %>% unlist(),
                        'MapMan_Category' = sapply(annot_full$MapMan_terms %>% as.character(), function(x) MapMan_split(x) %>% unlist()) %>% unlist())


MapMan_df %<>% mutate(MapMan_Category = sapply(MapMan_Category, function(x) gsub('^(mino)', 'amino', x))) %>% mutate(MapMan_Category = 
                                      sapply(MapMan_Category, function(x) gsub('^(ell)', 'cell',  x))) %>% mutate(MapMan_Category = sapply(MapMan_Category, function(x) gsub(',$', '', x)))

MapMan_df %<>% mutate('Zhewan1_10_mean' = sapply(c(1:nrow(chinese_splot)), function(x) rep(chinese_splot$Zhewan1_10_mean[x], sum(MapMan_df$transcript_id == chinese_splot$target_id[x]))) %>% unlist(),
                      'Zhewan1_25_mean' = sapply(c(1:nrow(chinese_splot)), function(x) rep(chinese_splot$Zhewan1_25_mean[x], sum(MapMan_df$transcript_id == chinese_splot$target_id[x]))) %>% unlist(),
                      'Zhongwan6_10_mean' = sapply(c(1:nrow(chinese_splot)), function(x) rep(chinese_splot$Zhongwan6_10_mean[x], sum(MapMan_df$transcript_id == chinese_splot$target_id[x]))) %>% unlist(),
                      'Zhongwan6_25_mean' = sapply(c(1:nrow(chinese_splot)), function(x) rep(chinese_splot$Zhongwan1_25_mean[x], sum(MapMan_df$transcript_id == chinese_splot$target_id[x]))) %>% unlist(),
                      'Sprint2_10_mean' = sapply(c(1:nrow(chinese_splot)), function(x) rep(chinese_splot$Sprint2_10_mean[x], sum(MapMan_df$transcript_id == chinese_splot$target_id[x]))) %>% unlist(),
                      'Sprint2_30_mean' = sapply(c(1:nrow(chinese_splot)), function(x) rep(chinese_splot$Sprint2_30_mean[x], sum(MapMan_df$transcript_id == chinese_splot$target_id[x]))) %>% unlist())

MapMan_df %<>% mutate(transcript_id = as.character(transcript_id))
MapMan_df %<>% mutate(MapMan_Category = sapply(MapMan_Category, function(x) gsub('aminor', 'minor', x)))

MapMan_df %<>% mutate(MapMan_Category1 = sapply(MapMan_df$MapMan_Category %>% as.character, function(x) x %>% strsplit('\\.') %>% unlist() %>% first()),
                      MapMan_Category2 = sapply(MapMan_df$MapMan_Category %>% as.character, function(x) x %>% strsplit('\\.') %>% unlist() %>% nth(2)),
                      MapMan_Category3 = sapply(MapMan_df$MapMan_Category %>% as.character, function(x) x %>% strsplit('\\.') %>% unlist() %>% nth(3)),
                      MapMan_Category4 = sapply(MapMan_df$MapMan_Category %>% as.character, function(x) x %>% strsplit('\\.') %>% unlist() %>% nth(4)))

MapMan_df %<>% mutate(ZheEarly_m = Zhewan1_10_mean  - Sprint2_10_mean, ZheLate_m = Zhewan1_25_mean - Sprint2_30_mean,
                      ZheEarly_p  = Zhewan1_10_mean + Sprint2_10_mean, ZheLate_p = Zhewan1_25_mean + Sprint2_30_mean,
                      ZhoEarly_m = Zhongwan6_10_mean - Sprint2_10_mean, ZhoLate_m = Zhongwan6_25_mean - Sprint2_30_mean,
                      ZhoEarly_p = Zhongwan6_10_mean + Sprint2_10_mean, ZhoLate_p = Zhongwan6_25_mean + Sprint2_30_mean)
MapMan_df %<>% mutate(normZheEarly = ZheEarly_m/ZheEarly_p, normZheLate = ZheLate_m/ZheLate_p,
                      normZhoEarly = ZhoEarly_m/ZhoEarly_p, normZhoLate = ZhoLate_m/ZhoLate_p)


MapMan_df$MapMan_Category1 %>% unique()

  ##### Now it gets personal: scatterplots

group_1 <- c('amino acid metabolism', 'N-metabolism', 'nucleotide metabolism', 'polyamine metabolism', 'protein')
group_2 <- c('Biodegradation of Xenobiotics', 'metal handling', 'stress')
group_3 <- c('C1-metabolism', 'CHO metabolism', 'fermentation', 'glycolysis', 'gluconeogenesis / glyoxylate cycle', 'OPP')
group_4 <- c('mitochondrial electron transport / ATP synthesis', 'PS', 'redox', 'TCA / org')
group_5 <- c('cell wall', 'development', 'signalling', 'hormone metabolism', 'transport')
group_6 <- c('DNA', 'RNA', 'micro RNA')
group_7 <- c('tetrapyrrole synthesis', 'secondary metabolism', 'S-assimilation', 'lipid metabolism', 'Co-factor and vitamine metabolism')
group_8 <- c('misc', 'not assigned', 'cell')

pea_scatter_plotter
pea_scatter_plotter <- function(category1, category2 = NA, category3 = NA, category4 = NA, uprmax =NA){
  if (!is.na(category4)) {
    x <- quo(MapMan_Category4 == category4)
    plot_title <- paste(category1, category4, sep = ':')
  } else if(!is.na(category3)){
    x <- quo(MapMan_Category3 == category3)
    plot_title <- paste(category1, category3, sep = ':')
  } else if(!is.na(category2)) {
    x <- quo(MapMan_Category2 == category2)
    plot_title <- paste(category1, category2, sep = ':')
  } else {
    x <- quo(MapMan_Category1 == category1)
    plot_title <- category1
  }
  MapMan_df_loc <- MapMan_df %>% filter(!!x)
  
  if (!is.na(uprmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.<=uprmax))
  }
  
  lm_Zhe <- lm(Zhewan1_25_mean ~ Zhewan1_10_mean, MapMan_df_loc)
  Zhewan_line <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
  Zhewan_lower <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
  Zhewan_upper <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()
  
  lm_Zho <- lm(Zhongwan6_25_mean ~ Zhongwan6_10_mean, MapMan_df_loc)
  Zhongwan_line <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
  Zhongwan_lower <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
  Zhongwan_upper <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()
  
  lm_Spr <- lm(Sprint2_30_mean ~ Sprint2_10_mean, MapMan_df_loc)
  Sprint2_line <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
  Sprint2_lower <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
  Sprint2_upper <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()
  
  ggplot(MapMan_df_loc, aes(x = value, y = value, alpha  = 0.1, size = I(3))) + 
    # geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan6_25_mean, shape = 'Zhongwan6'), alpha = 0.4, col = '#dfc9e2') +
    # geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, shape = 'Sprint2'), alpha = 0.4, col = '#228AD0') +
    # geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean, shape = 'Zhewan1'), alpha = 0.4, col = '#00b6c2') +
    xlab(label = 'Earlier condition') + ylab('Latter condition')+
    guides(alpha = 'none')+
  # geom_line(aes(Sprint2_10_mean, Sprint2_line), size = 0.8, col = '#228AD0', inherit.aes = F) + geom_ribbon(aes(
      # Sprint2_10_mean, ymin = Sprint2_lower, ymax = Sprint2_upper), fill = '#61CDFF', alpha = 0.3, inherit.aes = F) +
  # geom_line(aes(Zhewan1_10_mean, Zhewan_line), size = 0.8, col = '#00b6c2', inherit.aes = F) + geom_ribbon(aes(
  #     Zhewan1_10_mean, ymin = Zhewan_lower, ymax = Zhewan_upper), fill = '#ccffff', alpha = 0.3, inherit.aes = F)+
  # geom_line(aes(Zhongwan6_10_mean, Zhongwan_line), size = 0.8, col = '#dfc9e2', inherit.aes = F) + geom_ribbon(aes(
      # Zhongwan6_10_mean, ymin = Zhongwan_lower, ymax = Zhongwan_upper), fill = '#e3e8ff', alpha = 0.3, inherit.aes = F) +
    labs(title = paste('Transcript distribution (MapMan Category: ', plot_title, ')', sep = ''), shape = 'Line') + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 13),
                       legend.text = element_text(size = 11))
    
}

MapMan_df %>% select(MapMan_Category1) %>% unlist() %>% unique()

pea_scatter_plotter('protein', category2 = 'postranslational modification')

outliers <- data.frame()
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'protein') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 2000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'development') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                             Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'cell wall') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                 Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'major CHO metabolism') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                               Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'misc') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 2500)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 400)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'lipid metabolism') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'minor CHO metabolism') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 200)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'misc') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1500)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'not assigned') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'major CHO metabolism') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'signalling') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'stress') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'transport') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                                          Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 1000)))
outliers %<>% rbind(MapMan_df %>% filter(MapMan_Category1 == 'RNA', MapMan_Category2 == 'regulation of transcription') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                                               Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 500)))
MapMan_df %>% filter(MapMan_Category2 == 'biotic') %>% select(MapMan_Category3) %>% unlist() %>% unique()


outliers %<>% distinct(transcript_id, .keep_all = T)
local_dum <- annot_full %>% filter(transcript_id %in% outliers$transcript_id) %>% select(transcript_id, BLASTP_names, BLASTX_names, plants_BLASTP, plants_BLASTX, MapMan_terms)
local_dum <- local_dum[match(outliers$transcript_id, local_dum$transcript_id),]
outliers %<>% cbind(local_dum[,c(6)])
colnames(outliers)[29] <- 'MapMan_all_categories'
outliers %<>% select(c(1:8, 21:29))


MapMan_df %>% filter(MapMan_Category1 == 'not assigned') %>% nrow()
MapMan_df %>% filter(MapMan_Category1 == 'signalling') %>% filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean,
                                                                      Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), any_vars(.> 14000))
write.table(outliers, '~/Pea_transcriptomics/250719_pics/outliers.tsv', col.names = T, row.names = F, sep = '\t')

# lm_Zhe <- lm(Zhewan1_25_mean ~ Zhewan1_10_mean, MapMan_df %>% filter(MapMan_Category2 == 'storage proteins' | transcript_id %in% storage_proteins$target_id) %>% 
#   filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean, Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), all_vars(.<=25000)))
# Zhewan_line <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
# Zhewan_lower <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
# Zhewan_upper <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()
# 
# lm_Zho <- lm(Zhongwan6_25_mean ~ Zhongwan6_10_mean, MapMan_df %>% filter(MapMan_Category2 == 'storage proteins' | transcript_id %in% storage_proteins$target_id) %>% 
#                filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean, Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), all_vars(.<=25000)))
# Zhongwan_line <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
# Zhongwan_lower <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
# Zhongwan_upper <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()
# 
# lm_Spr <- lm(Sprint2_30_mean ~ Sprint2_10_mean, MapMan_df %>% filter(MapMan_Category2 == 'storage proteins' | transcript_id %in% storage_proteins$target_id) %>% 
#                filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean, Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), all_vars(.<=25000)))
# Sprint2_line <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% select(fit) %>% unlist()
# Sprint2_lower <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% select(lwr) %>% unlist()
# Sprint2_upper <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% select(upr) %>% unlist()
# 
# ggplot(MapMan_df %>% filter(MapMan_Category2 == 'storage proteins' | transcript_id %in% storage_proteins$target_id) %>% 
#          filter_at(vars(c(Zhewan1_10_mean, Zhewan1_25_mean, Zhongwan6_10_mean, Zhongwan6_25_mean, Sprint2_10_mean, Sprint2_30_mean)), all_vars(.<=25000)),
#        aes(x = value, y = value, alpha  = 0.1, size = I(3))) + 
#   geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan6_25_mean, shape = 'Zhongwan6'), alpha = 0.4, col = '#dfc9e2') +
#   geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, shape = 'Sprint2'), alpha = 0.4, col = '#228AD0') +
#   geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean, shape = 'Zhewan1'), alpha = 0.4, col = '#00b6c2') +
#   xlab(label = 'Earlier condition') + ylab('Latter condition') +
#   guides(alpha = 'none') +
#   scale_y_continuous(labels = scales::comma) +
#   guides(alpha = 'none')+
#   geom_line(aes(Sprint2_10_mean, Sprint2_line), size = 0.8, col = '#228AD0', inherit.aes = F) + geom_ribbon(aes(
#     Sprint2_10_mean, ymin = Sprint2_lower, ymax = Sprint2_upper), fill = '#61CDFF', alpha = 0.3, inherit.aes = F)+
#   geom_line(aes(Zhewan1_10_mean, Zhewan_line), size = 0.8, col = '#00b6c2', inherit.aes = F) + geom_ribbon(aes(
#     Zhewan1_10_mean, ymin = Zhewan_lower, ymax = Zhewan_upper), fill = '#ccffff', alpha = 0.3, inherit.aes = F)+
#   geom_line(aes(Zhongwan6_10_mean, Zhongwan_line), size = 0.8, col = '#dfc9e2', inherit.aes = F) + geom_ribbon(aes(
#     Zhongwan6_10_mean, ymin = Zhongwan_lower, ymax = Zhongwan_upper), fill = '#e3e8ff', alpha = 0.3, inherit.aes = F) +
#   labs(title = paste('Transcript distribution of storage proteins', sep = ''), shape = 'Line') + 
#   theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 13),
#                      legend.text = element_text(size = 11))


####### deltaplots
pea_deltaplots <- function(category1, category2 = NA, category3 = NA, category4 = NA, uprmax = NA, lwrmax = NA) {
  MapMan_df_loc <- MapMan_df %>% select(transcript_id, ZheEarly_m, ZheLate_m, ZhoEarly_m, ZhoLate_m, MapMan_Category, MapMan_Category1, MapMan_Category2,
                                        MapMan_Category3, MapMan_Category4)
  if (!is.na(category4)) {
    x <- quo(MapMan_Category4 == category4)
    plot_title <- paste(category1, category4, sep = ':')
  } else if(!is.na(category3)){
    x <- quo(MapMan_Category3 == category3)
    plot_title <- paste(category1, category3, sep = ':')
  } else if(!is.na(category2)) {
    x <- quo(MapMan_Category2 == category2)
    plot_title <- paste(category1, category2, sep = ':')
  } else {
    x <- quo(MapMan_Category1 == category1)
    plot_title <- category1
  }
  MapMan_df_loc %<>% filter(!!x)
  
  if (!is.na(uprmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.<=uprmax))
  }
  if (!is.na(lwrmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.>=lwrmax))
  }
  
  ggplot(MapMan_df_loc ,
         aes(x = value, y = value, size = I(3))) + 
    geom_point(aes(x = ZheEarly_m, y = ZheLate_m, shape = 'Zhewan1'), alpha = 0.4, size = I(2.5), col = '#00b6c2', inherit.aes = F) +
    geom_point(aes(x = ZhoEarly_m, y = ZhoLate_m, shape = 'Zhongwan6'), alpha = 0.4, size = I(2.5), col = '#dfc9e2', inherit.aes = F) +
    labs(title = paste('transcript Δtpm (MapMan Category: ', plot_title, ')', sep = ''), x = 'Early', y ='Late') +
    guides(alpha = 'none') +
    scale_y_continuous(labels = scales::comma) +  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 13),
    legend.text = element_text(size = 11))
}

pea_deltaplots('RNA', 'regulation of transcription', 'HSF,Heat-shock transcription factor family')
MapMan_df %>% filter(MapMan_Category1 == 'RNA', MapMan_Category2 == 'regulation of transcription') %>% select(MapMan_Category3) %>% unlist() %>% unique()

annot_full %>% filter(transcript_id == 'TRINITY_10_DN8344_c0_g1_i1') %>% select(BLASTP_names, BLASTX_names)


pea_normdeltaplots <- function(category1, category2 = NA, category3 = NA, category4 = NA, uprmax = NA, lwrmax = NA) {
  MapMan_df_loc <- MapMan_df %>% select(transcript_id, normZheEarly, normZheLate, normZhoEarly, normZhoLate, MapMan_Category, MapMan_Category1, MapMan_Category2,
                                        MapMan_Category3, MapMan_Category4)
  if (!is.na(category4)) {
    x <- quo(MapMan_Category4 == category4)
    plot_title <- paste(category1, category4, sep = ':')
  } else if(!is.na(category3)){
    x <- quo(MapMan_Category3 == category3)
    plot_title <- paste(category1, category3, sep = ':')
  } else if(!is.na(category2)) {
    x <- quo(MapMan_Category2 == category2)
    plot_title <- paste(category1, category2, sep = ':')
  } else {
    x <- quo(MapMan_Category1 == category1)
    plot_title <- category1
  }
  MapMan_df_loc %<>% filter(!!x)
  
  if (!is.na(uprmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.<=uprmax))
  }
  if (!is.na(lwrmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.>=lwrmax))
  }
  
  ggplot(MapMan_df_loc ,
         aes(x = value, y = value, size = I(3))) + 
    geom_point(aes(x = normZheEarly, y = normZheLate, shape = 'Zhewan1'), alpha = 0.3, size = I(2.5), col = '#00b6c2', inherit.aes = F) +
    geom_point(aes(x = normZhoEarly, y = normZhoLate, shape = 'Zhongwan6'), alpha = 0.3, size = I(2.5), col = '#dfc9e2', inherit.aes = F) +
    labs(title = paste('Normalized Δtpm (MapMan Category: ', plot_title, ')', sep = ''), x = 'Early', y ='Late') +
    guides(alpha = 'none') +
    scale_y_continuous(labels = scales::comma) +  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 13),
                                                          legend.text = element_text(size = 11))
}

pea_normdeltaplots('RNA', 'regulation of transcription', 'HSF,Heat-shock transcription factor family')

MapMan_df %>% filter(MapMan_Category1 == 'RNA', MapMan_Category2 == 'regulation of transcription') %>% select(MapMan_Category3) %>% unlist() %>% length()


############## Transcription factors
TF_df <- annot_full[grepl('transcription factor', annot_full$MapMan_terms),] %>% select(transcript_id, BLASTP_names, BLASTX_names, Pfam, MapMan_terms)
TF_df1 <- annot_full[grepl('transcription factor', annot_full$Pfam),] %>% filter(!(transcript_id %in% TF_df$transcript_id)) %>% select(transcript_id, BLASTP_names, BLASTX_names, Pfam, MapMan_terms)
TF_df %<>% rbind(TF_df1)
rm(TF_df1)
TF_df1 <- annot_full[grepl('(GeBP like)|(\\.ARR)|(CCAAT box binding factor family)|(C2C2\\(Zn\\) YABBY family)|(C2C2\\(Zn\\) DOF zinc finger family)|(ovate family OFP)|
                           (BBR/BPC)|(plant TF \\(pbf2\\))|(BSD domain containing family)|(AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family)|
                           (zf-HD)|(sigma like plant)|(ARF, Auxin Response Factor family)|(ARF, Auxin Response Factor family)|(C2C2\\(Zn\\) CO-like, Constans-like zinc finger family)|
                           (CCAAT box binding factor family, HAP2)|(NIN-like bZIP-related family)|(GRF zinc finger family)|(Alfin-like)|(ELF3)',
                           annot_full$MapMan_terms),] %>% filter(!(transcript_id %in% TF_df$transcript_id)) %>% select(transcript_id, BLASTP_names, BLASTX_names, Pfam, MapMan_terms)
TF_df %<>% rbind(TF_df1)
TF_df %>% nrow()
rm(TF_df1)

annot_full[grepl('\\.ARR', annot_full$MapMan_terms),] %>% select(MapMan_terms)

#### Graph vector summation
### Sprint2_10 -> Sprint2_20
S10S20 <- c(rep(0, nrow(annot_full)))
names(S10S20) <- annot_full$transcript_id

for (i in c(1:length(S10S20))) {
  if (names(S10S20)[i]  %in% upregulated_30_days$target_id) {
    S10S20[i] <- 1
  }
  else if (names(S10S20)[i] %in% downregulated_30_days$target_id) {
    S10S20[i] <- -1
  }
}

### Sprint2_20 -> Zhewan1_25
S20Zhe25 <- c(rep(0, nrow(annot_full)))
names(S20Zhe25) <- annot_full$transcript_id
for (i in c(1:length(S20Zhe25))) {
  if (names(S20Zhe25)[i]  %in% cross_4_down$target_id) {
    S20Zhe25[i] <- 1
  }
  else if (names(S20Zhe25)[i] %in% cross_4_up$target_id) {
    S20Zhe25[i] <- -1
  }
}

### Sprint2_20 -> Zhongwan6_25
S20Zho25 <- c(rep(0, nrow(annot_full)))
names(S20Zho25) <- annot_full$transcript_id
for (i in c(1:length(S20Zho25))) {
  if (names(S20Zho25)[i]  %in% cross_5_down$target_id) {
    S20Zho25[i] <- 1
  }
  else if (names(S20Zho25)[i] %in% cross_5_up$target_id) {
    S20Zho25[i] <- -1
  }
}


#### Sprint2_10 -> Zhewan1_25
S10Zhe25 <- c(rep(0, nrow(annot_full)))
names(S10Zhe25) <- annot_full$transcript_id
for (i in c(1:length(S10Zhe25))) {
  if (names(S10Zhe25)[i]  %in% cross_6_up$target_id) {
    S10Zhe25[i] <- 1
  }
  else if (names(S10Zhe25)[i] %in% cross_6_down$target_id) {
    S10Zhe25[i] <- -1
  }
}

#### Sprint2_10 -> Zhongwan6_25
S10Zho25 <- c(rep(0, nrow(annot_full)))
names(S10Zho25) <- annot_full$transcript_id
for (i in c(1:length(S10Zho25))) {
  if (names(S10Zho25)[i]  %in% cross_7_up$target_id) {
    S10Zho25[i] <- 1
  }
  else if (names(S10Zho25)[i] %in% cross_7_down$target_id) {
    S10Zho25[i] <- -1
  }
}


### Zhewan1_25 -> Sprint2_20 
Zhe25S20 <- c(rep(0, nrow(annot_full)))
names(Zhe25S20) <- annot_full$transcript_id
for (i in c(1:length(Zhe25S20))) {
  if (names(Zhe25S20)[i]  %in% cross_4_down$target_id) {
    Zhe25S20[i] <- -1
  }
  else if (names(Zhe25S20)[i] %in% cross_4_up$target_id) {
    Zhe25S20[i] <- 1
  }
}

### Sprint2_20 -> Zhongwan6_25
Zho25S20 <- c(rep(0, nrow(annot_full)))
names(Zho25S20) <- annot_full$transcript_id
for (i in c(1:length(Zho25S20))) {
  if (names(Zho25S20)[i]  %in% cross_5_down$target_id) {
    Zho25S20[i] <- -1
  }
  else if (names(Zho25S20)[i] %in% cross_5_up$target_id) {
    Zho25S20[i] <- 1
  }
}

###### Zhe10 -> Zhe25
Zhe10Zhe25 <- c(rep(0, nrow(annot_full)))
names(Zhe10Zhe25) <- annot_full$transcript_id
for (i in c(1:length(Zhe10Zhe25))) {
  if (names(Zhe10Zhe25)[i]  %in% zhe_down$target_id) {
    Zhe10Zhe25[i] <- -1
  }
  else if (names(Zhe10Zhe25)[i] %in% zhe_up$target_id) {
    Zhe10Zhe25[i] <- 1
  }
}

Zho10Zho25 <- c(rep(0, nrow(annot_full)))
names(Zho10Zho25) <- annot_full$transcript_id
for (i in c(1:length(Zho10Zho25))) {
  if (names(Zho10Zho25)[i]  %in% zhong_down$target_id) {
    Zho10Zho25[i] <- -1
  }
  else if (names(Zho10Zho25)[i] %in% zhong_up$target_id) {
    Zho10Zho25[i] <- 1
  }
}


one <- S10S20 + S20Zhe25
two <- S10S20 + S20Zho25
three <- S10Zhe25 + Zhe25S20
four <- S10Zho25 + Zho25S20
cosines <- data.frame('one' = c(cosine(one, one), cosine(one, two), cosine(one, three), cosine(one, four)),
                      'two' = c(cosine(two, one), cosine(two, two), cosine(two, three), cosine(two, four)),
                      'three' = c(cosine(three, one), cosine(three, two), cosine(three, three), cosine(three, four)),
                      'four' = c(cosine(four, one), cosine(four, two), cosine(four, three), cosine(four, four)),
                      row.names = c('one', 'two', 'three', 'four'))


data.frame('S-10->S-20' = c(cosine(S10S20, S10S20), cosine(S10S20, Zhe10Zhe25), cosine(S10S20, Zho10Zho25)),
           'C1-10->C1-25' = c(cosine(S10S20, Zhe10Zhe25), cosine(Zhe10Zhe25, Zhe10Zhe25), cosine(Zhe10Zhe25, Zho10Zho25)),
           'C2-10->C2-25' = c(cosine(Zho10Zho25, S10S20), cosine(Zhe10Zhe25, Zho10Zho25), cosine(Zho10Zho25, Zho10Zho25)),
           row.names = c('S-10->S-20', 'C1-10->C1-25', 'C2-10->C2-25'))



####### SOM
library(oposSOM)
env_test <- opossom.new(list(dataset.name="Tissues",
                          dim.1stLvlSom=20))
data(opossom.tissues)
str(opossom.tissues, vec.len=3)
head(opossom.tissues)

dummy <- MapMan_df %>% filter(MapMan_Category1 == 'protein', MapMan_Category2 == 'synthesis') %>% dplyr::select(c(1, 3:8)) %>% mutate(transcript_id = 
        as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE)
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
colnames(dummy) <- c('Zhewan1_10', 'Zhewan1_25', 'Zhongwan6_10', 'Zhongwan6_25', 'Sprint2_10', 'Sprint2_20')
env_synthesis <- opossom.new(list(dataset.name = 'protein:synthesis', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Zhewan 1', 2), rep('Zhongwan 6', 2), rep('Sprint 2', 2))
env_synthesis$group.colors <- c(rep('#00b6c2', 2), rep('#dfc9e2', 2), rep('#228AD0', 2))
opossom.run(env = env_synthesis)

alab_full %>% head()
raw_for_oposSOM <- alab_full %>% dplyr::select(-c("MapMan_Category", "MapMan_Category2", "MapMan_Category3"))
raw_for_oposSOM$MapMan <- annot_full[match(raw_for_oposSOM$target_id, annot_full$transcript_id),]$MapMan_terms
raw_for_oposSOM %<>% mutate(target_id = as.character(target_id))
raw_for_oposSOM_processed <- MapMan_df %>% dplyr::select(transcript_id, MapMan_Category, MapMan_Category1, MapMan_Category2, MapMan_Category3)
raw_for_oposSOM_processed %<>% cbind(raw_for_oposSOM[rep(seq_len(nrow(raw_for_oposSOM)), temporary_list), c(2:25)])
rm(raw_for_oposSOM)


dummy <- raw_for_oposSOM_processed %>% filter(transcript_id %in% annot_full[grepl('GO:0009738', annot_full$GOs),]$transcript_id) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                   as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(5:8)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4")
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'ABA', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_ABA/')
opossom.run(env = env_synthesis)

annot_full[grepl(tf_template, annot_full$MapMan_terms), "MapMan_terms"] %>% head()
dummy <- raw_for_oposSOM_processed %>% filter(transcript_id %in% annot_full[grepl(tf_template, annot_full$MapMan_terms),]$transcript_id) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                                                                                                          as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(5:8)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4")
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'TF', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_TF/')
opossom.run(env = env_synthesis)

dummy <- raw_for_oposSOM_processed %>% filter(transcript_id %in% annot_full[grepl('regulation of transcription', annot_full$MapMan_terms),]$transcript_id) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                                                                                                                                                   as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(5:8)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4")
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'Transcription', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_transcription_all/')
opossom.run(env = env_synthesis)

MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism') %>% dplyr::select(MapMan_Category2) %>% unlist() %>% unname() %>% unique()

dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 == 'hormone metabolism', MapMan_Category2 == 'salicylic acid') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
      as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(5:8)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4")
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'salicylic acid', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_salicylate/')
opossom.run(env = env_synthesis)

dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 == 'hormone metabolism') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
dummy %>% nrow()
colnames(dummy)[c(5:8)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4")
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'Hormones', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_hormones/')
opossom.run(env = env_synthesis)



raw_for_oposSOM <- alab_full %>% dplyr::select(-c("MapMan_Category", "MapMan_Category2", "MapMan_Category3"))
raw_for_oposSOM$MapMan <- annot_full[match(raw_for_oposSOM$target_id, annot_full$transcript_id),]$MapMan_terms
raw_for_oposSOM %<>% mutate(target_id = as.character(target_id))
raw_for_oposSOM_processed <- MapMan_df %>% dplyr::select(transcript_id, MapMan_Category, MapMan_Category1, MapMan_Category2, MapMan_Category3)
raw_for_oposSOM_processed %<>% cbind(raw_for_oposSOM[rep(seq_len(nrow(raw_for_oposSOM)), temporary_list), c(2:25)])
rm(raw_for_oposSOM)

dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 == 'protein', MapMan_Category2 == 'synthesis') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                                                                                                                        as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy %>% head()
env_synthesis <- opossom.new(list(dataset.name = 'protein:synthesis', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)

MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism') %>% dplyr::select(MapMan_Category3) %>% unlist() %>% unname() %>% unique()



##### Full protein
setwd('~/test_oposSOM_meaningful_clusterization')
dummy %>% nrow()
dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 == 'protein') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 

            as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy %>% head()
env_synthesis <- opossom.new(list(dataset.name = 'protein', dim.1stLvlSom = 30))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)

### same with kohonen package
dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 == 'protein')%>% mutate(transcript_id = 
 as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% `rownames<-`(.[,1])  %>%  
  dplyr::select(c(6:29)) %>% scale() %>% as.matrix()
dummy %>% head()
dummy_model <- som(dummy, somgrid(xdim=6, ydim=6, topo = 'hexagonal'), rlen=10000, alpha = c(0.05,0.01), keep.data = T)
plot(dummy_model, type="changes")
plot(dummy_model, type="count", main="Node Counts")
plot(dummy_model, type="dist.neighbours", main = "SOM neighbour distances")

pea_scatter_plotter('protein', category2 = 'postranslational modification', uprmax = 300)



dummy <- raw_for_oposSOM_processed %>% mutate(transcript_id = 
            as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% `rownames<-`(.[,1])  %>%  
  dplyr::select(c(6:29))  %>% scale() %>% as.matrix()

setwd('~/oposSOM_full_trlist/')
env_synthesis <- opossom.new(list(dataset.name = 'all', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)


setwd('~/Pea_transcriptomics/')


###### All transcript SOM - exclude 
group_A <- read.table('/home/reverend_casy/Pea_transcriptomics/Group_Overexpression_Spots_A.csv',
                      row.names = 1, header = T,  sep = '\t')
group_B <- read.table('/home/reverend_casy/Pea_transcriptomics/Group_Overexpression_Spots_B.csv',
                      row.names = 1, header = T,  sep = '\t') %>% dplyr::select(-c(Symbol))
group_C <- read.table('/home/reverend_casy/Pea_transcriptomics/Group_Overexpression_Spots_C.csv',
                      row.names = 1, header = T,  sep = '\t') %>% dplyr::select(-c(Symbol))
group_A %>% View()
group_B %>% View()
group_C %>% View()
colnames(group_A) <- colnames(group_C)
groups_all <- rbind(group_A, group_B, group_C)
groups_all$BLASTP <- annot_full %>% filter(transcript_id %in% all_groups) %>% dplyr::select(BLASTP_names) %>% unlist() %>% unname()
groups_all$BLASTP_code <- annot_full %>% filter(transcript_id %in% all_groups) %>% dplyr::select(plants_BLASTP) %>% 
  unlist() %>% unname() %>% sapply(function(x) strsplit(x, '\\^') %>% unlist() %>% unname() %>% dplyr::first()) %>% unname()


groups_all %>% View()
groups_all %<>% mutate(target_id = rownames(.))
group_A %<>% mutate(target_id = rownames(group_A))
group_B %<>% mutate(target_id = rownames(group_B))
group_C %<>% mutate(target_id = rownames(group_C))


sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, ont = x, cutoff = 0.01))

BP_group_C %>% View()
CC_group_B %>% View()
MF_group_C %>% View()

all_groups <- c(rownames(group_A), rownames(group_B), rownames(group_C)) %>% unique()


dummy <- raw_for_oposSOM_processed %>% mutate(transcript_id = 
                                                as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% 
  filter(!(transcript_id %in% all_groups)) %>% 
  `rownames<-`(.[,1])  %>%  
  dplyr::select(c(6:29))  %>% scale() %>% as.matrix()

  
setwd('~/oposSOM_120819//')
env_synthesis <- opossom.new(list(dataset.name = 'all', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)
write.table(groups_all, '~/oposSOM_120819/excluded_genes.tsv', col.names = T, row.names = T, sep = '\t', dec = '.')
setwd('~/Pea_transcriptomics/')

###### 3k
dummy <- raw_for_oposSOM_processed %>% filter_at(c(6:29), all_vars(.<300)) %>% 
  mutate(transcript_id = 
                                                as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% 
  filter(!(transcript_id %in% all_groups)) %>% 
  `rownames<-`(.[,1])  %>%  
  dplyr::select(c(6:29))  %>% scale() %>% as.matrix()


setwd('~/oposSOM_120819_0.3k')
env_synthesis <- opossom.new(list(dataset.name = 'all', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)
write.table(groups_all, '~/oposSOM_120819/excluded_genes.tsv', col.names = T, row.names = T, sep = '\t', dec = '.')
setwd('~/Pea_transcriptomics/')




sapply(ls()[grep('((BP)|(CC)|(MF))_group_(A|B|C)$', ls())], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'group_annotations',
                  'all_transcripts', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.'))

group_A2 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819/Group Overexpression Spots A.csv',
                      row.names = 1, header = T,  sep = ';')
group_B2 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819/Group Overexpression Spots B.csv',
                       row.names = 1, header = T,  sep = ';')
group_C2 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819/Group Overexpression Spots C.csv',
                       row.names = 1, header = T,  sep = ';')
group_D2 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819/Group Overexpression Spots D.csv',
                       row.names = 1, header = T,  sep = ';')


groups_all2 <- rbind(group_A2, group_B2, group_C2, group_D2)
all_groups2 <- c(rownames(group_A2), rownames(group_B2), rownames(group_C2), rownames(group_D2))

groups_all2$BLASTP <- annot_full %>% filter(transcript_id %in% all_groups2) %>% dplyr::select(BLASTP_names) %>% unlist() %>% unname()
groups_all2$BLASTP_code <- annot_full %>% filter(transcript_id %in% all_groups2) %>% dplyr::select(plants_BLASTP) %>% 
  unlist() %>% unname() %>% sapply(function(x) strsplit(x, '\\^') %>% unlist() %>% unname() %>% dplyr::first()) %>% unname()

group_A2 %<>% mutate(target_id = rownames(group_A2))
group_B2 %<>% mutate(target_id = rownames(group_B2))
group_C2 %<>% mutate(target_id = rownames(group_C2))
group_D2 %<>% mutate(target_id = rownames(group_D2))

sapply(ontologies, function(x) iterating_function_GO(group_A2, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_B2, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_C2, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_D2, ont = x, cutoff = 0.01))

sapply(ls()[grep('((BP)|(CC)|(MF))_group_(A|B|C|D)2', ls())], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'group_annotations',
                                                                                                '120819', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))

group_A3 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819_0.3k//Group Overexpression Spots A.csv',
                       row.names = 1, header = T,  sep = ';')
group_B3 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819_0.3k//Group Overexpression Spots B.csv',
                       row.names = 1, header = T,  sep = ';')
group_C3 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819_0.3k//Group Overexpression Spots C.csv',
                       row.names = 1, header = T,  sep = ';')
group_D3 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819_0.3k//Group Overexpression Spots D.csv',
                       row.names = 1, header = T,  sep = ';')
group_E3 <- read.table('/home/reverend_casy/Pea_transcriptomics/group_annotations/120819_0.3k//Group Overexpression Spots E.csv',
                       row.names = 1, header = T,  sep = ';')

groups_all3 <- rbind(group_A3, group_B3, group_C3, group_D3, group_E3)
all_groups3 <- c(rownames(group_A3), rownames(group_B3), rownames(group_C3), rownames(group_D3), rownames(group_E3))

groups_all3$BLASTP <- annot_full %>% filter(transcript_id %in% all_groups3) %>% dplyr::select(BLASTP_names) %>% unlist() %>% unname()
groups_all3$BLASTP_code <- annot_full %>% filter(transcript_id %in% all_groups3) %>% dplyr::select(plants_BLASTP) %>% 
  unlist() %>% unname() %>% sapply(function(x) strsplit(x, '\\^') %>% unlist() %>% unname() %>% dplyr::first()) %>% unname()

group_A3 %<>% mutate(target_id = rownames(group_A3))
group_B3 %<>% mutate(target_id = rownames(group_B3))
group_C3 %<>% mutate(target_id = rownames(group_C3))
group_D3 %<>% mutate(target_id = rownames(group_D3))
group_E3 %<>% mutate(target_id = rownames(group_E3))

sapply(ontologies, function(x) iterating_function_GO(group_A3, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_B3, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_C3, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_D3, ont = x, cutoff = 0.01))
sapply(ontologies, function(x) iterating_function_GO(group_E3, ont = x, cutoff = 0.01))

sapply(ls()[grep('((BP)|(CC)|(MF))_group_(A|B|C|D|E)3', ls())], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'group_annotations',
                                                                                           '120819_0.3k', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))

write.table(groups_all, '~/Pea_transcriptomics/group_annotations/all_transcripts/pivot.tsv', row.names = T, col.names = T, sep = '\t')
write.table(groups_all2, '~/Pea_transcriptomics/group_annotations/120819/pivot.tsv', row.names = T, col.names = T, sep = '\t')
write.table(groups_all3, '~/Pea_transcriptomics/group_annotations/120819_0.3k//pivot.tsv', row.names = T, col.names = T, sep = '\t')



 


##### SOM - no stress-related proteins

setwd('~/Pea_transcriptomics/oposSOM_no_stress/')
dummy %>% nrow()
dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 != 'stress') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                                                                                       as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))

rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy %>% head()
env_synthesis <- opossom.new(list(dataset.name = 'no_stress', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)

setwd('~/Pea_transcriptomics/oposSOM_no_stress_no_PH/')
dummy %>% nrow()
dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 != 'stress', MapMan_Category1 != 'PS') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                                                                                      as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
raw_for_oposSOM_processed %>% dplyr::select(MapMan_Category1) %>% unlist() %>% unique()
rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy %>% head()
env_synthesis <- opossom.new(list(dataset.name = 'no_stress_no_PH', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)

###### 
setwd('~/Pea_transcriptomics/oposSOMs/oposSOM_folding/')

raw_for_oposSOM_processed[!(grepl('stress', raw_for_oposSOM_processed$MapMan_Category)),] %>% filter(MapMan_Category1 == 'protein') %>%
  dplyr::select(MapMan_Category2) %>% unlist() %>% unique()

dummy %>% nrow()
dummy[dummy == Inf] %>% nrow()
dummy <- raw_for_oposSOM_processed %>%
  filter(MapMan_Category1 == 'protein', MapMan_Category2 == 'folding') %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id =
                      as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
dummy[dummy == Inf] <- 0
dummy[dummy == -Inf] <- 0
# dummy <- raw_for_oposSOM_processed %>%
#   filter(transcript_id %in% storage_proteins$target_id) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
#                                                                                                                    as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))

rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))

env_synthesis <- opossom.new(list(dataset.name = 'folding', dim.1stLvlSom = 20))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)

raw_for_oposSOM_processed %>% filter(transcript_id == 'TRINITY_10_DN23187_c0_g1_i1')
dummy %>% View()
dummy %>% filter_at(vars(colnames(dummy)), all_vars(. == 0))

setwd('~/Pea_transcriptomics/oposSOMs/oposSOM_no_stress_no_PH/')
dummy %>% nrow()
dummy <- raw_for_oposSOM_processed %>% filter(MapMan_Category1 != 'stress',MapMan_Category1 != 'PS', !(transcript_id %in% storage_proteins$target_id) ) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
 as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))

rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy %>% head()
env_synthesis <- opossom.new(list(dataset.name = 'no_stress_no_PS', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
opossom.run(env = env_synthesis)

#### No_stress - GO
setwd('~/Pea_transcriptomics/oposSOMs/oposSOM_no_stress')
group_B <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots B.csv',
                      row.names = 1, header = T,  sep = ';')
group_C <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots C.csv',
                      row.names = 1, header = T,  sep = ';')
group_D <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots D.csv',
                      row.names = 1, header = T,  sep = ';')
group_E <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots E.csv',
                      row.names = 1, header = T,  sep = ';')
group_F <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots F.csv',
                      row.names = 1, header = T,  sep = ';')
group_G <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots G.csv',
                      row.names = 1, header = T,  sep = ';')
group_H <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots H.csv',
                      row.names = 1, header = T,  sep = ';')
group_I <- read.table('no_stress - Results/CSV Sheets/Spot Lists/Group Overexpression Spots I.csv',
                      row.names = 1, header = T,  sep = ';')


BP_group_H %>% View()
colnames(group_E)[2] <- 'target_id'
colnames(group_F)[2] <- 'target_id'
colnames(group_G)[2] <- 'target_id'
colnames(group_H)[2] <- 'target_id'
colnames(group_I)[2] <- 'target_id'
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_I, ont = x, cutoff = 0.01))

sapply(ls()[grepl('((BP)|(CC)|(MF))_group', ls())], function(x) write.table(get(x), file.path('GO_tables', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))


#### No_stress_no_PS - GO
setwd('~/Pea_transcriptomics/oposSOMs/oposSOM_no_stress_no_PH/')
group_B <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots B.csv',
                      row.names = 1, header = T,  sep = ';')
group_C <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots C.csv',
                      row.names = 1, header = T,  sep = ';')
group_D <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots D.csv',
                      row.names = 1, header = T,  sep = ';')
group_E <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots E.csv',
                      row.names = 1, header = T,  sep = ';')
group_F <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots F.csv',
                      row.names = 1, header = T,  sep = ';')
group_G <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots G.csv',
                      row.names = 1, header = T,  sep = ';')
group_H <- read.table('no_stress_no_PS - Results/CSV Sheets/Spot Lists/Group Overexpression Spots H.csv',
                      row.names = 1, header = T,  sep = ';')

colnames(group_B)[2] <- 'target_id'
colnames(group_C)[2] <- 'target_id'
colnames(group_D)[2] <- 'target_id'
colnames(group_E)[2] <- 'target_id'
colnames(group_F)[2] <- 'target_id'
colnames(group_G)[2] <- 'target_id'
colnames(group_H)[2] <- 'target_id'
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_H, ont = x, cutoff = 0.01))
sapply(ls()[grepl('((BP)|(CC)|(MF))_group', ls())], function(x) write.table(get(x), file.path('GO_tables', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))


####### VARIANT CALLING
#### Gff filtering
gff <- read.table('~/Pea_transcriptomics/VC_annotation_database/assembly_proteins.gff3', sep = '\t', header = F)
gff %<>% mutate(V9 = as.character(V9))
head(gff)
gff$protein_name <- apply(gff %>% select(V9), 1, function(x) x %>% strsplit(';') %>% unlist() %>% unname() %>% first() %>% gsub('ID=', '', .) %>% 
                            gsub('TRINITY.*~~', '', ., perl = T) %>% gsub('(cds\\.)|(\\.((utr.*)|(exon.*)))', '', .))
pep_list <- read.table('~/Pea_transcriptomics/VC_annotation_database/pep_names_novel.txt', header = F) %>% unlist() %>% unname()
pep_list %>% length()
gff %>% nrow()
gff_new <- gff %>% filter(protein_name %in% pep_list) %>% select(-c(protein_name))
  gff_new %>% nrow()
write.table(gff_new, '~/Pea_transcriptomics/VC_annotation_database/novel_assembly_annotation.gff3', sep = '\t', row.names = F, col.names = F)
rm(list = c(gff, gff_new, pep_list))

gtf_exons <- read.table('~/Pea_transcriptomics/trinity_genes.gtf', header = F, sep = '\t') %>% mutate(V9 = as.character(V9))
gtf_exons$transcript_id <- sapply(gtf_exons$V9, function(x) x %>% strsplit(';') %>% unlist() %>% unname() %>% nth(2) %>% 
                                    gsub(' transcript_id ', '', .))
gtf_exons %<>% group_by(transcript_id)
gtf_exons %>% dplyr::select(transcript_id) %>% head()
gtf_exons %>% head()
write.table(gtf_exons, '~/Pea_transcriptomics/VC_annotation_database/exons.gtf', sep = '\t', row.names = F, col.names = F)
gtf_exons %>% filter(transcript_id == 'TRINITY_10_DN3332_c0_g1_i2')

snp_high_Zhewan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhewan_HIGH_no_err_snp.vcf', sep = '\t', header = F) %>% 
snp_high_Zhongwan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhongwan_HIGH_no_err_snp.vcf', sep = '\t', header = F)
annot_full %>% colnames()
annot_full %>% filter(transcript_id %in% as.character(snp_Zhewan$V1)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/Zhewan_snp_effects.tsv', sep = '\t', row.names = F)
annot_full %>% filter(transcript_id %in% as.character(snp_Zhongwan$V1)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/Zhongwan_snp_effects.tsv', sep = '\t')

snp_moderate_Zhewan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhewan_MODERATE_no_err_snp.vcf', header = F, sep = '\t')
snp_moderate_Zhongwan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhongwan_MODERATE_no_err_snp.vcf', header = F, sep = '\t')

annot_full %>% filter(transcript_id %in% as.character(snp_moderate_Zhewan$V1)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/Zhewan_MODERATE_snp_effects.tsv', sep = '\t', row.names = F)

annot_full %>% filter(transcript_id %in% as.character(snp_moderate_Zhongwan$V1)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/Zhongwan_MODERATE_snp_effects.tsv', sep = '\t', row.names = F)

library(UpSetR)
movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
movies %>% View()
mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")
mutations %>% View()

snp_moderate_Zhewan %<>% mutate(target_id = V1)
snp_moderate_Zhongwan %<>% mutate(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_moderate_Zhongwan, x, cutoff = 0.01))
iterating_function_GO(snp_moderate_Zhewan, 'BP', cutoff = 0.01)
`BP_%>%` %>% View()
`CC_%>%` %>% View()
`MF_%>%` %>% View()
BP_snp_moderate_Zhewan %>% View()
BP_snp_moderate_Zhongwan %>% View()

sapply(ls()[grepl('((BP)|(CC)|(MF))_snp*', ls())], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'VC_annotation_database', 'GO_moderate', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))
upregulated_30_days %>% filter(target_id == 'TRINITY_10_DN244_c0_g1_i7') %>% View()

annot_full %>% filter(transcript_id %in% snp_moderate_Zhewan$V1) %>% dplyr::select(MapMan_terms)  %>% unlist() %>% unique()

snp_high_Sprint <- read.table('~/Pea_transcriptomics/VC_annotation_database/Sprint2_snp_HIGH_no_warn.vcf', sep = '\t', header = F)
snp_moderate_Sprint <- read.table('~/Pea_transcriptomics/VC_annotation_database/Sprint2_snp_MODERATE_no_warn.vcf', sep = '\t', header = F)
annot_full %>% filter(transcript_id %in% as.character(snp_moderate_Sprint$V1)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/Sprint2_MODERATE_snp_effects.tsv', sep = '\t', row.names = F)
annot_full %>% filter(transcript_id %in% as.character(snp_high_Sprint$V1)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/Sprint2_HIGH_snp_effects.tsv', sep = '\t', row.names = F)

snp_moderate_Sprint %<>% mutate(target_id = V1)
snp_high_Sprint %<>% mutate(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_moderate_Sprint, x, cutoff = 0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_high_Sprint, x, cutoff = 0.01))
sapply(ls()[grepl('((BP)|(CC)|(MF))_(.*)_Sprint', ls(), perl = T)], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'VC_annotation_database', 'GO_Sprint', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))


### PSIA
psia <- read.table('~/Pea_transcriptomics/20190906_new_pis_sat.csv', sep = '\t', header = T)

library(AnnotationDbi)

ontologies <- Ontology(GOTERM)
terms <- Term(GOTERM)
filter_GO <- function(z, ont) {
  res <- sapply(z %>% strsplit(', ') %>% unlist(), function(x) ifelse(ontologies[names(ontologies) == x] == ont, x, ''))
  res <- res[res != '']
  return(paste(res, collapse = ', '))
}

filter_GO(psia$GO[1], 'BP')
psia$GO[1]
psia$GO <- annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),]$GOs
psia$Biological_Process <- sapply(psia$GO, function(x) filter_GO(x, 'BP'))
psia$Cellular_Component <- sapply(psia$GO, function(x) filter_GO(x, 'CC'))
psia$Molecular_Function <- sapply(psia$GO, function(x) filter_GO(x, 'MF'))
psia$Biological_Process <- sapply(psia$Biological_Process, function(y)
  paste(sapply(y %>% strsplit(', ') %>% unlist(), function(x) terms[names(terms) == x]) %>% unname(), collapse = '; '))
psia$Cellular_Component <-  sapply(psia$Cellular_Component, function(y)
  paste(sapply(y %>% strsplit(', ') %>% unlist(), function(x) terms[names(terms) == x]) %>% unname(), collapse = '; '))
psia$Molecular_Function <-  sapply(psia$Molecular_Function, function(y)
  paste(sapply(y %>% strsplit(', ') %>% unlist(), function(x) terms[names(terms) == x]) %>% unname(), collapse = '; '))
psia$NCBI_Accession <- sapply(annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),]$plants_BLASTP,
  function(x) x %>% strsplit(.,'\\^') %>% unlist() %>% unname() %>% dplyr::first())
psia$NCBI_Accession[is.na(psia$NCBI_Accession)] <- sapply(annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),][is.na(psia$NCBI_Accession),]$plants_BLASTX,
                                                          function(x) x %>% strsplit(.,'\\^') %>% unlist() %>% unname() %>% dplyr::first())
psia$NCBI_Name <- annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),]$BLASTP_names
psia$NCBI_Name[is.na(psia$NCBI_Name)] <- annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),][is.na(psia$NCBI_Name),]$BLASTX_names
psia$MapMan <- annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),]$MapMan_terms
annot_full %<>% mutate(prot_id = as.character(prot_id), GOs = as.character(GOs))
psia %>% dplyr::select(-c(GO, Biological_Process, Cellular_Component, Molecular_Function)) %>% View()
annot_full[annot_full$prot_id %in% psia$Accession,][order(match(annot_full[annot_full$prot_id %in% psia$Accession,]$prot_id, psia$Accession)),]$MapMan_terms[1]

annot_full[annot_full$prot_id %in% psia$Accession,]
annot_full$MapMan_terms[626]


colnames(annot_full)
psia %>% dplyr::select(-c(GO)) %>% View()
annot_full %>% filter(prot_id == 'TRINITY_10_DN907_c0_g1_i2.p1') %>% dplyr::select(prot_id, GOs, BLASTP_names, MapMan_terms)
psia$Biological_Process[2]
write.table(psia %>% dplyr::select(-c(GO)), '~/Pea_transcriptomics/psia_annotation_2.tsv', sep = '\t', row.names = F)

##### COMMON SNP ANNOTATION
snp_Zhe_Zho <- read.table('~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe-Zho.txt', header = F, sep = '\t')
snp_Zhe_Sprint <- read.table('~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe-Sprint.txt', header = F, sep = '\t')
snp_Zho_Sprint <- read.table('~/Pea_transcriptomics/VC_annotation_database/common_snps/Zho-Sprint.txt', header = F, sep = '\t')
snp_Zhe_Zho %<>% mutate(target_id = V1)
snp_Zhe_Sprint %<>% mutate(target_id = V1)
snp_Zho_Sprint %<>% mutate(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_Zhe_Sprint, x, cutoff = 0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_Zhe_Zho, x, cutoff = 0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_Zho_Sprint, x, cutoff = 0.01))
sapply(ls()[grepl('((BP)|(CC)|(MF))_snp_Zh[eo]_((Zho)|(Sprint))', ls(), perl = T)], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'VC_annotation_database', 'common_snps', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))

annot_full %>% filter(transcript_id %in% as.character(snp_Zhe_Sprint$target_id)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe_Sprint_effects.tsv', sep = '\t', row.names = F)
annot_full %>% filter(transcript_id %in% as.character(snp_Zhe_Zho$target_id)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe_Zho_effects.tsv', sep = '\t', row.names = F)
annot_full %>% filter(transcript_id %in% as.character(snp_Zho_Sprint$target_id)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/common_snps/Zho_Sprint_effects.tsv', sep = '\t', row.names = F)

snp_Zhe_uniq <- read.table('~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe_uniq.txt', header = F, sep = '\t')
snp_Zhe_uniq %<>% mutate(target_id = V1)
snp_Zhe_uniq %<>% filter(!(target_id %in% snp_Zhe_Sprint$target_id))
annot_full %>% filter(transcript_id %in% as.character(snp_Zhe_uniq$target_id)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe_uniq_effects.tsv', sep = '\t', row.names = F)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_Zhe_uniq, x, cutoff = 0.01))

snp_Zho_uniq <- read.table('~/Pea_transcriptomics/VC_annotation_database/common_snps/Zho_uniq.txt', header = F, sep = '\t')
snp_Zho_uniq %<>% mutate(target_id = V1)
snp_Zho_uniq %<>% filter(!(target_id %in% snp_Zho_Sprint$target_id))
annot_full %>% filter(transcript_id %in% as.character(snp_Zho_uniq$target_id)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/common_snps/Zho_uniq_effects.tsv', sep = '\t', row.names = F)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_Zho_uniq, x, cutoff = 0.01))
sapply(ls()[grepl('((BP)|(CC)|(MF))_snp_Zh[eo]_uniq', ls(), perl = T)], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'VC_annotation_database', 'common_snps', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))

snp_common_contigs <- read.table('~/Pea_transcriptomics/VC_annotation_database/common_snps/common_transcript.txt', header = F, sep = '\t')
snp_common_contigs %<>% mutate(target_id = V1)
snp_common_contigs %<>% filter(!(target_id %in% snp_Zho_Sprint$target_id), !(target_id %in% snp_Zhe_Sprint))
annot_full %>% filter(transcript_id %in% as.character(snp_common_contigs$target_id)) %>% 
  dplyr::select(c(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)) %>% write.table(., '~/Pea_transcriptomics/VC_annotation_database/common_snps/Zhe_Zho_common_contigs.tsv', sep = '\t', row.names = F)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_common_contigs, x, cutoff = 0.01))
sapply(ls()[grepl('((BP)|(CC)|(MF))_snp_common', ls(), perl = T)], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'VC_annotation_database', 'common_snps', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))

snp_Zhewan_all <- rbind(snp_moderate_Zhewan, snp_high_Zhewan %>% mutate(target_id = V1))
snp_Zhongwan_all <- rbind(snp_moderate_Zhongwan, snp_high_Zhongwan %>% mutate(target_id = V1))

to_remove <- MapMan_df %>% filter(MapMan_Category1 %in% c('PS', 'stress'))
nrow(to_remove)
snp_moderate_Sprint_2 <- snp_moderate_Sprint %>% filter(!(target_id %in% to_remove$transcript_id))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(snp_moderate_Sprint_2, x, cutoff = 0.01))
BP_snp_moderate_Sprint_2 %>% View()
unique(BP_snp_moderate_Sprint$Term, BP_snp_moderate_Sprint_2$Term)
sapply(ls()[grepl('((BP)|(CC)|(MF))_snp_moderate_Sprint_2', ls(), perl = T)], function(x) write.table(get(x), file.path('~', 'Pea_transcriptomics', 'VC_annotation_database', paste(x, 'tsv', sep = '.')), sep = '\t', col.names = T, dec = '.', row.names = F))

############# MapMan GSEA
library(lmPerm)
library(gtools)Оставить
library(coin)

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
GSEA_df <- MapMan_df %>% dplyr::select(transcript_id, normZheEarly, normZheLate, normZhoEarly, normZhoLate, MapMan_Category, MapMan_Category1, MapMan_Category2, MapMan_Category3, MapMan_Category4) %>% mutate_at(vars(-c(transcript_id)), all.vars(as.numeric()))
###### Заменить NaN 
GSEA_df[,c(2:5)] <- sapply(GSEA_df[,c(2:5)], function(x) replace(x, is.nan(x), 0)) 

GSEA_df %<>% mutate(sumEarly = normZheEarly^2 + normZhoEarly^2, sumLate = normZheLate^2 + normZhoLate^2)
GSEA_df %<>% mutate(TotalDelta = sumEarly + sumLate)
GSEA_df <- GSEA_df[order(GSEA_df$TotalDelta, decreasing = T),]

GSEA_df %>% head()
summary(lmp(TotalDelta ~ MapMan_Category1, data=GSEA_df))


snp_Sprint_all <- rbind(snp_moderate_Sprint, snp_high_Sprint)
##### SNP Euler Mapping 

nrow(snp_Zhe_uniq)
nrow(snp_moderate_Zhewan %>% filter(!(target_id %in% snp_moderate_Sprint$target_id), !(target_id %in% snp_moderate_Zhongwan$target_id),
                                    ))
nrow(snp_high_Zhewan)
install.packages('venneuler')
all_snps <- data.frame('transcript' = as.character(annot_full$transcript_id)) %>% mutate(
  'Zhewan1_moderate' = transcript %in% as.character(snp_moderate_Zhewan$target_id), 'Zhongwan6_moderate' = transcript %in% as.character(snp_moderate_Zhongwan$target_id),
  'Sprint2_moderate' = transcript %in% as.character(snp_moderate_Sprint$target_id), 'Zhewan1_high' = transcript %in% as.character(snp_high_Zhewan$V1), 
  'Zhongwan6_high' = transcript %in% as.character(snp_high_Zhongwan$V1), 'Sprint2_high' = transcript %in% as.character(snp_high_Sprint$target_id),
  'Common transcripts' = transcript %in% as.character(snp_common_contigs$target_id),
   'Zhewan1' = transcript %in% as.character(snp_Zhewan_all$V1),
  'Zhongwan6' = transcript %in% as.character(snp_Zhongwan_all$V1),
  'Sprint2' = transcript %in% as.character(snp_Sprint_all$V1)
) %>% mutate_all(as.character())

grid.newpage()
g <- draw.triple.venn(area1 = length(subset(all_snps, Zhewan1_moderate == TRUE)$transcript),
                 area2 = length(subset(all_snps, Zhongwan6_moderate == TRUE)$transcript),
                 area3 = length(subset(all_snps, Sprint2_moderate == TRUE)$transcript),
                 n12 = nrow(subset(all_snps, Zhewan1_moderate == TRUE & Zhongwan6_moderate == TRUE)),
                 n13 = nrow(subset(all_snps, Zhewan1_moderate == TRUE & Sprint2_moderate == TRUE)),
                 n23 = nrow(subset(all_snps, Zhongwan6_moderate == TRUE & Sprint2_moderate == TRUE)),
                 n123 = nrow(subset(all_snps, Zhongwan6_moderate == TRUE & Sprint2_moderate == TRUE)),
                 category = c('Zhewan-1', 'Zhongwan-6', 'Sprint-2'), lty = 'blank', fill = wesanderson::wes_palettes$Moonrise3[1:3],
                 euler.d=TRUE, scaled = T, cex = 2, cat.cex = 2)
require(gridExtra)
grid.arrange(gTree(children=g), top=textGrob("Distribution of moderate effect SNPs in three pea lines", gp=gpar(fontsize=25, x = 0)))

grid.newpage()
g <- draw.triple.venn(area1 = length(subset(all_snps, Zhewan1_high == TRUE)$transcript),
                      area2 = length(subset(all_snps, Zhongwan6_high == TRUE)$transcript),
                      area3 = length(subset(all_snps, Sprint2_high == TRUE)$transcript),
                      n12 = nrow(subset(all_snps, Zhewan1_high == TRUE & Zhongwan6_high == TRUE)),
                      n13 = nrow(subset(all_snps, Zhewan1_high == TRUE & Sprint2_high == TRUE)),
                      n23 = nrow(subset(all_snps, Zhongwan6_high == TRUE & Sprint2_high == TRUE)),
                      n123 = nrow(subset(all_snps, Zhongwan6_high == TRUE & Sprint2_high == TRUE)),
                      category = c('Zhewan-1', 'Zhongwan-6', 'Sprint-2'), lty = 'blank', fill = wesanderson::wes_palettes$Moonrise3[1:3],
                      euler.d=TRUE, scaled = T, cex = 2, cat.cex = 2)
require(gridExtra)
grid.arrange(gTree(children=g), top=textGrob("Distribution of high effect SNPs in three pea lines", gp=gpar(fontsize=25, x = 0)))

grid.newpage()
g <- draw.triple.venn(area1 = length(subset(all_snps, Zhewan1 == TRUE)$transcript),
                      area2 = length(subset(all_snps, Zhongwan6 == TRUE)$transcript),
                      area3 = length(subset(all_snps, Sprint2 == TRUE)$transcript),
                      n12 = nrow(subset(all_snps, Zhewan1 == TRUE & Zhongwan6 == TRUE)),
                      n13 = nrow(subset(all_snps, Zhewan1 == TRUE & Sprint2 == TRUE)),
                      n23 = nrow(subset(all_snps, Zhongwan6 == TRUE & Sprint2 == TRUE)),
                      n123 = nrow(subset(all_snps, Zhongwan6 == TRUE & Sprint2 == TRUE)),
                      category = c('Zhewan-1', 'Zhongwan-6', 'Sprint-2'), lty = 'blank', fill = wesanderson::wes_palettes$Moonrise3[1:3],
                      euler.d=TRUE, scaled = T, cex = 2, cat.cex = 2)
require(gridExtra)
grid.arrange(gTree(children=g), top=textGrob("Distribution of all SNPs in three pea lines", gp=gpar(fontsize=25, x = 0)))


library(plotrix)
library(wesanderson)

wes_palettes


(annot_full %>% filter(transcript_id %in% snp_moderate_Sprint$V1))[grepl('GO:0033095', (annot_full %>% filter(transcript_id %in% snp_moderate_Sprint$V1))$GOs),] %>% 
  dplyr::select(transcript_id, BLASTP_names, BLASTX_names) %>% write.table(., file.path('~/Pea_transcriptomics/VC_annotation_database/aleurone_Sprint2_Zhongwan.tsv'), sep = '\t')

BP_downregulated_30_days %>% nrow()
BP_Sprint_DEGs <- rbind(BP_upregulated_30_days[c(1:10),], BP_downregulated_30_days[c(1:10),]) %>% mutate(condition = c(rep('Upregulated', 10), rep('Downregulated', 10)))
BP_Sprint_DEGs$Term[1] <- 'organic cyclic compound biosynthetic process'
BP_Sprint_DEGs$Term[6] <- 'regulation of timing of transition from vegetative to reproductive phase'
BP_Sprint_DEGs$Term[7] <- 'positive regulation of transcription, DNA-templated'
BP_Sprint_DEGs$Term[10] <- 'abscisic acid-activated signaling pathway'
BP_Sprint_DEGs$Term[16] <- 'proteolysis involved in cellular protein catabolic process'
BP_Sprint_DEGs$Term[18] <- 'purine ribonucleotide biosynthetic process'
BP_Sprint_DEGs$Term[20] <- 'regulation of phenylpropanoid metabolic process'

CC_Sprint_DEGs <- rbind(CC_upregulated_30_days, CC_downregulated_30_days[1:10,]) %>% mutate(condition = c(rep('Upregulated', 8), rep('Downregulated', 10)))
CC_Sprint_DEGs %>% View()



sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(downregulated_30_days, x, cutoff = 0.001))
sapply(ls()[grepl('((BP)|(CC)|(MF))_.*regulated', ls())], function(x)
  x %>% get() %>% write.table(., file.path('~/Pea_transcriptomics/Sprint2_GO_final/', paste(x, 'tsv', sep = '\t')), col.names = T, row.names = F,
                              sep = '\t'))
MF_upregulated_30_days %>% View()

MF_Sprint_DEGs <- rbind(MF_upregulated_30_days, MF_downregulated_30_days) %>% mutate(condition = c(rep('Upregulated', 8), rep('Downregulated', 11)))
MF_Sprint_DEGs$Term[3] <- 'transferase activity, transferring hexosyl groups'
MF_Sprint_DEGs$Term[4] <- 'transmembrane receptor protein serine/threonine kinase binding'
MF_Sprint_DEGs$Term[6] <- 'quercetin 3-O-glucosyltransferase activity'
MF_Sprint_DEGs$Term[7] <- 'quercetin 7-O-glucosyltransferase activity'
MF_Sprint_DEGs$Term[13] <- 'hydrolase activity, hydrolyzing O-glycosyl compounds'
MF_Sprint_DEGs$Term[17] <- 'oxidoreductase activity, CH-OH group of donors, NAD/NADP as acceptor'

MF_Sprint_DEGs %>% View()


bubble <- ggplot(data = BP_Sprint_DEGs, aes(x=condition, y=Term, color=Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes, size = weight01Fisher)) +
  geom_point() + scale_size_continuous(name = 'qvalue', trans='exp',
  range = c(6, 16)) + scale_color_gradientn(name='% of annotated genes', colours = c('#52B365', '#ffff00')) + theme_light() + theme_linedraw() +
  theme(text = element_text(size=20, family = 'arial'), axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank())
colorRampPalette(rev(brewer.pal(9, "Greens")) )(255)[30:230]
c('#00441B', '#ffff00')
?colorRampPalette
?brewer.pal
brewer.pal.info

###### Genome assembly backup
unmapped <- read.table('~/Pea_transcriptomics/unmapped_contigs.txt', header = F)
to_remove <- annot_full %>% filter(transcript_id %in% unmapped$V1) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, Pfam)
annot_full %>% filter(transcript_id %in% to_remove$transcript_id, is.na(BLASTP_names), is.na(BLASTX_names)) %>% View()


ls()[grepl('(TF)|(tf)', ls())]


####### Transcription Factors
snp_Sprint_all %<>% mutate(MapMan = annot_full[annot_full$transcript_id %in% snp_Sprint_all$target_id,][match(snp_Sprint_all$target_id,
                  annot_full[annot_full$transcript_id %in% snp_Sprint_all$target_id,]$transcript_id),]$MapMan_terms)
snp_Sprint_all %>% dplyr::select(V1, target_id, MapMan) %>% tail()

tf_template <- '(GeBP like)|(\\.ARR)|(CCAAT box binding factor family)|(C2C2\\(Zn\\) YABBY family)|(C2C2\\(Zn\\) DOF zinc finger family)|(ovate family OFP)|
                           (BBR/BPC)|(plant TF \\(pbf2\\))|(BSD domain containing family)|(AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family)|
                           (zf-HD)|(sigma like plant)|(ARF, Auxin Response Factor family)|(ARF, Auxin Response Factor family)|(C2C2\\(Zn\\) CO-like, Constans-like zinc finger family)|
                           (CCAAT box binding factor family, HAP2)|(NIN-like bZIP-related family)|(GRF zinc finger family)|(Alfin-like)|(ELF3)|(transcription factor)'


Sprint2_TF_snp <- snp_Sprint_all[grepl(tf_template, snp_Sprint_all$MapMan),] %>% dplyr::select(target_id) %>% unlist() %>% unname() %>% unique()
Sprint2_TF_snp <- annot_full %>% filter(transcript_id %in% Sprint2_TF_snp) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)
Sprint2_TF_snp %<>% mutate(SNP_annotation = snp_Sprint_all[snp_Sprint_all$target_id %in% Sprint2_TF_snp$transcript_id,][
  match(Sprint2_TF_snp$transcript_id, snp_Sprint_all[snp_Sprint_all$target_id %in% Sprint2_TF_snp$transcript_id,]$target_id),]$V8)

snp_Zhewan_all %<>% mutate(MapMan = annot_full[annot_full$transcript_id %in% snp_Zhewan_all$target_id,][match(snp_Zhewan_all$target_id,
                                      annot_full[annot_full$transcript_id %in% snp_Zhewan_all$target_id,]$transcript_id),]$MapMan_terms)
Zhewan_TF_snp <- snp_Zhewan_all[grepl(tf_template, snp_Zhewan_all$MapMan),] %>% dplyr::select(target_id) %>% unlist() %>% unname() %>% unique()
Zhewan_TF_snp <- annot_full %>% filter(transcript_id %in% Zhewan_TF_snp) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)
Zhewan_TF_snp %<>% mutate(SNP_annotation = snp_Zhewan_all[snp_Zhewan_all$target_id %in% Zhewan_TF_snp$transcript_id,][
  match(Zhewan_TF_snp$transcript_id, snp_Zhewan_all[snp_Zhewan_all$target_id %in% Zhewan_TF_snp$transcript_id,]$target_id),]$V8)

snp_Zhongwan_all %<>% mutate(MapMan = annot_full[annot_full$transcript_id %in% snp_Zhongwan_all$target_id,][match(snp_Zhongwan_all$target_id,
                                                                                                              annot_full[annot_full$transcript_id %in% snp_Zhongwan_all$target_id,]$transcript_id),]$MapMan_terms)
Zhongwan_TF_snp <- snp_Zhongwan_all[grepl(tf_template, snp_Zhongwan_all$MapMan),] %>% dplyr::select(target_id) %>% unlist() %>% unname() %>% unique()
Zhongwan_TF_snp <- annot_full %>% filter(transcript_id %in% Zhongwan_TF_snp) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, GOs, MapMan_terms)
Zhongwan_TF_snp %<>% mutate(SNP_annotation = snp_Zhongwan_all[snp_Zhongwan_all$target_id %in% Zhongwan_TF_snp$transcript_id,][
  match(Zhongwan_TF_snp$transcript_id, snp_Zhongwan_all[snp_Zhongwan_all$target_id %in% Zhongwan_TF_snp$transcript_id,]$target_id),]$V8)

snp_Zhongwan_all %>% head(2)
Sprint2_TF_snp %>% View()
Zhewan_TF_snp[grepl('abscisic', Zhewan_TF_snp$MapMan_terms),]

snp_moderate_Zhewan %>% filter(V1 == 'TRINITY_10_DN3455_c1_g1_i1')
Zhewan_TF_snp %>% filter(transcript_id %in% Zhongwan_TF_snp$transcript_id) %>% nrow()

sapply(ls()[grepl('_TF_snp', ls())], function(x) get(x) %>% write.table(., file.path('~', 'Pea_transcriptomics', paste(x, 'tsv', sep = '.')),
                                                                        sep = '\t', col.names = T, row.names = F))

annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% upregulated_30_days$target_id)

unmapped <- read.table('~/Pea_transcriptomics/unmapped_contigs_splice.txt', header = F, sep = '\t')
unmapped <- annot_full %>% filter(transcript_id %in% unmapped$V1) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names) %>% filter(is.na(BLASTP_names) &
                                                                                                                                is.na(BLASTX_names))
unmapped %>% View()
snp_Zhewan_all$target_id %>% unique() %>% length()

###### Cosine similarity - new attempt
vector_populator <- function(x, up, down) {
  if (annot_full$transcript_id[x] %in% up$target_id) return(1)
  else if (annot_full$transcript_id[x] %in% down$target_id) return(-1)
  else return(0)
}

## V1: Sprint2-10 -> Sprint2-20
V1 <- rep(0, nrow(annot_full))
V1 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, upregulated_30_days, downregulated_30_days))
names(V1) <- annot_full$transcript_id

## V2: Zhewan1-10 -> Zhewan1-25
V2 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, zhe_up, zhe_down))
names(V2) <- annot_full$transcript_id

## V3: Zhongwan6-10 -> Zhongwan6-25
V3 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, zhong_up, zhong_down))
names(V3) <- annot_full$transcript_id

## V4: Sprint2-10 -> Zhewan1-10 (note the inverted order of up- and downregulation subsets) 
V4 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, cross_1_down, cross_1_up))
names(V4) <- annot_full$transcript_id

## V5: Sprint2-10 -> Zhongwan6-10 (note the inverted order of up- and downregulation subsets) 
V5 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, cross_2_down, cross_2_up))
names(V5) <- annot_full$transcript_id

## V6: Sprint2-20 -> Zhewan1-25 (note the inverted order of up- and downregulation subsets)
V6 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, cross_4_down, cross_4_up))
names(V6) <- annot_full$transcript_id

## V7: Sprint2-20 -> Zhongwan6-25 (note the inverted order of up- and downregulation subsets)
V7 <- sapply(c(1:nrow(annot_full)), function(x) vector_populator(x, cross_5_down, cross_5_up))
names(V7) <- annot_full$transcript_id

cosine(V1, V4)
cosine(V1, V5)
cosine(V4, V5)
cosine(V1, V6)
cosine(V1, V7)



#### heterologous SNP
hetero_Zhewan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhewan1_corr_het.ann.vcf', header = F, sep = '\t')
hetero_Zhongwan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhongwan6_corr_het.ann_entries.vcf', header = F, sep = '\t')
hetero_Zhewan %<>% mutate(target_id = V1)
hetero_Zhongwan %<>% mutate(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(hetero_Zhewan, ont = x, cutoff = 0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(hetero_Zhongwan, ont = x, cutoff = 0.01))


install.packages("UpSetR")
sapply(ls()[grepl('((BP)|(CC)|(MF))_hetero_', ls())], function(x) get(x) %>% write.table(., file.path('~', 'Pea_transcriptomics', 'VC_annotation_database',
                              'hetero_exclusive', paste(x, 'tsv', sep = '.')),
                                                                        sep = '\t', col.names = T, row.names = F))


pea_scatter_plotter <- function(MapMan_df_loc = MapMan_df, category1 = NA, category2 = NA, category3 = NA, category4 = NA, uprmax = NA,
                                plot_title = NA){
  if (!is.na(category4)) {
    x <- quo(MapMan_Category4 == category4)
    plot_title <- paste(category1, category4, sep = ':')
  } else if(!is.na(category3)){
    x <- quo(MapMan_Category3 == category3)
    plot_title <- paste(category1, category3, sep = ':')
  } else if(!is.na(category2)) {
    x <- quo(MapMan_Category2 == category2)
    plot_title <- paste(category1, category2, sep = ':')
  } else if(!is.na(category1)) {
    x <- quo(MapMan_Category1 == category1)
    plot_title <- category1
  } else {x <- 'None'}

  if (x != 'None') {MapMan_df_loc %<>% filter(!!x)}
  
  if (!is.na(uprmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.<=uprmax))
  }
  
  lm_Zhe <- lm(Zhewan1_25_mean ~ Zhewan1_10_mean, MapMan_df_loc)
  Zhewan_line <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(fit) %>% unlist()
  Zhewan_lower <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(lwr) %>% unlist()
  Zhewan_upper <- lm_Zhe %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(upr) %>% unlist()
  
  lm_Zho <- lm(Zhongwan6_25_mean ~ Zhongwan6_10_mean, MapMan_df_loc)
  Zhongwan_line <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(fit) %>% unlist()
  Zhongwan_lower <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(lwr) %>% unlist()
  Zhongwan_upper <- lm_Zho %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(upr) %>% unlist()
  
  lm_Spr <- lm(Sprint2_30_mean ~ Sprint2_10_mean, MapMan_df_loc)
  Sprint2_line <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(fit) %>% unlist()
  Sprint2_lower <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(lwr) %>% unlist()
  Sprint2_upper <- lm_Spr %>% predict(interval='confidence') %>% as.data.frame() %>% dplyr::select(upr) %>% unlist()
  
  ggplot(MapMan_df_loc, aes(x = value, y = value, alpha  = 0.1, size = I(3))) + 
    geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan6_25_mean, shape = 'Zhongwan6'), alpha = 0.4, col = '#dfc9e2') +
    geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean, shape = 'Sprint2'), alpha = 0.4, col = '#228AD0') +
    geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean, shape = 'Zhewan1'), alpha = 0.4, col = '#00b6c2') +
    xlab(label = 'Earlier condition') + ylab('Latter condition')+
    guides(alpha = 'none')+
    geom_line(aes(Sprint2_10_mean, Sprint2_line), size = 0.8, col = '#228AD0', inherit.aes = F) + geom_ribbon(aes(
    Sprint2_10_mean, ymin = Sprint2_lower, ymax = Sprint2_upper), fill = '#61CDFF', alpha = 0.3, inherit.aes = F) +
    geom_line(aes(Zhewan1_10_mean, Zhewan_line), size = 0.8, col = '#00b6c2', inherit.aes = F) + geom_ribbon(aes(
        Zhewan1_10_mean, ymin = Zhewan_lower, ymax = Zhewan_upper), fill = '#ccffff', alpha = 0.3, inherit.aes = F)+
    geom_line(aes(Zhongwan6_10_mean, Zhongwan_line), size = 0.8, col = '#dfc9e2', inherit.aes = F) + geom_ribbon(aes(
    Zhongwan6_10_mean, ymin = Zhongwan_lower, ymax = Zhongwan_upper), fill = '#e3e8ff', alpha = 0.3, inherit.aes = F) +
    labs(title = paste('Distribution of ', plot_title, '-related transcripts', sep = ''), shape = 'Line') + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 13),
                       legend.text = element_text(size = 11))
  
}

pea_scatter_plotter(MapMan_df %>% filter(transcript_id %in% annot_full[grepl(tf_template, annot_full$MapMan_terms),]$transcript_id), plot_title = 'TF')

MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism') %>% dplyr::select(MapMan_Category3) %>% unlist() %>% unname() %>% unique()
MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism', MapMan_Category2 == 'abscisic acid') %>% nrow() 

tf_mutants <- c('TRINITY_10_DN11300_c0_g1_i1', 'TRINITY_10_DN11300_c0_g1_i3', 'TRINITY_10_DN13275_c0_g1_i2',
                'TRINITY_10_DN3445_c0_g7_i1', 'TRINITY_10_DN3455_c1_g1_i1', 'TRINITY_10_DN1211_c0_g1_i2',
                'TRINITY_10_DN434_c0_g1_i1', 'TRINITY_10_DN27017_c0_g1_i1', 'TRINITY_10_DN11886_c0_g2_i1',
                'TRINITY_10_DN11860_c0_g1_i1', 'TRINITY_10_DN795_c0_g2_i1')

mutants_Zhewan <- rbind(zhe_up %>% filter(target_id %in% tf_mutants) %>% dplyr::select(-c(Pfam, gene_ontology_pfam, seed_eggNOG_ortholog, GOs)),
          zhe_down %>% filter(target_id %in% tf_mutants) %>% dplyr::select(-c(Pfam, gene_ontology_pfam, seed_eggNOG_ortholog, GOs))) %>% mutate(
            Condition = c(rep('up', 3), rep('down', 5))
          )

mutants_Zhongwan <- rbind(zhong_up %>% filter(target_id %in% tf_mutants) %>% dplyr::select(-c(Pfam, gene_ontology_pfam, seed_eggNOG_ortholog, GOs)),
                        zhong_down %>% filter(target_id %in% tf_mutants) %>% dplyr::select(-c(Pfam, gene_ontology_pfam, seed_eggNOG_ortholog, GOs))) %>% mutate(
                          Condition = c(rep('up', 2), rep('down', 2))
                        )

sapply(ls()[grepl('mutants_', ls())], function(x) get(x) %>% write.table(., file.path('~', 'Pea_transcriptomics',
                                      'mutant_tf_expression', paste(x, '.tsv', sep='')), col.names = T, row.names = F, sep = '\t'))

homo_Zhewan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhewan1_homo.corrected.entries.vcf', sep = '\t', header = F)
homo_Zhongwan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhongwan6_homo.corrected.entries.vcf', sep = '\t', header = F)
homo_Zhewan_uniq <- sqldf('SELECT * FROM homo_Zhewan[,c(1:5)] EXCEPT SELECT * FROM homo_Zhongwan[,c(1:5)]')
homo_Zhongwan_uniq <- sqldf('SELECT * FROM homo_Zhongwan[,c(1:5)] EXCEPT SELECT * FROM homo_Zhewan[,c(1:5)]')
homo_Zhongwan_uniq %<>% mutate(target_id = V1)
homo_Zhewan_uniq %<>% mutate(target_id = V1)

sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(homo_Zhewan_uniq, x, 0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(homo_Zhongwan_uniq, x, 0.01))
sapply(ls()[grepl('((BP)|(CC)|(MF))_homo_((Zhewan)|Zhongwan)_uniq', ls())], function(x) get(x) %>% write.table(., file.path('~', 'Pea_transcriptomics',
    'VC_annotation_database', 'unique_snps_GO', paste(x, '.tsv', sep='')), col.names = T, row.names = F, sep = '\t'))
library(sqldf)

BP_homo_Zhewan_uniq %>% filter(Term %in% c('seed maturation', 'starch catabolic process')) %>% View()
BP_homo_Zhongwan_uniq %>% filter(Term %in% c('seed maturation', 'starch catabolic process')) %>% View()

seed_development_snp <- annot_full[grepl('(GO:0010431)|(GO:0010162)|(GO:1990068)|(GO:2000034)', annot_full$GOs),] %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names,
      MapMan_terms) %>% mutate('Zhewan-1' = as.numeric(transcript_id %in% homo_Zhewan_uniq$V1),
                               'Zhongwan-6' = as.numeric(transcript_id %in% homo_Zhongwan_uniq$V1))
starch_catabolysm_snp <- annot_full[grepl('GO:0005983', annot_full$GOs),] %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names,
                                                                                                                                    MapMan_terms) %>% mutate('Zhewan-1' = as.numeric(transcript_id %in% homo_Zhewan_uniq$V1),
                                                                                                                                                             'Zhongwan-6' = as.numeric(transcript_id %in% homo_Zhongwan_uniq$V1))
starch_catabolysm_snp %>% View()
sapply(ls()[grepl('[^TF]_snp$', ls(), perl = T)], function(x) x %>% get() %>% 
         write.table(., file.path('~','Pea_transcriptomics', 'VC_annotation_database', 'mutant_snp_geneset', paste(x, 'tsv', sep = '.')), 
                     sep = '\t', col.names = T, row.names = F))

#### TF Overexpression Spots
group_A <- read.table('~/Pea_transcriptomics/TF_overexpression_spots290919/Group Overexpression Spots A.csv', sep = ';', header = T)
group_B <- read.table('~/Pea_transcriptomics/TF_overexpression_spots290919/Group Overexpression Spots B.csv', sep = ';', header = T)
group_C <- read.table('~/Pea_transcriptomics/TF_overexpression_spots290919/Group Overexpression Spots C.csv', sep = ';', header = T)
group_D <- read.table('~/Pea_transcriptomics/TF_overexpression_spots290919/Group Overexpression Spots D.csv', sep = ';', header = T)
group_E <- read.table('~/Pea_transcriptomics/TF_overexpression_spots290919/Group Overexpression Spots E.csv', sep = ';', header = T)
group_A <- annot_full %>% filter(transcript_id %in% group_A$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(target_id = transcript_id)
group_B <- annot_full %>% filter(transcript_id %in% group_B$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(target_id = transcript_id)
group_C <- annot_full %>% filter(transcript_id %in% group_C$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(target_id = transcript_id)
group_D <- annot_full %>% filter(transcript_id %in% group_D$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(target_id = transcript_id)
group_E <- annot_full %>% filter(transcript_id %in% group_E$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(target_id = transcript_id)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_B, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_C, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_D, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_E, ont=x, cutoff=0.01))
MF_group_A %>% View()
BP_group_C %>% View()
BP_group_B %>% View()
sapply(ls()[grepl('_group_(A|B|C|D|E)', ls())], function(x) x %>% get() %>%
         write.table(., file.path('~/Pea_transcriptomics/TF_overexpression_spots290919/reports/', paste(x, 'tsv', sep = '\t')),
                     sep = '\t', col.names = T, row.names = F))

#### Transcription-related Genes Overexpression Spots
group_A <- read.table('~/Pea_transcriptomics/transcription_all_overexpression_spots_290919/Group Overexpression Spots A.csv', sep = ';', header = T)
group_B <- read.table('~/Pea_transcriptomics/transcription_all_overexpression_spots_290919/Group Overexpression Spots B.csv', sep = ';', header = T)
group_C <- read.table('~/Pea_transcriptomics/transcriptio n_all_overexpression_spots_290919/Group Overexpression Spots C.csv', sep = ';', header = T)
group_D <- read.table('~/Pea_transcriptomics/transcription_all_overexpression_spots_290919/Group Overexpression Spots D.csv', sep = ';', header = T)
group_E <- read.table('~/Pea_transcriptomics/transcription_all_overexpression_spots_290919/Group Overexpression Spots E.csv', sep = ';', header = T)

group_A <- annot_full %>% filter(transcript_id %in% group_A$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_B <- annot_full %>% filter(transcript_id %in% group_B$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_C <- annot_full %>% filter(transcript_id %in% group_C$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_D <- annot_full %>% filter(transcript_id %in% group_D$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_E <- annot_full %>% filter(transcript_id %in% group_E$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)

sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_B, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_C, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_D, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_E, ont=x, cutoff=0.01))

sapply(ls()[grepl('_group_(A|B|C|D|E)', ls())], function(x) x %>% get() %>%
         write.table(., file.path('~/Pea_transcriptomics/transcription_all_overexpression_spots_290919/reports/', paste(x, 'tsv', sep = '.')),
                     sep = '\t', col.names = T, row.names = F))

#### Hormones Overexpression Spots
group_A <- read.table('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/Group Overexpression Spots A.csv', sep = '\t', header = T)
group_B <- read.table('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/Group Overexpression Spots B.csv', sep = '\t', header = T)
group_C <- read.table('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/Group Overexpression Spots C.csv', sep = '\t', header = T)
group_D <- read.table('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/Group Overexpression Spots D.csv', sep = '\t', header = T)

group_A %>% View()
group_A <- annot_full %>% filter(transcript_id %in% group_A$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_B <- annot_full %>% filter(transcript_id %in% group_B$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_C <- annot_full %>% filter(transcript_id %in% group_C$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)
group_D <- annot_full %>% filter(transcript_id %in% group_D$X) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% dplyr::rename(target_id = transcript_id)

sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_B, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_C, ont=x, cutoff=0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_D, ont=x, cutoff=0.01))


sapply(ls()[grepl('^group_(A|B|C|D)$', ls())], function(x) x %>% get() %>%
         write.table(., file.path('~/Pea_transcriptomics/hormones_oxerexpression_results_290919', paste(x, 'tsv', sep = '.')),
                     sep = '\t', col.names = T, row.names = F))

sapply(ls()[grepl('_group_(A|B|C|D)', ls())], function(x) x %>% get() %>%
         write.table(., file.path('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/reports/', paste(x, 'tsv', sep = '.')),
                     sep = '\t', col.names = T, row.names = F))

Sprint_tf_upregulated <- annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% upregulated_30_days$target_id) %>%
  dplyr::select(transcript_id, BLASTP_names, MapMan_terms) %>% mutate(MapMan_terms = as.character(MapMan_terms))
Sprint_tf_upregulated %<>% mutate(TF_class = sapply(Sprint_tf_upregulated$MapMan_terms, function(x) x %>% strsplit(.,'\\.') %>% unlist() %>% unname() %>% nth(3)))

Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN12709_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN1571_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN3872_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN41153_c0_g1_i1',]$TF_class <- 'ARR'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN3012_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN12091_c0_g1_i1',]$TF_class <- 'PHD finger transcription factor'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN3455_c1_g1_i1',]$TF_class <- 'ABI3/VP1-related B3-domain-containing transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN3455_c1_g1_i2',]$TF_class <- 'ABI3/VP1-related B3-domain-containing transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN1381_c0_g1_i2',]$TF_class <- 'MYB-related transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN7424_c1_g1_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN6525_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN17077_c1_g1_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN17077_c0_g1_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN1036_c1_g2_i2',]$TF_class <- 'PHD finger transcription factor'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN2730_c0_g1_i2',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN34830_c0_g1_i1',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN3035_c0_g1_i4',]$TF_class <- 'ABI3/VP1-related B3-domain-containing transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN3301_c0_g1_i6',]$TF_class <- 'GRF zinc finger family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN7257_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN2200_c0_g3_i1',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN5964_c0_g2_i4',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN15229_c0_g1_i2',]$TF_class <- 'HB,Homeobox transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN3507_c0_g1_i3',]$TF_class <- 'HB,Homeobox transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN6426_c0_g1_i1',]$TF_class <- 'NIN-like bZIP-related family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN4913_c0_g1_i1',]$TF_class <- 'GRAS transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN7359_c0_g4_i2',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN7359_c0_g2_i1',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN418_c0_g1_i1',]$TF_class <- 'TCP transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN418_c0_g2_i1',]$TF_class <- 'TCP transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN6579_c0_g1_i1',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN31308_c0_g1_i1',]$TF_class <- 'GRAS transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN10239_c0_g2_i1',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_10_DN2300_c0_g3_i1',]$TF_class <- 'ovate family OFP'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN27876_c0_g1_i1',]$TF_class <- 'ovate family OFP'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN7264_c0_g2_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_upregulated[Sprint_tf_upregulated$transcript_id == 'TRINITY_30_DN2565_c0_g1_i1',]$TF_class <- 'bZIP transcription factor family'

Sprint_tf_upregulated <- rbind(Sprint_tf_upregulated, data.frame(transcript_id = 'TRINITY_10_DN1381_c0_g1_i2',
                                                                     BLASTP_names = 'trihelix transcription factor GT-4 isoform X1 [Medicago', MapMan_terms = 'RNA.regulation of transcription.MYB-related transcription factor family,RNA.regulation of transcription.Trihelix, Triple-Helix transcription factor family',
                                                                     TF_class = 'Trihelix, Triple-Helix transcription factor family'))
Sprint_tf_upregulated <- rbind(Sprint_tf_upregulated, data.frame(transcript_id = 'TRINITY_30_DN34830_c0_g1_i1',
                                                                 BLASTP_names = 'probable transcription factor KAN4 [Medicago truncatula]', MapMan_terms = 'RNA.regulation of transcription.G2-like transcription factor family, GARP,RNA.regulation of transcription.MYB domain transcription factor family',
                                                                 TF_class = 'G2-like transcription factor family, GARP'))
Sprint_tf_upregulated <- rbind(Sprint_tf_upregulated, data.frame(transcript_id = 'TRINITY_30_DN2200_c0_g3_i1',
                                                                 BLASTP_names = 'myb-related protein 2 isoform X2 [Medicago truncatula]', MapMan_terms = 'RNA.regulation of transcription.G2-like transcription factor family, GARP,RNA.regulation of transcription.MYB domain transcription factor family',
                                                                 TF_class = 'G2-like transcription factor family, GARP'))

Sprint_tf_upregulated[grepl('transport', Sprint_tf_upregulated$TF_class),] %>% View()
Sprint_tf_upregulated[grepl('MYB', Sprint_tf_upregulated$TF_class),] %>% nrow()
Sprint_tf_upregulated[grepl('MYB', Sprint_tf_upregulated$MapMan_terms),][Sprint_tf_upregulated$TF_class != 'MYB domain transcription factor family',] %>% View()

Sprint_tf_upregulated_pivot <-  table(unlist(Sprint_tf_upregulated$TF_class)) %>% as.data.frame() %>% arrange(desc(Freq))




Sprint_tf_downregulated <- annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% downregulated_30_days$target_id) %>%
  dplyr::select(transcript_id, BLASTP_names, MapMan_terms) %>% mutate(MapMan_terms = as.character(MapMan_terms))

Sprint_tf_downregulated %<>% mutate(TF_class = sapply(Sprint_tf_downregulated$MapMan_terms, function(x) x %>% strsplit(., '\\.') %>% unlist() %>% unname() %>% nth(3)))
Sprint_tf_downregulated[grepl('heat', Sprint_tf_downregulated$MapMan_terms),]$TF_class <- 'HSF,Heat-shock transcription factor'
Sprint_tf_downregulated[grepl('Psudo', Sprint_tf_downregulated$MapMan_terms),]$TF_class <- 'Pseudo ARR transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN2567_c0_g3_i2',]$TF_class <- 'ABI3/VP1-related B3-domain-containing transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN41917_c0_g1_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN11715_c0_g1_i1',]$TF_class <- 'ARR'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN7892_c0_g1_i2',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_30_DN9382_c0_g1_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN46484_c0_g1_i1',]$TF_class <- 'HB,Homeobox transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN5964_c0_g3_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN24994_c0_g2_i1',]$TF_class <- 'ovate family OFP'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN3507_c0_g1_i2',]$TF_class <- 'HB,Homeobox transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN3529_c0_g1_i1',]$TF_class <- 'MADS box transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN40686_c0_g1_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN11553_c2_g1_i3',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN9418_c0_g1_i1',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN46024_c0_g1_i1',]$TF_class <- 'NIN-like bZIP-related family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN5126_c0_g1_i5',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN7359_c0_g1_i3',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN13396_c0_g1_i1',]$TF_class <- 'HB,Homeobox transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN1700_c0_g2_i1',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN1643_c0_g1_i4',]$TF_class <- 'MADS box transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN1643_c0_g1_i12',]$TF_class <- 'MADS box transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN1643_c0_g1_i13',]$TF_class <- 'MADS box transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN3905_c0_g1_i2',]$TF_class <- 'NAC domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN347_c0_g2_i2',]$TF_class <- 'HB,Homeobox transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN9972_c0_g2_i1',]$TF_class <- 'GRAS transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN9996_c0_g2_i1',]$TF_class <- 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN3427_c0_g1_i2',]$TF_class <- 'bZIP transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN4578_c0_g1_i4',]$TF_class <- 'WRKY domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN4578_c0_g1_i1',]$TF_class <- 'WRKY domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN10515_c0_g2_i3',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN10515_c0_g2_i2',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN395_c0_g1_i1',]$TF_class <- 'C2C2(Zn) DOF zinc finger family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN395_c0_g3_i1',]$TF_class <- 'C2C2(Zn) DOF zinc finger family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN10969_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN21536_c0_g1_i1',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN1926_c0_g1_i2',]$TF_class <- 'MYB-related transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN11302_c0_g1_i2',]$TF_class <- 'PHD finger transcription factor'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN1801_c0_g1_i7',]$TF_class <- 'PHD finger transcription factor'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN2729_c0_g1_i3',]$TF_class <- 'bZIP transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN3666_c4_g1_i1',]$TF_class <- 'ARF, Auxin Response Factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN27478_c0_g1_i1',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN14028_c0_g2_i1',]$TF_class <- 'putative transcription regulator'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN2482_c1_g1_i3',]$TF_class <- 'MYB domain transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN7695_c0_g1_i1',]$TF_class <- 'bZIP transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN10471_c0_g1_i4',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN5676_c0_g1_i1',]$TF_class <- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN6881_c0_g1_i9',]$TF_class<- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN6881_c0_g1_i12',]$TF_class<- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN10619_c0_g2_i1',]$TF_class<- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN3012_c0_g1_i2',]$TF_class<- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN5160_c0_g1_i1',]$TF_class<- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN10969_c0_g1_i1',]$TF_class<- 'HSF,Heat-shock transcription factor family'
Sprint_tf_downregulated[Sprint_tf_downregulated$transcript_id == 'TRINITY_10_DN9922_c0_g1_i1',]$TF_class<- 'HSF,Heat-shock transcription factor family'

Sprint_tf_downregulated <- rbind(Sprint_tf_downregulated, data.frame(transcript_id = 'TRINITY_10_DN27478_c0_g1_i1',
    BLASTP_names = 'protein PHR1-LIKE 3 [Medicago truncatula]', MapMan_terms = 'RNA.regulation of transcription.G2-like transcription factor family, GARP,RNA.regulation of transcription.MYB domain transcription factor family',
    TF_class = 'G2-like transcription factor family, GARP'))
Sprint_tf_downregulated <- rbind(Sprint_tf_downregulated, data.frame(transcript_id = 'TRINITY_10_DN2482_c1_g1_i3',
                                                                     BLASTP_names = 'protein PHR1-LIKE 3 [Medicago truncatula]', MapMan_terms = 'RNA.regulation of transcription.G2-like transcription factor family, GARP,RNA.regulation of transcription.MYB domain transcription factor family',
                                                                     TF_class = 'G2-like transcription factor family, GARP'))

Sprint_tf_downregulated_pivot <- table(unlist(Sprint_tf_downregulated$TF_class)) %>% as.data.frame() %>% arrange(desc(Freq))
Sprint_tf_downregulated_pivot %>% View()
Sprint_tf_downregulated[grepl('NIN-like bZIP-related family', Sprint_tf_downregulated$TF_class),]


annot_full[annot_full$transcript_id %in% Sprint_tf_downregulated[grepl('signal transduction', Sprint_tf_downregulated$MapMan_terms),]$transcript_id,] %>% dplyr::select(transcript_id, MapMan_terms)



annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% upregulated_30_days$target_id) %>%
  dplyr::select(transcript_id, BLASTP_names, MapMan_terms) %>% View()

sapply(ls()[grep('Sprint_tf_((up)|(down))', ls())], function(x) x %>% get() %>% 
         write.table(., file.path('~/Pea_transcriptomics/TF_lists/', paste(x, 'tsv', sep = '.')), col.names = T, row.names = F, sep = '\t'))

MapMan_df %>% filter(MapMan_Category2 == 'regulation of transcription') %>% dplyr::select(MapMan_Category3) %>% unlist() %>% unname() %>% unique()
annot_full[grepl('GO:0045735', annot_full$GOs),]  %>%
  dplyr::select(transcript_id, BLASTP_names, MapMan_terms) %>% View()

annot_full[grepl('regulation of transcription', annot_full$MapMan_terms) & grepl('abscisic', annot_full$MapMan_terms),] %>% filter(transcript_id %in% downregulated_30_days$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)

annot_full[grepl('(heat shock)|(chaperone)|(hsp)', annot_full$BLASTP_names),] %>% 
  dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% filter(transcript_id %in% upregulated_30_days$target_id) %>% View()

group_B <- read.table('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/Group Overexpression Spots B.csv', sep = '\t', header = T) %>% dplyr::rename(target_id = X)
group_C <- read.table('~/Pea_transcriptomics/hormones_oxerexpression_results_290919/Group Overexpression Spots C.csv', sep = '\t', header = T) %>% dplyr::rename(target_id = X)

abscisic_B <- annot_full[grepl('(GO:0009688)', annot_full$GOs),] %>% filter(transcript_id %in% group_B$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)
gibberellin_B <- annot_full[grepl('GO:0009686', annot_full$GOs),] %>% filter(transcript_id %in% group_B$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)
gibberellin_B %>% View()
abscisic_B %>% View()

gibberellin_C <- annot_full[grepl('(GO:0009685)|(GO:0009740)|(GO:0009686)', annot_full$GOs),] %>% filter(transcript_id %in% group_C$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)
abscisic_C <- annot_full[grepl('(GO:0046345)|(GO:0009737)', annot_full$GOs),] %>% filter(transcript_id %in% group_C$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)
abscisic_C %>% View()
gibberellin_C %>% View()
sapply(ls()[grepl('^zh((e)|(ong))_((up)|(down))', ls())], function(x) x %>% get() %>% write.table(
   ., file.path('~/Pea_transcriptomics/chinese_sleuth/', paste(x, 'tsv', sep = '.')), sep = '\t', row.names = F, col.names = T))

sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_1_up, x, cutoff = 0.001, make_table = T))
sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_1_down, x, cutoff = 0.001, make_table = T))
sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_2_up, x, cutoff = 0.001, make_table = T))
sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_2_down, x, cutoff = 0.001, make_table = T))

sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_4_down, x, cutoff = 0.001, make_table = T))
sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_4_up, x, cutoff = 0.001, make_table = T))
sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_5_down, x, cutoff = 0.001, make_table = T))
sapply(c('BP', 'CC', "MF"), function(x) iterating_function_GO(cross_5_up, x, cutoff = 0.001, make_table = T))

sapply(ls()[grepl('((BP)|(CC)|(MF))_cross_[12]_((up)|(down))', ls())], function(x) x %>% 
         get() %>% write.table(., file.path('~/Pea_transcriptomics/10vs10_new/', paste(x, 'tsv', sep = '.')), 
                               col.names = T, row.names = F, sep = '\t'))

sapply(ls()[grepl('((BP)|(CC)|(MF))_cross_[45]_((up)|(down))', ls())], function(x) x %>% 
         get() %>% write.table(., file.path('~/Pea_transcriptomics/20vs25_new/', paste(x, 'tsv', sep = '.')), 
                               col.names = T, row.names = F, sep = '\t'))

annot_full[grepl('abscisic', annot_full$MapMan_terms),] %>% filter(transcript_id %in% cross_4_down$target_id) %>% dplyr::select(transcript_id, MapMan_terms, BLASTP_names) %>% View()

BP_cross_2_up %>% View()
cross_4_up %>% filter(target_id %in% cross_5_up$target_id) %>% nrow()
models(sleuth_cross_1)

### Find all gibberellin signalling absent in Sprint-2_10
MapMan_df %>% filter(MapMan_Category2 == 'gibberelin') %>% dplyr::select(MapMan_Category3) %>% unlist() %>% unname() %>% unique()
MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism') %>% dplyr::select(MapMan_Category2) %>% unlist() %>% unique()

annot_full[grepl('GO:0004407', annot_full$GOs),] %>% filter(transcript_id %in% cross_5_up$target_id) %>% dplyr::select(transcript_id, MapMan_terms, BLASTP_names) %>% View()
cross_4_up %>% nrow()
cross_5_up %>% nrow()

cross_4_down_TF <- MapMan_df %>% filter(transcript_id %in% cross_4_down$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_4_down_TF %>% View()
cross_4_down_TF_pivot <- table(unlist(cross_4_down_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_4_up_TF <- MapMan_df %>% filter(transcript_id %in% cross_4_up$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_4_up_TF_pivot <- table(unlist(cross_4_up_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_5_down_TF <- MapMan_df %>% filter(transcript_id %in% cross_5_down$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_5_down_TF_pivot <- table(unlist(cross_5_down_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_5_up_TF <- MapMan_df %>% filter(transcript_id %in% cross_5_up$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_5_up_TF_pivot <- table(unlist(cross_5_up_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_1_up_TF <- MapMan_df %>% filter(transcript_id %in% cross_1_up$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_1_up_TF_pivot <- table(unlist(cross_1_up_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_1_down_TF <- MapMan_df %>% filter(transcript_id %in% cross_1_down$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_1_down_TF_pivot <- table(unlist(cross_1_down_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_2_up_TF <- MapMan_df %>% filter(transcript_id %in% cross_2_up$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_2_up_TF_pivot <- table(unlist(cross_2_up_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

cross_2_down_TF <- MapMan_df %>% filter(transcript_id %in% cross_2_down$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
cross_2_down_TF_pivot <- table(unlist(cross_2_down_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

Sprint2_down_TF <- MapMan_df %>% filter(transcript_id %in% downregulated_30_days$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
Sprint2_down_TF_pivot <- table(unlist(Sprint2_down_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

Sprint2_up_TF <- MapMan_df %>% filter(transcript_id %in% upregulated_30_days$target_id, MapMan_Category2 == 'regulation of transcription') %>%
  dplyr::select(MapMan_Category2, MapMan_Category3)
Sprint2_up_TF_pivot <- table(unlist(Sprint2_up_TF$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq))

sapply(ls()[grepl('TF_pivot', ls())], function(x) x %>% get() %>% write.table(.,
     file.path('~/Pea_transcriptomics/TF_lists_new_061019/', paste(x, 'tsv', sep = '.')), col.names = T, row.names = F, sep = '\t'))

annot_full[grepl('regulation of transcription\\.unclassified', annot_full$MapMan_terms),]%>% filter(transcript_id %in% cross_4_down$target_id) %>% dplyr::select(transcript_id, MapMan_terms, BLASTP_names) %>% View()


table(unlist((MapMan_df %>% filter(transcript_id %in% upregulated_30_days$target_id, MapMan_Category2 == 'regulation of transcription') %>%
                dplyr::select(MapMan_Category2, MapMan_Category3))$MapMan_Category3)) %>% as.data.frame() %>% arrange(desc(Freq)) %>% View()
models(sleuth_cross_4)

annot_full %>% colnames()
annot_full[grepl('GO:0009686', annot_full$GOs) | grepl('gibberelin.synthesis-degradat', annot_full$MapMan_terms),] %>% filter(transcript_id %in% upregulated_30_days$target_id) %>% dplyr::select(transcript_id, MapMan_terms, BLASTP_names, KEGG_Pathway) %>% View()
annot_full[grepl('GO:0009686', annot_full$GOs),] %>% filter(transcript_id %in% upregulated_30_days$target_id) %>% dplyr::select(transcript_id, MapMan_terms, BLASTP_names) %>% View()

annot_full[grepl('regulation of transcription\\.unclassified', annot_full$MapMan_terms),]%>% filter(transcript_id %in% cross_4_down$target_id) %>% dplyr::select(transcript_id, MapMan_terms, BLASTP_names) %>% View()

MF_Sprint_DEGs %>% str()
ggplot(MF_Sprint_DEGs, aes(Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes)) + coord_flip() +
  geom_bar(aes(x=MF_Sprint_DEGs[MF_Sprint_DEGs$condition == 'up',]),stat = 'identity', position = 'identity')

GO_barplotter <- function(dataset) {
  ggplot(dataset, aes(y=Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes,
                                                               x=reorder(Term, Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes), fill = condition)) + 
    geom_bar(stat = 'identity', position = 'identity', width = 0.8) + coord_flip() +
    theme(text = element_text(size=28, family = 'arial', color='black'), axis.text.x = element_text(color = 'black', family = 'arial'),
          axis.title.x = element_text(color = 'black', family = 'arial', size = 20),
          axis.text.y = element_text(color = 'black', family = 'arial', hjust = 1),
          axis.title.y = element_blank(),
          axis.line.x = element_line(color = 'black'),
          axis.ticks.length.y =  unit(0, 'lines'),
          panel.background = element_blank(), legend.text = element_text(size = 22), legend.title = 
            element_text(size = 28)) + 
    labs(y = '% of annotated genes') + 
    scale_fill_manual('condition', values = c('Upregulated' = '#068013', 'Downregulated' = '#9c1718'))
  }


ggplot(MF_Sprint_DEGs, aes(y=Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes,
  x=reorder(Term, Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes), fill = condition)) + 
  geom_bar(stat = 'identity', position = 'identity', width = 0.9) + coord_flip() +
  theme(text = element_text(size=30, family = 'arial', color='black'), axis.text.x = element_text(color = 'black', family = 'arial'),
        axis.text.y = element_text(color = 'black', family = 'arial', hjust = 1),
        axis.title.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.length.y =  unit(0, 'lines'),
        panel.background = element_blank(), legend.text = element_text(size = 20), legend.title = 
          element_text(size = 22)) + 
  labs(y = '% of annotated genes') + 
  scale_fill_manual('condition', values = c('up' = '#4db850', 'down' = '#910f27'))

##### бины поуже, текст побольше, приблизить текст к бинам
BP_Sprint_DEGs %>% View()
ggplot(BP_Sprint_DEGs, aes(y=Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes,
                           x=reorder(Term, Percentage_of_Differentially_Expressed_Genes_from_Number_Of_Annotated_Genes), fill = condition)) + 
  geom_bar(stat = 'identity', position = 'identity', width = 0.9) + coord_flip() +
  theme(text = element_text(size=30, family = 'arial', color='black'), axis.text.x = element_text(color = 'black', family = 'arial'),
        axis.text.y = element_text(color = 'black', family = 'arial', hjust = 1),
        axis.title.y = element_blank(),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.length.y =  unit(0, 'lines'),
        panel.background = element_blank(), legend.text = element_text(size = 20), legend.title = 
          element_text(size = 22)) + 
  labs(y = '% of annotated genes') + 
  scale_fill_manual('condition', values = c('Upregulated' = '#4db850', 'Downregulated' = '#910f27'))

GO_barplotter(CC_Sprint_DEGs)
GO_barplotter(BP_Sprint_DEGs)
GO_barplotter(MF_Sprint_DEGs)
sapply(ls()[grepl('^zh(e|(ong))_((up)|(down))$', ls())], function(x) x %>% get() %>% 
         write.table(., file.path('~/Pea_transcriptomics/cross_071019/', paste(x, 'tsv', sep = '.')), row.names = F, col.names = T, sep = '\t'))

sapply(ls()[grepl('^cross_[45]_((up)|(down))$', ls())], function(x) x %>% get() %>% 
         write.table(., file.path('~/Pea_transcriptomics/cross_071019/', paste(x, 'tsv', sep = '.')), row.names = F, col.names = T, sep = '\t'))


sapply(ls()[grepl('^((BP)|(CC)|(MF))_zh', ls())], function(x) x %>% get() %>% 
         write.table(., file.path('~/Pea_transcriptomics/cross_071019/', paste(x, 'tsv', sep = '.')), row.names = F, col.names = T, sep = '\t'))

BP_cross_1_up %>% View()
BP_cross_1_down %>% View()

BP_cross_1_DEGs <- rbind(BP_cross_1_up, BP_cross_1_down[c(1,3,4,5,7,8, 9, 10, 14, 17),]) %>% mutate(condition = c(rep('Downregulated', 11), rep('Upregulated', 10)))
BP_cross_1_DEGs %>% View()
BP_cross_1_DEGs$Term[2] <- 'photosynthesis, light harvesting in photosystem I'
BP_cross_1_DEGs$Term[5] <- 'cellular response to abscisic acid stimulus'
BP_cross_1_DEGs$Term[7] <- 'mitochondrial ATP synthesis coupled electron transport'
BP_cross_1_DEGs$Term[11] <- 'negative regulation of catalytic activity'
BP_cross_1_DEGs$Term[15] <- 'cellular response to abscisic acid stimulus'
BP_cross_1_DEGs$Term[17] <- 'abscisic acid-activated signaling pathway'
BP_cross_1_DEGs$Term[19] <- 'positive regulation of response to water deprivation'

GO_barplotter(BP_cross_1_DEGs)


CC_cross_1_up %>% View()
CC_cross_1_down %>% View()

CC_cross_1_DEGs <- rbind(CC_cross_1_up[c(1,4,5,8,9,12,13,20, 22, 23),], CC_cross_1_down) %>% mutate(condition = c(rep('Downregulated', 10), rep('Upregulated', 4)))
CC_cross_1_DEGs$Term[4] <- 'inner mitochondrial membrane protein complex'
GO_barplotter(CC_cross_1_DEGs)

MF_cross_1_up %>% nrow()

MF_cross_1_DEGs <- rbind(MF_cross_1_up, MF_cross_1_down) %>% mutate(condition = c(rep('Downregulated', 8), rep('Upregulated', 6)))
MF_cross_1_DEGs$Term[5] <- 'oxidoreductase activity, CH-OH group of donors, NAD/NADP as acceptor'
MF_cross_1_DEGs$Term[7] <- 'oxidoreductase activity, oxidizing metal ions'
MF_cross_1_DEGs$Term[11] <- ' DNA-binding transcription factor activity'
GO_barplotter(MF_cross_1_DEGs)

BP_cross_2_down %>% View()
BP_cross_2_up %>% View()

BP_cross_2_DEGs <- rbind(BP_cross_2_up, BP_cross_2_down[c(2:7,11, 13, 15),]) %>% mutate(condition = c(rep('Downregulated', 6), rep('Upregulated', 9)))
BP_cross_2_DEGs %>% View()
BP_cross_2_DEGs$Term[2] <- 'nucleic acid phosphodiester bond hydrolysis'
BP_cross_2_DEGs$Term[4] <- 'regulation of stomatal complex development'
BP_cross_2_DEGs$Term[10] <- 'cellular response to abscisic acid stimulus'
BP_cross_2_DEGs$Term[15] <- 'abscisic acid-activated signaling pathway'

GO_barplotter(BP_cross_2_DEGs)

CC_cross_2_down %>% nrow()
CC_cross_2_up %>% nrow()

CC_cross_2_DEGs <- rbind(CC_cross_2_up, CC_cross_2_down) %>% mutate(condition = c(rep('Downregulated', 9), rep('Upregulated', 6)))
CC_cross_2_DEGs$Term[4] <- 'NAD(P)H dehydrogenase complex (plastoquinone)'
GO_barplotter(CC_cross_2_DEGs)

MF_cross_2_down %>% nrow()
MF_cross_2_up %>% nrow()

MF_cross_2_DEGs <- rbind(MF_cross_2_up, MF_cross_2_down) %>% mutate(condition = c(rep('Downregulated', 3), rep('Upregulated', 5)))
MF_cross_2_DEGs$Term[4] <- 'hydrolase activity, acting on ester bonds'
MF_cross_2_DEGs$Term[7] <- 'DNA-binding transcription factor activity'
GO_barplotter(MF_cross_2_DEGs)

BP_cross_4_down %>% View()
BP_cross_4_up %>% View()

BP_cross_4_DEGs <- rbind(BP_cross_4_up, BP_cross_4_down[c(1,3:7,10,12, 14,17),]) %>% mutate(condition = c(rep('Downregulated', 6), rep('Upregulated', 10)))
BP_cross_4_DEGs %>% View()
BP_cross_4_DEGs$Term[10] <- 'maturation of SSU-rRNA from tricistronic rRNA transcript'
BP_cross_4_DEGs$Term[11] <- 'tRNA aminoacylation for protein translation'
GO_barplotter(BP_cross_4_DEGs)

CC_cross_4_down %>% View()
CC_cross_4_up %>% nrow()

CC_cross_4_DEGs <- rbind(CC_cross_4_up, CC_cross_4_down[c(2, 3, 4, 5, 7, 8, 9, 19, 21, 24),]) %>% mutate(condition = c(rep('Downregulated', 3), rep('Upregulated', 10)))
CC_cross_4_DEGs %>% View()
GO_barplotter(CC_cross_4_DEGs)

MF_cross_4_down %>% View()
MF_cross_4_up %>% View()

MF_cross_4_DEGs <- rbind(MF_cross_4_up, MF_cross_4_down) %>% mutate(condition = c(rep('Downregulated', 3), rep('Upregulated', 10)))
MF_cross_4_DEGs %>% View()
MF_cross_4_DEGs$Term[10] <- 'protein folding chaperone'
GO_barplotter(MF_cross_4_DEGs)

BP_cross_5_down %>% View()
BP_cross_5_up %>% nrow()

BP_cross_5_DEGs <- rbind(BP_cross_5_up[-c(2,10),], BP_cross_5_down[-c(2, 3, 6, 8, 14, 15),]) %>% mutate(condition = c(rep('down', 9), rep('up', 10)))
BP_cross_5_DEGs %>% View()
BP_cross_5_DEGs$Term[2] <- 'double-strand break repair via homologous recombination'
BP_cross_5_DEGs$Term[3] <- 'positive regulation of stem cell population maintenance'
BP_cross_5_DEGs$Term[5] <- 'cellular response to abscisic acid stimulus'
BP_cross_5_DEGs$Term[7] <- 'regulation of multicellular organismal development'
BP_cross_5_DEGs$Term[15] <- 'maturation of LSU-rRNA from tricistronic rRNA transcript'
GO_barplotter(BP_cross_5_DEGs)

CC_cross_5_down %>% View()
CC_cross_5_up %>% nrow()

CC_cross_5_DEGs <- rbind(CC_cross_5_up, CC_cross_5_down[-c(3, 9, 10, 11),]) %>% mutate(condition = c(rep('down', 3), rep('up', 10)))
CC_cross_5_DEGs %>% View()
GO_barplotter(CC_cross_5_DEGs)

MF_cross_5_down %>% View()
MF_cross_5_up %>% View()

MF_cross_5_DEGs <- rbind(MF_cross_5_up, MF_cross_5_down) %>% mutate(condition = c(rep('down', 2), rep('up', 7)))
MF_cross_5_DEGs %>% View()
MF_cross_5_DEGs$Term[7] <- 'phosphoenolpyruvate carboxykinase activity'
GO_barplotter(MF_cross_5_DEGs)
lrt %>% mutate(b = wt[match(lrt$target_id, wt$target_id),]$b) %>% write.table('~/Pea_transcriptomics/Sprint2_diffexprt.tsv',
                                                                              col.names = T, row.names = F, sep = '\t')
zhe_lrt %>% mutate(b = zhe_wt[match(zhe_lrt$target_id, zhe_wt$target_id),]$b) %>% write.table('~/Pea_transcriptomics/Zhewan1_diffexpr.tsv',
                                                                              col.names = T, row.names = F, sep = '\t')
zhong_lrt %>% mutate(b = zhong_wt[match(zhong_lrt$target_id, zhong_wt$target_id),]$b) %>% write.table('~/Pea_transcriptomics/Zhongwan6_diffexpr.tsv',
                                                                                              col.names = T, row.names = F, sep = '\t')



##### Hormone-restricted GO terms: expansion
hormone_GOs <- c('GO:0016707', 'GO:0010336', 'GO:001047', 'GO:0045487', 'GO:0009739', 'GO:0009686', 'GO:0009685',
                 'GO:0033470', 'GO:0033469', 'GO:0071370', 'GO:0010371', 'GO:0010372', 'GO:0010373', 'GO:0010331',
                 'GO:0045544', 'GO:0102713', 'GO:1905201', 'GO:0047927', 'GO:0102111', 'GO:0102119', 'GO:0102118',
                 'GO:0102117', 'GO:0102122', 'GO:0045543', 'GO:0102714', 'GO:0102716', 'GO:0102715', 'GO:0102712',
                 'GO:0102711', 'GO:0016707', 'GO:0102653', 'GO:0102652', 'GO:0102663', 'GO:0102924', 'GO:0051779', 'GO:0047928',
                 'GO:0010341', 'GO:0102123', 'GO:0102125', 'GO:0102124', 'GO:0102738', 'GO:0102739', 'GO:0103057', 'GO:0103056',
                 'GO:0052634', 'GO:0052635', 'GO:0103010', 'GO:0102972', 'GO:0103054', 'GO:0102755',
                 'GO:0080168', 'GO:1902265', 'GO:0010427', 'GO:0009687', 'GO:0010294', 'GO:0010293',
                 'GO:1902266', 'GO:0046345', 'GO:0009737', 'GO:0009688', 'GO:0090440', 'GO:0010295', 'GO:0009724', 'GO:0009738',
                 'GO:0071215', 'GO:0010115', 'GO:0090359', 'GO:0051993', 'GO:1902418', 'GO:0010116', 'GO:0009787', 'GO:1902417',
                 'GO:0009789', 'GO:0009788', 'GO:0075343', 'GO:1901527', 'GO:1990218',
                 'GO:0060918', 'GO:0060919', 'GO:0010011', 'GO:0010252', 'GO:0010315', 'GO:0009850', 'GO:0038198', 'GO:0080162',
                 'GO:0009733', 'GO:0009852', 'GO:0009851', 'GO:0010541', 'GO:0010540', 'GO:0009926', 'GO:0010249', 'GO:0010328',
                 'GO:0080161', 'GO:0009721', 'GO:0009734', 'GO:0009921', 'GO:0090354', 'GO:2000012', 'GO:0071365', 'GO:0010329',
                 'GO:0010600', 'GO:0090355', 'GO:0090356', 'GO:0010601', 'GO:0010928', 'GO:1901703', 'GO:0010930', 'GO:0010929',
                 'GO:0060774', 'GO:0090015', 'GO:0044032', 'GO:1990210', 'GO:0044846', 'GO:0009672',
                 'GO:0010184', 'GO:0044373', 'GO:0009690', 'GO:0009735', 'GO:0009691', 'GO:0009823', 'GO:0009884', 'GO:0019139',
                 'GO:0009722', 'GO:0009736', 'GO:0080062', 'GO:1903856', 'GO:0047807', 'GO:0071368', 'GO:0080036', 'GO:1903857',
                 'GO:0001647', 'GO:0009885', 'GO:0080037', 'GO:0080038', 'GO:1990223', 'GO:0009824', 'GO:0050447',
                 'GO:0051740', 'GO:0042457', 'GO:0038199', 'GO:0009723', 'GO:0009693', 'GO:0009692', 'GO:0009727', 'GO:0009873',
                 'GO:0038200', 'GO:0071369', 'GO:0010364', 'GO:0075022', 'GO:0010365', 'GO:0010366', 'GO:0010104', 'GO:0102276',
                 'GO:0010105', 'GO:0009861', 'GO:0009866', 'GO:0052021', 'GO:0052084', 'GO:1990212', 'GO:0052076',
                 'GO:0052005', 'GO:0043272', 'GO:0009871', 'GO:0009868', 'GO:0032260', 'GO:0052441', 'GO:0009815', 'GO:0004085',
                 'GO:0052276',
                 'GO:0030795', 'GO:0080123', 'GO:0102078', 'GO:0080032', 'GO:0009867', 'GO:0120091', 'GO:0009694', 'GO:0009753',
                 'GO:0009695', 'GO:0009754', 'GO:0009754', 'GO:0080140', 'GO:0080141', 'GO:0071395', 'GO:2000022', 'GO:0009861',
                 'GO:0009864', 'GO:0052022', 'GO:0052088', 'GO:1990211', 'GO:0052075', 'GO:0052068', 'GO:0043272', 'GO:0009871',
                 'GO:0009868', 'GO:0052073', 'GO:0032260', 'GO:0080032', 'GO:0052480', 'GO:0030795', 'GO:0080123',
                 'GO:0090411', 'GO:0010268', 'GO:0016131', 'GO:0016132', 'GO:0016133', 'GO:0080118', 'GO:0009741', 'GO:0009742',
                 'GO:0009729', 'GO:0071367', 'GO:0010422', 'GO:1900457', 'GO:2000488', 'GO:0010423', 'GO:1900458', 'GO:1900459',
                 'GO:1901149', 'GO:0009696', 'GO:0046244', 'GO:0009751', 'GO:0009697', 'GO:0009752', 'GO:0009863', 'GO:0010337',
                 'GO:0009627', 'GO:0080142', 'GO:0052639', 'GO:0052640', 'GO:0071446', 'GO:2000031', 'GO:0080151', 'GO:0009862',
                 'GO:0052023', 'GO:0010679', 'GO:0052089', 'GO:1990213', 'GO:0052074', 'GO:0052004', 'GO:0052081', 'GO:0052072',
                 'GO:0052003', 'GO:0052284', 'GO:0052253', 'GO:0052272', 'GO:0080031', 'GO:0018658', 'GO:0052624', 'GO:1901672',
                 'GO:1902347', 'GO:1902348', 'GO:1901601', 'GO:1901600')

annot_full[grepl('hormone', annot_full$MapMan_terms),] %>% nrow()

hormone_annotations <-annot_full[grepl(paste0('(', paste(hormone_GOs, collapse = ')|('), ')', sep = ''), annot_full$GOs) | 
             grepl('hormone', annot_full$MapMan_terms),]
MapMan_df %>% filter(MapMan_Category1 == 'RNA') %>% dplyr::select(MapMan_Category2) %>% unique()

# Now perform a SOM clusterization
dummy %>% colnames()
dummy %>% head()
dummy <- raw_for_oposSOM_processed %>% filter(transcript_id %in% hormone_annotations$transcript_id) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                          as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(6:25)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4",
                             'Zhewan1_10_1', 'Zhewan1_10_2', 'Zhewan1_10_3', 'Zhewan1_10_4',
                             'Zhewan1_25_1','Zhewan1_25_2','Zhewan1_25_3','Zhewan1_25_4',
                             'Zhongwan6_10_1', 'Zhongwan6_10_2', 'Zhongwan6_10_3', 'Zhongwan6_10_4',
                             'Zhongwan6_20_1', 'Zhongwan6_25_2', 'Zhongwan6_25_3', 'Zhongwan6_25_4')

rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'Hormones', dim.1stLvlSom = 30))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_all_hormones_131019//')
opossom.run(env = env_synthesis)

group_A %>% View()
group_A <- read.table('/home/reverend_casy/Pea_transcriptomics/hormones_overexpression_results_131019/Group Overexpression Spots A.csv',
                      header = T,  sep = '\t', row.names = 1) %>% mutate(target_id = rownames(.))
group_B <- read.table('/home/reverend_casy/Pea_transcriptomics/hormones_overexpression_results_131019/Group Overexpression Spots B.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_C <- read.table('/home/reverend_casy/Pea_transcriptomics/hormones_overexpression_results_131019/Group Overexpression Spots C.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_D <- read.table('/home/reverend_casy/Pea_transcriptomics/hormones_overexpression_results_131019/Group Overexpression Spots D.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_E <- read.table('/home/reverend_casy/Pea_transcriptomics/hormones_overexpression_results_131019/Group Overexpression Spots E.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_F <- read.table('/home/reverend_casy/Pea_transcriptomics/hormones_overexpression_results_131019/Group Overexpression Spots F.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
annot_full[grepl(paste('(', paste(GA, collapse = ')|('), ')', sep = ''), annot_full$GOs) | grepl('gibber', annot_full$MapMan_terms),] %>% filter(
  transcript_id %in% group_F$target_id) %>% nrow()

annot_full[grepl('GO:0009739', annot_full$GOs),] %>% filter(transcript_id %in% group_F$target_id) %>% 
  dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% write.table('~/Pea_transcriptomics/hormones_overexpression_results_131019/spot_f_GA_response.tsv',
                                                                           col.names = T, row.names = F, sep = '\t')
annot_full %>% filter(transcript_id == 'TRINITY_30_DN6673_c0_g1_i1') %>% dplyr::select(GOs, MapMan_terms)
annot_full %>% View()
group_F %>% View()
upregulated_30_days %>% filter(target_id)

sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_B, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_C, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_D, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_E, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_F, x, cutoff = 0.01, make_table = T))

sapply(ls()[grepl('((BP)|(CC)|(MF))_group', ls())], function(x) get(x) %>% write.table(
  file.path('~/Pea_transcriptomics/hormones_overexpression_results_131019/GO_reports/', paste(x, 'tsv', sep= '.')),
col.names = T, row.names = F, sep = '\t'))

library(UpSetR)


#### All transcripts - a final attempt
dummy <- raw_for_oposSOM_processed  %>% mutate(transcript_id = 
      as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(6:25)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4",
                              'Zhewan1_10_1', 'Zhewan1_10_2', 'Zhewan1_10_3', 'Zhewan1_10_4',
                              'Zhewan1_25_1','Zhewan1_25_2','Zhewan1_25_3','Zhewan1_25_4',
                              'Zhongwan6_10_1', 'Zhongwan6_10_2', 'Zhongwan6_10_3', 'Zhongwan6_10_4',
                              'Zhongwan6_25_1', 'Zhongwan6_25_2', 'Zhongwan6_25_3', 'Zhongwan6_25_4')

rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'Global Transcriptomic Ladnscape', dim.1stLvlSom = 40))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_global_141019//')
opossom.run(env = env_synthesis)


group_A <- read.table('/home/reverend_casy/Pea_transcriptomics/global_spots_141019/Group Overexpression Spots A.csv',
                      header = T,  sep = '\t', row.names = 1) %>% mutate(target_id = rownames(.))
group_B <- read.table('/home/reverend_casy/Pea_transcriptomics/global_spots_141019/Group Overexpression Spots B.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_C <- read.table('/home/reverend_casy/Pea_transcriptomics/global_spots_141019/Group Overexpression Spots C.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_D <- read.table('/home/reverend_casy/Pea_transcriptomics/global_spots_141019/Group Overexpression Spots D.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_E <- read.table('/home/reverend_casy/Pea_transcriptomics/global_spots_141019/Group Overexpression Spots E.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))


sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_B, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_C, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_D, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_E, x, cutoff = 0.01, make_table = T))
sapply(ls()[grepl('((BP)|(CC)|(MF))_group_[ABCDE]', ls())], function(x) x %>% get() %>% 
         write.table(., file.path('~/Pea_transcriptomics/global_spots_141019/GO_reports/', paste(x, 'tsv', sep = '.')),
                     col.names = T, row.names = F, sep = '\t'))

### TF SOM - a final attempt
TF_GO <- c('GO:0003700', 'GO:0001217', 'GO:0001216', 'GO:0051091', 'GO:0043433', 'GO:0051090', 'GO:0098531',
           'GO:0034246', 'GO:0005667', 'GO:0000981', 'GO:0003712', 'GO:0003714', 'GO:0003713',
           'GO:0006325', 'GO:1990700', 'GO:0006338', 'GO:0035092', 'GO:0070827', 'GO:0016569',
           'GO:0070828', 'GO:0090202', 'GO:0006333', 'GO:1905269', 'GO:1902275', 'GO:1905268',
           'GO:0034401', 'GO:0001301', 'GO:0034728', 'GO:0039525', 'GO:0061641')

TF_annotation <- annot_full[grepl(paste0('(', paste(TF_GO, collapse = ')|('), ')', sep = ''), annot_full$GOs) | 
                                   grepl('regulation of transcription', annot_full$MapMan_terms),]

dummy <- raw_for_oposSOM_processed %>% filter(transcript_id %in% TF_annotation$transcript_id) %>% dplyr::select(c(1,6:29)) %>% mutate(transcript_id = 
                                                                                                                                              as.character(transcript_id)) %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate_at(vars(-transcript_id), all_vars(log(.)))
colnames(dummy)[c(6:25)] <- c("Sprint2_20_1", "Sprint2_20_2", "Sprint2_20_3", "Sprint2_20_4",
                              'Zhewan1_10_1', 'Zhewan1_10_2', 'Zhewan1_10_3', 'Zhewan1_10_4',
                              'Zhewan1_25_1','Zhewan1_25_2','Zhewan1_25_3','Zhewan1_25_4',
                              'Zhongwan6_10_1', 'Zhongwan6_10_2', 'Zhongwan6_10_3', 'Zhongwan6_10_4',
                              'Zhongwan6_20_1', 'Zhongwan6_25_2', 'Zhongwan6_25_3', 'Zhongwan6_25_4')

rownames(dummy) <- dummy$transcript_id
dummy %<>% dplyr::select(-c(transcript_id))
dummy[dummy == -Inf] <- 0
dummy[dummy == Inf] <- 0
dummy <- dummy[!(apply(dummy, 1, function(x) all(x == 0))),]
env_synthesis <- opossom.new(list(dataset.name = 'Transcription-related proteins', dim.1stLvlSom = 30))
env_synthesis$indata <- dummy %>% as.matrix()
env_synthesis$group.labels <- c(rep('Sprint2_10', 4), rep('Sprint2_20', 4), rep('Zhewan1_10', 4), rep('Zhewan1_25', 4), rep('Zhongwan6_10', 4), rep('Zhongwan6_25', 4))
env_synthesis$group.colors <- c(rep('#00b6c2', 4),rep('#0006c2', 4), rep('#dfc9e2', 4),rep('#dfe9e2', 4), rep('#228AD0', 4), rep('#228AC0', 4))
setwd('~/Pea_transcriptomics/oposSOM_TF_new_201019/')
opossom.run(env = env_synthesis)



### GA-related transcripts
GA_annotations <- read.table('~/Pea_transcriptomics/GA_refs/annotated_GA_proteins.tsv', sep = '\t', header = T, stringsAsFactors = F)
GA_annotations %<>% mutate(Protein = sapply(GA_annotations$Protein, function(x) substring(x, 1, nchar(x) - 3) %>% unname()))

GA_annotations %>% filter(Protein %in% upregulated_30_days$target_id)
GA_annotations %>% filter(Protein %in% downregulated_30_days$target_id)

GA_annotations %>% filter(Protein %in% zhe_up$target_id)
GA_annotations %>% filter(Protein %in% zhe_down$target_id)

GA_annotations %>% filter(Protein %in% zhong_up$target_id)
GA_annotations %>% filter(Protein %in% zhong_down$target_id)

GA_annotations %>% filter(Protein %in% cross_1_up$target_id)
GA_annotations %>% filter(Protein %in% cross_1_down$target_id)

GA_annotations %>% filter(Protein %in% cross_2_up$target_id)
GA_annotations %>% filter(Protein %in% cross_2_down$target_id)


### ABA-related transcripts



### TF SOM Spot GO Annotation
group_A <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots A.csv',
                      header = T,  sep = '\t', row.names = 1) %>% mutate(target_id = rownames(.))
group_B <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots B.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_C <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots C.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_D <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots D.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_E <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots E.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_F <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots F.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_G <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots G.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))
group_H <- read.table('/home/reverend_casy/Pea_transcriptomics/TFs_spots_201019/Group Overexpression Spots H.csv',
                      row.names = 1, header = T,  sep = '\t') %>% mutate(target_id = rownames(.))

sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_A, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_B, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_C, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_D, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_E, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_F, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_G, x, cutoff = 0.01, make_table = T))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(group_H, x, cutoff = 0.01, make_table = T))

sapply(ls()[grepl('((BP)|(CC)|(MF))_group_[ABCDEFGH]', ls())], function(x) x %>% get() %>% 
         write.table(., file.path('~/Pea_transcriptomics/TFs_spots_201019/', paste(x, 'tsv', sep = '.')),
                     col.names = T, row.names = F, sep = '\t'))




##### new scatters
lipids <- c('GO:0005811', 'GO:0055089', 'GO:0015908', 'GO:0006631', 'GO:0006633', 'GO:0030497', 'GO:0019915')
storage_proteins <- c('GO:0030497', 'GO:0034422', 'GO:0032578', 'GO:0045735', 'GO:0000326', 'GO:1990019')
### GO:0043556 - regulation of translation in response to oxidative stress; GO:0032939 - positive, GO:0032938 - negative
### GO:1990497 - regulation to cytoplasmic translation in response to stress
### GO:0032057 - negative modulation of translational initiation in response to stress
translation <- c('GO:0006412', 'GO:0002181', 'GO:0002188', 'GO:0032543', 'GO:0032544', 'GO:0030371', 'GO:0006417',
                 'GO:0070992', 'GO:0070993', 'GO:0008494', 'GO:0045182', 'GO:1903502', 'GO:0070126', 'GO:0070124',
                 'GO:0070125', 'GO:0017148', 'GO:0018444', 'GO:0031369', 'GO:0070129', 'GO:0045727', 'GO:0008079',
                 'GO:2000765', 'GO:0003743', 'GO:0003746', 'GO:0003747', 'GO:0097167', 'GO:0061770', 'GO:0044207',
                 'GO:0006418', 'GO:0045974', 'GO:0001731', 'GO:0070130', 'GO:0070131', 'GO:0008135', 'GO:0035278',
                 'GO:2000767', 'GO:2000766', 'GO:0046011', 'GO:1905082', 'GO:0006414', 'GO:0070132', 'GO:0001677',
                 'GO:0040033', 'GO:0005851', 'GO:0005850', 'GO:0005853', 'GO:0005852', 'GO:0090079', 'GO:0016281',
                 'GO:0016150', 'GO:0016149', 'GO:0045975', 'GO:0001732', 'GO:0070127', 'GO:0043143', 'GO:0043614',
                 'GO:0070133', 'GO:0070134', 'GO:0000900', 'GO:0043555', 'GO:0039606', 'GO:1905173', 'GO:0097010',
                 'GO:1905143', 'GO:1905616', 'GO:0004694', 'GO:0045183', 'GO:0070196', 'GO:1901193', 'GO:0000901',
                 'GO:0071541', 'GO:0071540', 'GO:0006413', 'GO:1901190', 'GO:1990500', 'GO:0043556', 'GO:0043557',
                 'GO:0070549', 'GO:1905618', 'GO:1905617', 'GO:0070321', 'GO:1901194', 'GO:1901195', 'GO:1990497',
                 'GO:0032055', 'GO:0032056', 'GO:0071891', 'GO:0043558', 'GO:0080149', 'GO:0006415', 'GO:0036490',
                 'GO:0032939', 'GO:0032938', 'GO:1905535', 'GO:0070322', 'GO:0070323', 'GO:1901192', 'GO:1901191',
                 'GO:1904803', 'GO:0032061', 'GO:0032062', 'GO:0010892', 'GO:0036493', 'GO:0036491', 'GO:1905537',
                 'GO:1905536', 'GO:1902010', 'GO:0036494', 'GO:0036495', 'GO:0006451', 'GO:0006452', 'GO:0009386',
                 'GO:0036497', 'GO:0036492', 'GO:0002182', 'GO:0002183', 'GO:0002184', 'GO:0002190', 'GO:0002191',
                 'GO:1990145', 'GO:0075522', 'GO:0006446', 'GO:0006450', 'GO:0006448', 'GO:0006449', 'GO:2001124',
                 'GO:0032057', 'GO:0045948', 'GO:0045947', 'GO:0045904', 'GO:0045905', 'GO:0045902', 'GO:0045903',
                 'GO:0045900', 'GO:0045901', 'GO:1903674', 'GO:1903677', 'GO:1904688', 'GO:1990580', 'GO:1900247',
                 'GO:2001126', 'GO:2001125', 'GO:0140018', 'GO:0032063', 'GO:1903916', 'GO:0032058', 'GO:1905084',
                 'GO:1905083', 'GO:0002192', 'GO:0006447', 'GO:0106074', 'GO:0110019', 'GO:0110017', 'GO:0110018',
                 'GO:1903675', 'GO:1903678', 'GO:1903679', 'GO:1903676', 'GO:1904690', 'GO:1904689', 'GO:1900249',
                 'GO:1900248', 'GO:0097622', 'GO:0071263', 'GO:0032064', 'GO:0045993', 'GO:0045994', 'GO:0036499',
                 'GO:0071264', 'GO:1903270', 'GO:0036496', 'GO:0071262', 'GO:0010998', 'GO:0043561', 'GO:1903272',
                 'GO:1903271', 'GO:1990611', 'GO:1990625', 'GO:1903917', 'GO:1903912', 'GO:0007571', 'GO:0010998')

repair <- c('GO:0097551', 'GO:0036297', 'GO:000630', 'GO:0006281', 'GO:0006298', 'GO:0000719', 'GO:0000725',
                'GO:0006284', 'GO:0006289', 'GO:0006290', 'GO:0043504', 'GO:0001778', 'GO:1990391', 'GO:1990516',
                'GO:0010206', 'GO:0009380', 'GO:0000710', 'GO:0000726', 'GO:0032300', 'GO:1990710', 'GO:0010213',
                'GO:0006302', 'GO:0006282', 'GO:0090735', 'GO:0070914', 'GO:0000012', 'GO:0000109', 'GO:0000711',
                'GO:0032423', 'GO:0032404', 'GO:0051103', 'GO:0097552', 'GO:0006283', 'GO:0006288', 'GO:0006287',
                'GO:0033683', 'GO:0042275', 'GO:1905051', 'GO:0070911', 'GO:1903823', 'GO:1905684', 'GO:0045739',
                'GO:0045738', 'GO:2000819', 'GO:0097520', 'GO:0036298', 'GO:0032425', 'GO:0032424', 'GO:0006294',
                'GO:0003684', 'GO:1990414', 'GO:0000716', 'GO:0000715', 'GO:0006307', 'GO:0006285', 'GO:0006293',
                'GO:0006297', 'GO:1905053', 'GO:1905052', 'GO:0000111', 'GO:0000110', 'GO:0000113', 'GO:0000112',
                'GO:1905686', 'GO:1905685', 'GO:1903516', 'GO:1902113', 'GO:2000779', 'GO:1990598', 'GO:0097698',
                'GO:1990731', 'GO:0036299', 'GO:0000731', 'GO:0000718', 'GO:0000717', 'GO:0140274', 'GO:1990249',
                'GO:0098783', 'GO:0090262', 'GO:0098504', 'GO:1903517', 'GO:1903518', 'GO:1990396', 'GO:2000780',
                'GO:2000781', 'GO:0000724', 'GO:0000720', 'GO:0140273', 'GO:0006295', 'GO:0006296', 'GO:0006303',
                'GO:0006286', 'GO:1904161', 'GO:0090699', 'GO:1903824', 'GO:0045002', 'GO:1990918', 'GO:0000727',
                'GO:0010776', 'GO:0010777', 'GO:1990250', 'GO:1903110', 'GO:1901255', 'GO:0045003', 'GO:0097681',
                'GO:0097680', 'GO:0061306', 'GO:0000734', 'GO:0010779', 'GO:0010778', 'GO:0010569', 'GO:0043150',
                'GO:1904162', 'GO:1903111', 'GO:1903112', 'GO:1905168', 'GO:1901591', 'GO:2000042', 'GO:2001032',
                'GO:0097510', 'GO:0043151', 'GO:0042276', 'GO:1901291', 'GO:1901592', 'GO:2001033', 'GO:2001034',
                'GO:1902346', 'GO:0043765', 'GO:0043739', 'GO:0070716', 'GO:0100026', 'GO:0061674', 'GO:0000736',
                'GO:0010792')
###TRINITY_10_DN2672_c0_g1_i3 ubiquitin-conjugating enzyme E2 variant 1D [Cucurbita moschata] - how in the world could it be conncted to repair?
starch <- c('GO:0043036', 'GO:2001070', 'GO:0005982', 'GO:0062052', 'GO:0005983', 'GO:0009011', 'GO:0009569',
            'GO:0009568', 'GO:0019252', 'GO:0044654', 'GO:0004373', 'GO:0102218', 'GO:0044570', 'GO:2000904',
            'GO:0102502', 'GO:0033840', 'GO:2000881', 'GO:0010581', 'GO:0044574', 'GO:2000905', 'GO:2000906',
            'GO:1904160', 'GO:2000465', 'GO:2000882', 'GO:2000883', 'GO:0102222', 'GO:2000467', 'GO:2000466',
            'GO:1900512', 'GO:1900513', 'GO:1900514', 'GO:0035956', 'GO:0050521', 'GO:0003844')

desiccation <- c('GO:0009269', 'GO:0071465', 'GO:0097439', 'GO:0048700', 'GO:0072515')

postmod <- c('GO:0043687', 'GO:1901873', 'GO:0036211', 'GO:0042160', 'GO:0031179', 'GO:0030047', 'GO:0018208', 'GO:0018207',
             'GO:0018209', 'GO:0018211', 'GO:0018210', 'GO:0018213', 'GO:0018212', 'GO:0018202', 'GO:0018201', 'GO:0018204',
             'GO:0018203', 'GO:0018206', 'GO:0018205', 'GO:0018194', 'GO:0018196', 'GO:0018195', 'GO:0018198', 'GO:0018199',
             'GO:0050844', 'GO:1901875', 'GO:1901874', 'GO:0072580', 'GO:0043686', 'GO:0006464', 'GO:0018200', 'GO:0018193',
             'GO:0018197', 'GO:0071587', 'GO:0140030', 'GO:0035610', 'GO:0036210', 'GO:0034421', 'GO:0006497', 'GO:0018410',
             'GO:0031400', 'GO:0031401', 'GO:0031365', 'GO:0032446', 'GO:0140035', 'GO:0070647', 'GO:0098823', 'GO:1903059',
             'GO:1903320', 'GO:1903322', 'GO:1903321', 'GO:1903061', 'GO:1903060', 'GO:0031289', 'GO:0006468', 'GO:0007258',
             'GO:0035404', 'GO:0035405', 'GO:0035406', 'GO:1990164', 'GO:0018109', 'GO:0018108', 'GO:0018105', 'GO:0018107',
             'GO:0018106', 'GO:0018218', 'GO:1990245', 'GO:0043989', 'GO:0043988', 'GO:0043987', 'GO:0043990', 'GO:0043991',
             'GO:0043538', 'GO:0035408', 'GO:0035409', 'GO:0035407', 'GO:0035978', 'GO:0033127', 'GO:0018217', 'GO:0072370',
             'GO:0072355', 'GO:0001932', 'GO:0022828', 'GO:0140031', 'GO:0042501', 'GO:0050730', 'GO:0033135', 'GO:0033129',
             'GO:0033128', 'GO:0001933', 'GO:0001934', 'GO:0023015', 'GO:0023014', 'GO:0023016', 'GO:1990853', 'GO:0010799',
             'GO:0002030', 'GO:0050731', 'GO:0050732', 'GO:0033137', 'GO:0033138', 'GO:2000281', 'GO:0010801', 'GO:0010800',
             'GO:2000775', 'GO:0036492', 'GO:0010998', 'GO:1904324', 'GO:1904325', 'GO:0100002', 'GO:0060734', 'GO:0060733',
             'GO:1903125', 'GO:1903912', 'GO:1901408', 'GO:1901409', 'GO:0071619', 'GO:0071620', 'GO:0044387', 'GO:2000751',
             'GO:2001163', 'GO:2000817', 'GO:2001164', 'GO:2001165', 'GO:1903654', 'GO:1903655', 'GO:1900018', 'GO:0019901',
             'GO:0004672', 'GO:1902911', 'GO:0032147', 'GO:0019887', 'GO:0006486', 'GO:0018280', 'GO:0033576', 'GO:0033578',
             'GO:0042076', 'GO:0018103', 'GO:0006487', 'GO:0006493', 'GO:0060049', 'GO:0140032', 'GO:0033577', 'GO:0033575',
             'GO:0060050', 'GO:0060051', 'GO:0042543', 'GO:0018317', 'GO:0018279', 'GO:0018258', 'GO:0018240', 'GO:0018242',
             'GO:0018241', 'GO:0018244', 'GO:0018243', 'GO:0018245', 'GO:0042077', 'GO:0090283', 'GO:1904098', 'GO:1904100',
             'GO:0090285', 'GO:0090284', 'GO:1904099', 'GO:0035629', 'GO:0018406', 'GO:0032147', 'GO:0043549', 'GO:1990782')

degradation <- c('GO:0016574', 'GO:0016567', 'GO:0007014', 'GO:0033522', 'GO:0033523', 'GO:0140372', 'GO:0070534', 'GO:0085020', 'GO:0070979',
                 'GO:0031396', 'GO:0070936', 'GO:0033182', 'GO:0035519', 'GO:1990390', 'GO:0036351', 'GO:0036352', 'GO:0140373', 'GO:0044314',
                 'GO:0070535', 'GO:0031398', 'GO:0031397', 'GO:0033184', 'GO:0033183', 'GO:2001166', 'GO:1900044', 'GO:0061945', 'GO:2001167',
                 'GO:2001168', 'GO:0140035', 'GO:0039648', 'GO:1901314', 'GO:1902523', 'GO:1902524', 'GO:1900045', 'GO:0061944', 'GO:0071894',
                 'GO:0075346', 'GO:1901316', 'GO:1901315', 'GO:1990756', 'GO:2001173', 'GO:2001174', 'GO:2001175', 'GO:0039542', 'GO:0051865',
                 'GO:1904666', 'GO:1902498', 'GO:0061630', 'GO:2000060', 'GO:2000059', 'GO:0006511', 'GO:2000058', 'GO:1904668', 'GO:1905524',
                 'GO:1904667', 'GO:1902499', 'GO:0043248', 'GO:0070628', 'GO:0031144', 'GO:0000502', 'GO:0031597', 'GO:0031595', 'GO:0005839',
                 'GO:0005838', 'GO:1903009', 'GO:0008537', 'GO:0022624', 'GO:0034515', 'GO:0080129', 'GO:0031603', 'GO:0031600', 'GO:0031601',
                 'GO:0031598', 'GO:0070682', 'GO:0090364', 'GO:0036402', 'GO:1902906', 'GO:1902907', 'GO:1990919', 'GO:1904855', 'GO:1904854',
                 'GO:0022623', 'GO:0008540', 'GO:0008541', 'GO:1990920', 'GO:0031615', 'GO:0031613', 'GO:0031612', 'GO:0031610', 'GO:1904327',
                 'GO:1990236', 'GO:0019773', 'GO:0019774', 'GO:0090363', 'GO:1902885', 'GO:0031604', 'GO:0031609', 'GO:0031606', 'GO:0031607',
                 'GO:0043161', 'GO:1902886', 'GO:1902887', 'GO:0039588', 'GO:0004175', 'GO:1990237', 'GO:0071630', 'GO:0071629', 'GO:0010498',
                 'GO:1904379', 'GO:0036369', 'GO:0004298', 'GO:1901483', 'GO:0006515', 'GO:1901799', 'GO:1901800', 'GO:1901484', 'GO:1901485',
                 'GO:0031593', 'GO:0071795', 'GO:0071796', 'GO:0070530', 'GO:0036435', 'GO:0071795', 'GO:0071796')

replication <- c('GO:0005657', 'GO:0006260', 'GO:0046809', 'GO:0006269', 'GO:0043596', 'GO:0006264', 'GO:0006270', 'GO:0006274', 'GO:0031297',
                 'GO:0000076', 'GO:0033260', 'GO:0033259', 'GO:0043111', 'GO:0098689', 'GO:0045004', 'GO:1990078', 'GO:0048478', 'GO:1902969',
                 'GO:0019045', 'GO:0071932', 'GO:0044786', 'GO:0031634', 'GO:0006261', 'GO:0003688', 'GO:0006275', 'GO:0031261', 'GO:0031298',
                 'GO:0097047', 'GO:0033314', 'GO:1902975', 'GO:1902979', 'GO:1902317', 'GO:0005662', 'GO:0005663', 'GO:0006336', 'GO:0006335',
                 'GO:0033567', 'GO:0090296', 'GO:0008156', 'GO:0045740', 'GO:2000621', 'GO:1990506', 'GO:1902294', 'GO:1902292', 'GO:1902973',
                 'GO:1902333', 'GO:0071170', 'GO:0071163', 'GO:0071946', 'GO:0036387', 'GO:0043598', 'GO:0043599', 'GO:0031582', 'GO:0005664',
                 'GO:0043599', 'GO:0006268', 'GO:0070517', 'GO:0090329', 'GO:0090592', 'GO:0051104', 'GO:1903459', 'GO:1903460', 'GO:1903466',
                 'GO:1903468', 'GO:0090001', 'GO:0097046', 'GO:0090298', 'GO:0090297', 'GO:0072438', 'GO:0043110', 'GO:0043137', 'GO:1902595',
                 'GO:1990505', 'GO:1990426', 'GO:1902297', 'GO:1902291', 'GO:1902971', 'GO:1902977', 'GO:1902315', 'GO:1902320', 'GO:0000809',
                 'GO:0032201', 'GO:0101017', 'GO:1990414', 'GO:0000808', 'GO:0030174', 'GO:0006271', 'GO:1903463', 'GO:1903467', 'GO:1903211',
                 'GO:0110035', 'GO:0033262', 'GO:0072437', 'GO:0072444', 'GO:0072441', 'GO:0098673', 'GO:0011000', 'GO:0045005', 'GO:2000105',
                 'GO:2000104', 'GO:1902596', 'GO:1902597', 'GO:1902298', 'GO:1902990', 'GO:1990943', 'GO:0101018', 'GO:0036388', 'GO:1903461',
                 'GO:1903465', 'GO:1903464', 'GO:0110025', 'GO:0072436', 'GO:0072443', 'GO:1902576', 'GO:1902681', 'GO:1902982', 'GO:1902983',
                 'GO:1904860', 'GO:0071171', 'GO:0010571', 'GO:0071807', 'GO:0032213', 'GO:0032297', 'GO:0032298', 'GO:0030894', 'GO:1903469',
                 'GO:1903221', 'GO:0072442', 'GO:1902296', 'GO:1902981', 'GO:0032214', 'GO:0032215', 'GO:0110026', 'GO:1902295', 'GO:1902299',
                 'GO:1902319', 'GO:0006267', 'GO:0110027', 'GO:1902318', 'GO:1990099', 'GO:1902561', 'GO:1902985')

folding <- c('GO:0006457', 'GO:0044183', 'GO:0006458', 'GO:0061077', 'GO:1903332', 'GO:1990727', 'GO:0051086', 'GO:0051084', 'GO:0051083',
             'GO:1903333', 'GO:1903334', 'GO:0007023', 'GO:0022417', 'GO:0034975', 'GO:1903644', 'GO:1990507', 'GO:0061992', 'GO:0060904',
             'GO:1903646', 'GO:1903645', 'GO:0051085', 'GO:0051087', 'GO:0101031', 'GO:0016531', 'GO:0033254', 'GO:1990565', 'GO:0034663',
             'GO:0016532', 'GO:0051131', 'GO:1902694', 'GO:1990507', 'GO:0061992', 'GO:0090034', 'GO:0072323', 'GO:1903645', 'GO:0090035',
             'GO:0051082')

chromatin <- c('GO:0000785', 'GO:0031497', 'GO:0031498', 'GO:0006325', 'GO:0006342', 'GO:0006338', 'GO:0030874', 'GO:0003682', 'GO:0070827',
               'GO:0001739', 'GO:0000789', 'GO:0000790', 'GO:0016569', 'GO:0031933', 'GO:0031490', 'GO:0005677', 'GO:0035327', 'GO:0035328',
               'GO:1990700', 'GO:0061793', 'GO:0006333', 'GO:0005521', 'GO:0061641', 'GO:0031935', 'GO:0090223', 'GO:0000183', 'GO:0099114',
               'GO:0006343', 'GO:0006344', 'GO:0006348', 'GO:0006346', 'GO:0030527', 'GO:1990280', 'GO:0031055', 'GO:0043035', 'GO:0043044',
               'GO:0035561', 'GO:0030702', 'GO:0061638', 'GO:1902275', 'GO:1990841', 'GO:0010848', 'GO:0010847', 'GO:0071168', 'GO:0090579',
               'GO:0031011', 'GO:0035101', 'GO:1905269', 'GO:1905268', 'GO:0048096', 'GO:0031936', 'GO:0031937', 'GO:0031048', 'GO:0035562',
               'GO:0035563', 'GO:1904834', 'GO:1990141', 'GO:0008623', 'GO:0033186', 'GO:0001672', 'GO:0090308', 'GO:0031938', 'GO:0090052',
               'GO:0061187', 'GO:1905634', 'GO:0070869', 'GO:0070868', 'GO:0070870', 'GO:0098578', 'GO:0035392', 'GO:0035390', 'GO:0071169',
               'GO:0097240', 'GO:0090309', 'GO:0031940', 'GO:0031939', 'GO:0090053', 'GO:0090310', 'GO:0061188', 'GO:0043156', 'GO:0045798',
               'GO:0045799', 'GO:2000749', 'GO:1904499', 'GO:0061644', 'GO:0120186', 'GO:0120187', 'GO:0034401', 'GO:0010964', 'GO:0060906',
               'GO:0001301', 'GO:1904500', 'GO:1904501', 'GO:0097549', 'GO:0016590', 'GO:0090115', 'GO:0070924', 'GO:0035391', 'GO:0001305',
               'GO:0001304', 'GO:0070923', 'GO:0070919', 'GO:0001308', 'GO:0070921', 'GO:0032116', 'GO:1904497', 'GO:1902368', 'GO:1902368',
               'GO:1902360', 'GO:1905635', 'GO:0071959', 'GO:0099404', 'GO:1905309', 'GO:0071922', 'GO:1905644', 'GO:2000715', 'GO:0000408',
               'GO:0071923', 'GO:2000716', 'GO:2000717', 'GO:0099403', 'GO:0071960', 'GO:0070603', 'GO:2000718', 'GO:1904907', 'GO:2000719',
               'GO:2000720', 'GO:1905646', 'GO:1905645', 'GO:1904908', 'GO:1904909', 'GO:0036414', 'GO:0016570', 'GO:0016571', 'GO:0016572',
               'GO:0016573', 'GO:0016574', 'GO:0016575', 'GO:0016576', 'GO:0016577', 'GO:0016578', 'GO:0043486', 'GO:0042393', 'GO:0106077',
               'GO:0001207', 'GO:0010390', 'GO:0071110', 'GO:0042054', 'GO:0035404', 'GO:0035405', 'GO:0035406', 'GO:0035035', 'GO:0000123', 
               'GO:0000412', 'GO:0061649', 'GO:0140372', 'GO:0043967', 'GO:0043966', 'GO:0043969', 'GO:0043968', 'GO:0070933', 'GO:0070932',
               'GO:0031493', 'GO:0004402', 'GO:0004407', 'GO:0000118', 'GO:1990258', 'GO:1990226', 'GO:1990164', 'GO:0033522', 'GO:0033523',
               'GO:0035064', 'GO:0106229', 'GO:0035363', 'GO:0106153', 'GO:0035097', 'GO:0035518', 'GO:0035521', 'GO:0070076', 'GO:0070077',
               'GO:0035173', 'GO:0106078', 'GO:0042826', 'GO:0061922', 'GO:0036205', 'GO:0032452', 'GO:0140069', 'GO:0140068', 'GO:0034969',
               'GO:0034968', 'GO:0006334', 'GO:0070577', 'GO:0018110', 'GO:0035400', 'GO:0035034', 'GO:0035184', 'GO:0035174', 'GO:0099077',
               'GO:0043989', 'GO:0043988', 'GO:0043981', 'GO:0043980', 'GO:0043983', 'GO:0043982', 'GO:0043985', 'GO:0043984', 'GO:0043987',
               'GO:0043990', 'GO:0043991', 'GO:0043998', 'GO:0043978', 'GO:0043977', 'GO:0043979', 'GO:0043970', 'GO:0043972', 'GO:0043971',
               'GO:0043974', 'GO:0043973', 'GO:0043976', 'GO:0043975', 'GO:0036123', 'GO:0036124', 'GO:0097043', 'GO:0070775', 'GO:1990245',
               'GO:1990259', 'GO:0080182', 'GO:0070734', 'GO:0070544', 'GO:0044648', 'GO:0031060', 'GO:0031063', 'GO:0031056', 'GO:0031059',
               'GO:0033182', 'GO:0033169', 'GO:0106185', 'GO:0035408', 'GO:0035409', 'GO:0035407', 'GO:0008334', 'GO:0035978', 'GO:0033100',
               'GO:0033127', 'GO:0098532', 'GO:0035033', 'GO:0035065', 'GO:0072370', 'GO:0072355', 'GO:0051567', 'GO:0051568', 'GO:1903584',
               'GO:0035522', 'GO:0035574', 'GO:0070078', 'GO:0070079', 'GO:0035267', 'GO:1902562', 'GO:0036410', 'GO:1990619', 'GO:1990596',
               'GO:0062072', 'GO:1990467', 'GO:1990468', 'GO:0061647', 'GO:0097692', 'GO:0071044', 'GO:0097676', 'GO:0061628', 'GO:1900049',
               'GO:0036413', 'GO:1990889', 'GO:0097198', 'GO:0036353', 'GO:0036351', 'GO:0036352', 'GO:1990678', 'GO:1990679', 'GO:0044013',
               'GO:0034729', 'GO:0034720', 'GO:0034773', 'GO:0034772', 'GO:0034771', 'GO:0034770', 'GO:0097725', 'GO:0044154', 'GO:0010452',
               'GO:0010484', 'GO:0010485', 'GO:0071572', 'GO:0071557', 'GO:0140373', 'GO:0034972', 'GO:0034971', 'GO:0034970', 'GO:0046811',
               'GO:0006398', 'GO:1990104', 'GO:0090239', 'GO:0043189', 'GO:0070776', 'GO:1901725', 'GO:0070537', 'GO:0070535', 'GO:0017136',
               'GO:0018024', 'GO:0031061', 'GO:0031064', 'GO:0031065', 'GO:0031062', 'GO:0031057', 'GO:0031058', 'GO:0033184', 'GO:0033183',
               'GO:0008469', 'GO:0045129', 'GO:0033129', 'GO:0033128', 'GO:0035066', 'GO:0035067', 'GO:1903586', 'GO:1903585', 'GO:0032777',
               'GO:0001208', 'GO:0036409', 'GO:0062060', 'GO:1990483', 'GO:0071045', 'GO:0071208', 'GO:0036206', 'GO:1990853', 'GO:2001166',
               'GO:1900050', 'GO:1900051', 'GO:0035979', 'GO:0035093', 'GO:0071204', 'GO:0043999', 'GO:0043992', 'GO:0043994', 'GO:0043993',
               'GO:0043996', 'GO:0043995', 'GO:0043997', 'GO:0061085', 'GO:0090241', 'GO:0090240', 'GO:0042800', 'GO:1990244', 'GO:1905471',
               'GO:1905435', 'GO:1990162', 'GO:1901727', 'GO:1901726', 'GO:0033749', 'GO:0033746', 'GO:0031151', 'GO:0070510', 'GO:0070611',
               'GO:0070612', 'GO:0031078', 'GO:0031066', 'GO:0032931', 'GO:0035401', 'GO:0035402', 'GO:0035403', 'GO:2000281', 'GO:0072354',
               'GO:0072371', 'GO:0051569', 'GO:0051570', 'GO:1901674', 'GO:0035575', 'GO:0035642', 'GO:0035175', 'GO:0042799', 'GO:0000414',
               'GO:0036408', 'GO:1902028', 'GO:0062122', 'GO:0062073', 'GO:1900112', 'GO:1900109', 'GO:0071207', 'GO:2001253', 'GO:1902464',
               'GO:0051864', 'GO:0036207', 'GO:0036208', 'GO:2000618', 'GO:2000615', 'GO:0044012', 'GO:0044014', 'GO:0044016', 'GO:0044015',
               'GO:0044018', 'GO:0044017', 'GO:0044019', 'GO:0044020', 'GO:0044023', 'GO:0044022', 'GO:0044025', 'GO:0044024', 'GO:0034739',
               'GO:0032129', 'GO:2001160', 'GO:2001167', 'GO:2001168', 'GO:0071440', 'GO:0071558', 'GO:0032454', 'GO:0032453', 'GO:0046974',
               'GO:0046975', 'GO:0046972', 'GO:0046976', 'GO:0034649', 'GO:0034648', 'GO:0034647', 'GO:0034721', 'GO:0061086', 'GO:0061087',
               'GO:1905473', 'GO:1905472', 'GO:1905437', 'GO:1905436', 'GO:0070511', 'GO:0070512', 'GO:0031068', 'GO:0031067', 'GO:1901314',
               'GO:0051574', 'GO:0051572', 'GO:0051573', 'GO:0051571', 'GO:1901675', 'GO:1901676', 'GO:0035616', 'GO:0000416', 'GO:0000415',
               'GO:2000620', 'GO:1902030', 'GO:1902029', 'GO:2001254', 'GO:2001255', 'GO:1900111', 'GO:1900110', 'GO:1900113', 'GO:1900114',
               'GO:0061866', 'GO:1902465', 'GO:1902466', 'GO:2000617', 'GO:2000619', 'GO:2000616', 'GO:1902649', 'GO:2001162', 'GO:2001161',
               'GO:0071442', 'GO:0071441', 'GO:0071894', 'GO:0035042', 'GO:1904173', 'GO:1901316', 'GO:1901315', 'GO:2000775', 'GO:0097372',
               'GO:0032041', 'GO:1902651', 'GO:1902650', 'GO:0046969', 'GO:0046970', 'GO:0032221', 'GO:0097044', 'GO:1904175', 'GO:1904174',
               'GO:2001173', 'GO:0033698', 'GO:0000125', 'GO:2000776', 'GO:2001174', 'GO:2001175', 'GO:0016581', 'GO:1990121', 'GO:2000751',
               'GO:2000873', 'GO:2000817', 'GO:0034080', 'GO:1903097', 'GO:0097549', 'GO:0034401', 'GO:1903098', 'GO:1903099')

GA <- c('GO:0010336', 'GO:0010331', 'GO:0010476', 'GO:0045487', 'GO:0009739', 'GO:0009686', 'GO:0009685', 'GO:0102713', 'GO:1905201', 'GO:0045544',
        'GO:0047927', 'GO:0033470', 'GO:0033469', 'GO:0102714', 'GO:0102716', 'GO:0102715', 'GO:0102712', 'GO:0102711', 'GO:0016707', 'GO:0102653',
        'GO:0102652', 'GO:0102663', 'GO:0102924', 'GO:0051779', 'GO:0102111', 'GO:0102119', 'GO:0102118', 'GO:0102117', 'GO:0102122', 'GO:0045543',
        'GO:0047928', 'GO:0071370', 'GO:0010341', 'GO:0010371', 'GO:0102738', 'GO:0102739', 'GO:0102123', 'GO:0102123', 'GO:0102125', 'GO:0102124',
        'GO:0010372', 'GO:0010373', 'GO:0103057', 'GO:0103056', 'GO:0052634', 'GO:0052635', 'GO:0103010', 'GO:0102972', 'GO:0103054', 'GO:0102755')

ABA <- c('GO:0080168', 'GO:1902265', 'GO:0010427', 'GO:0009687', 'GO:1902266', 'GO:0010294', 'GO:0010293', 'GO:0046345', 'GO:0009737', 'GO:0009688',
         'GO:0090440', 'GO:0010295', 'GO:0009724', 'GO:0009738', 'GO:0071215', 'GO:0010115', 'GO:0090359', 'GO:0051993', 'GO:1902418', 'GO:0009787',
         'GO:0010116', 'GO:1902417', 'GO:0009789', 'GO:0009788', 'GO:0075343', 'GO:1901527', 'GO:1990218')

splicing <- c('GO:0008380', 'GO:0034247', 'GO:0000394', 'GO:0000374', 'GO:0000373', 'GO:0000372', 'GO:1990935', 'GO:0006388', 'GO:0070054', 'GO:0016539',
              'GO:0043484', 'GO:0072669', 'GO:0000366', 'GO:0000398', 'GO:0033120', 'GO:0033119', 'GO:0045292', 'GO:0045291', 'GO:0002564', 'GO:0000380',
              'GO:0000365', 'GO:0000375', 'GO:0030909', 'GO:0048024', 'GO:0019801', 'GO:0019802', 'GO:0070274', 'GO:0000214', 'GO:0048026', 'GO:0048025',
              'GO:1905744', 'GO:0035048', 'GO:0000381', 'GO:1905746', 'GO:1905745', 'GO:0002563', 'GO:0071030', 'GO:0000376', 'GO:0000376', 'GO:0016607',
              'GO:0000213')

cell_wall <- c('GO:0005618', 'GO:0042545', 'GO:0042546', 'GO:0070726', 'GO:0052386', 'GO:0009530', 'GO:0009531', 'GO:0071555', 'GO:0044277',
               'GO:0009505', 'GO:0010384', 'GO:0071668', 'GO:0044347', 'GO:0099613', 'GO:0031506', 'GO:0031504', 'GO:0000032', 'GO:0000196',
               'GO:0016998', 'GO:0052543', 'GO:0052546', 'GO:0070592', 'GO:0005199', 'GO:0010339', 'GO:0010383', 'GO:0052325', 'GO:0009664',
               'GO:0044036', 'GO:0044038', 'GO:0009832', 'GO:0009828', 'GO:0009827', 'GO:0071554', 'GO:0009830', 'GO:0009829', 'GO:0034406',
               'GO:1990394', 'GO:2000652', 'GO:0052482', 'GO:0009834', 'GO:0009833', 'GO:0034410', 'GO:0044568', 'GO:0044567', 'GO:2000966',
               'GO:0010404', 'GO:0044348', 'GO:1903338', 'GO:0052541', 'GO:1901347', 'GO:1901348', 'GO:1902066', 'GO:0052324', 'GO:0071669',
               'GO:0010981', 'GO:2000967', 'GO:2000968', 'GO:1903339', 'GO:1903340', 'GO:0042547', 'GO:0052544', 'GO:0070597', 'GO:0070598',
               'GO:0060870', 'GO:1990076', 'GO:1902088', 'GO:2000939', 'GO:0090379', 'GO:1905588', 'GO:0080157', 'GO:2001009', 'GO:2000940',
               'GO:2000941', 'GO:1902089', 'GO:0009831', 'GO:0009920')

tubulin <- c('GO:0090042', 'GO:0015631', 'GO:0045298', 'GO:0042903', 'GO:0043015', 'GO:0043014', 'GO:0048487', 'GO:0000930', 'GO:0007021',
             'GO:0071929', 'GO:0070738', 'GO:0090043', 'GO:0004835', 'GO:0033566', 'GO:0019799', 'GO:0008275', 'GO:0008274', 'GO:0070463',
             'GO:1990727', 'GO:0000931', 'GO:1902481', 'GO:0034741', 'GO:0090044', 'GO:0070740', 'GO:0150067', 'GO:1904428', 'GO:0000927',
             'GO:0000924', 'GO:0007023', 'GO:0150068', 'GO:0150069', 'GO:0008568', 'GO:0055032', 'GO:0000928', 'GO:0055031', 'GO:0055033',
             'GO:0061494', 'GO:0061495', 'GO:1990735', 'GO:0110121', 'GO:0110120', 'GO:0005874', 'GO:0005828', 'GO:0005827', 'GO:0051013',
             'GO:0051012', 'GO:0015630', 'GO:0055028', 'GO:0008017', 'GO:0005880', 'GO:0005881', 'GO:0005876', 'GO:0005879', 'GO:0097427',
             'GO:1990752', 'GO:0007019', 'GO:0000235', 'GO:0007020', 'GO:0046785', 'GO:0034453', 'GO:0001578', 'GO:0099070', 'GO:0099071',
             'GO:0099609', 'GO:0099111', 'GO:0043622', 'GO:0031121', 'GO:0031122', 'GO:0005815', 'GO:0030981', 'GO:0030954', 'GO:0030953',
             'GO:0080175', 'GO:0060404', 'GO:0035371', 'GO:1905720', 'GO:0060172', 'GO:0005875', 'GO:1990498', 'GO:0003777', 'GO:0036449',
             'GO:0000226', 'GO:1990644', 'GO:0007018', 'GO:0007018', 'GO:0010938', 'GO:0010970', 'GO:0032932', 'GO:0099098', 'GO:0090222',
             'GO:0090064', 'GO:0099118', 'GO:0030473', 'GO:1990295', 'GO:0051011', 'GO:0051010', 'GO:0031113', 'GO:0031114', 'GO:0031109',
             'GO:0031023', 'GO:0031021', 'GO:0035372', 'GO:0047496', 'GO:1905759', 'GO:0098840', 'GO:1903754', 'GO:0106006', 'GO:1990537',
             'GO:1904511', 'GO:1904526', 'GO:0061673', 'GO:1990755', 'GO:0000923', 'GO:0061842', 'GO:0061863', 'GO:1990941', 'GO:0032117',
             'GO:0032118', 'GO:0034454', 'GO:0010968', 'GO:0036078', 'GO:0090221', 'GO:0031887', 'GO:0072384', 'GO:0036250', 'GO:0090221', 
             'GO:0090063', 'GO:0031534', 'GO:0031535', 'GO:0099112', 'GO:1904186', 'GO:1904185', 'GO:0031116', 'GO:0031117', 'GO:0031115',
             'GO:0070507', 'GO:0032886', 'GO:0031024', 'GO:0031025', 'GO:0070462', 'GO:1905833', 'GO:1905721', 'GO:1905725', 'GO:1905755',
             'GO:0072698', 'GO:0060632', 'GO:0072699', 'GO:1902513', 'GO:2000576', 'GO:2000575', 'GO:1904519', 'GO:0010005', 'GO:1902838',
             'GO:1902850', 'GO:1902816', 'GO:1990933', 'GO:1904825', 'GO:0007027', 'GO:0034643', 'GO:0034631', 'GO:1902619', 'GO:0140274',
             'GO:0034994', 'GO:0140210', 'GO:1905121', 'GO:0099606', 'GO:1903032', 'GO:1903033', 'GO:0099117', 'GO:1905509', 'GO:0030951',
             'GO:0031112', 'GO:0031111', 'GO:1903696', 'GO:0072385', 'GO:0072383', 'GO:0072382', 'GO:0072386', 'GO:0051415', 'GO:0075519',
             'GO:1901610', 'GO:1901609', 'GO:0098863', 'GO:1904518', 'GO:1902840', 'GO:1902839', 'GO:1902817', 'GO:1990810', 'GO:1904759',
             'GO:0034640', 'GO:1902620', 'GO:0140273', 'GO:1905185', 'GO:0090172', 'GO:0008569', 'GO:0008574', 'GO:0055031', 'GO:0055033',
             'GO:0106007', 'GO:1990734', 'GO:1990852', 'GO:0090176', 'GO:0090226', 'GO:0099110', 'GO:1903562', 'GO:1903563', 'GO:1990976',
             'GO:2000580', 'GO:2000577', 'GO:0140024', 'GO:0007099', 'GO:2000581', 'GO:2000582', 'GO:2000578', 'GO:2000579', 'GO:0008608',
             'GO:0051315', 'GO:0099607', 'GO:0051987', 'GO:0051986', 'GO:1902423', 'GO:1905115', 'GO:1902424', 'GO:1902425', 'GO:1905116',
             'GO:0051455', 'GO:1904967', 'GO:1904968', 'GO:0047497')

### GO:0010970 - transport along microtubule ; GO:0099118 - microtubule-based protein transport
motility <- c('GO:0030286', 'GO:0070840', 'GO:0005868', 'GO:0045504', 'GO:0045505', 'GO:0045503', 'GO:0051959', 'GO:0003777', 'GO:2000574',
              'GO:2000576', 'GO:2000575', 'GO:0019894', 'GO:0005871', 'GO:0016938', 'GO:0016939', 'GO:0005873', 'GO:0005872', 'GO:2000577',
              'GO:2000578', 'GO:2000579', 'GO:2000580', 'GO:2000581', 'GO:2000582', 'GO:2000576', 'GO:2000575', 'GO:0008569', 'GO:0008574')

pea_scatter_plotter <- function(MapMan_df_loc, uprmax = NA,
                                plot_title = NA, x_breaks = c(), y_breaks = c()){

  if (!is.na(uprmax)){
    MapMan_df_loc %<>% filter_at(vars(c(colnames(MapMan_df_loc[,unlist(lapply(MapMan_df_loc, is.numeric))]))), all_vars(.<=uprmax))
  }
  
  ggplot(MapMan_df_loc, aes(x = value, y = value, size = I(2.8))) + 
    geom_point(aes(x = Sprint2_10_mean, y = Sprint2_30_mean), alpha = 0.4, col = '#068013') +
    geom_point(aes(x = Zhewan1_10_mean, y = Zhewan1_25_mean), alpha = 0.4, col = '#f1cb20') +
    geom_point(aes(x = Zhongwan6_10_mean, y = Zhongwan6_25_mean), alpha = 0.4, col = '#b30a07') +
    xlab(label = 'Earlier condition') + ylab('Latter condition')+
    guides(alpha = 'none')+
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 20),
                       legend.text = element_text(size = 16), axis.title = element_text(size=14),
                       axis.text = element_text(size=12, color = 'black')) +
    scale_x_continuous(labels = scales::comma) + scale_y_continuous(labels = scales::comma) + scale_color_manual(name = 'cultivar',
          values = c('Sprint-2'='#068013', 'Zhewan-1'='#f1cb20', 'Zhongwan-6'='#b30a07'))
  
}


pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(storage_proteins, collapse = ')|('), ')', sep = ''),
                annot_full$GOs) | grepl('storage proteins', annot_full$MapMan_terms),]$transcript_id,], uprmax = 50000)
MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(storage_proteins, collapse = ')|('), ')', sep = ''),
                        annot_full$GOs) | grepl('storage proteins', annot_full$MapMan_terms),]$transcript_id,] %>% filter(Zhewan1_10_mean > 200000)
annot_full %>% filter(transcript_id == 'TRINITY_10_DN3139_c0_g1_i3') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names)


pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(lipids, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,], uprmax = 2000)
MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(lipids, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Zhewan1_10_mean > 2000)
annot_full %>% filter(transcript_id == 'TRINITY_30_DN977_c0_g1_i2') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names)

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(translation, collapse = ')|('), ')', sep = ''),
                    annot_full$GOs) | grepl('protein.synthesis', annot_full$MapMan_terms),]$transcript_id,], uprmax = 1600)


MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(translation, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs) | grepl('protein.synthesis', annot_full$MapMan_terms),]$transcript_id,] %>% filter(Sprint2_10_mean > 4000)


pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(repair, collapse = ')|('), ')', sep = ''),
                                                    annot_full$GOs) | grepl('DNA.repair', annot_full$MapMan_terms),]$transcript_id,])

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(translation, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Sprint2_10_mean > 1000)

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(repair, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Sprint2_10_mean > 250)
annot_full %>% filter(transcript_id == 'TRINITY_10_DN2672_c0_g1_i3') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names)

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(starch, collapse = ')|('), ')', sep = ''),
                                              annot_full$GOs) | grepl('starch', annot_full$MapMan_terms),]$transcript_id,], uprmax = 1000)
MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(starch, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Zhongwan6_25_mean > 2000)
annot_full %>% filter(transcript_id == 'TRINITY_10_DN4889_c0_g1_i3') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names)

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(desiccation, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(desiccation, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Zhewan1_25_mean > 400)
annot_full %>% filter(transcript_id == 'TRINITY_10_DN25_c0_g2_i1') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names)

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(postmod, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(degradation, collapse = ')|('), ')', sep = ''),
                    annot_full$GOs) | grepl('protein.degradation', annot_full$MapMan_terms),]$transcript_id,], uprmax = 1000)

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(replication, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(TF_GO, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs) | grepl('regulation of transcription', annot_full$MapMan_terms),]$transcript_id,])

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(folding, collapse = ')|('), ')', sep = ''),
                    annot_full$GOs)| grepl('folding', annot_full$MapMan_terms),]$transcript_id,], uprmax = 1000)

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(folding, collapse = ')|('), ')', sep = ''),
                      annot_full$GOs)| grepl('folding', annot_full$MapMan_terms),]$transcript_id,] %>% filter(Zhewan1_25_mean > 2000)


pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(chromatin, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(degradation, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Sprint2_30_mean > 4000)

annot_full %>% filter(transcript_id == 'TRINITY_10_DN1016_c0_g1_i1') %>% dplyr::select(BLASTP_names, BLASTX_names)


pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(GA, collapse = ')|('), ')', sep = ''),
                    annot_full$GOs) | grepl('gibberelin', annot_full$MapMan_terms),]$transcript_id,])

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(ABA, collapse = ')|('), ')', sep = ''),
                            annot_full$GOs) | grepl('abscisic acid', annot_full$MapMan_terms),]$transcript_id,] %>% filter(MapMan_Category2 != 'storage proteins'))

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(ABA, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs) | grepl('abscisic acid', annot_full$MapMan_terms),]$transcript_id,] %>% filter(MapMan_Category2 != 'storage proteins')

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(splicing, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(cell_wall, collapse = ')|('), ')', sep = ''),
                            annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms) | grepl('cell wall', annot_full$MapMan_terms),]$transcript_id,])

MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(cell_wall, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs) | grepl('cell wall', annot_full$MapMan_terms),]$transcript_id,] %>% filter(Zhewan1_25_mean > 90000)

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(tubulin, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

pea_scatter_plotter(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(motility, collapse = ')|('), ')', sep = ''),
                                                                            annot_full$GOs),]$transcript_id,])

annot_full %>% filter(transcript_id == 'TRINITY_30_DN244_c0_g1_i7') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)

cross_5_down %>% filter(target_id == 'TRINITY_10_DN4509_c0_g1_i2')

MapMan_df %>% filter(MapMan_Category1 == 'hormone metabolism') %>% dplyr::select(MapMan_Category2) %>% unique()
MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(starch, collapse = ')|('), ')', sep = ''),
                                                        annot_full$GOs),]$transcript_id,] %>% filter(Sprint2_30_mean > 1000)


annot_full %>% filter(transcript_id == 'TRINITY_10_DN2919_c0_g1_i1') %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms, plants_BLASTP)

MapMan_df %>% filter(transcript_id == 'TRINITY_10_DN47076_c0_g1_i1')
upregulated_30_days %>% filter(target_id == 'TRINITY_30_DN34986_c0_g1_i1')
zhong_down %>% filter(target_id == 'TRINITY_30_DN34986_c0_g1_i1')

models(so_pea)
cross_5_up %>% filter(target_id == 'TRINITY_10_DN1645_c0_g1_i10')

##### outlier detection
find_outliers <- function(df, x, y, additional_categories = NA, group, category) {
  cd <- cooks.distance(lm(get(y)~get(x), data = df))
  plot(cd,pch=1, cex = 1)
  abline(h = 4*mean(cd,na.rm = T), col = "blue")
  text(x=1:length(cd)+2, y=cd, labels=ifelse(cd>4*mean(cd, na.rm=T),names(cd),""), col="blue")
  outliers <- cd[cd>4*mean(cd, na.rm=T)] %>% names()
  df <- df[rownames(df) %in% outliers, c(x,y, additional_categories)]
  colnames(df)[c(1,2)] <- c('Early', 'Late')
  df$group <- rep(group, nrow(df)) %>% as.character()
  df$category <- rep(category, nrow(df)) %>% as.character()
  return(df)
}

### now track down the outliers  

# ABA
outlier_list <- find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(ABA, collapse = ')|('), ')', sep = ''),
      annot_full$GOs) | grepl('abscisic acid', annot_full$MapMan_terms),]$transcript_id,] %>% filter(MapMan_Category2 != 'storage proteins') %>% filter_at(vars(-c(transcript_id, MapMan_Category, MapMan_Category1, MapMan_Category2, MapMan_Category3, MapMan_Category4)), all_vars(.<=5000)), 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'ABA metabolism and signaling')

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(ABA, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('abscisic acid', annot_full$MapMan_terms),]$transcript_id,] %>% filter(MapMan_Category2 != 'storage proteins') %>% filter_at(vars(-c(transcript_id, MapMan_Category, MapMan_Category1, MapMan_Category2, MapMan_Category3, MapMan_Category4)), all_vars(.<=5000)), 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'ABA metabolism and signaling'))
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(ABA, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('abscisic acid', annot_full$MapMan_terms),]$transcript_id,] %>% filter(MapMan_Category2 != 'storage proteins') %>% filter_at(vars(-c(transcript_id, MapMan_Category, MapMan_Category1, MapMan_Category2, MapMan_Category3, MapMan_Category4)), all_vars(.<=5000)), 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'ABA metabolism and signaling'))
# chromatin

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(chromatin, collapse = ')|('), ')', sep = ''),
                annot_full$GOs),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Chromatin constituent and remodeling'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(chromatin, collapse = ')|('), ')', sep = ''),
        annot_full$GOs),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Chromatin constituent and remodeling'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(chromatin, collapse = ')|('), ')', sep = ''),
                                               annot_full$GOs),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'Chromatin constituent and remodeling'))

# degradation

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(degradation, collapse = ')|('), ')', sep = ''),
                annot_full$GOs) | grepl('protein.degradation', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Protein degradation'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(degradation, collapse = ')|('), ')', sep = ''),
      annot_full$GOs) | grepl('protein.degradation', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Protein degradation'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(degradation, collapse = ')|('), ')', sep = ''),
          annot_full$GOs) | grepl('protein.degradation', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'Protein degradation'))

# folding

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(folding, collapse = ')|('), ')', sep = ''),
                annot_full$GOs)| grepl('folding', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Protein folding'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(folding, collapse = ')|('), ')', sep = ''),
                            annot_full$GOs)| grepl('folding', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Protein folding'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(folding, collapse = ')|('), ')', sep = ''),
          annot_full$GOs)| grepl('folding', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'Protein folding'))

# GA
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(GA, collapse = ')|('), ')', sep = ''),
  annot_full$GOs) | grepl('gibberelin', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'GA metabolism and signaling'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(GA, collapse = ')|('), ')', sep = ''),
    annot_full$GOs) | grepl('gibberelin', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'GA metabolism and signaling'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(GA, collapse = ')|('), ')', sep = ''),
            annot_full$GOs) | grepl('gibberelin', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'GA metabolism and signaling'))

# lipids
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(lipids, collapse = ')|('), ')', sep = ''),
                        annot_full$GOs),]$transcript_id,],  'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Lipid storage'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(lipids, collapse = ')|('), ')', sep = ''),
             annot_full$GOs),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Lipid storage'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(lipids, collapse = ')|('), ')', sep = ''),
                annot_full$GOs),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'Lipid storage'))

# post-translational modifications
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(postmod, collapse = ')|('), ')', sep = ''),
      annot_full$GOs),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Post-translational modification'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(postmod, collapse = ')|('), ')', sep = ''),
            annot_full$GOs),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Post-translational modification'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(postmod, collapse = ')|('), ')', sep = ''),
            annot_full$GOs),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'Post-translational modification'))

# repair
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(repair, collapse = ')|('), ')', sep = ''),
  annot_full$GOs) | grepl('DNA.repair', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'DNA repair'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(repair, collapse = ')|('), ')', sep = ''),
      annot_full$GOs) | grepl('DNA.repair', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'DNA reair'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(repair, collapse = ')|('), ')', sep = ''),
      annot_full$GOs) | grepl('DNA.repair', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhongwan-6', 'DNA repair'))

# replication
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(replication, collapse = ')|('), ')', sep = ''),
      annot_full$GOs),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'DNA replication'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(replication, collapse = ')|('), ')', sep = ''),
 annot_full$GOs),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'DNA replication'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(replication, collapse = ')|('), ')', sep = ''),
        annot_full$GOs),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'DNA replication'))

# splicing
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(splicing, collapse = ')|('), ')', sep = ''),
    annot_full$GOs),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'RNA spicing'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(splicing, collapse = ')|('), ')', sep = ''),
            annot_full$GOs),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'RNA spicing'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(splicing, collapse = ')|('), ')', sep = ''),
  annot_full$GOs),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'RNA spicing'))

# starch 
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(starch, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('starch', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Starch metabolism'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(starch, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('starch', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Starch metabolism'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(starch, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('starch', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Starch metabolism'))

# storage proteins

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(storage_proteins, collapse = ')|('), ')', sep = ''),
                                                                                             annot_full$GOs) | grepl('storage proteins', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Protein storage'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(storage_proteins, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('storage proteins', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Protein storage'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(storage_proteins, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('storage proteins', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Protein storage'))


# transcription
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(TF_GO, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('regulation of transcription', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Transcription regulation'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(TF_GO, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('regulation of transcription', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Transcription regulation'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(TF_GO, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) | grepl('regulation of transcription', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Transcription regulation'))

# translation
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(translation, collapse = ')|('), ')', sep = ''),
  annot_full$GOs) | grepl('protein.synthesis', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Translation'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(translation, collapse = ')|('), ')', sep = ''),
        annot_full$GOs) | grepl('protein.synthesis', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Translation'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(translation, collapse = ')|('), ')', sep = ''),
    annot_full$GOs) | grepl('protein.synthesis', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Translation'))


# cell wall
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(cell_wall, collapse = ')|('), ')', sep = ''),
          annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms) | grepl('cell wall', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Cell wall'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(cell_wall, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms) | grepl('cell wall', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Cell wall'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(cell_wall, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms) | grepl('cell wall', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Cell wall'))

# tubulin
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(tubulin, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Tubulin'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(tubulin, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Tubulin'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(tubulim, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Tubulin'))

# tubulin-associated motility
outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(motility, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms),]$transcript_id,], 'Sprint2_10_mean', 'Sprint2_30_mean', c('MapMan_Category', 'transcript_id'), 'Sprint-2', 'Tubulin-associated motility'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(motility, collapse = ')|('), ')', sep = ''),
annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms),]$transcript_id,], 'Zhewan1_10_mean', 'Zhewan1_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-1', 'Tubulin-associated motility'))

outlier_list %<>% rbind(find_outliers(MapMan_df[MapMan_df$transcript_id %in% annot_full[grepl(paste0('(', paste(motility, collapse = ')|('), ')', sep = ''),
                                                                                              annot_full$GOs) & !grepl('storage', annot_full$MapMan_terms),]$transcript_id,], 'Zhongwan6_10_mean', 'Zhongwan6_25_mean', c('MapMan_Category', 'transcript_id'), 'Zhewan-6', 'Tubulin-associated motility'))

outlier_list$'BLASTP name' <- sapply(outlier_list$transcript_id, function(x) annot_full[annot_full$transcript_id == x,]$BLASTP_names)
outlier_list$'MapMan category' <- sapply(outlier_list$transcript_id, function(x) annot_full[annot_full$transcript_id == x,]$MapMan_terms)
outlier_list <- outlier_list[,c(4,5,6,1,2,3,7)]
outlier_list %>% View()
write.table(outlier_list, '~/Pea_transcriptomics/publication_pics/Table_S13_extended.tsv', col.names = T, row.names = F, sep = '\t')

##### Gotta plot some TFs!
library(OmicCircos)
Sprint2_all_tfs <- Sprint_tf_downregulated %>% rbind(Sprint_tf_upregulated)
Sprint2_all_tfs %<>% mutate(b = sapply(Sprint2_all_tfs$transcript_id, function(x) joint_results[joint_results$target_id == x,]$b))
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'AP2/EREBP, APETALA2/Ethylene-responsive element binding protein family',]$TF_class <- 'AP2/EREPB'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'MYB domain transcription factor family',]$TF_class <- 'MYB'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'GeBP like',]$TF_class <- 'GeBP'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'WRKY domain transcription factor family',]$TF_class <- 'WRKY'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'PHD finger transcription factor',]$TF_class <- 'PHD'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'bZIP transcription factor family',]$TF_class <- 'bZIP'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'GRAS transcription factor family',]$TF_class <- 'GRAS'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'CCAAT box binding factor family, HAP3',]$TF_class <- 'CCAAT box (HAP3)'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'CCAAT box binding factor family, HAP5',]$TF_class<- 'CCAAT box (HAP5)'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'C2C2(Zn) CO-like, Constans-like zinc finger family',]$TF_class <- 'C2C2(Zn), CO-like'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'ovate family OFP',]$TF_class <- 'OFP'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'C2C2(Zn) DOF zinc finger family',]$TF_class <- 'C2C2(Zn), DOF'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'C2C2(Zn) YABBY family',]$TF_class <- 'C2C2(Zn), YABBY' 
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'C2C2(Zn) GATA transcription factor family',]$TF_class <- 'C2C2(Zn), GATA' 
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'C2C2(Zn) GATA transcription factor family,DNA',]$TF_class <- 'C2C2(Zn), GATA'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'Global transcription factor group',]$TF_class <- 'Global'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'EIN3-like(EIL) transcription factor family',]$TF_class <- 'EIL'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'CPP(Zn),CPP1-related transcription factor family',]$TF_class <- 'CPP1-related'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'E2F/DP transcription factor family',]$TF_class <-'E2F/DP' 
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'FHA transcription factor',]$TF_class <- 'FHA'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'CCAAT box binding factor family, DR1',]$TF_class <- 'CCAAT box (DR1)'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'MYB domain transcription factor family,DNA',]$TF_class <- 'MYB'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'GRF zinc finger family',]$TF_class <- 'GRF'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'HB,Homeobox transcription factor family',]$TF_class <- 'Homeobox'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'G2-like transcription factor family, GARP',]$TF_class <- 'GARP'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'HSF,Heat-shock transcription factor family',]$TF_class <- 'Heat-shock factors'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'ARF, Auxin Response Factor family',]$TF_class <- 'ARF'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'sigma like plants',]$TF_class <- 'sigma like'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'MADS box transcription factor family',]$TF_class <- 'MADS'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'MYB-related transcription factor family',]$TF_class <- 'MYB-related'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'ABI3/VP1-related B3-domain-containing transcription factor family',]$TF_class <- 'ABI3/VP1'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'B3 transcription factor family',]$TF_class <- 'B3'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'TCP transcription factor family',]$TF_class <- 'TCP'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'NAC domain transcription factor family',]$TF_class <- 'NAC'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'NIN-like bZIP-related family',]$TF_class <- 'NIN-like'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'Pseudo ARR transcription factor family',]$TF_class <- 'Pseudo ARR'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'Psudo ARR transcription factor family',]$TF_class <- 'Pseudo ARR'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'CCAAT box binding factor family, HAP2',]$TF_class <- 'CCAAT box (HAP2)'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'Trihelix, Triple-Helix transcription factor family',]$TF_class <- 'Trihelix'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'putative transcription regulator',]$TF_class <- 'putative'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'AT-rich interaction domain containing transcription factor family',]$TF_class <- 'AT-rich interaction'
Sprint2_all_tfs[Sprint2_all_tfs$TF_class == 'transcription factor jumonji (jmjC) domain-containing protein',]$TF_class <- 'jmjC'


Sprint2_all_tfs %>% View()

chromosome_grouper <- function(df) {
  df$start <- rep(0, nrow(df))
  df$end <- rep(0, nrow(df))
  local_pivot <- table(unlist(df$TF_class)) %>% as.data.frame() %>% arrange(desc(Freq)) %>% mutate(Var1 = as.character(Var1), Freq = as.numeric(Freq))
  for (i in c(1:nrow(local_pivot))) {
    start_counter <- 0
    end_counter <- 1
    for (j in c(1:nrow(df))) {
      if (local_pivot[i,1] == df[j, 'TF_class']) {
        df[j, 'start'] <- df[j, 'start'] + start_counter
        df[j, 'end'] <- df[j, 'end'] + end_counter
        start_counter <- start_counter + 1
        end_counter <- end_counter + 1
      }
    }
  }
  return(df)
}

Sprint2_all_tfs %<>% chromosome_grouper()
Sprint2_segment <- Sprint2_all_tfs[,c('TF_class', 'start', 'end', 'transcript_id')]
Sprint2_mapping <- Sprint2_all_tfs[,c('TF_class', 'start', 'transcript_id', 'b', 'BLASTP_names')]

Sprint2_segment %>% View()
Sprint2_mapping %>% View()
Sprint2_all_tfs %>% View()

#### test dummy
seg.num <- 10
ind.num <- 20
seg.po <- c(20:50)
link.num <- 10
link.pg.num <- 10

sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, link.pg=link.pg.num)
names(sim.out)
head(sim.out$seg.frame[,c(1:3)])

set.seed(1234)
seg.f <- sim.out$seg.frame
seg.v <-  sim.out$seg.mapping
link.v  <-  sim.out$seg.link
link.pg.v <- sim.out$seg.link.pg
seg.num <- length(unique(seg.f[,1]))

seg.name <- paste("chr", 1:seg.num, sep="")
db <- segAnglePo(seg.f, seg=seg.name)
colors <- rainbow(seg.num, alpha=0.5)

par(mar=c (2, 2, 2, 2))
plot(c(1, 800) , c(1, 800) , type="n", axes=FALSE, xlab="", ylab="", main="" )
circos(R=400, cir=db, type="chr",  col=colors, print.chr.lab=TRUE, W=4, scale=TRUE)
circos(R=360, cir=db , W=40, mapping=seg.v, col.v =3, type="l", B=TRUE, col=colors[1], lwd=2, scale=TRUE)
circos(R=320, cir=db, W=40, mapping=seg.v, col.v=3, type="ls", B=FALSE, col=colors[9], lwd=2, scale=TRUE)
circos(R=280, cir=db , W=40, mapping=seg.v, col.v =3, type="lh", B=TRUE, col=colors[7], lwd=2, scale=TRUE)
circos (R=240, cir=db , W=40, mapping=seg.v, col.v=19, type="ml" , B=FALSE, col=colors, lwd=2, scale=TRUE)
circos(R=200, cir=db, W=40, mapping=seg.v, col.v =19, type="ml2" , B=TRUE, col=colors, lwd=2)
circos(R=160, cir=db , W=40, mapping=seg.v, col.v=19, type="ml3" , B=FALSE, cutoff=5, lwd=2)
# 11 c i r c o s (R=150, c i r=db , W=40, mapping=l i n k . v , type=" l ink" , lwd=2, co l=co lors [ c ( 1 , 7 ) ] ) ;
# 12 c i r c o s (R=150, c i r=db , W=40, mapping=l i n k . p g . v , type=" l i n k . p g " , lwd=2, co l=sample(
#   co lors , li n k. p g. n um ) )

### Not so test dummy: Sprint-2

Sprint2_segment %>% View()
Sprint2_mapping %>% View()
Sprint2_all_tfs %>% View()

par(mar=c (2, 2, 2, 2))
plot(c(1, 800) , c(1, 800) , type="n", axes=FALSE, xlab="", ylab="", main="" )

sim.out %>% head()
# sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, link.pg=link.pg.num)
sim.out <- sim.circos(Sprint2_segment, Sprint2_mapping)
nrow(table(unlist(Sprint2_all_tfs$TF_class)))

names(Sprint2_mapping)
db <- segAnglePo(seg.f, seg=seg.name)
par(mar=c (2, 2, 2, 2))
plot(c(1, 800) , c(1, 800) , type="n", axes=FALSE, xlab="", ylab="", main="" )

circos(R=300, cir=segAnglePo(Sprint2_mapping, seg = Sprint2_mapping$TF_class %>% unique()),
       type="chr",  col=rainbow(nrow(table(unlist(Sprint2_all_tfs$TF_class))), alpha = 0.5),
       print.chr.lab=TRUE, W=4, scale=TRUE)
circos(R=300, cir=segAnglePo(Sprint2_mapping, seg = Sprint2_mapping$TF_class %>% unique()),
       W=100, mapping=Sprint2_mapping %>% dplyr::select(b), col.v=1, type="heatmap2" ,
              cluster=TRUE, col.bar=TRUE, lwd=0.1, col="blue") 
length(unique(Sprint2_mapping[,1]))
movies %>% View()

##### SNP UpSetR representation - a new attempt
het_Zhewan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhewan1_corr_het.ann.vcf', header = F, sep = '\t', colClasses = 'character')
het_Zhewan$Moderate <- ifelse(grepl('MODERATE', het_Zhewan$V8), 1, 0)
het_Zhewan$High <- ifelse(grepl('HIGH', het_Zhewan$V8), 1, 0)
het_Zhewan$identifier <- sapply(c(1:nrow(het_Zhewan)), function(x) paste(as.character(het_Zhewan[x,c(1,2,4,5)]), collapse = ''))

homo_Zhewan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhewan1_corr_homo.ann.vcf', header = F, sep = '\t', colClasses = 'character')
homo_Zhewan$Moderate <- ifelse(grepl('MODERATE', homo_Zhewan$V8), 1, 0)
homo_Zhewan$High <- ifelse(grepl('HIGH', homo_Zhewan$V8), 1, 0)
homo_Zhewan$identifier <- sapply(c(1:nrow(homo_Zhewan)), function(x) paste(as.character(homo_Zhewan[x,c(1,2,4,5)]), collapse = ''))

het_Zhongwan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhongwan6_corr_het.ann_entries.vcf', header = F, sep = '\t', colClasses = 'character')
het_Zhongwan$Moderate <- ifelse(grepl('MODERATE', het_Zhongwan$V8), 1, 0)
het_Zhongwan$High <- ifelse(grepl('HIGH', het_Zhongwan$V8), 1, 0)
het_Zhongwan$identifier <- sapply(c(1:nrow(het_Zhongwan)), function(x) paste(as.character(het_Zhongwan[x,c(1,2,4,5)]), collapse = ''))
homo_Zhongwan[grepl('HIGH', homo_Zhongwan$V8),] %>% nrow()

homo_Zhongwan <- read.table('~/Pea_transcriptomics/VC_annotation_database/Zhongwan6_corr_homo.ann.vcf', header = F, sep = '\t', colClasses = 'character')
homo_Zhongwan$Moderate <- ifelse(grepl('MODERATE', homo_Zhongwan$V8), 1, 0)
homo_Zhongwan$High <- ifelse(grepl('HIGH', homo_Zhongwan$V8), 1, 0)
homo_Zhongwan$identifier <- sapply(c(1:nrow(homo_Zhongwan)), function(x) paste(as.character(homo_Zhongwan[x,c(1,2,4,5)]), collapse = ''))

het_Sprint <- read.table('~/Pea_transcriptomics/VC_annotation_database/Sprint2_corr_het.ann.vcf', header = F, sep = '\t', colClasses = 'character')
het_Sprint$Moderate <- ifelse(grepl('MODERATE', het_Sprint$V8), 1, 0)
het_Sprint$High <- ifelse(grepl('HIGH', het_Sprint$V8), 1, 0)
het_Sprint$identifier <- sapply(c(1:nrow(het_Sprint)), function(x) paste(as.character(het_Sprint[x,c(1,2,4,5)]), collapse = ''))
het_Sprint$High %>% sum()


all_hets <- rbind(het_Zhewan, het_Zhongwan, het_Sprint)
all_hets <- all_hets[!duplicated(all_hets$identifier),]
all_hets$'Zhewan-1' <- ifelse(all_hets$identifier %in% het_Zhewan$identifier, 1, 0)
all_hets$'Zhongwan-6' <- ifelse(all_hets$identifier %in% het_Zhongwan$identifier, 1, 0)
all_hets$'Sprint-2' <- ifelse(all_hets$identifier %in% het_Sprint$identifier, 1, 0)
all_hets <- all_hets[,c(1,11,12,14,15,16)]
colnames(all_hets)[1] <- 'Transcript'

all_hets %>% View()
mutations %>% View()
upset(all_hets, sets = c('Zhewan-1', 'Zhongwan-6', 'Sprint-2', 'Moderate', 'High'), sets.bar.color = '#6dc585',
      order.by = 'freq', text.scale = c(2.1, 2, 1.8, 1.5, 2, 2), mainbar.y.label = 'Number of SNPs')

all_homos <- rbind(homo_Zhewan, homo_Zhongwan, het_Sprint)
all_homos <- all_homos[!duplicated(all_homos$identifier),]
all_homos$'Zhewan-1' <- ifelse(all_homos$identifier %in% homo_Zhewan$identifier, 1, 0)
all_homos$'Zhongwan-6' <- ifelse(all_homos$identifier %in% homo_Zhongwan$identifier, 1, 0)
all_homos$'Sprint-2' <- ifelse(all_homos$identifier %in% het_Sprint$identifier, 1, 0)
all_homos <- all_homos[,c(1,11,12,14,15,16)]
colnames(all_homos)[1] <- 'Transcript'

upset(all_homos, sets = c('Zhewan-1', 'Zhongwan-6', 'Sprint-2', 'Moderate', 'High'), sets.bar.color = '#6dc585',
      order.by = 'freq', text.scale = c(2.1, 2, 1.8, 1.5, 2, 2), mainbar.y.label = 'Number of SNPs')

zhezho_only <- rbind(homo_Zhewan, homo_Zhongwan, het_Zhewan, het_Zhongwan)
zhezho_only <- zhezho_only[!duplicated(zhezho_only$identifier),]
zhezho_only$'Zhewan-1' <- ifelse(zhezho_only$identifier %in% homo_Zhewan$identifier | zhezho_only$identifier %in% het_Zhewan$identifier , 1, 0)
zhezho_only$'Zhongwan-6' <- ifelse(zhezho_only$identifier %in% homo_Zhongwan$identifier | zhezho_only$identifier %in% het_Zhongwan$identifier, 1, 0)
zhezho_only$Heterozygous <- ifelse(zhezho_only$identifier %in% het_Zhewan$identifier | zhezho_only$identifier %in% het_Zhongwan$identifier, 1, 0)
zhezho_only$Homozygous <- ifelse(zhezho_only$identifier %in% homo_Zhongwan$identifier | zhezho_only$identifier %in% homo_Zhewan$identifier, 1, 0)

zhezho_only <- zhezho_only[,c(1,11,12,14,15, 16, 17)]
colnames(zhezho_only)[1] <- 'Transcript'
upset(zhezho_only, sets = c('Zhewan-1', 'Zhongwan-6', 'Moderate', 'High', 'Homozygous', 'Heterozygous'), sets.bar.color = '#6dc585',
      order.by = 'freq', text.scale = c(2.1, 2, 1.8, 1.5, 2, 2), mainbar.y.label = 'Number of SNPs')

homo_Zhewan %>% colnames()
homo_Zhewan %<>% dplyr::rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(homo_Zhewan, ont = x, cutoff = 0.01))
homo_Zhongwan %<>% dplyr::rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(homo_Zhongwan, ont = x, cutoff = 0.01))
het_Zhongwan %<>% dplyr::rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(het_Zhongwan, ont = x, cutoff = 0.01))
het_Zhewan %<>% dplyr::rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(het_Zhewan, ont = x, cutoff = 0.01))
het_Sprint %<>% dplyr::rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(het_Sprint, ont = x, cutoff = 0.01))

homo_chinese <- homo_Zhewan[homo_Zhewan$identifier %in% homo_Zhongwan$identifier,]
homo_chinese %<>% rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(homo_chinese, ont = x, cutoff = 0.01))
het_chinese <- het_Zhewan[het_Zhewan$identifier %in% het_Zhongwan$identifier,]
het_chinese %<>% rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(het_chinese, ont = x, cutoff = 0.01))
Sprint_het_Zhewan <- het_Sprint[het_Sprint$identifier %in% het_Zhewan$identifier,]
Sprint_het_Zhewan %<>% rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(Sprint_het_Zhewan, ont = x, cutoff = 0.01))
Sprint_het_Zhongwan <- het_Sprint[het_Sprint$identifier %in% het_Zhewan$identifier,]
Sprint_het_Zhongwan %<>% rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(Sprint_het_Zhongwan, ont = x, cutoff = 0.01))

het_Sprint %<>% rename(target_id = V1)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(Sprint_het_Zhongwan, ont = x, cutoff = 0.01))

sapply(ls()[grepl('het_Sprint', ls())],
       function(x) get(x) %>% write.table(., file.path('~/Pea_transcriptomics/VC_annotation_database/GO_novel_annotations/', paste(x, 'tsv', sep='.')),
                                          sep = '\t', col.names = T, row.names = F))
annot_full[grepl('GO:0046686', annot_full$GOs),] %>% filter(transcript_id %in% het_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names)

annot_full[grepl(paste0('(', paste(TF_GO, collapse = ')|('), ')', sep = ''),
  annot_full$GOs) | grepl('regulation of transcription', annot_full$MapMan_terms),] %>% filter(transcript_id %in% het_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names) %>% nrow()

homo_Zhongwan %>% filter(target_id %in% annot_full[grepl('GO:0071215', annot_full$GOs),]$transcript_id) %>% filter(!(identifier %in% homo_Zhewan$identifier))
annot_full[grepl('GO:0071215', annot_full$GOs),] %>% filter(transcript_id %in% het_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms)


homo_Zhewan %>% filter(target_id == 'TRINITY_30_DN1752_c0_g1_i1') %>% filter(identifier)

TF_Zhongwan_snps_het <- annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% het_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(MapMan_terms = as.character(MapMan_terms)) %>% mutate(TF_class = sapply(MapMan_terms, function(x) x %>% strsplit(., '\\.') %>% unlist() %>% unname() %>% nth(3)))
TF_Zhongwan_snps_homo <- annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% homo_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(MapMan_terms = as.character(MapMan_terms)) %>% mutate(TF_class = sapply(MapMan_terms, function(x) x %>% strsplit(., '\\.') %>% unlist() %>% unname() %>% nth(3)))

TF_Zhewan_snps_homo <- annot_full[grepl(tf_template, annot_full$MapMan_terms),] %>% filter(transcript_id %in% homo_Zhewan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, MapMan_terms) %>% mutate(MapMan_terms = as.character(MapMan_terms)) %>% mutate(TF_class = sapply(MapMan_terms, function(x) x %>% strsplit(., '\\.') %>% unlist() %>% unname() %>% nth(3)))
TF_Zhewan_snps_homo %>% View()
TF_Zhewan_snps_homo %<>% mutate(TF_class = as.character(TF_class))
TF_Zhewan_snps_homo[TF_Zhewan_snps_homo$TF_class == 'HB,Homeobox transcription factor family,development',]$TF_class <- 'HB'
TF_Zhewan_snps_homo[grepl('AP2/EREBP, APETALA2/Ethylene-responsive elemen', TF_Zhewan_snps_homo$TF_class),]$TF_class <- 'AP2/EREBP'
TF_Zhewan_snps_homo[grepl('Trihelix, Triple-Helix transcription factor family', TF_Zhewan_snps_homo$TF_class),]$TF_class <- 'Trihelix'


homo_Zhewan %>% filter(target_id %in% annot_full[grepl(tf_template, annot_full$MapMan_terms),]$transcript_id) %>% filter(identifier %in% homo_Zhongwan$identifier) %>% View()


TF_Zhewan_snps_homo %>% View()

annot_full %>% filter(transcript_id == 'TRINITY_10_DN13275_c0_g1_i2') %>% dplyr::select(BLASTP_names, BLASTX_names)

annot_full %>% filter(transcript_id == 'TRINITY_10_DN11860_c0_g1_i1') %>% dplyr::select(BLASTP_names, BLASTX_names, plants_BLASTP, MapMan_terms)

zhezho_only %>% colnames()

common_chinese_het <- zhezho_only %>% filter(`Zhewan-1` == 1, `Zhongwan-6` == 1, Heterozygous == 1, Homozygous == 0)
common_chinese_het %<>% mutate(target_id = Transcript)
common_chinese_homo <- zhezho_only %>% filter(`Zhewan-1` == 1, `Zhongwan-6` == 1, Heterozygous == 0, Homozygous == 1)
common_chinese_homo %<>% mutate(target_id = Transcript)
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(common_chinese_het, ont = x, cutoff = 0.01))
sapply(c('BP', 'CC', 'MF'), function(x) iterating_function_GO(common_chinese_homo, ont = x, cutoff = 0.01))

sapply(ls()[grepl('_common_chinese', ls())],
       function(x) get(x) %>% write.table(., file.path('~/Pea_transcriptomics/VC_annotation_database/GO_novel_annotations/', paste(x, 'tsv', sep='.')),
                                          sep='\t', col.names = T, row.names = F))

annot_full[grepl('GO:0009686', annot_full$GOs),] %>% filter(transcript_id %in% het_Zhewan$target_id, transcript_id %in% het_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, plants_BLASTP)
annot_full[grepl('GO:0009686', annot_full$GOs),] %>% filter(transcript_id %in% homo_Zhewan$target_id, transcript_id %in% homo_Zhongwan$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names, plants_BLASTP)
homo_Zhewan %>% filter(identifier %in% homo_Zhongwan$identifier, target_id %in% annot_full[grepl('GO:0009686', annot_full$GOs),]$transcript_id) %>% View()
annot_full[grepl('GO:0003774', annot_full$GOs),] %>% filter(transcript_id %in% het_Sprint$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names) %>% View()
annot_full[grepl('GO:0009615', annot_full$GOs),] %>% filter(transcript_id %in% het_Sprint$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names) %>% View()
annot_full[grepl('GO:0009294', annot_full$GOs),] %>% filter(transcript_id %in% het_Sprint$target_id) %>% dplyr::select(transcript_id, BLASTP_names, BLASTX_names) %>% View()



###### Interpoint dispersion
install.packages('cluster')
install.packages('factoextra')
library(cluster)
library(factoextra)

scatter_anova <- function(subset){
#### here be function relying on   
}

colnames(alab_full)
alab_full %>% ncol()
rownames(alab_full) <- alab_full$target_id
alab_dummy <- alab_full[,-c(1,26:28)]
alab_dummy <- sapply(c(1:4,9:12,17:20), function(x) sqrt(alab_dummy[,x]^2 + alab_dummy[,x+4]^2))
colnames(alab_dummy) <- c('Sprint-2,1', 'Sprint-2,2', 'Sprint-2,3', 'Sprint-2,4',
                          'Zhewan-1,1', 'Zhewan-1,2', 'Zhewan-1,3', 'Zhewan-1,4',
                          'Zhongwan-6,1', 'Zhongwan-6,2', 'Zhongwan-6,3', 'Zhongwan-6,4')
rownames(alab_dummy) <- rownames(alab_full)
alab_dummy %>% head()
alab_dummy %<>% scale()
fviz_dist(get_dist(alab_dummy[c(1:1000),], stand = T, method = 'kendall'), 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_nbclust(alab_dummy[c(1:10000),], kmeans, method = "gap_stat")

t(data.frame(split(unname(unlist(alab_dummy[1,])), ceiling(seq_along(unname(unlist(alab_dummy[1,]))) / 4))))
data.frame(V1 = unname(unlist(alab_dummy[1,])), V2 = sapply(seq(0,2), function(x) rep(x,4)) %>% as.vector())


dotwise_aov <- function(xrow, nor, nog) {
  local_array <- data.frame(V1 = unname(unlist(xrow)), V2 = sapply(seq(0,nog-1), function(x) rep(x,nor)) %>% as.vector()) %>% mutate(V2 = as.factor(V2))
return(TukeyHSD(aov(V1 ~ V2, local_array), meth))
  }
?TukeyHSD
dotwise_aov(alab_dummy[1,], 4, 3)
