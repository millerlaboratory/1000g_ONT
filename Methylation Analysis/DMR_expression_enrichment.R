library(rtracklayer)
library(data.table)
library(tidyverse)
library(plyranges)
library(magrittr)

gff <- readGFFAsGRanges("./ReferenceFiles/GENCODEv35_GRCh38.p13/gencode.v35.primary_assembly.annotation.gtf")
genes <- gff[gff$type == "gene",]

dmrs <- fread("./DMRs_d1.5_p01_1KGp_samples.tsv")
dmrs
dmrs %<>% separate_rows(samples, sep=",")
dmr_gr <- makeGRangesFromDataFrame(dmrs, keep.extra.columns =T)
# flank DMR by 10kb 
start(dmr_gr) <- start(dmr_gr) - 10000
end(dmr_gr) <- end(dmr_gr) + 10000

dmr_genes <- dmr_gr %>% join_overlap_left(genes) %>% select(cpg,gene,samples,type.x,gene_id,gene_type,gene_name) %>% filter(gene_type=="protein_coding")
dmr_genes$samples <- gsub(dmr_genes$samples, pattern="GM", replacement="NA")

e_zscores <- fread("./residuals.eqtls.protein_coding.lincRNA.cov1.pcs30.gpc8.tsv")
e_zscores

e_zscores_filt <- e_zscores[e_zscores$SubjectID %in% dmr_genes$samples,]
dmr_genes_filt <- dmr_genes[dmr_genes$samples %in% e_zscores_filt$SubjectID,]
dmr_genes_filt$gene_id <- gsub(dmr_genes_filt$gene_id, pattern="(ENSG\\d+)\\.\\d+", replacement="\\1")
e_zscores_filt <- e_zscores[e_zscores$GeneName %in% dmr_genes_filt$gene_id,]

dmr_genes_ezscores <- as.data.frame(dmr_genes_filt) %>% rename("SubjectID" = samples, "GeneName" = gene_id) %>% full_join(e_zscores_filt)
dmr_genes_ezscores

dmr_genes_ezscores$DMR <- !is.na(dmr_genes_ezscores$cpg)
table(dmr_genes_ezscores$DMR, abs(dmr_genes_ezscores$eOutliers) > 4)


set.seed(1)
dmr_genes_ezscores$permuted <- sample(dmr_genes_ezscores$eOutliers)
enrich.result <- data.frame(Reduce(rbind,lapply(c(.1,seq(1, 4, 1)), function(x) {
    test <- data.frame(dmr_genes_ezscores$DMR, outlier=abs(dmr_genes_ezscores$eOutliers) > x)
    model <- glm(outlier ~ ., family=binomial(link="logit"), data=test)
    c(ezscore_threshold=x, summary(model)$coefficients[2,])
})))
enrich.result
enrich.result$Estimate %<>% as.numeric()
enrich.result$`Std..Error` %<>% as.numeric()
enrich.result$higher <- enrich.result$Estimate + 1.96*enrich.result$`Std..Error`
enrich.result$lower <- enrich.result$Estimate - 1.96*enrich.result$`Std..Error`

ggplot(enrich.result, aes(ezscore_threshold, Estimate, ymin=lower, ymax=higher)) + geom_line() + 
	geom_pointrange(shape=21, fatten=10, fill="cyan") + 
	geom_hline(yintercept=0,linetype="dashed") + 
	ylab("log(Odds Ratio)") + 
	xlab("Expression Zscore Threshold")
ggsave('eOutlier_enrichment_near_DMRS.pdf')
enrich.result

zscore_annotated_dmrs <- dmr_genes_ezscores %>% filter(DMR) %>% rename("e_zscore" = eOutliers)
write.table(zscore_annotated_dmrs, "./DMRs_d1.5_p01_1KGp_samples.AFGR_expression_zscore_annotated.tsv", quote=F, sep="\t", col.names=T, row.names=F)
enrich.result

length(unique(zscore_annotated_dmrs$GeneName))
length(unique(zscore_annotated_dmrs$SubjectID))
zscore_annotated_dmrs

length(unique(dmr_genes$samples))

zscore_annotated_dmrs %>% filter(abs(e_zscore) > 2) %>% ggplot(aes(type.x, e_zscore)) + 
	geom_jitter()
ggsave('DMR_outlier_direction.ezscore_jitter.pdf')
