#script to check how many OMIM genes are in segdups

segDup <- read.table('../../../../Desktop/LRS/IGVtracks/segDups.hg38.bed') %>% select(V1, V2, V3)
segDup = unique(segDup)

names(segDup) <- c("chr", "start", "end")
#editing annotation based on 10 kb margin
new_segdup_annotation = segDup %>% mutate(start = start - 10000,
                                         end = end + 10000) 

new_segdup_annotation %>% arrange(chr, start)

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                 "chr20", "chr21", "chr22")

omim_anno = read.csv("../hg38_annotation_by_chromosome/OMIM_gene_phen_hg38.txt", sep = "\t")
colnames(omim_anno) <- c("chr", "start","end","gene_name", "ensembl_id")
omim_anno = omim_anno %>% filter(chr %in% chromosomes) %>% select(-ensembl_id)


#checking intersection
gr = GenomicRanges::GRanges(seqnames = new_segdup_annotation$chr,
                           ranges = IRanges::IRanges(start = new_segdup_annotation$start, end = new_segdup_annotation$end))

compressed_gr = GenomicRanges::reduce(gr)

new_segdup_annotation = as.data.frame(compressed_gr) %>% select(-strand, -width) 


omim_gr = GenomicRanges::GRanges(seqnames = omim_anno$chr,
                            ranges = IRanges::IRanges(start = omim_anno$start, end = omim_anno$end))

overlaps = findOverlaps(compressed_gr, omim_gr)
overlaps = as.data.frame(overlaps)

overlaps$segDup_chr = new_segdup_annotation$seqnames[overlaps$queryHits]
overlaps$segDup_start = new_segdup_annotation$start[overlaps$queryHits]
overlaps$segDup_end = new_segdup_annotation$end[overlaps$queryHits]
overlaps$chr = omim_anno$chr[overlaps$subjectHits]
overlaps$gene_name = omim_anno$gene_name[overlaps$subjectHits]
overlaps$segDup = "yes"

overlaps = overlaps %>% select(chr, gene_name, segDup)
overlaps = unique(overlaps)
overlaps = left_join(omim_anno, overlaps, by = c("chr","gene_name")) %>% 
  mutate(segDup = ifelse(is.na(segDup), "no", segDup)) %>%
  arrange(gene_name)

write.table(overlaps, "overlaps_segdup_omim.txt", col.names  = TRUE, row.names = FALSE, quote = FALSE)
