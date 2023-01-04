library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
RNA_cluster_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/cluster/RADICL_iPSC_RNA_chr22.bed_cluster.txt",
      col_names = F,delim = "\t")

cluster_summary_tbl<-RNA_cluster_tbl %>% 
  group_by(X7) %>% 
  summarise(start=min(X2),
            end=max(X3),
            n=n())
cluster_summary_tbl %>% 
  filter(n > 1) %>% 
  ggplot(.,aes(n))+
  geom_density()+
  scale_x_log10()


cluster_summary_tbl %>% 
  filter(n > 1) %>% 
  mutate(d=abs(end-start)) %>% 
  ggplot(.,aes(d))+
  geom_density()+
  scale_x_log10()

read_GRanges<-GRanges(seqnames=RNA_cluster_tbl$X1,
                      ranges = IRanges(start=RNA_cluster_tbl$X2,
                                       end=RNA_cluster_tbl$X3
                      ),
                      strand=RNA_cluster_tbl$X6
                      )
#-------------------------------------------------------------
# enforce column type based on column type listed in https://genome.ucsc.edu/cgi-bin/hgTables
gene_anno<-vroom("~/Documents/FANTOM6/data/annotation/refseq_all_hg38.tsv",
                 col_names = T,delim = "\t",col_types = list("i",'c','c','c','i','i','i','i','i','c','c','i','c','c','c'))

anno_GRanges<-GRanges(seqnames=gene_anno$chrom,
                      ranges = IRanges(start=gene_anno$txStart,
                                       end=gene_anno$txEnd
                      ),
                      strand=gene_anno$strand)
mcols(anno_GRanges)<-tibble(name=gene_anno$name)

exon_tbl<-gene_anno %>% 
  dplyr::select(name,chrom,txStart,txEnd,strand,exonCount,exonStarts,exonEnds) %>% 
  mutate(exonStarts=map(exonStarts,function(x){
    as.integer(unlist(str_split(x,',')))
  }),
  exonEnds=map(exonEnds,function(x){
    as.integer(unlist(str_split(x,',')))
  })) %>% 
  unnest(cols=c(exonStarts,exonEnds)) %>% 
  filter(!(is.na(exonStarts))|!(is.na(exonEnds)))

exon_GRanges<-GRanges(seqnames=exon_tbl$chrom,
                      IRanges(
                                       start=exon_tbl$exonStarts,
                                       end=exon_tbl$exonEnds,
                      ),
                      strand=exon_tbl$strand)
          
mcols(exon_GRanges)<-tibble(name=exon_tbl$name)

length(unique(queryHits(findOverlaps(read_GRanges,anno_GRanges))))/length(read_GRanges)
length(unique(queryHits(findOverlaps(read_GRanges,exon_GRanges))))/length(read_GRanges)
gene_anno %>% 
  dplyr::select(name,chrom, strand, txStart,  txEnd,exonCount) %>% 
  mutate(read.count=countOverlaps(anno_GRanges,read_GRanges)) %>% 
  filter(read.count>0) %>% 
  mutate(rd=read.count/(txEnd-txStart)) %>% 
  ggplot(.,aes(rd))+
  geom_density()+
  scale_x_log10()

exon_tbl %>% 
  mutate(read.count=countOverlaps(exon_GRanges,read_GRanges)) %>% 
  filter(read.count>0) %>% 
  ggplot(.,aes(read.count))+
  geom_density()+
  scale_x_log10()
