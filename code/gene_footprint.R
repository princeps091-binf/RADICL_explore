# Examine specific read distribution of individual RefSeq genes
library(tidyverse)
library(readr)
library(vroom)
library(GenomicRanges)
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

chromo<-'chr19'
chr_exon<-exon_tbl %>% 
  filter(chrom == chromo)
chr_gene<-gene_anno %>% 
  filter(chrom == chromo)
chr_anno_GRanges<-anno_GRanges[seqnames(anno_GRanges)==chromo]
chr_exon_GRanges<-exon_GRanges[seqnames(exon_GRanges)==chromo]
# intron
## setdiff for each gene 
plan(multisession,workers=4)
intron_tbl<-future_map_dfr(unique(chr_gene$name),function(i){
  tmp_gene_GRange<-chr_anno_GRanges[chr_anno_GRanges@elementMetadata$name == i]
  tmp_exon_GRanges<-chr_exon_GRanges[chr_exon_GRanges@elementMetadata$name == i]
  tmp_intron<-setdiff(tmp_gene_GRange,tmp_exon_GRanges)
  tibble(as.data.frame(tmp_intron)) %>% 
    mutate(name=i)
  
})
plan(sequential)

exon_intron_tbl<-chr_exon %>% 
  dplyr::select(name,chrom,strand,exonStarts,exonEnds) %>% 
  dplyr::rename(start=exonStarts,end=exonEnds) %>% 
  dplyr::mutate(type='exon') %>% 
  bind_rows(.,
            intron_tbl %>% 
              dplyr::select(name,seqnames,strand,start,end) %>% 
              dplyr::rename(chrom=seqnames) %>% 
              dplyr::mutate(type='intron'))

exon_intron_GRanges<-GRanges(seqnames=exon_intron_tbl$chrom,
                             IRanges(
                               start=exon_intron_tbl$start,
                               end=exon_intron_tbl$end,
                             ),
                             strand=exon_intron_tbl$strand)

mcols(exon_intron_GRanges)<-tibble(name=exon_intron_tbl$name,type=exon_intron_tbl$type)


RNA_cluster_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/cluster/RADICL_iPSC_RNA_chr19.bed_cluster.txt",
                       col_names = F,delim = "\t")
read_GRanges<-GRanges(seqnames=RNA_cluster_tbl$X1,
                      ranges = IRanges(start=RNA_cluster_tbl$X2,
                                       end=RNA_cluster_tbl$X3
                      ),
                      strand=RNA_cluster_tbl$X6
)
mcols(read_GRanges)<-tibble(ID=RNA_cluster_tbl$X4)
i<-'NR_038237.1'

tmp_refseq_GRange<-exon_intron_GRanges[exon_intron_GRanges@elementMetadata$name==i]
as_tibble(tmp_refseq_GRange) %>% 
  mutate(count=countOverlaps(tmp_refseq_GRange,read_GRanges),
         rate=count/width) %>% 
  mutate(p.val=pmap_dbl(list(count,width),function(count,width){
    poisson.test(count,width,r=gene_rate,alternative='greater')$p.value
  }))

tmp_over<-as_tibble(findOverlaps(tmp_refseq_GRange,RNA_read_GRanges)) %>% 
  mutate(seqnames=as.vector(seqnames(tmp_refseq_GRange))[queryHits],
         start=start(tmp_refseq_GRange)[queryHits],
         end=end(tmp_refseq_GRange)[queryHits],
         strand=as.vector(strand(tmp_refseq_GRange))[queryHits],
         width=abs(start(tmp_refseq_GRange)[queryHits] - end(tmp_refseq_GRange)[queryHits]),
         read.ID=mcols(RNA_read_GRanges)$ID[subjectHits]) 
# Read DNA side with read ID filter from RNA side