# Examine specific read distribution of individual RefSeq genes
library(tidyverse)
library(readr)
library(vroom)
library(GenomicRanges)
library(furrr)
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


RNA_cluster_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/RADICL_iPSC_RNA_chr19.bed",
                       col_names = F,delim = "\t")
read_GRanges<-GRanges(seqnames=RNA_cluster_tbl$X1,
                      ranges = IRanges(start=RNA_cluster_tbl$X2,
                                       end=RNA_cluster_tbl$X3
                      ),
                      strand=RNA_cluster_tbl$X6
)
mcols(read_GRanges)<-tibble(ID=RNA_cluster_tbl$X4)



# Read DNA side with read ID filter from RNA side
i<-'XM_011527205.3'
tmp_refseq_GRange<-exon_intron_GRanges[exon_intron_GRanges@elementMetadata$name==i]
tmp_over<-as_tibble(findOverlaps(tmp_refseq_GRange,read_GRanges)) %>% 
  mutate(seqnames=as.vector(seqnames(tmp_refseq_GRange))[queryHits],
         start=start(tmp_refseq_GRange)[queryHits],
         end=end(tmp_refseq_GRange)[queryHits],
         strand=as.vector(strand(tmp_refseq_GRange))[queryHits],
         width=abs(start(tmp_refseq_GRange)[queryHits] - end(tmp_refseq_GRange)[queryHits]),
         read.ID=mcols(read_GRanges)$ID[subjectHits])
read_set<-tmp_over$read.ID
f <- function(x, pos){
  x %>% 
    filter(X4 %in% read_set)
  
}
DNA_side_tbl<-readr::read_delim_chunked('~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_iPSC_DNA.bed',
                          delim = "\t",
                          chunk_size = 1e4,
                          col_names = F,
                          DataFrameCallback$new(f))
tmp_over %>% 
  left_join(.,DNA_side_tbl,by=c('read.ID'='X4'))%>% 
  ggplot(.,aes(X2,1,color=as.factor(queryHits)))+
  geom_point(size=0.5)+
  facet_grid(X1~.)+
  theme(legend.position="none")
gene_chunk_set<-unique(tmp_over$queryHits)
chunk_combo_tbl<-expand_grid(a=gene_chunk_set,b=gene_chunk_set) %>% 
  filter(a != b)
plan(multisession,workers=4)
#chunk_combo_tbl<-
test_chunk<-  chunk_combo_tbl %>%
  dplyr::slice(1:20) %>% 
  mutate(inter.n=future_pmap_dfr(list(a,b),function(a,b){
    tmp_a<-tmp_over %>% 
      filter(queryHits == a) %>% 
      left_join(.,DNA_side_tbl,by=c('read.ID'='X4'))
    tmp_b<-tmp_over %>% 
      filter(queryHits == b) %>% 
      left_join(.,DNA_side_tbl,by=c('read.ID'='X4'))
    a_GRanges<-GRanges(seqnames=tmp_a$X1,
                       ranges = IRanges(start=tmp_a$X2,
                                        end=tmp_a$X3
                       ))
    b_GRanges<-GRanges(seqnames=tmp_b$X1,
                       ranges = IRanges(start=tmp_b$X2,
                                        end=tmp_b$X3
                       ))
    
    a_chr<-as.vector(seqnames(a_GRanges))
    b_chr<-as.vector(seqnames(b_GRanges))
    # change metric to instead reflect proportion of reads sharing chromosome with at least one other read from the other set
    chr_sim<-expand_grid(a_chr,b_chr) %>% 
      mutate(same=a_chr == b_chr) %>% 
      summarise(sim=sum(same)/n()) %>% 
      unlist()
    n_inter<-length(findOverlaps(a_GRanges,b_GRanges))
    dist_tbl<-as_tibble(distanceToNearest(a_GRanges,b_GRanges)) %>% 
      bind_rows(.,as_tibble(distanceToNearest(b_GRanges,a_GRanges)))
    return(tibble(chr.sim=chr_sim,n.inter=n_inter,dist=list(dist_tbl$distance)))
  }))
plan(sequential)
j<-1
g<-31

as_tibble(distanceToNearest(a_GRanges,b_GRanges)) %>% 
  bind_rows(.,as_tibble(distanceToNearest(b_GRanges,a_GRanges))) %>% 
  ggplot(.,aes(distance))+
  geom_density()
