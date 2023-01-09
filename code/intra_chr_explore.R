# Distance decay
library(tidyverse)
library(vroom)

RADICL_intra_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_iPSC_intra.bed",
                       col_names = F,delim = "\t")


RADICL_intra_tbl %>% 
  mutate(DNA.mid=X2 + (X3-X2)/2,
         RNA.mid=X8 + (X8-X9)/2) %>% 
  ggplot(.,aes(abs(DNA.mid - RNA.mid)))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(X1~.)

DNA_read_GRanges<-GRanges(seqnames=RADICL_intra_tbl$X7,
                      ranges = IRanges(start=RADICL_intra_tbl$X8,
                                       end=RADICL_intra_tbl$X9
                      ),
                      strand=RADICL_intra_tbl$X12
)
mcols(DNA_read_GRanges)<-tibble(ID=RADICL_intra_tbl$X10)

RNA_read_GRanges<-GRanges(seqnames=RADICL_intra_tbl$X1,
                          ranges = IRanges(start=RADICL_intra_tbl$X2,
                                           end=RADICL_intra_tbl$X3
                          ),
                          strand=RADICL_intra_tbl$X6
)
mcols(RNA_read_GRanges)<-tibble(ID=RADICL_intra_tbl$X4)

#-------------------------------------------------------
# Per gene examination
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

i<-'NM_003422.3'
tmp_refseq_GRange<-exon_intron_GRanges#[exon_intron_GRanges@elementMetadata$name==i]
tmp_over<-as_tibble(findOverlaps(tmp_refseq_GRange,RNA_read_GRanges)) %>% 
  mutate(seqnames=as.vector(seqnames(tmp_refseq_GRange))[queryHits],
         start=start(tmp_refseq_GRange)[queryHits],
         end=end(tmp_refseq_GRange)[queryHits],
         strand=as.vector(strand(tmp_refseq_GRange))[queryHits],
         width=abs(start(tmp_refseq_GRange)[queryHits] - end(tmp_refseq_GRange)[queryHits]),
         type=mcols(tmp_refseq_GRange)$type[queryHits],
         read.ID=mcols(RNA_read_GRanges)$ID[subjectHits]) 
tmp_over<-tmp_over%>% 
  left_join(.,as_tibble(DNA_read_GRanges[DNA_read_GRanges@elementMetadata$ID %in% tmp_over$read.ID]),
            by=c('read.ID'='ID'))

tmp_over %>% 
  mutate(ref.mid=start.x + width.x/2,
         DNA.mid=start.y + width.y/2) %>% 
  mutate(d= abs(ref.mid - DNA.mid)) %>% 
  ggplot(.,aes(d,color=type))+
  geom_density()+
  scale_x_log10()


as_tibble(tmp_refseq_GRange) %>% 
  mutate(count=countOverlaps(tmp_refseq_GRange,RNA_read_GRanges))

anno_GRanges

tmp_refseq_GRange<-exon_GRanges#[exon_intron_GRanges@elementMetadata$name==i]
tmp_over<-as_tibble(findOverlaps(tmp_refseq_GRange,RNA_read_GRanges)) %>% 
  mutate(seqnames=as.vector(seqnames(tmp_refseq_GRange))[queryHits],
         start=start(tmp_refseq_GRange)[queryHits],
         end=end(tmp_refseq_GRange)[queryHits],
         strand=as.vector(strand(tmp_refseq_GRange))[queryHits],
         width=abs(start(tmp_refseq_GRange)[queryHits] - end(tmp_refseq_GRange)[queryHits]),
         read.ID=mcols(RNA_read_GRanges)$ID[subjectHits]) 
tmp_over<-tmp_over%>% 
  left_join(.,as_tibble(DNA_read_GRanges[DNA_read_GRanges@elementMetadata$ID %in% tmp_over$read.ID]),
            by=c('read.ID'='ID'))

tmp_over %>% 
  mutate(ref.mid=start.x + width.x/2,
         DNA.mid=start.y + width.y/2) %>% 
  mutate(d= abs(ref.mid - DNA.mid)) %>% 
  ggplot(.,aes(1,d))+
  geom_boxplot()+
  scale_y_log10()
