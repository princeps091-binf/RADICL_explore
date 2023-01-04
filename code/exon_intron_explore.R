# exon-intron explore
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
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

chr_exon<-exon_tbl %>% 
  filter(chrom == 'chr22')
chr_gene<-gene_anno %>% 
  filter(chrom == 'chr22')
chr_anno_GRanges<-anno_GRanges[seqnames(anno_GRanges)=='chr22']
chr_exon_GRanges<-exon_GRanges[seqnames(exon_GRanges)=='chr22']
# intron
## gaps for each gene
plan(multisession,workers=4)
intron_tbl<-future_map_dfr(unique(chr_gene$name),function(i){
  tmp_gene_GRange<-chr_anno_GRanges[chr_anno_GRanges@elementMetadata$name == i]
  tmp_exon_GRanges<-chr_exon_GRanges[chr_exon_GRanges@elementMetadata$name == i]
  tmp_intron<-gaps(tmp_exon_GRanges,start=start(tmp_gene_GRange))
  tibble(as.data.frame()) %>% 
    mutate(name=i)
  
})
plan(sequential)

exon_intron_ratio<-chr_exon %>% 
  dplyr::select(name,chrom,strand,exonStarts,exonEnds) %>% 
  dplyr::rename(start=exonStarts,end=exonEnds) %>% 
  dplyr::mutate(type='exon') %>% 
  bind_rows(.,
            intron_tbl %>% 
              dplyr::select(name,seqnames,strand,start,end) %>% 
              dplyr::rename(chrom=seqnames) %>% 
              dplyr::mutate(type='intron')) %>% 
  group_by(name,type) %>% 
  summarise(w=sum(abs(end-start)),
            n=n()) 
chr_gene_size<-chr_gene %>% 
  group_by(name) %>% 
  summarise(g.size=sum(abs(txStart - txEnd))) 
# check exon + intron cover whole the gene
exon_intron_ratio %>% 
  ungroup() %>% 
  group_by(name) %>% 
  summarise(sw=sum(w),
            nc=sum(n)-1) %>% 
  full_join(.,chr_gene_size) %>% 
  mutate(diffei=g.size-(sw+nc)) %>% 
  ggplot(.,aes(diffei))+
  geom_density()

corr_gene_size<-exon_intron_ratio %>% 
  ungroup() %>% 
  group_by(name) %>% 
  summarise(nc=sum(n)-1) %>% 
  full_join(.,chr_gene_size) %>% 
  mutate(gc.size=g.size-nc) %>% 
  dplyr::select(name,gc.size)

order_name<-exon_intron_ratio %>% 
  full_join(.,corr_gene_size) %>% 
  mutate(rw=w/gc.size) %>% 
  filter(type=='exon') %>% 
  arrange(desc(rw))
  
exon_intron_ratio%>% 
  mutate(name=factor(name,levels=order_name$name)) %>% 
  ggplot(.,aes(name,w,fill=type))+
  geom_bar(stat = 'identity',position='fill')
