# exon-intron explore
library(tidyverse)
library(vroom)
library(GenomicRanges)
library(furrr)
#-------------------------------------------------------------
transcript_annotation_file<-"~/Documents/FANTOM6/data/annotation/FANTOM_CAT.lv3_robust.bed"

transcript_tbl<-vroom(transcript_annotation_file,
                      col_types = list("c",'i','i','c','d','c','i','i','c','d','c','c'),
                      col_names = F,
                      delim="\t")

transctipt_GRanges<-GRanges(seqnames=transcript_tbl$X1,
                            ranges = IRanges(start=transcript_tbl$X2,
                                             end=transcript_tbl$X3
                            ),
                            strand=transcript_tbl$X6
)
mcols(transctipt_GRanges)<-tibble(ID=transcript_tbl$X4)

exon_tbl<-transcript_tbl %>% 
  dplyr::select(X4,X1,X2,X3,X6,X10,X11,X12) %>% 
  mutate(RelexonStarts=map(X12,function(x){
    as.integer(unlist(str_split(x,',')))
  }),
  exonLength=map(X11,function(x){
    as.integer(unlist(str_split(x,',')))
  })) %>% 
  unnest(cols=c(RelexonStarts,exonLength)) %>% 
  mutate(exonStart=X2+RelexonStarts) %>% 
  mutate(exonEnd=exonStart+exonLength) %>% 
  dplyr::select(X4,X1,exonStart,exonEnd,X6,X2,X3)

exon_GRanges<-GRanges(seqnames=exon_tbl$X1,
                      IRanges(
                        start=exon_tbl$exonStart,
                        end=exon_tbl$exonEnd,
                      ),
                      strand=exon_tbl$X6)

mcols(exon_GRanges)<-tibble(name=exon_tbl$X4)
#------------------------------------------------------
chromo<-'chr19'
chr_exon<-exon_tbl %>% 
  filter(X1 == chromo)
chr_gene<-transcript_tbl %>% 
  filter(X1 == chromo)
chr_transcript_GRanges<-transctipt_GRanges[seqnames(transctipt_GRanges)==chromo]
chr_exon_GRanges<-exon_GRanges[seqnames(exon_GRanges)==chromo]
# intron
## setdiff for each gene 
plan(multisession,workers=4)
intron_tbl<-future_map_dfr(unique(chr_gene$X4),function(i){
  tmp_gene_GRange<-chr_transcript_GRanges[chr_transcript_GRanges@elementMetadata$ID == i]
  tmp_exon_GRanges<-chr_exon_GRanges[chr_exon_GRanges@elementMetadata$name == i]
  tmp_intron<-setdiff(tmp_gene_GRange,tmp_exon_GRanges)
  tibble(as.data.frame(tmp_intron)) %>% 
    mutate(name=i)
  
})
plan(sequential)
# exon vs intro gene coverage
exon_intron_ratio<-chr_exon %>% 
  dplyr::select(X4,X1,X6,exonStart,exonEnd) %>% 
  dplyr::rename(seqnames=X1,strand=X6,start=exonStart,end=exonEnd,name=X4) %>% 
  dplyr::mutate(type='exon') %>% 
  bind_rows(.,
            intron_tbl %>% 
              dplyr::select(name,seqnames,strand,start,end) %>% 
              dplyr::mutate(type='intron')) %>% 
  group_by(name,type) %>% 
  summarise(w=sum(abs(end-start)),
            n=n()) 
chr_gene_size<-chr_gene %>% 
  group_by(X4) %>% 
  summarise(g.size=sum(abs(X3 - X2))) %>% 
  dplyr::rename(name=X4)
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
  ungroup %>% 
  full_join(.,corr_gene_size) %>% 
  mutate(rw=w/gc.size) %>% 
  filter(type=='exon') %>% 
  arrange(desc(rw))

exon_intron_ratio%>% 
  ungroup %>% 
  mutate(name_o=factor(name,levels=order_name$name)) %>% 
  ggplot(.,aes(name_o,w,fill=type))+
  geom_bar(stat = 'identity',position='fill')+
  scale_fill_brewer(palette='Dark2')+
  xlab('Genes')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
