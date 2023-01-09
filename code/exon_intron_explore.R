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
# exon vs intro gene coverage
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
  geom_bar(stat = 'identity',position='fill')+
  scale_fill_brewer(palette='Dark2')+
  xlab('Genes')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
#---------------------------------------------------
RNA_cluster_tbl<-vroom("~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/cluster/RADICL_iPSC_RNA_chr19.bed_cluster.txt",
                       col_names = F,delim = "\t")
read_GRanges<-GRanges(seqnames=RNA_cluster_tbl$X1,
                      ranges = IRanges(start=RNA_cluster_tbl$X2,
                                       end=RNA_cluster_tbl$X3
                      ),
                      strand=RNA_cluster_tbl$X6
)

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
length(unique(queryHits(findOverlaps(read_GRanges,exon_intron_GRanges))))/length(read_GRanges)

exon_intron_tbl %>% 
  mutate(read.count=countOverlaps(exon_intron_GRanges,read_GRanges)) %>% 
  filter(read.count>0) %>% 
  mutate(rd=read.count/(end-start)) %>% 
  ggplot(.,aes(rd,color=type))+
  geom_density()+
  scale_x_log10()
# within gene gini coeff and deviation from gene-wide rate
i<-'NR_038237.1'
gini_cpu_fn<-function(exon_intron_GRanges,read_GRanges,i){
  rate_vec<-countOverlaps(exon_intron_GRanges[exon_intron_GRanges@elementMetadata$name==i],read_GRanges)/width(exon_intron_GRanges[exon_intron_GRanges@elementMetadata$name==i])
  if(length(rate_vec)>1){
    combo_vec<-t(combn(1:length(rate_vec),2))
    return(sum(abs(rate_vec[combo_vec[,1]]-rate_vec[combo_vec[,2]]))/(length(rate_vec)**2*mean(rate_vec)))
  } else{
    return(1)
  }
}


gene_set<-unique(exon_intron_GRanges@elementMetadata$name)
plan(multisession,workers=3)
gene_gini<-future_map_dbl(gene_set,function(i){
  gini_cpu_fn(exon_intron_GRanges,read_GRanges,i)
})
plan(sequential)
gene_rate<-countOverlaps(anno_GRanges[anno_GRanges@elementMetadata$name==i],read_GRanges)/width(anno_GRanges[anno_GRanges@elementMetadata$name==i])

tmp_refseq_GRange<-exon_intron_GRanges[exon_intron_GRanges@elementMetadata$name==i]
as_tibble(tmp_refseq_GRange) %>% 
  mutate(count=countOverlaps(tmp_refseq_GRange,read_GRanges),
         rate=count/width) %>% 
  mutate(p.val=pmap_dbl(list(count,width),function(count,width){
    poisson.test(count,width,r=gene_rate,alternative='greater')$p.value
  }))

