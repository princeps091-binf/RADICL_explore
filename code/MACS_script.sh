# select relevant DNA columns in RADICL bed file
cut -d$'\t' -f7-11 15.AGTTCC.2.bed > Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA.bed
awk '$7 == "chr19"' 15.AGTTCC.2.bed > RADICL_iPSC_chr19.bed
awk '$7 == $1' 15.AGTTCC.2.bed > RADICL_iPSC_intra.bed
cut -d$'\t' -f7-11 RADICL_iPSC_chr19.bed
cut -d$'\t' -f1-6 RADICL_iPSC_chr19.bed

# 
macs3 callpeak  -f BED --nomodel --shift -37 --extsize 75 -t Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA.bed > test.bed
# Size of full nucleosome
macs3 callpeak  -f BED --nomodel --shift -72 --extsize 147 -t Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA_chr19.bed -B -n chr19_nucleosome

macs3 callpeak  -f BED --nomodel --shift -72 --extsize 147 -t Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA.bed -n naive_broad --broad

macs3 callpeak -t Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA_chr19.bed -g hs -f BED --nomodel --shift -72 --extsize 147 --broad -n chr19_broad


macs3 pileup -f BED -i Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA.bed --extsize 25 -o RADICL_DNA_track.bdg

macs3 pileup -f BED -i Set18-Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA_chr19.bed --extsize 1000000 -o RADICL_DNA_track.bdg

macs3 pileup -f BED -i BCL3_GM12878_ENCFF629LEU_sorted_chr19.bed --extsize 500000  -B -o BCL3_GM12878_ENCFF629LEU_sorted_chr19.bdg
## Smmothing signal
macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/RADICL_iPSC_RNA_chr19.bed --extsize 12500 -B -o RADICL_RNA_smooth_chr19.bdg

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 12500 -B -o RADICL_DNA_smooth_chr19.bdg

#-------------------------------------
macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/chr19_DNA_singleton.bed --extsize 5000 -B -o RADICL_DNA_singleton_bg_chr19.bdg

macs3 bdgcmp -t RADICL_DNA_smooth_chr19.bdg -c RADICL_DNA_singleton_bg_chr19.bdg -m qpois -o RADICL_DNA_qval_chr19.bdg
## Smmothing signal
# 147/1e3 = 0.147
macs3 pileup -f BED -i Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA.bed --extsize 500 -o smooth_RADICL_DNA_track.bdg
macs3 bdgopt -i smooth_RADICL_DNA_track.bdg -m multiply -p 0.147 -o smooth_RADICL_DNA_track_norm.bdg

###############################################
samtools sort -o ff_iPSC_sorted.bam AThi10063_190416_D00670_0187_BCDLTDANXX_ACCACTGT_L007_R1_001.F6-002-DNA-033-00018-CE13.NoBarcode.31000.bam 
samtools index ff_iPSC_sorted.bam 
samtools view -b ff_iPSC_sorted.bam chr19 > ff_iPSC_ATAC_sorted_chr19.bam
~/bedtools2.3 bamtobed -bed12 -i ff_iPSC_ATAC_sorted_chr19.bam > ff_iPSC_ATAC_sorted_chr19.bed
macs3 pileup -f BED -i ff_iPSC_ATAC_sorted_chr19.bed --extsize 25000 -o ff_iPSC_ATAC_sorted_chr19.bdg
macs3 pileup -f BED -i Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001_DNA_chr19.bed --extsize 25000 -o Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/RADICL_DNA_track_chr19.bdg

macs3 callpeak -t ff_iPSC_ATAC_sorted_chr22.bed -g hs -f BED --nomodel --shift -72 --extsize 147 --broad -n chr22_ATAC

macs3 predictd -i ff_iPSC_ATAC_sorted_chr19.bed -g hs -m 5 50
###############################################
# Modified script to call peaks with custom background
macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 147 -o RADICL_DNA_smooth_chr19.bdg

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/chr19_DNA_singleton.bed --extsize 147 -o RADICL_DNA_singleton_bg_chr19.bdg

macs3 bdgcmp -t RADICL_DNA_smooth_chr19.bdg -c RADICL_DNA_singleton_bg_chr19.bdg -m qpois -p 1 -o RADICL_DNA_qval_chr19.bdg

# to adjust until removal of centromeric region?
# examine qval track and focus on centromeric region qval distribution!
macs3 bdgpeakcall -i RADICL_DNA_qval_chr19.bdg -c 5 -l 75 -g 147 -o RADICL_DNA_chr19_peaks.bed

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 12500 -B -o RADICL_DNA_smooth_chr19.bdg
