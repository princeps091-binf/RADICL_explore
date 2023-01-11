# 

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 147 -o RADICL_DNA_smooth_chr19.bdg

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/chr19_DNA_singleton.bed --extsize 147 -o RADICL_DNA_singleton_bg_chr19.bdg

macs3 bdgcmp -t RADICL_DNA_smooth_chr19.bdg -c RADICL_DNA_singleton_bg_chr19.bdg -m qpois -p 1 -o RADICL_DNA_qval_chr19.bdg

# to adjust until removal of centromeric region?
# examine qval track and focus on centromeric region qval distribution!
macs3 bdgpeakcall -i RADICL_DNA_qval_chr19.bdg -c 5 -l 75 -g 147 -o RADICL_DNA_chr19_peaks.bed

for i in {1..22}
do
chrom="chr$i"
echo $chrom

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_${chrom}.bed --extsize 147 -o RADICL_DNA_smooth_${chrom}.bdg

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/${chrom}_DNA_singleton.bed --extsize 147 -o RADICL_DNA_singleton_bg_${chrom}.bdg

macs3 bdgcmp -t RADICL_DNA_smooth_${chrom}.bdg -c RADICL_DNA_singleton_bg_${chrom}.bdg -m qpois -p 1 -o RADICL_DNA_qval_${chrom}.bdg

# to adjust until removal of centromeric region?
# examine qval track and focus on centromeric region qval distribution!
macs3 bdgpeakcall -i RADICL_DNA_qval_${chrom}.bdg -c 5 -l 75 -g 147 -o RADICL_DNA_${chrom}_peaks.bed
rm RADICL_DNA_singleton_bg_${chrom}.bdg RADICL_DNA_smooth_${chrom}.bdg
done
