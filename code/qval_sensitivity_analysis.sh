macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 147 -o RADICL_DNA_smooth_chr19.bdg

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/chr19_DNA_singleton.bed --extsize 147 -o RADICL_DNA_singleton_bg_chr19.bdg

macs3 bdgcmp -t RADICL_DNA_smooth_chr19.bdg -c RADICL_DNA_singleton_bg_chr19.bdg -m qpois -p 1 -o RADICL_DNA_qval_chr19.bdg


for thresh in {1,2,3,4,5,6,10}
do
echo ${thresh}
macs3 bdgpeakcall -i RADICL_DNA_qval_chr19.bdg -c ${thresh}  -l 75 -g 147 -o RADICL_DNA_${thresh}_chr22_peaks.bed
done

