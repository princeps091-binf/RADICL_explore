# script to explore reads patterning
cut -d$'\t' -f7-11 15.AGTTCC.2.bed > RADICL_iPSC_DNA.bed
cut -d$'\t' -f1-6 15.AGTTCC.2.bed > RADICL_iPSC_RNA.bed
## Subset chromosome specific RADICL data
awk '$1 == "chr19"' RADICL_iPSC_DNA.bed > RADICL_iPSC_DNA_chr19.bed
awk '$1 == "chr19"' RADICL_iPSC_RNA.bed > RADICL_iPSC_RNA_chr19.bed

## Subset chromosome specific RADICL data
~/bedtools2.3 sort -i RADICL_iPSC_chr19_DNA.bed > RADICL_iPSC_chr19_DNA_sorted.bed
~/bedtools2.3 cluster -i RADICL_iPSC_chr19_DNA_sorted.bed > RADICL_iPSC_chr19_DNA_cluster.txt
~/bedtools2.3 merge -i RADICL_iPSC_chr19_DNA_sorted.bed > RADICL_iPSC_chr19_DNA_merged.bed
~/bedtools2.3 closest -io -t "all" -d -a RADICL_iPSC_chr19_DNA_merged.bed -b RADICL_iPSC_chr19_DNA_merged.bed > RADICL_iPSC_chr19_DNA_cluster_neighbour.tsv
## ATAC
~/bedtools2.3 cluster -i ff_iPSC_ATAC_sorted_chr19.bed > ff_iPSC_ATAC_sorted_chr19.txt

~/bedtools2.3 merge -i RADICL_iPSC_chr19_DNA_sorted.bed > RADICL_iPSC_chr19_DNA_merged.bed
~/bedtools2.3 closest -io -t "all" -d -a RADICL_iPSC_chr19_DNA_merged.bed -b RADICL_iPSC_chr19_DNA_merged.bed > RADICL_iPSC_chr19_DNA_cluster_neighbour.tsv

########################
for i in {1..22}
do
chrom="chr$i"
echo $chrom
awk -v chr=$chrom '$1 == chr' RADICL_iPSC_RNA.bed > ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RNA/RADICL_iPSC_RNA_$chrom.bed
done

for i in {1..22}
do
chrom="chr$i"
echo $chrom
awk -v chr=$chrom '$1 == chr' RADICL_iPSC_DNA.bed > ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_$chrom.bed
done

for i in {1..22}
do
chrom="chr$i"
echo $chrom
awk -v chr=$chrom '$7 == chr' 15.AGTTCC.2.bed > ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/RADICL_iPSC_$chrom.bed
done

# Ran in DNA and RNA folders
for f in ./*
do
echo $f
echo sorting
~/bedtools2.3 sort -i $f > "$f"_sorted.bed
echo cluster
~/bedtools2.3 cluster -i "$f"_sorted.bed > ./cluster/"$f"_cluster.txt
rm "$f"_sorted.bed
done
#--------------------
# Produce peak with custom background
macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 75 -B -o RADICL_DNA_smooth_chr19.bdg

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/cluster/singleton_data/chr19_DNA_singleton.bed --extsize 75 -B -o RADICL_DNA_singleton_bg_chr19.bdg

macs3 bdgcmp -t RADICL_DNA_smooth_chr19.bdg -c RADICL_DNA_singleton_bg_chr19.bdg -m qpois -p 1 -o RADICL_DNA_qval_chr19.bdg

# to adjust until removal of centromeric region?
# examine qval track and focus on centromeric region qval distribution!
macs3 bdgpeakcall -i RADICL_DNA_qval_chr19.bdg -c 5 -l 75 -g 147 -o RADICL_DNA_chr19_peaks.bed

macs3 pileup -f BED -i ~/Documents/FANTOM6/data/RADICL/Set18-7_iPSC_rep2_AGTTCC_S0_L002_R1_001.bed/chr_data/DNA/RADICL_iPSC_DNA_chr19.bed --extsize 12500 -B -o RADICL_DNA_smooth_chr19.bdg
