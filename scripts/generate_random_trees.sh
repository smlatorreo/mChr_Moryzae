VCF=$1
BED=$2
out=$(basename $BED | sed 's/\.bed//g')

bcftools view -R $BED $VCF | bgzip > $out.tmp.vcf.gz

plink --aec --vcf $out.tmp.vcf.gz --recode transpose --missing-genotype "?" --out $out.tmp
bash tped2fasta.sh $out.tmp > $out.tmp.fasta

rm $out.tmp.log $out.tmp.nosex $out.tmp.tfam $out.tmp.tped $out.tmp.vcf.gz

iqtree -s $out.tmp.fasta -m GTR -fast
