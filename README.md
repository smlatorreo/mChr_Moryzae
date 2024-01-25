# Data and code for the manuscript: "Multiple Horizontal Mini-Chromosome Transfers Drive Genome Evolution of Clonal Blast Fungus Lineages"
DOI: XX.XX.XX

## Software requirements
Program                  | Location
------------------------ | ----------------------------
*AdapterRemoval2 v.2*    | (https://github.com/mikkelschubert/adapterremoval)
*Bwa-mem2 v.2.1*         | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*        | (https://github.com/samtools/samtools)
*sambamba v0.8.0*        | (https://github.com/biod/sambamba)
*GATK v.4.2*             | (https://github.com/broadinstitute/gatk/releases)
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)
*IQ-Tree v.2*            | (https://github.com/iqtree/iqtree2)
*TreeTime*               | (https://github.com/neherlab/treetime)
*Dstat.py*               | (https://github.com/smlatorreo/Dstats)
*tped2fasta.sh*          | (https://github.com/smlatorreo/misc_tools/blob/main/tped2fasta.sh)
*popstats*               | (https://github.com/pontussk/popstats)
*BEAST2 v.2.7.6*         | (https://github.com/CompEvol/beast2)

## Preprocessing and mapping of short reads to the rice-infecting *M. oryzae* reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the rice-infecting *Magnaporthe oryzae* AG006 assembly as the reference genome and indexed this genome using *Bwa-mem2*.
```bash
bwa index AG006.fa
```

*BWA mem2* was used to map the trimmed reads to the reference genome, *samtools* to discard non-mapped reads, and *sambamba* to sort and mark PCR optical duplicates.
```bash
bwa-mem2 mem -R "@RG\tID:$sample\tSM:$sample" AG006.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R2.fastq.gz > sample1.sam
samtools view -SbhF 4 > sample1_mapped.bam
sambamba sort -o sample1_mapped_sorted.bam sample1_mapped.bam
sambamba markdup sample1_mapped_sorted.bam sample1_mapped_sorted.dd.bam
```

Depth and coverage statistics were computed using samtools coverage
```bash
# Average depth
samtools coverage sample1_mapped_sorted.dd.bam | awk '{sum += $3}END{print sum / NR}'
```
A summary of these statistics was organized in a [table](/data/preprocess_and_SNPs/summary_coverage.tsv)


## Variant calling
We used the *HaplotypeCaller* from *GATK* to generate genomic haplotype calls per individual using the duplicate-marked BAM file as input.
```bash
gatk HaplotypeCaller -R AG006.fa -I sample1_mapped_sorted.dd.bam -O sample1.g.vcf.gz
```

We used *CombineGVCFs*, *GenotypeGVCFs* and *SelectVariants* from *GATK* to combine the individual genomic VCFs, call genotypes and filter SNPs, respectively.
```bash
gatk CombineGVCFs -R AG006.fa -V sample1.g.vcf.gz -V sample2.g.vcf.gz -V sampleN.g.vcf.gz -O ALL.g.vcf.gz
gatk GenotypeGVCFs -R AG006.fa -ploidy 1 -V ALL.g.vcf.gz -O ALL.raw.vcf.gz
gatk SelectVariants -select-type SNP -V ALL.raw.vcf.gz -O ALL.raw.snps.vcf.gz
```

We extracted all Quality-by-Depth (QD) values
```bash
bcftools view -H ALL.raw.snps.vcf.gz | cut -f8 | \
awk -F "QD=" '{print $2}' | cut -f1 -d ";" | gzip >  ALL.raw.snps.QD.gz
```

Based on the distribution of Quality-by-Depth values, we set filters of one standard deviation around the median value.
```python
# Python
import pandas as pd
QD = pd.read_csv('ALL.raw.snps.QD.gz', header = None, compression = 'gzip')
med = QD.median()
lower = med - QD.std()
upper = med + QD.std()
print(lower, upper)
```

Finally, using the above-mentioned scheme, we filtered SNPs using *GATK VariantFiltration* and created a new VCF file, keeping non-missing positions, using *bcftools*.
```bash
gatk VariantFiltration --filter-name "QD" \
--filter-expression "QD <= $lower || QD >= $upper" \
-V wheat-blast.raw.snps.QD.gz \
-O wheat-blast.snps.filter.vcf.gz

bcftools view -g ^miss wALL.snps.filter.vcf.gz | bgzip > ALL.snps.filtered.vcf.gz
```

The filtered [VCF file can be found here](/data/preprocess_and_SNPs/ALL.snps.filtered.vcf.gz)


## Bayesian-based dated phyloenetic reconstruction
We submitted 6 independent chains to [CIPRES](doi: 10.1109/GCE.2010.5676129).  
Here, we provide the [submitted configuration .xlm file](/data/BEAST2/Rice_Outgroup.Fullinfo.HYK.LogNormPriorAdj.20M.xml.gz)  
As we use whole genomic SNPs (N = 95044), we rescaled the calculations using the number of invariant A/C/G/T sites using the following tag inside xml file:
`<data id='clonal_rice' spec='FilteredAlignment' filter='-' data='@clonal_rice_Original' constantSiteWeights='8166700 8817938 8831984 8177939'/>`  

The resulting combined [phylogenetic tree can be found here](/data/BEAST2/BEAST2_Rice_Outgroup.Fullinfo.HYK.LogNormPrior.COMBINED6X.MC.tree)

## Using a *mugration* model to ascertain the phylogentic acquisition of the mChrA
We used the GTR-based implementation of [TreeTime](https://doi.org/10.1093/ve/vex042) to ascertain the phylogenetic acquisition of the mChrA in different isolates of rice-infecting *Magnaporthe oryzae*. Here, we provide input and resulting files of this analysis:

File | Description
---- | -----------
[Input tree](/data/mugration_analyses/clonal_rice.withSetariaout.snps.filtered.maxmiss10.rooted.contree) | Input ML-tree generated with IQTree, using a GTR substitution model
[GTR](/data/mugration_analyses/GTR.txt) | Parameters of the GTR-based mugration model
[Annotated tree](/data/mugration_analyses/annotated_tree.nexus) | Resulting annotated tree with mChrA presence/absence annotation
[Confidence](/data/mugration_analyses/confidence.csv) | Per-branch mChrA ascertainment confidence

## Genetic distances between Rice-infecting and Eleusine-infecting isolates
We measured the genetic Hamming distances between rice-infecting isolates and the Eleusine-infecting isolate Br62 at the [whole core chromosome](/data/distances/AG006_coreChr.bed) and in [100 randomly selected regions](/data/distances/ALL.random_regions.bed) as follows:

```bash
vcf=$1
sample=$2

core_distance = $(bcftools view -R AG006_mChrA.bed -a -s Br62,$sample $vcf | bcftools view -m2 -M2 -i "(AN == 2)&&(AC != AN)" -H | awk '$5 != "*"' | wc -l)
echo -e "Core_distance\t$core_distance"

for i in {1..100}; do
    d = $(bcftools view -R random_regions_bed/r$i.bed -a -s Br62,$sample $vcf | bcftools view -m2 -M2 -i "(AN == 2)&&(AC != AN)" -H | awk '$5 != "*"' | wc -l)
    echo -e "Radom_region_$i_distance\t$d"
done

```

[A summary table with the distances can be found here](/data/distances/summary_distances_to_Br62.tsv)

## Comparing mChrA tree phylogeny with random selected tree topologies
We subset two VCF files with all the putative carriers of the mChrA. [The first VCF file contains SNPs belonging to the core chromosome regions](/data/random_trees/mChrA_samples.coreChr.snps.filtered.maxmiss10.vcf.gz), while [the second VCF contains SNPs only from the mChrA region](data/random_trees/mChrA_samples.mChrA.snps.filtered.maxmiss10.vcf.gz)

Using plink and the bash script [tped2fasta.sh](https://github.com/smlatorreo/misc_tools/blob/main/tped2fasta.sh), we transformed the VCF files into fasta-like format files.  

We then computed NJ trees for the two subsets using IQ-TREE
```bash
iqtree -s SNPs.fasta -m GTR -fast
```

Output trees from the [core chromosome](/data/random_trees/mChrA_samples.coreChr.snps.filtered.maxmiss10.fasta.treefile) and the [mChrA](/data/random_trees/mChrA_samples.mChrA.snps.filtered.maxmiss10.fasta.bionj) are provided.  

We then selected [100 random regions](/data/distances/ALL.random_regions.bed) from the core chromosome region and calculated NJ trees as decribed. We used [this bash script](/scripts/generate_random_trees.sh) to automate the process.  

Finally, we assessed the monophyly of the host-infecting isolates to count for the number of instances where the topology of the estimated trees resemble the one of the whole core chromosome. This step was done using the custom python script [check_monophyl.py](/scripts/check_monophyl.py)

## Differentiating between sexual mating and horizontal genetic transfer through *D* statistics
We subset genomic SNP information from all the mChrA putative carrier isolates, together with random selected isolates as controls. We also included the *Digitaria*-infecting isolate Dig41 to be used as an outgroup. [The resulting VCF file is provided here](/data/dstats/subset.coreChr.snps.filtered.vcf.gz)  

We used popsats as well as a [*Dstat.py*](https://github.com/smlatorreo/Dstats/) to compute the Patterson's D statistic using the following configurations:  

[\(Dig41, Br62 ; Isolate_without_mChrA, Isolate_with_mChrA\)](/data/dstats/Dtest_configurations.txt)  
[\(Dig41, Br62 ; Isolate_without_mChrA, Isolate_without_mChrA\)](/data/dstats/Dtest_configurations.control_without_mChrA.txt)  
[\(Dig41, Br62 ; Isolate_with_mChrA, Isolate_with_mChrA\)](/data/dstats/Dtest_configurations.control_with_mChrA.txt)



