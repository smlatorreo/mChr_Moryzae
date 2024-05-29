[![DOI](https://zenodo.org/badge/747116432.svg)](https://zenodo.org/doi/10.5281/zenodo.10628812)

# Data and code for the manuscript: "Multiple horizontal mini-chromosome transfers drive genome evolution of clonal blast fungus lineages"

## Software requirements
Program                     | Location
--------------------------- | ----------------------------
*AdapterRemoval2 v.2*       | (https://github.com/mikkelschubert/adapterremoval)
*Bwa-mem2 v.2.1*            | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*           | (https://github.com/samtools/samtools)
*sambamba v0.8.0*           | (https://github.com/biod/sambamba)
*GATK v.4.2*                | (https://github.com/broadinstitute/gatk/releases)
*bcftools v.1.11*           | (https://github.com/samtools/bcftools)
*IQ-Tree v.2*               | (https://github.com/iqtree/iqtree2)
*K2P_from_VCF.py*           | (https://github.com/smlatorreo/mChr_Moryzae/blob/main/scripts/K2P_from_VCF.py)
*TreeTime*                  | (https://github.com/neherlab/treetime)
*Dstat.py*                  | (https://github.com/smlatorreo/Dstats)
*tped2fasta.sh*             | (https://github.com/smlatorreo/misc_tools/blob/main/tped2fasta.sh)
*popstats*                  | (https://github.com/pontussk/popstats)
*BEAST2 v.2.7.6*            | (https://github.com/CompEvol/beast2)
*from_VCF_to_bin_fasta.py*  | (https://github.com/smlatorreo/mChr_Moryzae/blob/main/scripts/simulations/from_VCF_to_bin_fasta.py)
*simulate_D.py*             | (https://github.com/smlatorreo/mChr_Moryzae/blob/main/scripts/simulations/simulate_D.py)
*simulate_D_mate_choice.py* | (https://github.com/smlatorreo/mChr_Moryzae/blob/main/scripts/simulations/simulate_D_mate_choice.py)

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
We submitted 6 independent chains to [CIPRES](https://www.phylo.org/).  
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
We measured Kimura two-parameter (K2P) distances in [contiguous, non-overlapping 100bp-sized windows between mChrA](/data/distances/mChrA_regions_100Kb.bed) in all [Oryza-infecting isolates carrying this sequence (n=32) and Eleusine-infecting isolate Br62](/data/distances/samples_mChrA.list). We used a the python script [K2P_from_VCF.py](/scripts/K2P_from_VCF.py) to compute the distances:

```bash
# Iterate through the bed regions
while read region; do
	# Fields from a bed-format file
	contig=$(echo $region | cut -f1 -d " ")
	start=$(echo $region | cut -f2 -d " ")
	end=$(echo $region | cut -f3 -d " ")

	# Subset and index new temporary VCF for every bed region
	bcftools view -S samples_mChrA.list -r $contig\:$start\-$end ALL.snps.filtered.vcf.gz | bgzip > $contig.$start.$end.tmp.vcf.gz
	tabix -p vcf $contig.$start.$end.tmp.vcf.gz

	# Compute Kimura two-parameter distances
	python K2P_from_VCF.py $contig.$start.$end.tmp.vcf.gz samples_mChrA.list > $contig.$start.$end.K2P.dist

	# Remove temporary VCF
	rm $contig.$start.$end.tmp.vcf.gz*
done < mChrA_regions_100Kb.bed
```
A [summary of the distances can be found here](/data/distances/summary_distances_mChrA_windows100kb_kimura-2-param_to_Br62.tsv)  

Similarly, we also measured K2P distances across [100 non-overlapping and randomly sampled 100kb core chromosomal regions](/data/distances/coreChr_random_regions_100Kb.bed) between each isolate and Br62. The output [summary can be found here](/data/distances/summary_distances_coreChr_windows100kb_kimura-2-param_to_Br62.tsv)

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


## Simulations to estimate the D-statistic detection power under two scenarios
Our *D* statistic analyses aimed to investigate the likely origin of mChrA. We hypothesized that within a 4-taxa configuration `Dig41, (Br62, (Isolate_without_mChrA, Isolate_with_mChrA))`, a Patterson *D* value  > 0 would indicate a likely case of introgression as the cause of mChrA acquisition. Conversely, a *D* statistic of 0 would suggest a scenario of horizontal gene transfer.  

Howeverm, to further understand the implications of ancestral introgression followed by multiple generations of backcrossing, we examined the detection power of our hypothesis considering: i) The probability of sexual reproduction per generation; ii) The number of generations; and iii) Linkage Disequilibrium.  

We simulated two scenarios, modeling a single pulse of introgression followed by multiple generations of backcrosses, and measured the resulting *D* statistic. Below is a brief description of our process.

Based on the genome-wide SNP data contained in the VCf file, we used [*from_VCF_to_bin_fasta.py*](/scripts/simulations/from_VCF_to_bin_fasta.py) to codify invariable and variable positions for each isolate in respect to the reference genome. As a result, we generated fasta-like files where reference-type alleles are codified as 0's and alternative-type alleles are codified as 1's. The files can be found here: [Dig41](/data/simulations/Dig41.fa.gz), [Br62](/data/simulations/Br62.fa.gz), [FJ12JN-084-3](/data/simulations/FJ12JN-084-3.fa.gz), [658](/data/simulations/658.fa.gz)

We measured a Patterson's *D* statistic of 0 (Z-score = -0.22) in the following 4-taxa configuration:
`Dig41, (Br62, (FJ12JN-084-3, 658))`

We then carried out simulations under two scenarios:

### Scenario 1. Backcrosses

We then recreated a scenario of a single pulse of introgression from the Eleusine-infecting isolate Br62 into a rice-infecting isolate (658). We again meassured a Patterson's D statistic of ~1 (Z-score >> 3) in this hypothetical F1 individual.  

We then recreated a back-cross between the F1 individual and the original 658 isolate. The backcross was modelled as the binomial probability (*P*) for the cross to happen. After each generation, a poisson distribution with a *lambda* of 7E-8  * genome_size was used to add new random mutations.

We logged Patterson's *D* statistics after each generation and estimated the number of generations at which D is statistically equal to 0.

```bash
LD50=35000
jackknife_block=5000000
generations=1000

for sex_prob in 0.0{1..9} 0.{1..9} 1.0; do
    python simulate_D.py 658.fa FJ12JN-084-3.fa Br62.fa Dig41.fa $sex_prob $LD50 $jackknife_block $generations
done
```

All results can be found at [simulations_output.tar.gz](/data/simulations/simulations_output.tar.gz)


### Scenario 2. Allowing mate choice

In Scenario 1, each generation is allowed to backcross with the parental rice-infecting isolate with a probability *P*. Here, we allow a process of choosing a mate. This process is governed by a binomial probability *B*. Specifically, if an individual does not backcross with the parental isolate, it will instead mate with another individual from the same generation.

Again, we logged Patterson's *D* statistics after each generation and estimated the number of generations at which D is statistically equal to 0.

```bash
LD50=35000
jackknife_block=5000000
generations=1000

for sex_prob in 0.0{1..9} 0.{1..9} 1.0; do
    python simulate_D_mate_choice.py 658.fa FJ12JN-084-3.fa Br62.fa Dig41.fa $sex_prob $LD50 $jackknife_block $generations
done
```

All results can be found at [simulations_with_mate_choice_output.tar.gz](/data/simulations/simulations_with_mate_choice_output.tar.gz)

