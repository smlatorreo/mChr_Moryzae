# Some scripts and data for the manuscript: "Multiple Horizontal Mini-Chromosome Transfers Drive Genome Evolution of Clonal Blast Fungus Lineages"
DOI: XX.XX.XX

## Bayesian-based dated phyloenetic reconstruction
We submitted 6 independent chains to [CIPRES](doi: 10.1109/GCE.2010.5676129).  
Here, we provide the [submitted configuration .xlm file](/data/BEAST2/Rice_Outgroup.Fullinfo.HYK.LogNormPriorAdj.20M.xml.gz)  
As we use whole genomic SNPs (N = 95044), we rescaled the calculations using the number of invariant A/C/G/T sites using the following tag inside xml file:
`<data id='clonal_rice' spec='FilteredAlignment' filter='-' data='@clonal_rice_Original' constantSiteWeights='8166700 8817938 8831984 8177939'/>`  

The resulting combined [phylogenetic tree can be found here](/data/BEAST2/BEAST2_Rice_Outgroup.Fullinfo.HYK.LogNormPrior.COMBINED6X.MC.tree)

## Using a mugration model to ascertain the phylogentic acquisition of the mChrA
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

Output trees from the [core chromosome](/data/random_trees/mChrA_samples.new.coreChr.snps.filtered.maxmiss10.fasta.treefile) and the [mChrA](/data/random_trees/mChrA_samples.new.mChrA.snps.filtered.maxmiss10.fasta.bionj) are provided.  

We then selected [100 random regions](/data/distances/ALL.random_regions.bed) from the core chromosome region and calculated NJ trees as decribed. We used [this bash script](/scripts/generate_random_trees.sh) to automate the process.  

Finally, we assessed the monophyly of the host-infecting isolates to count for the number of instances where the topology of the estimated trees resemble the one of the whole core chromosome. This step was done using the custom python script [check_monophyl.py](/scripts/check_monophyl.py)

## Differentiating between sexual mating and horizontal genetic transfer through *D* statistics
We subset genomic SNP information from all the mChrA putative carrier isolates, together with random selected isolates as controls. We also included the *Digitaria*-infecting isolate Dig41 to be used as an outgroup. [The resulting VCF file is provided here](/data/dstats/subset.coreChr.snps.filtered.vcf.gz)  

We used popsats as well as a [*Dstat.py*](/scripts/Dstats/Dstat.py) to compute the Patterson's D statistic using the following configurations:  

[\(Dig41, Br62 ; Isolate_without_mChrA, Isolate_with_mChrA\)](/data/dstats/Dtest_configurations.txt)
[\(Dig41, Br62 ; Isolate_without_mChrA, Isolate_without_mChrA\)](/data/dstats/Dtest_configurations.control_without_mChrA.txt)
[\(Dig41, Br62 ; Isolate_with_mChrA, Isolate_with_mChrA\)](/data/dstats/Dtest_configurations.control_with_mChrA.txt)



