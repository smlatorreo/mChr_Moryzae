import sequencetools # https://github.com/smlatorreo/modules_python/blob/master/sequencetools.py
import vcfpytools_class # https://github.com/smlatorreo/modules_python/blob/master/vcfpytools_class.py
from itertools import combinations
from sys import argv, stderr

vcf = argv[1]
list_samples = argv[2]
OUT = 'Br62'

VCF = vcfpytools_class.vcf_object(vcf)

samples = [i.strip() for i in open(list_samples, 'r').readlines()]
genotypes = {id:[] for id in samples}
for record in VCF.get_genotypes_hap(samples=samples, filtered=False):
    alleles = record['Genotypes']
    if alleles.count('.') <= (len(samples) - 2):
        [genotypes[sample_id].append(alleles[i]) for i, sample_id in enumerate(genotypes)]

print('Distances to ', OUT, file = stderr)
for s in samples:
    if s != OUT:
        d = sequencetools.kimura2parameter(genotypes[OUT], genotypes[s])
        print(OUT, s, d, sep = '\t')
