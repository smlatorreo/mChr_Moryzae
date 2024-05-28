from sys import argv
import vcfpytools_class
from Bio import SeqIO

vcf = argv[1]
VCF = vcfpytools_class.vcf_object(vcf)
ref = argv[2]

for sample in range(6):
    REF = {}
    for record in SeqIO.parse(ref, 'fasta'):
        REF[record.id] = list(record.seq)

    for record in VCF.get_genotypes_hap(binary=True):
        CHR = record['Position'][0]
        POS = record['Position'][1] - 1 # Python
        Allele = record['Genotypes'][sample]
        REF[CHR][POS] = Allele

    with open('sample_{}.fa'.format(sample), 'w') as f:
        for chr in REF:
            print('>', chr, sep = '', file = f)
            seq = ''.join(REF[chr]).upper().replace('A', '0').replace('C','0').replace('G','0').replace('T','0')
            print(seq, file = f)

