# USAGE:
# python simulate <P1.fa> <P2.fa> <Test.fa> <Out.fa> <sexP> <LD50> <jackknife_block_length> <generations>

#----------------------------- LIBRARIES ---------------------------------------------------------------
import numpy as np
import random
from Bio import SeqIO
from sys import stderr, argv


#------------------------------ FUNCTIONS --------------------------------------------------------------

def recomb(P1, P2, LD50):
    '''
    Function to simulate recombination after a mating event between two parental
    organisms (P1) and (P2), given a (LD50) value.
    The LD decay is modelled as an exponential decay distribution.
    The LD50 value corresponds to the median of the exponential distribution.
    Hence, the scale value (mean of the distribution) is equivalent to LD50 / log(2)
    '''
    scale = LD50 / np.log(2)
    # Convert the seuquences (str) into lists for better handling
    P1 = list(P1)
    P2 = list(P2)

    # Genome length
    genome_length = len(P1)

    # Create an exponential distribution with a big enough size to samples from it
    size_dist = int((genome_length*10) / scale)
    random_exp_dist = [int(i) for i in np.random.exponential(scale = scale, size = size_dist)]
    # Estimate the amount and length of required blocks to achieve half of the genome size
    cumul = np.cumsum(random_exp_dist)
    stop = len(cumul[cumul <= (genome_length / 2)])
    block_lengths = random_exp_dist[:stop]

    # Create a mosaic of non-recombining and recombining blocks with the sampled sizes
    no_recomb = [None] * (genome_length - sum(block_lengths))
    T = block_lengths + no_recomb
    random.shuffle(T)

    # Extract start and end coordinates for the recombining blocks
    recomb_blocks = []
    count_index = 0
    for i in range(len(T)):
        if T[i] == None:
            count_index += 1
        else:
            recomb_blocks.append(range(count_index, count_index  + T[i]))
            count_index += T[i]
    # Fill the recombining blocks with P2 information
    recombinant = P1
    for b in recomb_blocks:
        recombinant[b.start:b.stop] = P2[b.start:b.stop]
    return ''.join(recombinant)

def drift(seq, rate=7E-8):
    '''
    Function to add new mutations into an existing seuquence.
    The number of new mutations are sampled from a Poisson distribution with a lambda equal to a mutation rate.
    The units of the rate are expressed as mutations per nucleotide base per generation (assuming 1 generation per year).
    '''
    # Dictionary to flip the allele states
    mutdict = {'1':'0', '0':'1', '.':'.'}
    # Genome length
    genome_size = len(seq)
    # Estimate the number of mutations per generation
    N_mutations = np.random.poisson(rate * genome_size)
    # Randomly apply the mutations
    mutations = np.random.randint(0,genome_size, N_mutations)
    seq = list(seq)
    for m in mutations:
        try:
            seq[m] = mutdict[seq[m]]
        except KeyError: # In case of different characters
            continue
    return ''.join(seq)
    
    
def Dstat(ABBA, BABA):
    '''
    Paterson's D statistic
    '''
    try:
        D = (ABBA - BABA) / (ABBA + BABA)
    except ZeroDivisionError:
        D = np.nan
    return D


def Dstats(A, B, X, OUT, jack_block = 100000):
    '''
    Function to compute the Patterson's D statistic given sequences in a 4-taxa configuration:
    D: ((A,B), X), OUT
    '''
    # Dictionary with possible ABBA / BABA configurations
    configs = {'1001':'ABBA','0110':'ABBA',
               '1010':'BABA','0101':'BABA'}

    # Store the number and location of ABBA/BABA sites
    counts = {'ABBA':[], 'BABA':[]}
    total = 0
    for i in range(len(Out)):
        alleles = [s[i] for s in (A, B, X, OUT)]
        if '.' not in alleles: # Only sites with full information
            total += 1
            alleles = ''.join(alleles)
            if alleles in configs:
                counts[configs[alleles]].append(i) # Log the position of the ABBA / BABA site
    # Convert lists into numpy arrays for faster and better handling
    counts = {'ABBA':np.array(counts['ABBA']), 'BABA':np.array(counts['BABA'])}
    # Estimate global D
    D = Dstat(len(counts['ABBA']), len(counts['BABA']))

    # Jackknife
    # Number of jackknife blocks
    T = len(range(0,len(Out), jack_block))
    # Compute pseudovalues of D by removing one block every time
    pseudoDs = []
    for bl in range(0,len(Out), jack_block):
        A = len(counts['ABBA'][(counts['ABBA'] <= bl)|(counts['ABBA'] >= (bl + jack_block))])
        B = len(counts['BABA'][(counts['BABA'] <= bl)|(counts['BABA'] >= (bl + jack_block))])
        pseudoD = (D * T) - ((Dstat(A,B) * (T-1)))
        pseudoDs.append(pseudoD)
    # Compute standard error and Z scores
    D_err = np.std(pseudoDs) / np.sqrt(len(pseudoDs))
    D_Z = D / D_err
    # Output dictionary
    out = {'ABBA':len(counts['ABBA']),'BABA':len(counts['BABA']),'Total':total,
           'D':D, 'N_blocks':len(pseudoDs), 'Std_err':D_err, 'Z':D_Z}
    return out
    
#---------------------------- PARAMETERS --------------------------------------------------------
# Load individual sequences as single contigs
# The paths from every individual are taken from the command line positions
P1 = []
for record in SeqIO.parse(argv[1], 'fasta'):
    P1 = P1 + list(record.seq)
P2 = []
for record in SeqIO.parse(argv[2], 'fasta'):   
    P2 = P2 + list(record.seq)
Test = []
for record in SeqIO.parse(argv[3], 'fasta'):   
    Test = Test + list(record.seq)
Out = []
for record in SeqIO.parse(argv[4], 'fasta'):   
    Out = Out + list(record.seq)

# Catch parameters from the command line positions
sexP = float(argv[5])
LD50 = int(argv[6])
block_len = int(argv[7])
generations = int(argv[8])

#------------------------------ SIMULATION ------------------------------------------------------
print('State\tABBA\tBABA\tTotal\tD\tN_blocks\tStd_Err\tZ', flush=True)

# Compute D at the original state
print('State0', *Dstats(P2, P1, Test, Out, block_len).values(), sep = '\t', flush=True)

# Single pulse of introgression and compute D
F1 = recomb(Test, P1, LD50)
F1P = recomb(Test, P1, LD50)
print('F1', *Dstats(P2, F1, Test, Out, block_len).values(), sep = '\t', flush=True)

# Simulate recurrent crosses between offscript and either a the parental individual (backcross) or another mate from the same offspring (offspringP) for N generations given a sex probability P
# Compute D after each generation
offspring = F1[:]
offspringP = F1P[:]
for gen in range(generations):
    print('#Generation', gen, file = stderr)
    sex = bool(np.random.binomial(n = 1, p = sexP))
    if sex == True:
        mate = random.choice(['P1','offspringP'])
        offspring = recomb(offspring, globals()[mate], LD50)
        if mate == 'offspringP':
            offspringP = recomb(offspring, offspringP, LD50)
    offspring = drift(offspring)
    offspringP = drift(offspringP)
    print('Back-cross{}'.format(gen+1), *Dstats(P2, offspring, Test, Out, block_len).values(), sep='\t', flush=True)

print('#sexP={};LD50={};block_len={}'.format(sexP,LD50,block_len), file=stderr)
print('#Last state. Generation:', gen+1, file=stderr)
print('>offspring{}\n{}'.format(gen+1, offspring), file=stderr)
print('>offspringP{}\n{}'.format(gen+1, offspringP), file=stderr)
