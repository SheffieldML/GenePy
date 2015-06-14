import collections
from pandas import DataFrame
import numpy as np
from pandas import concat
import RNA
import random
#from Bio.Graphics.GenomeDiagram._FeatureSet import FeatureSet
#from gi.overrides.GObject import features



## -----------------
## FUNCTIONS TO PRE-PROCESS THE FASTA FILES
## -----------------

def process_fasta(name_origin,destination_name):
    file_read = open(name_origin,'r')
    file_write = open(destination_name,'w')
    file_write.write('ID gene_sequence')
    line_r = file_read.readline()
    line = line_r.replace("\n", "")

    while line!='':
        if list(line)[0]=='>':
            file_write.write('\n')
            file_write.write(line[20:38]+ ' ')
   
        if list(line)[0]!='>':	
            line = line.replace("Sequence unavailable","NA")
            file_write.write(line)
  
        line_r = file_read.readline()
        line = line_r.replace("\n", "")

    file_read.close()
    file_write.close()

## -----------------
## FUNCTIONS TO EXTRACT THE FEATURES
## -----------------

## define the table of codons/aminoacids
def create_codon_table():
   bases = ['t', 'c', 'a', 'g']
   codons = [a+b+c for a in bases for b in bases for c in bases]
   amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
   codon_table = dict(zip(codons, amino_acids))
   return(codon_table)

## define the table of aminoacids
def create_aminoacids_table():
  codon_dict = create_codon_table()
  amin_dict = {}
  for k, v in codon_dict.iteritems():
    amin_dict.setdefault(v, []).append(k)
  return(amin_dict)

## generates a random coherent sequence with another one given
def generate_random_sequence(seq):
  amin_dict = create_aminoacids_table()
  codon_dict = create_codon_table()
  trans_seq = translate(seq, codon_dict)
  Namin = len(trans_seq)
  new_seq = ''
  for k in range(Namin):
        redundant_codons = amin_dict[trans_seq[k]]
        Range = range(len(redundant_codons))
        np.random.shuffle(Range)
        new_seq=new_seq +redundant_codons[Range[0]]   
  return(new_seq)

## define the table of codons/aminoacids pairs
def create_codon_pairs_table():
   bases = ['t', 'c', 'a', 'g']
   codons = [a+b+c for a in bases for b in bases for c in bases]
   amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
   codon_pairs = [a+b for a in codons for b in codons]
   amino_acids_pairs = [a+b for a in amino_acids for b in amino_acids]
   codon_pairs_table = dict(zip(codon_pairs, amino_acids_pairs))
   return(codon_pairs_table)

## define the table of codons/aminoacids triplets
def create_codon_triplets_table():
   bases = ['t', 'c', 'a', 'g']
   codons = [a+b+c for a in bases for b in bases for c in bases]
   amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
   codon_triplets = [a+b+c for a in codons for b in codons for c in codons]
   amino_acids_triplets = [a+b+c for a in amino_acids for b in amino_acids for c in amino_acids]
   codon_triplets_table = dict(zip(codon_triplets, amino_acids_triplets))
   return(codon_triplets_table)

## translate gene seqence to aminoacids sequence
def translate(seq, code):
   return "".join((code[seq[i:i+3]] for i in range(0, len(seq)-len(seq)%3, 3)))

## function to counts from a sequence
def count_codons(seq):
  codon_table = create_codon_table()
  counts = collections.defaultdict(int)
  for codon in codon_table.keys(): counts[codon] = 0
  for i in range(0,len(seq)-len(seq)%3,3): codon = seq[i:i+3] ; counts[codon] += 1
  return(counts)

## calculates the number of codons of a sequence
def number_codons(seq):
   number = len(seq)-len(seq)%3
   return number

## fucntion to find lists within other lists
find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]

def geomean(num_list):
   return sum(num_list) ** (1.0/len(num_list))

## codon adaptation index
def codon_adaptation_index(seq):     
    codon_table = create_codon_table()
    trans_seq = translate(seq, codon_table)
    codon_counts = count_codons(seq)
    CU = [(codon, codon_table[codon], float(codon_counts[codon])) for codon in codon_counts.keys()]
    headers_CU = ['codon','aminoacid','n_codons',]
    CU = DataFrame(CU, columns=headers_CU)
    W = [CU['n_codons'][k]/float(max(CU['n_codons'][CU['aminoacid']==CU['aminoacid'][k]])) for k in range(0,64)]  
    W = np.array(W)
    W[np.isnan(W)]=0
    cai = geomean(W)
    return cai

## relative codon usage bias following Jesse M. Fox and Ivan Erill*
def relative_codon_bias_sequence(seq):
    codon_table = create_codon_table()
    trans_seq = translate(seq, codon_table)
    codon_counts = count_codons(seq)
    CU = [(codon, codon_table[codon], float(codon_counts[codon])) for codon in codon_counts.keys()]
    headers_CU = ['codon','aminoacid','n_codons',]
    CU = DataFrame(CU, columns=headers_CU)

    # sequences on certain codon positions
    seq1 = seq[0::3] 
    seq2 = seq[1::3]   
    seq3 = seq[2::3]

    # weights in eq 4 of the paper
    Wxyz = [CU['n_codons'][k] / (seq1.count(CU['codon'][k][0]) + seq2.count(CU['codon'][k][1]) + seq3.count(CU['codon'][k][2])) for k in range(0,64)]  
    Wxyz = np.array(Wxyz)
    Wxyz[np.isnan(Wxyz)]=0

    # eq/ (4)
    rscb = geomean(Wxyz)-1
    return rscb 

## Relative synonymous codon usage (RSCU) is defined as the ratio of the observed frequency of codons to the expected 
## frequency given that all the synonymous codons for the same amino acids are used equally.
def rs_codon_usage(seq):     
    codon_table = create_codon_table()
    trans_seq = translate(seq, codon_table)
    codon_counts = count_codons(seq)
    CU = [(codon, codon_table[codon], float(codon_counts[codon])) for codon in codon_counts.keys()]
    headers_CU = ['feature','aminoacid','value',]
    CFP = DataFrame(CU, columns=headers_CU)
    W = [CFP['value'][k]/float(np.mean(CFP['value'][CFP['aminoacid']==CFP['aminoacid'][k]])) for k in range(0,64)]  
    CFP['value'] = W
    CFP = CFP.sort_index(by=['feature'], ascending=[True])
    X_CFP = CFP.ix[:,[0,2]]
    df_CFP = X_CFP.copy()
    df_CFP['feature'] = X_CFP['feature'] + '_rs' 
    return df_CFP

def early_codon_adaptation_index(seq):
     if len(seq)<149:
         early_cai = codon_adaptation_index(seq)
     else:
         early_cai = codon_adaptation_index(seq[0:149]) ## codon usage in the      
     return early_cai

## function to count codon pairs from a sequence
def count_codon_pairs(seq):
  codon_pairs_table = create_codon_pairs_table()
  counts = collections.defaultdict(int)
  for codon in codon_pairs_table.keys(): counts[codon] = 0
  for i in range(0,len(seq)-len(seq)%6,6): codon = seq[i:i+6] ; counts[codon] += 1
  return(counts)

## function to count codon triplets from a sequence
def count_codon_triplets(seq):
  codon_triplets_table = create_codon_triplets_table()
  counts = collections.defaultdict(int)
  for codon in codon_triplets_table.keys(): counts[codon] = 0
  for i in range(0,len(seq)-len(seq)%9,9): codon = seq[i:i+9] ; counts[codon] += 1
  return(counts)

## calculates the number of codons of a sequence
def count_number_codons(seq):
   number = len(seq)/3-len(seq)%3
   return number

## counts the number of basis
def count_number_basis(seq):
   number = len(seq)
   return number

## calculates the GC content
def gc_content(seq,loc=0): 
    ## loc represents the location of the bases position
    if loc==1: seq = seq[0::3] 
    if loc==2: seq = seq[1::3]   
    if loc==3: seq = seq[2::3]
    content = (seq.count('g')+seq.count('c'))/float(count_number_basis(seq))
    return(content) 

## calculates the GC in the first 50 codons
def early_gc_content(seq):
     if len(seq)<149:
         early_gc_content = gc_content(seq)
     else:
         early_gc_content = gc_content(seq[0:149]) ## codon usage in the 
     return early_gc_content

## calculates the GC ratio
def gc_ratio(seq): 
   content = (seq.count('g')+seq.count('c'))/float((seq.count('a')+seq.count('t')))
   return(content) 

## calculates the AT content
def at_content(seq): 
   content = (seq.count('a')+seq.count('t'))/float(count_number_basis(seq))*100
   return(content) 

## calculates the A content (for the 5utr)
def a_content(seq): 
   content = (seq.count('a'))/float(count_number_basis(seq))*100
   return(content) 

## calculates the ATG content (for 5utr)
def atg_frequency(seq): 
   content = seq.count('atg')/float(count_number_basis(seq))*100
   return(content) 

## polyadenylation motif AATAAA
def polyadenylation(seq): 
   content = seq.count('aataaa')/float(count_number_basis(seq))*100
   return(content) 

## calculates the AT ratio
def at_ratio(seq): 
   content = (seq.count('a')+seq.count('t'))/float((seq.count('g')+seq.count('c')))
   return(content) 

## codon usage
def codon_usage(seq):
    codon_counts = count_codons(seq)
    codon_table = create_codon_table()
    CU = [(codon, codon_table[codon], float(codon_counts[codon])/sum(codon_counts.values())) for codon in codon_counts.keys()]
    df_CU = DataFrame(CU, columns=['feature','aminoacid','value']).sort_index(by=['feature'], ascending=[True])
    X_CU = df_CU.ix[:,[0,2]]
    return X_CU

## codon pairs usage
def codon_pairs_usage(seq):
    codon_pairs_counts = count_codon_pairs(seq)
    codon_pairs_table = create_codon_pairs_table()
    CU = [(codon, codon_pairs_table[codon], float(codon_pairs_counts[codon])/sum(codon_pairs_counts.values())) for codon in codon_pairs_counts.keys()]
    df_CU = DataFrame(CU, columns=['feature','aminoacid','value']).sort_index(by=['feature'], ascending=[True])
    X_CU = df_CU.ix[:,[0,2]]
    return X_CU

## function to extract the aminoacid frequencies from the translated region of a gene
def aminoacid_frequency(seq):
    codon_counts = count_codons(seq)
    codon_table = create_codon_table()
    CU = [(codon, codon_table[codon], float(codon_counts[codon])/sum(codon_counts.values())) for codon in codon_counts.keys()]
    df_CU = DataFrame(CU, columns=['feature','aminoacid','value']).sort_index(by=['aminoacid'], ascending=[True])  
    df_grouped = df_CU.groupby('aminoacid') 
    df_amin = DataFrame([df_CU['aminoacid'].unique(),df_grouped.value.sum()]).T   
    df_amin.columns = ['feature','value']
    return df_amin
 
## calculates the codon usage in the first 50 codons
def early_codon_usage(seq):
     if len(seq)<149:
         early_CU = codon_usage(seq)
     else:
         early_CU = codon_usage(seq[0:149]) ## codon usage in the 
     
     early_CU['feature'] = early_CU['feature']+'_early50' 
     return early_CU

def TISmotif(seq): 
    counter = 0 
    for y in ['t','c']:
        for m in ['a','c']:
            for r in ['a','g']:
                for v in ['a','c','g']:
                    motif =  r + y + m + r + m + v + 'atggc'
                    counter += seq.count(motif)
    return counter

## free fold energy of the whole sequence
def free_energy(seq):
    return RNA.fold(seq)[1]

## free fold energy of the first 50 bases    
def free_energy_early40(seq):
    return RNA.fold(seq[0:40])[1]
    
def free_energy_late50(seq):
    return RNA.fold(seq[0:40])[1]
    

def local_secondary_structure(seq,length,type='min'):
    def chunkseq(sequence, length): ## splits the original sequence in parts 
        return (sequence[0+i:length+i] for i in range(0, len(sequence), length))
       
    if len(seq)< 2*length:       # no standard deviation otherwise
        return free_energy(seq) 
    
    else:
        max_length = len(seq)/length * length    
        seq_chunks = list(chunkseq(seq[:max_length], length))
        n_chunks   = len(seq_chunks)
        chunks_energy = np.zeros(n_chunks)
    
        for k in range(len(seq_chunks)):
            chunks_energy[k] = free_energy(seq_chunks[k])  ## free energy of the chunks

        if type == 'min':  ## chunk the the minimum free energy
            return chunks_energy.min()                     

        elif type == 'normalized': ## chunk the the minimum free energy normalized by the chunks
            return ((chunks_energy - np.mean(chunks_energy))/np.std(chunks_energy)).min()
            
        else:
            n_replicates = 500 ## permute the sequences in the chunks trying to improve their minimum energy 
            randomized_chunks_energy = np.zeros((n_chunks,n_replicates))
            for i in range(n_chunks):
                for j in range(n_replicates):
                    new_seq_chunk = ''.join(random.sample(seq_chunks[i],length))            
                    randomized_chunks_energy[i,j] = free_energy(new_seq_chunk) 
        
            best_random_chunks_energy = randomized_chunks_energy.min(1)
            return ((chunks_energy - np.mean(best_random_chunks_energy))/np.std(best_random_chunks_energy)).min()

# Terminal oligopyrimidine repeats, Motif tasken from 'The ever-evolving role of mTOR in translation'
def TOP_iniatiation(seq):
  motifs =  ['ctttt',
            'ctctttcc',
            'ctctcc',
            'cctttc',
            'ctctttt',
            'ctttcc',
            'cttctctctc',
            'ctcttttcct',
            'cttttcttt',
            'cttctttctc',
            'ctctttcttcct',
            'ctcttcct',
            'ctctttcc',
            'cttttc',
            'ctttt']
  counter = 0 
  for motif in motifs:
    counter += seq.count(motif)
  return counter

def TOP_elongation(seq):
  motifs =  ['ctttttc',
            'ctttttcctctcttc'
            'ccctttc'
            'ctttcttt'
            'ctcttcc']

  counter = 0 
  for motif in motifs:
    counter += seq.count(motif)
  return counter

def TOP_rnabinding(seq):
  motifs =  ['ctctctttc',
            'ctcctttct',
            'cccttctcccc',
            'ctctctcc',
            'ctttt']
  counter = 0 
  for motif in motifs:
    counter += seq.count(motif)
  return counter

# only those in refererence [131] of the paper
def TOP_ribosomal_proteins(seq):
  motifs = ['ccctttt',
            'ctctttcctt',
            'cctttcc',
            'ctctttt',
            'cttttt',
            'cttccttttc',
            'cctttt',
            'ctctttctc',
            'cctct',
            'cttttcct',
            'cttcctctttttcc',
            'cttcttccttctc',
            'cctttcc',
            'ctctttccct',
            'cctttttc',
            'cttcctttcc',
            'cttcc']
  counter = 0 
  for motif in motifs:
    counter += seq.count(motif)
  return counter

# AREs and variants
def AREs(seq):
  counter = 0 
  at = ['a','t']
  for i in at:
    for j in at:
      for k in at:
        for s in at:
          motif =  k + s + j + 'tatttattt' + i
          counter += seq.count(motif)
  return counter

# Kozak consensus sequence and variants
def kozak_sequence(seq):
  counter = 0
  for R in ['a','g']:
    motif = 'gcc' + R + 'ccatgg'
    counter += seq.count(motif)
  return counter

# Glycosylation consensus sequence and variants
def glycosylation_consensous_sequence(seq):
  code = create_codon_table()
  translated_seq = translate(seq,code)
  yy = ['S','T','C']
  xx = [ 'A', 'C','D','E','F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'Q', 'S', 'V', 'W', 'Y'] 
  counter = 0
  for x in xx:
    for y in yy:  
      motif = 'N' + x + y  
      counter += translated_seq.count(motif)
  return counter

## -----------------
## FUNCTIONS TO WRAP THE FEATURES UP
## -----------------

## Features for 3'utr
def extract_features_3utr(seq):
    names_other = ['number_codons_3utr',
                   'polyadenylation_3utr',
                   'free_energy_3utr',
                   'best_local_sec_structure40_3utr',
                   'best_normalized_sec_structure40_3utr',
                   'best_randomized_sec_structure40_3utr',
                   'best_local_sec_structure60_3utr',
                   'best_normalized_sec_structure60_3utr',
                   'best_randomized_sec_structure60_3utr',
                   'AREs_3utr']
 
    other    =  [count_number_codons(seq),
                 polyadenylation(seq),
                 free_energy(seq),
                 local_secondary_structure(seq,40,type='min'),
                 local_secondary_structure(seq,40,type='normalized'),
                 local_secondary_structure(seq,40,type='randomized'),
                 local_secondary_structure(seq,60,type='min'),
                 local_secondary_structure(seq,60,type='normalized'),
                 local_secondary_structure(seq,60,type='randomized'),
                 AREs(seq)
                 ]
   
    df_others = DataFrame([names_other,other]).T
    df_others.columns = ['feature','value']   
    return df_others
    
## Features for 5'utr
def extract_features_5utr(seq):
    names_other = ['number_codons_5utr',
                   'aug_frequency_5utr',
                   'a_content_5utr',    
                   'free_energy_5utr',
                   'free_energy_late50_5utr',
                   'best_local_sec_structure40_5utr',
                   'best_normalized_sec_structure40_5utr',
                   'best_randomized_sec_structure40_5utr',
                   'best_local_sec_structure60_5utr',
                   'best_normalized_sec_structure60_5utr',
                   'best_randomized_sec_structure60_5utr',
                   'TOP_iniatiation_5utr',
                   'TOP_elongation_5utr',
                   'TOP_rnabinding_5utr',
                   'TOP_ribosomal_proteins',
                   'AREs_5utr']
    
    other    =  [count_number_codons(seq),
                 atg_frequency(seq),
                 a_content(seq),
                 free_energy(seq),
                 free_energy_late50(seq),
                 local_secondary_structure(seq,40,type='min'),
                 local_secondary_structure(seq,40,type='normalized'),
                 local_secondary_structure(seq,40,type='randomized'),
                 local_secondary_structure(seq,60,type='min'),
                 local_secondary_structure(seq,60,type='normalized'),
                 local_secondary_structure(seq,60,type='randomized'),
                 TOP_iniatiation(seq),
                 TOP_elongation(seq),
                 TOP_rnabinding(seq),
                 TOP_ribosomal_proteins(seq),
                 AREs(seq)]
   
    df_others = DataFrame([names_other,other]).T
    df_others.columns = ['feature','value']   
    return df_others

## Features for translated region    
def extract_features_tregion(seq):   
    ## codon usage
    df_codon_usage = codon_usage(seq)

    ## early codon usage
    df_early_codon_usage = early_codon_usage(seq)
    
    ## codon pairs usage
    df_codon_pairs_usage = codon_pairs_usage(seq)

    ## codon frequency patterns
    df_rs_codon_usage = rs_codon_usage(seq)

    ## aminoacids frequency
    df_amin_frequency = aminoacid_frequency(seq)
        
    ## other features
    names_other = ['cai', 
                   'cai_early50', 
                   'rcbs',
                   'number_codons',
                   'gc_content',
                   'gc_content_early50',
                   'gc_content_cp1',
                   'gc_content_cp2',
                   'gc_content_cp3',
                   'gc_ratio',
                   'at_content',
                   'at_ratio',
                   'TIS_motif',
                   'free_energy',
                   'free_energy_early40',
                   'kozak_sequence',
                   'glycosylation_consensous_sequence']
    
    other   = [codon_adaptation_index(seq),
                early_codon_adaptation_index(seq),
                relative_codon_bias_sequence(seq),
                count_number_codons(seq),
                gc_content(seq),
                early_gc_content(seq),
                gc_content(seq,1),
                gc_content(seq,2),
                gc_content(seq,3),
                gc_ratio(seq),
                at_content(seq),
                at_ratio(seq),
                TISmotif(seq),
                free_energy(seq),
                free_energy_early40(seq),
                kozak_sequence(seq),
                glycosylation_consensous_sequence(seq)]
    
    
    df_others = DataFrame([names_other,other]).T
    df_others.columns = ['feature','value']
    
    ## join all features in a single data frame
    df_features = concat([df_codon_usage, 
                          df_early_codon_usage,
                          df_codon_pairs_usage,
                          df_rs_codon_usage,
                          df_amin_frequency,
                          df_others
                          ], ignore_index=True)
 
    return df_features

 
## Features from a whole gene sequence 3utr+tregion+4utr
def extract_features_total(seq_tregion,seq_3utr,seq_5utr):
    seq = seq_3utr + seq_tregion + seq_5utr
    
    names_other = ['gc_content_total']    
    other       =  [gc_content(seq)]
    df_others = DataFrame([names_other,other]).T
    df_others.columns = ['feature','value']
    return df_others
 
## Process the whole gene 
def process_gene(seq_tregion,seq_3utr,seq_5utr):
    features_tregion = extract_features_tregion(seq_tregion)
    features_3utr    = extract_features_3utr(seq_3utr)
    features_5utr    = extract_features_5utr(seq_5utr)
    features_total   = extract_features_total(seq_tregion,seq_3utr,seq_5utr)
    
    ## merge features in a vector
    df_features = concat([features_tregion, 
                          features_3utr,
                          features_5utr,
                          features_total], ignore_index=True)
 
    return df_features











