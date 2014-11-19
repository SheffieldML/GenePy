import collections
from pandas import DataFrame
import numpy as np
from pandas import concat
from scipy.special import erfc


## FUNCTIONS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

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
def gc_content(seq): 
   content = (seq.count('g')+seq.count('c'))/float(count_number_basis(seq))*100
   return(content) 

## calculates the GC ratio
def gc_ratio(seq): 
   content = (seq.count('g')+seq.count('c'))/float((seq.count('a')+seq.count('t')))
   return(content) 

## calculates the AT content
def at_content(seq): 
   content = (seq.count('a')+seq.count('t'))/float(count_number_basis(seq))*100
   return(content) 

## calculates the AT ratio
def at_ratio(seq): 
   content = (seq.count('a')+seq.count('t'))/float((seq.count('g')+seq.count('c')))
   return(content) 

## function to extract features from a gene sequence (more features can be easily added ) 
def extract_features(seq,gene_id):
 codon_counts = count_codons(seq)
 codon_table = create_codon_table()
 CU = [(codon, codon_table[codon], float(codon_counts[codon])/sum(codon_counts.values())) for codon in codon_counts.keys()]
 headers_CU = ['feature','aminoacid',gene_id,]
 df_CU = DataFrame(CU, columns=headers_CU)
 df_CU = df_CU.sort_index(by=['feature'], ascending=[True]) ## codons are sorted in aphabetical order
 X_CU = df_CU.ix[:,[0,2]]
 names_other = ['cdi','number_codons','gc_content','gc_ratio','at_content','at_ratio'] ## add more features here
 other = [codon_adaptation_index(seq),count_number_codons(seq),gc_content(seq),gc_ratio(seq),at_content(seq),at_ratio(seq)]
 X_others = DataFrame([names_other,other]).T
 X_others.columns = ['feature',gene_id]
 df_features = concat([X_CU, X_others], ignore_index=True)
 return(df_features)











