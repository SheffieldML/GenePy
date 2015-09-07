import numpy as np
from pandas import concat
import pandas as pd

## import functions
from functions import *

## pre-process fasta files
print '--------------------------'
print 'Pre-processing fasta files'
print '--------------------------'
process_fasta('./original/3utr-fasta.txt','./original/3utr_processed.txt')
process_fasta('./original/5utr-fasta.txt','./original/5utr_processed.txt')
process_fasta('./original/tregion-fasta.txt','./original/tregion_processed.txt')

## Load the sequences from the genes (removing non available sequences)
seqs_3utr    = pd.read_csv('./original/3utr_processed.txt', sep=' ').dropna().reset_index() 
seqs_5utr    = pd.read_csv('./original/5utr_processed.txt', sep=' ').dropna().reset_index() 
seqs_tregion = pd.read_csv('./original/tregion_processed.txt', sep=' ').dropna().reset_index()  
ids          = pd.read_csv('./original/ids.txt', sep='\t').dropna().reset_index().ID

## genes in the dataser of Schwanhausser 2011 from which the 3utr, 5utr the coding region and the rates are available
common_ids = list(set(seqs_3utr.ID) & set(seqs_5utr.ID) & set(seqs_tregion.ID) & set(ids)) 

print '--------------------------'
print 'Genes processing started'
print '--------------------------'
## We process each gene
for k in range(len(common_ids)):
    print 'Processing gene ... %s' % k
    
    # ID of the gene to process
    gene_id = common_ids[k]

    # Location of this gene in the different data bases
    indx_3utr       = list(seqs_3utr.ID).index(gene_id)
    indx_5utr       = list(seqs_5utr.ID).index(gene_id)
    indx_tregion    = list(seqs_tregion.ID).index(gene_id)

    # Extract the 3utr, 5utr and the translated region sequences of the gene
    sequence_3utr = seqs_3utr.ix[indx_3utr,'gene_sequence'].lower() 
    sequence_5utr = seqs_5utr.ix[indx_5utr,'gene_sequence'].lower()
    sequence_tregion = seqs_tregion.ix[indx_tregion,'gene_sequence'].lower() 

    # Process the gene
    features_gene = process_gene(sequence_tregion, sequence_3utr, sequence_5utr)
    features_gene.columns=['gene_ID',gene_id] 
    
    # Data frame for the results of this iteration
    if k == 0:
        ##--- features
        df_features = features_gene

        ##---- sequences
        df_sequences = DataFrame([sequence_3utr,sequence_tregion,sequence_5utr], ['3UTR','Translated_region','5UTR'])
        df_sequences.columns = [gene_id]
        
    else:
        ##--- features
        df_features = df_features.join(features_gene.ix[:,gene_id]) 
    
        ##----sequences
        df_sequence_k = DataFrame([sequence_3utr,sequence_tregion,sequence_5utr], ['3UTR','Translated_region','5UTR'])
        df_sequence_k.columns = [gene_id]
        df_sequences = df_sequences.join(df_sequence_k)
    
## save data frames
df_features.T.to_csv('./processed/features_genes0.csv',sep=' ',header=True)
df_sequences.T.to_csv('./processed/sequences_genes0.csv',sep=' ',header=True) 

print '--------------------------'
print 'Genes processing finished'
print '--------------------------'
    
    
    
    
    
    
