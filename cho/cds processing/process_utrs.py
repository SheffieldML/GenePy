
from functions import *
import pandas as pd

input_file  = 'CDS.fasta'
aux_file    = 'pre_processed.fasta'
output_file = 'merged.fasta'


## We process the utrs and we load it in a pandas df
process_utrs(input_file,aux_file)
df_sequences = pd.read_csv('pre_processed.fasta', sep=' ')
df_sequences_merged = df_sequences.drop_duplicates('Gene_ID').reset_index(drop=True)

## compute the unique utrs
unique_utrs = df_sequences['Gene_ID'].unique()

for utr in unique_utrs:
	print utr
	seq_utr = ''.join(list(df_sequences['gene_sequence'][df_sequences['Gene_ID']==utr]))
	idx = df_sequences_merged['Gene_ID']==utr
	df_sequences_merged['gene_sequence'][idx] = seq_utr

# save the sequences in fasta format
df2fasta(df_sequences_merged,output_file)

