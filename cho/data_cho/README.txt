# --- August 2015

This folder contains the main datafiles from the 'CHO recombiant project'. The notebook 'process_data.ipnb', also is this folder, contains the code to merge this files into a unique
database useful for later modelling purposes.


# --- Raw datafiles description

- mapping_table.csv: contains the different ids that are used across the files to indetify a gene/protein.

- p2_h2l.csv: contains the protein degradation rates (kdp) and mRNA half-life (hl) calculated as schwanhausser (and using the robust method). cds ids are used.

- 2pep_h2l: checking file. The ids of this file are the proteins identified with more that 2 peptides.

- Galaxy116: transciptomics data. The id is reported as tracking_id and the mRNA abundancy as FPKM

- 5utr.fasta, cds.fasta and 3utr.fasta contain the gene sequences in CHO. 

- The venn diagram in venn_datasets is computed using sequences un which the utrs and cdss are available.

- genes in the file genes_with_full_data_id contains the ids of the genes for which data are available (in all the datasets)





===================== Notes from jopseph

2pep_h2l

This file is the output from Persues And represnets the proteins identified in our total serach. Each is grouped into a protien family with one or more cds's being called the majority protein id.
The presence of the cds in this file indicates that the file has 2 or more unique peptides identified



p2_h2l

This csv is my rates output comparing the different ways of processing computing the rates. The 2nd column has the cds # . If there is a line with only NA in all the rates it is not worthwhile 
adding as nothing can be computed



Galaxy116-[Cufflinks_on_data_38,_data_10,_and_data_114__transcript_expression] (3)

This is the RNASEQ output the first column contains the tracking id an RNA number or a gene number. The abundance is reported as FPKM. If this is 0 then it should not be inclded as 
we can't calculated anything from it.


mapping table
