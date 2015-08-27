import pandas as pd

# load the mapping
df_mapping  = pd.read_csv('mapping_table.csv', sep=',').dropna()
df_mapping['cds'] = 'cds' + df_mapping['cds'].astype(str)
df_mapping['rna'] = 'rna' + df_mapping['rna'].astype(str)

file_read = open('cds_original.fasta','r')
file_write = open('cds.fasta','w')
line = file_read.readline()

while line!='':
    if list(line)[0]=='>':
        print line
        line = line.replace("\n", "").replace(">", "")
        df_map = df_mapping[df_mapping['cds']== line].reset_index()
        if not df_map.empty:
            rna_code = df_map['rna'][0]
            file_write.write('>'+rna_code+'\n')
            line = file_read.readline() 
        else:
            while True:
                line = file_read.readline()
                if list(line)[0]=='>':
                    break
    else:
        file_write.write(line) 
        line = file_read.readline() 
    
file_read.close()
file_write.close()
