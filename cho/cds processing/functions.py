
## Reads a fasta file and puts each sequence in a row
def process_utrs(input_file,output_file):
    file_read = open(input_file,'r')
    file_write = open(output_file,'w')
    file_write.write('Gene_ID gene_sequence')
    line_r = file_read.readline()
    line = line_r.replace("\n", "")

    while line!='':
        if list(line)[0]=='>':
            file_write.write('\n')
            file_write.write(line.partition(' | ')[0]+ ' ')
   
        if list(line)[0]!='>':	
            line = line.replace("Sequence unavailable","NA")
            file_write.write(line)
  
        line_r = file_read.readline()
        line = line_r.replace("\n", "")
    file_read.close()
    file_write.close()


def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]


def df2fasta(df,output_file):
    file_write = open(output_file,'w')

    for k in range(df.shape[0]):
        seq_chunks = split_len(df['gene_sequence'][k],60)  #60 is for the standard fasta format

        # write gene id
        file_write.write(df['Gene_ID'][k]+'\n')

        # write sequence
        for j in range(len(seq_chunks)):
            file_write.write(seq_chunks[j]+'\n')

    file_write.close()


