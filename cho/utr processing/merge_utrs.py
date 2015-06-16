
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
            file_write.write(line[20:38]+ ' ')
   
        if list(line)[0]!='>':	
            line = line.replace("Sequence unavailable","NA")
            file_write.write(line)
  
        line_r = file_read.readline()
        line = line_r.replace("\n", "")
    file_read.close()
    file_write.close()

## From a processed fasta, with sequences in a row, merges sequences with the same id
#def merge_utrs(input_file,output_file):
#    file_read = open(input_file,'r')
#    file_write = open(output_file,'w')
#    file_write.write('ID gene_sequence')
#    line_r = file_read.readline()

#    while line!='':
#        if list(line)[0]=='>':
#            id = line[1:]

## From a file with sequences in a row, creates a compatible fasta file
#def create_fasta(input_file,output_file)
#    file_read = open(input_file,'r')
#    file_write = open(output_file,'w')





