## define the table of codons/aminoacids
def create_codon_table():
   bases = ['t', 'c', 'a', 'g']
   codons = [a+b+c for a in bases for b in bases for c in bases]
   amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
   codon_table = dict(zip(codons, amino_acids))
   return(codon_table)

def create_aminoacids_table():
  codon_dict = create_codon_table()
  amin_dict = {}
  for k, v in codon_dict.iteritems():
    amin_dict.setdefault(v, []).append(k)
  return(amin_dict)

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
