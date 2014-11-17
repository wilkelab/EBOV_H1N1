# SJS. Prints out the id sets for nuc or aa seqs. See commented out lines and swap for nuc or amino acid.
# So, we'll end up, for each gene, with a list of ids which have the same sequence in th 2014 outbreak.


from Bio import AlignIO
from Bio import Seq
for gene in ['np', 'l', 'gp', 'vp24', 'vp30', 'vp35', 'vp40']:
    print "\n\n========== " + gene + " =========="
    aln = AlignIO.read('EBOV_all_outbreaks/' + gene + '_nuc.aln', 'fasta')
    seqs = []
    for rec in aln:
        if "2014" in str(rec.id):
            seqs.append(str(rec.seq))
            #seqs.append(str(rec.seq.translate()))
    seqs = list(set(seqs))
    final = {}
    for rec in aln:
        if "2014" in str(rec.id):
            index = seqs.index( str(rec.seq) )
            #index = seqs.index( str(rec.seq.translate()) )
            if index not in final:
                final[index] = [str(rec.id)]
            else:
                final[index].append(str(rec.id))
    for entry in final:
        print "\n".join( final[entry] )
        print

        
                    
            
