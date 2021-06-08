from Bio import motifs
from Bio.Seq import Seq
import os


map_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\mapping_data\uniprot_genename.txt'
fasta_dir = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\tmp_data'
taxa_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\mapping_data\taxid_name_assembly.tsv'

# motif_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\ts_motifs\hsa-miR-197-3p.xml'
motif_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\ts_motifs\hsa-miR-769-5p.xml'
target_prot = 'MAPK1'

ts_seq = 'TCCCCACTGCTCTCTGTGGTGAG'

# read motif
with open(motif_path, 'r') as mof_fh:
    mof = motifs.read(mof_fh, 'MEME')
print(mof.pssm.consensus)

genename_2_uniprot = {}
uniprot_2_genename = {}
# read protein mapping
with open(map_path, 'r') as fh:
    for line in fh:
        uni, name = line.strip().split()
        genename_2_uniprot[name] = uni
        uniprot_2_genename[uni] = name

# read taxon mapping
tax_dict = {}
with open(taxa_path, 'r') as fh:
    for line in fh:
        taxid, name, assembly = line.strip().split('\t')
        tax_dict[taxid] = name

# read multifasta
uni_target = genename_2_uniprot[target_prot]
mf_file = f'{fasta_dir}{os.sep}{uni_target}.fa'
print(mf_file)
with open(mf_file, 'r') as fh:
    for line in fh:
        if line.startswith('>'):
            taxid = line.strip().replace('>ncbi', '')
            taxon = tax_dict[taxid]
        else:
            seq = line.strip()
            if ts_seq in seq:
                print('')
                print(taxon)
                print('')
            count = 0
            for pos, score in mof.pssm.search(seq):
                # exlude reverse matches
                if pos < 0:
                    continue
                # print(seq[pos:pos+8], score)
                count += 1
            print(f'{taxon}, count: {count}')

