from Bio.Seq import Seq
from Bio import AlignIO

matmir = 'hsa-miR-197-3p'
align_path = r'C:\Users\felix\project\UTRanalysis\target_analysis\{}.aln'.format(matmir)
mat_path = r'C:\Users\felix\project\UTRanalysis\target_analysis\all_mats.fa'

mat_dict = {}
# parse mature sequences
with open(mat_path, 'r') as fh:
    for line in fh:
        if line.startswith('>'):
            mat_id = line.strip().split()[0].replace('>', '')
        else:
            seq = line.strip()
            mat_dict[mat_id] = Seq(seq)
print(mat_dict)

mat_seq = mat_dict[matmir]
rc_mat = mat_seq.reverse_complement()
align = AlignIO.read(align_path, "fasta")
align.add_sequence(f'rc_{matmir}', str(rc_mat).replace('U', 'T'))
print(align)



