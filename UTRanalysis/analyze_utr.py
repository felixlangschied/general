import glob
import json
import os

# human_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\GCF_000001405.39_GRCh38.p13_UTRs.json'
# mouse_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\GCF_000001635.27_GRCm39_UTRs.json'
utr_data_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data'
outdir = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\prot_align'
target_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'
tax_list = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\taxid_name_assembly.tsv'
# utr_data_path = '/share/project/felixl/ncOrtho/data/UTR/all_iso_all_genes'
# target_path = '/home/felixl/project/general/UTRanalysis/list_of_connected_targetprots.txt'
# tax_list = '/share/project/felixl/ncOrtho/data/UTR/taxid_name_assembly.tsv'
# outdir = '/home/felixl/project/general/UTRanalysis/utr_alignments/longest_utr'


longest_only = True

if not os.path.isdir(outdir):
    os.mkdir(outdir)

protein_list = []
with open(target_path, 'r') as fh:
    protnames = fh.read().split('\n')
    empty = protnames.pop()
    if empty:
      print('# Removed {} from protein list. Please double check'.format(empty))
for tname in protnames:
    mousename = tname.capitalize()
    protein_list.append((tname, mousename))

assembly_2_name = {}
name_2_assembly = {}
with open(tax_list, 'r') as fh:
    for line in fh:
        taxid, name, assembly = line.strip().split('\t')
        if assembly == 'None':
            continue
        assembly_2_name[assembly] = name
        name_2_assembly[name] = assembly

# with open(human_path, 'r') as fh:
#     human_data = json.load(fh)
# with open(mouse_path, 'r') as fh:
#     mouse_data = json.load(fh)

for protein in protein_list:
    print('# Starting UTR extraction for {}'.format(protein[0]))
    with open(f'{outdir}{os.sep}{protein[0]}_utr.fa', 'w') as of:
        for assem in assembly_2_name:
        # for assem in list(assembly_2_name)[0:1]:
            # load utr data
            ass_path = f'{utr_data_path}/{assem}_UTRs.json'
            with open(ass_path, 'r') as fh:
                utr_dict = json.load(fh)
            # try to extract UTR sequence with each type of identifier (e.g MAPK1 or Mapk1)
            # print(utr_dict.keys())
            longest_length = 0
            try:
                for c, prot_id in enumerate(utr_dict[protein[0]]):
                    # header = '>{}_{}_{}\n'.format(assembly_2_name[assem], protein[0], c)
                    seq = utr_dict[protein[0]][prot_id]
                    length = len(seq)
                    if longest_only and length > longest_length:
                        longest_length = length
                        out_seq = seq
                    else:
                        header = '>{}_{}_{}\n'.format(assembly_2_name[assem], prot_id, protein[0])
                        of.write(header)
                        of.write(seq)
                        of.write('\n')

            except KeyError:
                try:
                    for c, prot_id in enumerate(utr_dict[protein[1]]):
                        # header = '>{}_{}_{}\n'.format(assembly_2_name[assem], protein[0], c)
                        seq = utr_dict[protein[1]][prot_id]
                        length = len(seq)
                        if longest_only and length > longest_length:
                            longest_length = length
                            out_seq = seq
                        else:
                            header = '>{}_{}_{}\n'.format(assembly_2_name[assem], prot_id, protein[0])
                            of.write(header)
                            of.write(seq)
                            of.write('\n')
                except KeyError:
                    print('# No UTR found for {} in {}'.format(protein[0], assembly_2_name[assem]))
                    continue
            if longest_only:
                header = '>{}_{}\n'.format(assembly_2_name[assem], protein[0])
                of.write(header)
                of.write(out_seq)
                of.write('\n')
