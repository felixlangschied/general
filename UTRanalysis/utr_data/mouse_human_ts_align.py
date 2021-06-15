import json
import os

utr_data_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data'
outdir = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\prot_align\all_iso_no_X'
# target_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'
tax_list = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\taxid_name_assembly.tsv'
interaction_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\interactions_per_candidate_miRNA.csv'
mirtar_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\mirtarbase_targetsites.txt'

miRNA = 'hsa-miR-197-3p'
# miRNA = 'hsa-miR-769-5p'
print(miRNA)

longest_only = True

if not os.path.isdir(outdir):
    os.mkdir(outdir)

protein_list = []
with open(interaction_path, 'r') as fh:
    # skip header
    name_index = next(fh).strip().split(',').index('"shared name"')
    for line in fh:
        data = line.strip().replace('"', '').split(',')
        # print(data[name_index].split(' (interacts with) '))
        mirna, human_prot = data[name_index].split(' (interacts with) ')

        if mirna == miRNA:
            mouse_prot = human_prot.capitalize()
            protein_list.append((human_prot, mouse_prot))
print(protein_list)

protein_list.append(('MTHFD1L', 'Mthfd1l'))
protein_list.append(('CTNND1', 'Ctnnd1'))
protein_list.append(('CYLD', 'Cyld'))
protein_list.append(('IGFBP5', 'Igfbp5'))
# protein_list = [('MTHFD1L', 'Mthfd1l')]

assembly_2_name = {}
name_2_assembly = {}
with open(tax_list, 'r') as fh:
    for line in fh:
        taxid, name, assembly = line.strip().split('\t')
        if assembly == 'None':
            continue
        # only use mouse and human data
        if name not in ['Homo_sapiens', 'Mus_musculus']:
            continue
        assembly_2_name[assembly] = name
        name_2_assembly[name] = assembly

target_dict = {}
count = 1
with open(mirtar_path, 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')
        if data[1] == miRNA:
            seq = data[6]
            prot = data[3]
            if not prot in target_dict:
                target_dict[prot] = []
                header = f'>{miRNA}_ts1'
            else:
                count = int(target_dict[prot][-2][-1]) + 1
                header = f'>{miRNA}_ts{count}'
            target_dict[prot].append(header)
            target_dict[prot].append(seq)

print(target_dict)

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
                    if prot_id.startswith('X'):
                        continue
                    # header = '>{}_{}_{}\n'.format(assembly_2_name[assem], protein[0], c)
                    seq = utr_dict[protein[0]][prot_id]
                    length = len(seq)
                    if longest_only:
                        if length > longest_length:
                            longest_length = length
                            out_seq = seq
                            out_prot = prot_id
                    else:
                        header = '>{}_{}_{}\n'.format(assembly_2_name[assem], prot_id, protein[0])
                        of.write(header)
                        of.write(seq)
                        of.write('\n')

            except KeyError:
                try:
                    for c, prot_id in enumerate(utr_dict[protein[1]]):
                        if prot_id.startswith('X'):
                            continue
                        # header = '>{}_{}_{}\n'.format(assembly_2_name[assem], protein[0], c)
                        seq = utr_dict[protein[1]][prot_id]
                        # print(prot_id)
                        length = len(seq)
                        if longest_only:
                            if length > longest_length:
                                longest_length = length
                                out_prot = prot_id
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
                header = '>{}_{}_{}\n'.format(assembly_2_name[assem], out_prot, protein[0])
                of.write(header)
                of.write(out_seq)
                of.write('\n')
                print(header.strip())
                print(out_seq)
        try:
            out_ts = '\n'.join(target_dict[protein[0]])
            print(out_ts)
            of.write(out_ts)
            of.write('\n')
        except KeyError:
            pass

            # try:
            #     print('\n'.join(target_dict[protein[0]]))
            # except KeyError:
            #     pass
            # print(header)
            # print(out_seq)
