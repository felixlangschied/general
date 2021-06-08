import glob
import json
import os

human_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\GCF_000001405.39_GRCh38.p13_UTRs.json'
mouse_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data\GCF_000001635.27_GRCm39_UTRs.json'
outdir = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\utr_data'
target_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'

protein_list = []
with open(target_path, 'r') as fh:
    protnames = fh.read().split('\n')
    empty = protnames.pop()
    if empty:
      print('# Removed {} from protein list. Please double check'.format(empty))
for tname in protnames:
    mousename = tname.capitalize()
    protein_list.append((tname, mousename))

# # always human first
# protein_list = [
#     ('MAPK1', 'Mapk1'),
#     ('CDK2', 'Cdk2'),
#     ('BMPR1A', 'Bmpr1a'),
#     ('IGFBP5', 'Igfbp5'),
#     ('RAN', 'Ran'),
#
# ]

with open(human_path, 'r') as fh:
    human_data = json.load(fh)
with open(mouse_path, 'r') as fh:
    mouse_data = json.load(fh)

for protein in protein_list:
    with open(f'{outdir}{os.sep}{protein[0]}_utr.fa', 'w') as of:
        for prot_id in human_data[protein[0]]:
            seq = human_data[protein[0]][prot_id]
            print(len(seq))
            header = f'>human_{prot_id}\n'
            of.write(header)
            of.write(seq)
            of.write('\n')
        print('')
        for prot_id in mouse_data[protein[1]]:
            seq = mouse_data[protein[1]][prot_id]
            print(len(seq))
            header = f'>mouse_{prot_id}\n'
            of.write(header)
            of.write(seq)
            of.write('\n')