import requests, sys
import json
from collections import Counter

# extract go annotation for target proteins

genepath = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'
name_out = 'gene_GOname_counter.json'
go_out = 'gene_GOterm_counter.json'
# searching in mouse or human
# taxid = 10090
taxid = 9606

genenames = []
with open(genepath, 'r') as infile:
    for line in infile:
        genenames.append(line.strip())
print(len(genenames))

# TODO: remove break
# genenames = ['Abca1', 'Atf1']
ensembl_id = []
all_GOs = []
all_names = []
for gene in genenames:
    # print(f'### {gene} ###')
    fst_URL = f'https://www.ebi.ac.uk/QuickGO/services/geneproduct/search?taxonId={taxid}&query={gene}'
    fst_r = requests.get(fst_URL, headers={"Accept": "application/json"})
    hit_dict = json.loads(fst_r.text)['results']
    for hit in hit_dict:
        if hit['symbol'] == gene and hit['databaseSubset'] == 'Swiss-Prot':
            ensembl_id = hit['id']
            # print(ensembl_id)
            break
    # collect GOs
    GOs = []
    go_names = []
    snd_URL = f'https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId={ensembl_id}'
    snd_r = requests.get(snd_URL, headers={"Accept": "application/json"})
    anno_dict = json.loads(snd_r.text)['results']
    for entry in anno_dict:
        GOs.append(entry['goId'])
    all_GOs.extend(GOs)
    # extract names of GO terms
    golist_string = '%2C'.join(GOs).replace(':', '%3A')
    trd_URL = f'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{golist_string}'
    trd_r = requests.get(trd_URL, headers={"Accept": "application/json"})
    name_dict = json.loads(trd_r.text)['results']
    names = [entry['name'] for entry in name_dict]
    all_names.extend(names)
print('# done')

name_C = Counter(all_names)
go_C = Counter(all_GOs)

print(name_C)
with open(name_out, 'w') as of:
    json.dump(name_C, of)