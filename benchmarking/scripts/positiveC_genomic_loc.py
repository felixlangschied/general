import json

all_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/ALL.gff'
overview_path = '/home/felixl/PycharmProjects/general/benchmarking/data/05overview.txt'
spec_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/species_list.txt'
outpath = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb_locations.json'
res_out = '/home/felixl/PycharmProjects/general/benchmarking/data/result_locations.json'
pres_path = '/home/felixl/PycharmProjects/general/benchmarking/data/present.txt'

with open(pres_path, 'r') as fh:
    present = fh.read().split('\n')
    present.pop()


# read mirgenedb entries
spec_map = {}
with open(spec_path, 'r') as fh:
    for line in fh:
        abr, name, taxid = line.strip().split('|')
        spec_map[abr.capitalize()] = name
df_dict = {'species': [], 'mirid': [], 'mirfam': [], 'chrom': [], 'start': [], 'end': [], 'strand': []}
with open(all_path, 'r') as fh:
    for c, line in enumerate(fh):
        if c == 0 or c % 1000 == 0:
            print(f'progress: {c}')
        if not line.startswith('#'):
            data = line.strip().split('\t')
            if data[2] != 'pre_miRNA':
                continue
            mirid = data[-1].replace('ID=', '').split(';')[0]
            mirfam = '-'.join(mirid.split('-')[1:3]).replace('_pre', '')
            species = spec_map[mirid.split('-')[0]]
            if species not in present:
                continue
            df_dict['species'].append(species)
            df_dict['chrom'].append(data[0])
            df_dict['start'].append(data[3])
            df_dict['end'].append(data[4])
            df_dict['strand'].append(data[6])
            df_dict['mirid'].append(mirid)
            df_dict['mirfam'].append(mirfam)
with open(outpath, 'w') as of:
    json.dump(df_dict, of)

# read results
path_list = []
with open(overview_path, 'r') as fh:
    for line in fh:
        path_list.append(line.split('\t')[0])
res_dict = {'species': [], 'mirid': [], 'mirfam': [], 'chrom': [], 'start': [], 'end': [], 'strand': []}
for c, file in enumerate(path_list):
    if c == 0 or c % 50 == 0:
        print(f'Progress: {c}')
    species = file.split('/')[7]
    if species not in present:
        continue
    file = file.replace('mirGeneDB', 'mirgenedb')
    mirna = file.split('/')[-1].split('_')[0]
    mirfam = '-'.join(mirna.split('-')[1:3])
    with open(file, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                data = line.strip().split('|')
                res_dict['species'].append(species)
                res_dict['mirid'].append(mirna)
                res_dict['mirfam'].append(mirfam)
                res_dict['chrom'].append(data[1])
                res_dict['start'].append(data[2])
                res_dict['end'].append(data[3])
                res_dict['strand'].append(data[4])
with open(res_out, 'w') as of:
    json.dump(res_dict, of)

