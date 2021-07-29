import json
import pandas as pd
import glob

loc_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb_locations.json'
res_path = '/home/felixl/PycharmProjects/general/benchmarking/data/result_locations.json'
map_dir = '/home/felixl/PycharmProjects/general/benchmarking/data/chrom_map'
avail_path = '/home/felixl/PycharmProjects/general/benchmarking/data/available_specs.txt'




with open(loc_path, 'r') as fh:
    loc_dict = json.load(fh)
loc_df = pd.DataFrame.from_dict(loc_dict)
loc_df['start'] = pd.to_numeric(loc_df['start'])
loc_df['end'] = pd.to_numeric(loc_df['end'])


with open(res_path, 'r') as fh:
    res_dict = json.load(fh)
res_df = pd.DataFrame.from_dict(res_dict)
res_df['start'] = pd.to_numeric(res_df['start'])
res_df['end'] = pd.to_numeric(res_df['end'])
specs = set(res_df['species'])
print(res_df.columns)

map_dict = {}
with open(avail_path, 'r') as fh:
    for line in fh:
        taxid, name, gcf, ass, state = line.strip().split('|')
        name = name.replace(' ', '_')
        if name not in specs:
            continue
        spec_dict = {'gcf': gcf, 'assembly': ass}
        map_dict[name] = spec_dict

for species in map_dict:
    gcf = map_dict[species]['gcf']
    file = glob.glob(f'{map_dir}/{gcf}*')[0]
    # chrom_dict = {}
    with open(file, 'r') as fh:
        for line in fh:
            if not line.startswith('#'):
                data = line.strip().split()
                if data[3] == 'Chromosome':
                    nc = data[6]
                    # print(nc)
                    chr = data[2]
                    res_df.loc[nc, 'chrom'] = chr

# print(res_df.head())
# print(loc_df.head())
mirfams = sorted(list(set(loc_df['mirfam'])))
# print(mirfams)

for spec in specs:
    for fam in mirfams:
        spec_res = res_df[res_df['species'] == spec]
        spec_res = spec_res[spec_res['mirfam'] == fam].sort_values(by='start')
        if not spec_res.empty:
            print(spec_res)
            loc_res = loc_df[loc_df['species'] == spec]
            loc_res = loc_res[loc_res['mirfam'] == fam].sort_values(by='start')
            print(loc_res)
            # check if hits overlap
            for row in spec_res.iterrows():
                print(row)
            # if (
            #     loc_res['strand'].values[0] == spec_res['strand'].values[0]
            # ):
            #     print('yep')


            exit()



