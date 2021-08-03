import json
import pandas as pd
import glob

# loc_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb_locations.json'
# res_path = '/home/felixl/PycharmProjects/general/benchmarking/data/result_locations.json'
# map_dir = '/home/felixl/PycharmProjects/general/benchmarking/data/chrom_map'
# avail_path = '/home/felixl/PycharmProjects/general/benchmarking/data/available_specs.txt'

loc_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb_locations.json'
res_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\result_locations.json'
map_dir = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\chrom_map'
avail_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\available_specs.txt'
map_file = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mapped_result_locations.csv'
tp_out = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\tp_results.csv'


with open(loc_path, 'r') as fh:
    loc_dict = json.load(fh)
loc_df = pd.DataFrame.from_dict(loc_dict)
loc_df['start'] = pd.to_numeric(loc_df['start'], downcast='integer')
loc_df['end'] = pd.to_numeric(loc_df['end'], downcast='integer')
specs = set(loc_df['species'])


with open(res_path, 'r') as fh:
    res_dict = json.load(fh)
res_df = pd.DataFrame.from_dict(res_dict)
res_df['start'] = pd.to_numeric(res_df['start'], downcast='integer')
res_df['end'] = pd.to_numeric(res_df['end'], downcast='integer' )
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


# for species in map_dict:
#     print(species)
#     # platypus is given with refseq ids in mirGeneDB
#     if species in ['Ornithorhynchus_anatinus']:
#         continue
#     gcf = map_dict[species]['gcf']
#     file = glob.glob(f'{map_dir}/{gcf}*')[0]
#     # chrom_dict = {}
#     with open(file, 'r') as fh:
#         for line in fh:
#             if not line.startswith('#'):
#                 data = line.strip().split()
#                 # if data[3] == 'Chromosome':
#                 nc = data[6]
#                 # print(nc)
#                 chr = data[0]
#                 if nc in res_df['chrom'].values:
#                     res_df['chrom'] = res_df['chrom'].replace(nc, chr)
# print(res_df)
# res_df.to_csv(map_file)
#
res_df = pd.read_csv(map_file)
specs = set(res_df['species'])


# print(res_df.head())
# print(loc_df.head())
mirfams = sorted(list(set(loc_df['mirfam'])))
# print(mirfams)


### looks for each result in database
# for spec in specs:
#     print(spec)
#     # if spec != 'Macaca_mulatta':
#     #     continue
#     for fam in mirfams:
#         spec_res = res_df.copy()
#         spec_res = spec_res[spec_res['species'] == spec]
#         spec_res = spec_res[spec_res['mirfam'] == fam].sort_values(by='start')
#         results = spec_res.filter(items=['chrom', 'start', 'end', 'strand'])
#         results['chrom'] = results['chrom'].str.replace('chr', '')
#         results['chrom'] = results['chrom'].str.replace('Chr', '')
#         if not spec_res.empty:
#             loc_res = loc_df[loc_df['species'] == spec]
#             loc_res = loc_res[loc_res['mirfam'] == fam].sort_values(by='start')
#             db_loc = loc_res.filter(items=['chrom', 'start', 'end', 'strand'])
#             db_loc['chrom'] = db_loc['chrom'].str.replace('chr', '')
#             db_loc['chrom'] = db_loc['chrom'].str.replace('Chr', '')
#             # print(results)
#             # print(db_loc)
#             # print(db_loc['chrom'])
#
#             for rownum in range(results['start'].size):
#                 index = results.iloc[rownum].name
#                 # print(index)
#                 chrom, start, end, strand = results.iloc[rownum, :]
#                 # check if there are entries of this miRNA on the same chromosome in the database
#                 if chrom in db_loc['chrom'].values:
#                     db_chrom = db_loc[db_loc['chrom'] == chrom]
#                 else:
#                     # print('chromosomes do not match')
#                     res_df.loc[index, 'True_positive'] = 'False_chromo'
#                     continue
#                 # allow for fuzzy hits
#                 startrange = range(start - 3, start + 4)
#                 endrange = range(end - 3, end + 4)
#                 # check for overlap
#                 st_check = [start for start in db_chrom['start'].values if start in startrange]
#                 end_check = [end for end in db_chrom['end'].values if end in endrange]
#                 if st_check and end_check:
#                     # check strand
#                     hitted = db_chrom[(db_chrom['start'] == st_check[-1]) & (db_chrom['end'] == end_check[-1])]
#                     try:
#                         if hitted.strand.values[0] == strand:
#                             # print('ACCEPTED')
#                             res_df.loc[index, 'True_positive'] = 'True'
#                         else:
#                             res_df.loc[index, 'True_positive'] = 'False_strand'
#                     except IndexError:
#                         print(fam)
#                         print(st_check, end_check)
#                 else:
#                     res_df.loc[index, 'True_positive'] = 'False_location'
#
# # print(res_df['True_positive'])
# res_df.to_csv(tp_out)

### looks for each database entry in results
for spec in specs:
    print(spec)
    # if spec != 'Macaca_mulatta':
    #     continue
    for fam in mirfams:
        #load results
        spec_res = res_df.copy()
        spec_res = spec_res[spec_res['species'] == spec]
        spec_res = spec_res[spec_res['mirfam'] == fam].sort_values(by='start')
        results = spec_res.filter(items=['chrom', 'start', 'end', 'strand'])
        results['chrom'] = results['chrom'].str.replace('chr', '')
        results['chrom'] = results['chrom'].str.replace('Chr', '')

        # load database entries
        loc_res = loc_df.copy()
        loc_res = loc_res[loc_res['species'] == spec]
        loc_res = loc_res[loc_res['mirfam'] == fam].sort_values(by='start')
        db_loc = loc_res.filter(items=['chrom', 'start', 'end', 'strand'])
        db_loc['chrom'] = db_loc['chrom'].str.replace('chr', '')
        db_loc['chrom'] = db_loc['chrom'].str.replace('Chr', '')

        for rownum in range(db_loc['start'].size):
            index = db_loc.iloc[rownum].name
            # print(index)
            chrom, start, end, strand = db_loc.iloc[rownum, :]
            # check if there are entries of this miRNA on the same chromosome in the database
            if chrom in results['chrom'].values:
                res_chrom = results[results['chrom'] == chrom]
            else:
                # print('chromosomes do not match')
                loc_df.loc[index, 'True_positive'] = 'False_chromo'
                continue
            # allow for fuzzy hits
            startrange = range(start - 3, start + 4)
            endrange = range(end - 3, end + 4)
            # check for overlap
            st_check = [start for start in res_chrom['start'].values if start in startrange]
            end_check = [end for end in res_chrom['end'].values if end in endrange]
            if st_check and end_check:
                # check strand
                hitted = res_chrom[(res_chrom['start'] == st_check[-1]) & (res_chrom['end'] == end_check[-1])]
                try:
                    if hitted.strand.values[0] == strand:
                        # print('ACCEPTED')
                        loc_df.loc[index, 'True_positive'] = 'True'
                    else:
                        loc_df.loc[index, 'True_positive'] = 'False_strand'
                except IndexError:
                    print(fam)
                    print(st_check, end_check)
            else:
                loc_df.loc[index, 'True_positive'] = 'False_location'

print(loc_df)
loc_df.to_csv(tp_out)







