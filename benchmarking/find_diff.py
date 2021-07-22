import pandas as pd
from collections import Counter

five_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\05overview.txt'
six_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\06overview.txt'

def read_file(path):
    out_dict = {'org': [], 'mirna': [], 'family': [], 'coorth': []}
    with open(path, 'r') as fh:
        for line in fh:
            path, coorths = line.strip().split('\t')
            out_dict['coorth'].append(coorths)
            data = path.split('/')
            out_dict['org'].append(data[7])
            mirna = data[-1].replace('_pre_orthologs.fa', '')
            out_dict['mirna'].append(mirna)
            family = '-'.join(mirna.split('-')[1:3])
            out_dict['family'].append(family)
        return out_dict, pd.DataFrame.from_dict(out_dict)


five_dict, five_df = read_file(five_path)
six_dict, six_df = read_file(six_path)

five_count = Counter(five_dict['family'])
six_count = Counter(six_dict['family'])
five_df = pd.DataFrame.from_dict(five_count, orient='index', columns=['five_count'])
six_df = pd.DataFrame.from_dict(six_count, orient='index', columns=['six_count'])

com = five_df.join(six_df)
com = com.assign(diff=com['five_count'] - com['six_count'])
com = com.sort_values(by='diff', ascending=False)
print(com.head())

top_30_diffs = list(com.head(30).index)
with open('data/top_30_diff.txt', 'w') as of:
    of.write(','.join(top_30_diffs))
    of.write('\n')


# five_fam = five_df.filter(['org', 'family']).drop_duplicates().set_index('family')
# six_fam = six_df.filter(['org', 'family']).drop_duplicates().set_index('family')
# print(five_fam.head())
# print(six_fam.head())
# com = five_fam.join(six_fam, lsuffix='_five', rsuffix='_six', how='outer')
# print(com.head())