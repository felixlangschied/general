import json
import pandas as pd
import subprocess as sp
import sys

ncortho_path = '/home/felixl/PycharmProjects/general/def_benchmark/data/ncortho.json'

# load results
with open(ncortho_path, 'r') as np:
    nc_dict = json.load(np)
    nc_dict['species'] = [spec.replace('>', '') for spec in nc_dict['species']]
    ncdf = pd.DataFrame.from_dict(nc_dict)
# extract reference secondary structures
refdf = ncdf[ncdf['species'] == 'Homo_sapiens']
refdf = refdf.filter(['mirna', 'scheme'])
# join
comdf = ncdf.join(refdf.set_index('mirna'), on='mirna', lsuffix='', rsuffix='_ref')
print(comdf)
print(comdf.columns)


# run RNAdist
def run_rnadist(df, index):
    refi = df.columns.get_loc('scheme_ref')
    orthi = df.columns.get_loc('scheme')

    call = sp.Popen('RNAdistance', shell=True, stdin=sp.PIPE,
            stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8')
    instring = f'{df.iloc[index, refi]}\n{df.iloc[index, orthi]}'
    res, err = call.communicate(instring)
    # print(instring)
    if err:
        print(f'ERROR: {err}')
        sys.exit()
    else:
        # print(res)
        if res:
            distance = res.split(' ')[1]
        else:
            distance = 'NaN'
    return distance

scores = []
for rownum in range(comdf['species'].size):
    scores.append(run_rnadist(comdf, rownum))
comdf['distance'] = scores

comdf.to_csv('/home/felixl/PycharmProjects/general/def_benchmark/data/ncortho_wDistance.tsv', sep='\t')
print('# Done')