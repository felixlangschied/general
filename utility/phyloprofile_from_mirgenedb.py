import pandas as pd

gff_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/ALL.gff'
org_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/organisms.txt'

# out = '/share/project/felixl/ncOrtho/data/mirGeneDB/phyloprofile.long'
out = '/home/felixl/Desktop/tmp/mirgenedb.long'


org_dict = {}
with open(org_path, 'r') as fh:
    for line in fh:
        if not line.startswith('#'):
            data = line.strip().split('\t')
            abb = data[0]
            taxid = data[-1]
            org_dict[abb] = taxid

# manual entries (not in miRBase)
org_dict['sto'] = 75743  # Scyliorhinus torazame
org_dict['cgi'] = 29159  # Crassostrea gigas
org_dict['efe'] = 6396  # Eisenia fetida
org_dict['bge'] = 6973  # Blattella germanica
org_dict['pfl'] = 63121  # Ptychodera flava
org_dict['lan'] = 7574  # Lingula anatina



mirid_list = []
mirid_dict = {'taxid': [], 'miRNA': [], 'mir_fam': []}
with open(gff_path, 'r') as fh:
    for line in fh:
        if not line.startswith('#'):
            data = line.strip().split('\t')
            if data[2] == 'pre_miRNA':
                mirid = data[-1].split('ID=')[1].split(';')[0].replace('_pre', '')
                mirid_list.append(mirid)
                spec = mirid.split('-')[0].lower()
                taxid = int(org_dict[spec])
                mirna = '-'.join(mirid.split('-')[1:])
                mir_fam = '-'.join(mirid.split('-')[1:3])

                # print(taxid, mirna, mir_fam)

                mirid_dict['taxid'].append(taxid)
                mirid_dict['miRNA'].append(mirna)
                mirid_dict['mir_fam'].append(mir_fam)


df = pd.DataFrame.from_dict(mirid_dict)
df = df.sort_values(by='taxid')
print(df.tail())


with open(out, 'w') as of:
    header = 'geneID\tncbiID\torthoID\n'
    of.write(header)
    print(header.strip())
    for ind in range(1, df.shape[0]):
        taxid = df.iloc[ind]['taxid']
        mir_fam = df.iloc[ind]['mir_fam']
        mirna = df.iloc[ind]['miRNA']
        outstring = f'{mir_fam}\tncbi{taxid}\t{mirna}\n'
        print(outstring.strip())
        of.write(outstring)









