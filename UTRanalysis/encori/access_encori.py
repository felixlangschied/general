import requests

mirna = 'hsa-miR-197-3p'
tarprot = 'all'
organism = 'human'
# organism = 'mouse'

payload = {
    'assembly': '',
    'geneType': 'mRNA',
    'miRNA': mirna,
    'clipExpNum': '1',
    'degraExpNum': '0',
    'pancancerNum': '0',
    'programNum': '1',
    'program': 'None',
    'target': tarprot,
    'cellType': 'all'
}
if organism == 'human':
    payload['assembly'] = 'hg19'
elif organism == 'mouse':
    payload['assembly'] = 'mm10'

base_url = 'http://starbase.sysu.edu.cn/api/miRNATarget'

r = requests.get(base_url, params=payload)

encori_target_prots = []
for line in r.text.split('\n'):
    if line and not line.startswith('#') :
        data = line.split('\t')
        print(data)
        if data[3] == 'geneName':
            continue
        encori_target_prots.append(data[3])

# print(encori_target_prots)
print('# Total target proteins: {}'.format(len(encori_target_prots)))
print('# Unique target proteins: {}'.format(len(set(encori_target_prots))))
