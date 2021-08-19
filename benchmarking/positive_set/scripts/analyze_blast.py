import glob
from collections import Counter
import os

# blast_dir = r'C:\Users\felix\PycharmProjects\general\benchmarking\positive_set\data\TP_blast'
# blast_dir = '/home/felixl/project/ncOrtho/benchmark/ncortho/analysis/TP_blast'
blast_dir = '/home/felixl/project/ncOrtho/benchmark/filtered_mirbh/analysis/TP_blast'

blast_files = glob.glob(f'{blast_dir}/*')

best_hit = set()
accepted = {}
rejected = {}
for file in blast_files:
    species = file.split(os.sep)[-1].replace('_blastout.txt', '')
    accepted[species] = []
    rejected[species] = []
    with open(file, 'r') as fh:
        for line in fh:
            # keep only best hits
            data = line.strip().split('\t')
            if not data[0] in best_hit:
                best_hit.add(data[0])
                db_fam = '-'.join(data[0].split('-')[1:3]).replace('_pre', '')
                res_fam = '-'.join(data[1].split('-')[1:3]).split('_')[0]

                identity = float(data[2])
                coverage = float(data[-1])

                # check that hit is from same family
                if db_fam == res_fam:
                    # check for identity and coverage
                    if identity >= 95 and coverage >= 90:
                        accepted[species].append(db_fam)
                    else:
                        rejected[species].append(db_fam)
                        print(f'{species}: {data[0]}')
                else:
                    rejected[species].append(db_fam)
                    print(f'{species}: {data[0]}')
    accepted[species] = dict(Counter(accepted[species]))
    rejected[species] = dict(Counter(rejected[species]))


# print(accepted)
counter = 0
for species in accepted:
    counter += len(accepted[species].values())
    # print(len(accepted[species].values()))
print(counter)

# print(rejected)
# for species in rejected:
#     if rejected[species]:
#         print(f'{species}: {rejected[species]}')