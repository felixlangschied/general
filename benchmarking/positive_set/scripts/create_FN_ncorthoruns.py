import json
import os
import subprocess as sp
import sys

mis_fam = '/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/missing_fams_per_species.json'
mirna_path = '/share/project2/felix/ncOrtho/mirgenedb/data/mirgenedb.tsv'
outdir = '/home/felixl/project/ncOrtho/debugging'

with open(mis_fam, 'r') as fh:
    missing = json.load(fh)

# # make miRNA input files for each species
# mirna_dict = {}
# with open(mirna_path) as fh:
#     for line in fh:
#         fam = '-'.join(line.split()[0].split('-')[1:3]).replace('_pre', '')
#         if fam not in mirna_dict:
#             mirna_dict[fam] = []
#         mirna_dict[fam].append(line)
#
out_mir_dir = f'{outdir}/mirna_data'
# if not os.path.isdir(out_mir_dir):
#     os.makedirs(out_mir_dir)
# for species in missing:
#     if missing[species]:
#         with open(f'{out_mir_dir}/{species}.tsv', 'w') as of:
#             for missing_fam in missing[species]:
#                 for line in mirna_dict[missing_fam]:
#                     of.write(line)

with open('/home/felixl/project/ncOrtho/debugging/seed.txt', 'w') as seed:
    for species in missing:
        call = (
            'python /home/felixl/PycharmProjects/ncOrtho/ncOrtho/ncortho.py '
            '-m /share/project/felixl/ncOrtho/data/mirGeneDB/coreset/output/CMs '
            f'-n {out_mir_dir}/{species}.tsv '
            f'-o /home/felixl/project/ncOrtho/debugging/results '
            f'-q /share/project2/felix/ncOrtho/mirgenedb/05cutoff/{species}/data/{species}.fa '
            f'-r /share/project2/felix/ncOrtho/mirgenedb/05cutoff/Homo_sapiens/data/Homo_sapiens.fa '
            f'--queryname {species} '
            f'--queryblast /share/project2/felix/ncOrtho/mirgenedb/05cutoff/{species}/data/{species}.fa '
            f'--cpu 4 '
            f'--minlength 0.5 '
            f'--heuristic False '
            # f'--cm_cutoff 0 '
            # f'--heur_blast_evalue 10 '
            # f'--heur_blast_length 0'
            '\n'
        )
        seed.write(call)

        # res = sp.run(call, shell=True, stdout=sys.stdout)

    # print(call)
    # sp.check_output(call, stdout=sys.stdout, stderr=sp.STDOUT)
    # print(res)

    # with sp.Popen(call, stdout=sp.PIPE, bufsize=1, universal_newlines=True) as p:
    #     for line in p.stdout:
    #         print(line, end='')  # process line here
    #
    # if p.returncode != 0:
    #     raise sp.CalledProcessError(p.returncode, p.args)

