import subprocess as sp
import sys

mirna = 'Hsa-Mir-154-P4b'
species = 'Oryctolagus_cuniculus'


call = (
    'python /home/felixl/PycharmProjects/ncOrtho/ncOrtho/ncortho.py '
    '-m /share/project/felixl/ncOrtho/data/mirGeneDB/coreset/output/CMs '
    f'-n /home/felixl/project/ncOrtho/debugging/large_folding/test.tsv '
    f'-o /home/felixl/project/ncOrtho/debugging/large_folding '
    f'-q /share/project2/felix/ncOrtho/mirgenedb/05cutoff/{species}/data/{species}.fa '
    f'-r /share/project2/felix/ncOrtho/mirgenedb/05cutoff/Homo_sapiens/data/Homo_sapiens.fa '
    f'--queryname {species} '
    f'--queryblast /share/project2/felix/ncOrtho/mirgenedb/05cutoff/{species}/data/{species}.fa '
    f'--cpu 3 '
    f'--cleanup False '
    # f'--minlength 0.5 '
    # f'--heuristic False '
    # f'--cm_cutoff 0 '
    # f'--heur_blast_evalue 10 '
    # f'--heur_blast_length 0'
    '\n'
    )

sp.run(call, shell=True, stdout=sys.stdout)