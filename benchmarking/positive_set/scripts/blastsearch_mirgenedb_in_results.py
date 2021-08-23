import glob
import subprocess as sp
import os

spec_dir = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/spec_specific/'
mkblast_raw = 'makeblastdb -in {0} -out {0} -dbtype nucl'
blast_raw = 'blastn -task megablast -query {0} -db {1} -outfmt "6 std qcovhsp" -perc_identity 90 -qcov_hsp_perc 80 -out {2}'

# ncortho = True
# res_dir = '/home/felixl/project/ncOrtho/benchmark/ncortho'
# out_dir = '/home/felixl/project/ncOrtho/benchmark/ncortho/analysis/TP_blast'

ncortho = False
res_dir = '/home/felixl/project/ncOrtho/benchmark/filtered_mirbh'
out_dir = '/home/felixl/project/ncOrtho/benchmark/filtered_mirbh/analysis/TP_blast'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

spec_files = glob.glob(f'{spec_dir}/*')
print(spec_files)

for file in spec_files:
    name = '_'.join(file.split('/')[-1].split('_')[:2]).replace('.fa', '')
    print(name)
    if ncortho:
        results = f'{res_dir}/{name}/{name}_orthologs.fa'
    else:
        results = f'{res_dir}/{name}_orthologs.fa'
    # print(results)
    mkblast_cmd = mkblast_raw.format(results)
    mk_res = sp.run(mkblast_cmd, shell=True, capture_output=True)

    if mk_res.returncode != 0:
        print(mk_res.stderr.decode('utf-8'))
    outfile = f'{out_dir}/{name}_blastout.txt'
    blast_cmd = blast_raw.format(file, results, outfile)
    res = sp.run(blast_cmd, shell=True, capture_output=True)
    if res.returncode != 0:
        print(res.stderr.decode('utf-8'))
