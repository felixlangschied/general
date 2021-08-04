import glob
import subprocess as sp

spec_dir = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/spec_specific/'
res_dir = '/share/project2/felix/ncOrtho/mirgenedb/05cutoff/pooled_results'
out_dir = '/home/felixl/PycharmProjects/general/benchmarking/data/TP_blast/{}_blastout.txt'


mkblast_raw = 'makeblastdb -in {0}.fa -out {0} -dbtype nucl'
blast_raw = 'blastn -task megablast -query {0} -db {1} -outfmt "6 std qcovhsp" -perc_identity 90 -qcov_hsp_perc 80 -out {2}'

spec_files = glob.glob(f'{spec_dir}/*')
print(spec_files)

for file in spec_files:
    name = '_'.join(file.split('/')[-1].split('_')[:2])
    print(name)
    results = f'{res_dir}/{name}'
    # print(results)
    mkblast_cmd = mkblast_raw.format(results)
    mk_res = sp.run(mkblast_cmd, shell=True, capture_output=True)

    if mk_res.returncode != 0:
        print(mk_res.stderr.decode('utf-8'))

    blast_cmd = blast_raw.format(file, results, out_dir.format(name))
    res = sp.run(blast_cmd, shell=True, capture_output=True)
    if res.returncode != 0:
        print(res.stderr.decode('utf-8'))