import glob
import subprocess as sp
import os
import json
import pandas as pd

spec_dir = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/spec_specific/'

ncortho = True
res_dir = '/home/felixl/project/ncOrtho/benchmark/ncortho'
outf = '/home/felixl/project/ncOrtho/benchmark/filtered_mirbh/analysis/blast_results_in_mirgenedb.json'



def run_blast(mirf):
    mkblast_raw = 'makeblastdb -in {0} -out {1} -dbtype nucl'
    # blast_raw = 'blastn -task megablast -query {0} -db {1} -outfmt "6 std qcovhsp" -perc_identity 90 -qcov_hsp_perc 80'
    blast_raw = 'blastn -task megablast -query {0} -dust "no" -db {1} -outfmt "6 std qcovhsp"'

    blast_db_path = mirf.replace('.fa', '')
    species = blast_db_path.split(os.sep)[-1]
    # make BLAST database of specieswise mirgenedb entries
    if not os.path.isfile(f'{blast_db_path}.nhr'):
        mkblast_cmd = mkblast_raw.format(mirf, blast_db_path)
        mk_res = sp.run(mkblast_cmd, shell=True, capture_output=True)
        if mk_res.returncode != 0:
            print(mk_res.stderr.decode('utf-8'))
    # run BLAST
    if ncortho:
        results = f'{res_dir}/{species}/{species}_orthologs.fa'
    else:
        results = f'{res_dir}/{species}_orthologs.fa'
    # count number of results
    res_count = 0
    with open(results, 'r') as rh:
        for line in rh:
            if line.startswith('>'):
                res_count += 1

    blast_cmd = blast_raw.format(results, blast_db_path)
    # read and return
    res = sp.run(blast_cmd, shell=True, capture_output=True)
    if res.returncode != 0:
        print(res.stderr.decode('utf-8'))
    else:
        result = res.stdout.decode('utf-8').split('\n')
        return result, res_count


def check_best_hit(results, ic, cc):
    ortho_check = set()
    reject_count = 0
    return_dict = {'species': [], 'miRNA': [], 'identity': [], 'coverage': [], 'confirmed': []}
    for line in results:
        if line:  # skip emtpy lines
            data = line.strip().split('\t')
            if data[0] not in ortho_check:  # only check best hit per result
                ortho_check.add(data[0])  # e.g. Aedes_aegypti|Hsa-Mir-137-P1-v2_pre_c1|NC_035108.1|1896120...
                species = data[0].split('|')[0]
                res_mirid = data[0].split('|')[1].split('_')[0]
                res_fam = '-'.join(res_mirid.split('-')[1:3])
                db_fam = '-'.join(data[1].replace('_pre', '').split('-')[1:3])
                ident = float(data[2])
                coverage = float(data[-1])
                if res_fam == db_fam:
                    # check identity and coverage
                    if ident >= ic and coverage >= cc:
                        return_dict['species'].append(species)
                        return_dict['miRNA'].append(res_mirid)
                        return_dict['coverage'].append(coverage)
                        return_dict['identity'].append(ident)
                        return_dict['confirmed'].append(True)
                    else:
                        print('WARNING: rejecting hit because of failed cutoffs')
                        print(data[0])
                        print(f'Identity: {ident}')
                        print(f'Coverage: {coverage}')
                        print(f'Result Chromosome: {data[0].split("|")[2]}')
                        print(f'Query miRNA: {data[1]}')
                        print('\n')
                        reject_count += 1

                        return_dict['species'].append(species)
                        return_dict['miRNA'].append(res_mirid)
                        return_dict['coverage'].append(coverage)
                        return_dict['identity'].append(ident)
                        return_dict['confirmed'].append(False)
                else:
                    print('WARNING: Not same family detected')
                    print('\n')
                    return_dict['species'].append(species)
                    return_dict['miRNA'].append(res_mirid)
                    return_dict['coverage'].append(coverage)
                    return_dict['identity'].append(ident)
                    return_dict['confirmed'].append(False)
                    reject_count += 1
    return return_dict, reject_count

            # print(line)


mirgene_files = glob.glob(f'{spec_dir}/*.fa')
# df = pd.DataFrame()
df_dict = {'species': [], 'miRNA': [], 'identity': [], 'coverage': [], 'confirmed': []}
total_res = 0
total_rejected = 0
for mirg_f in mirgene_files:
    species = mirg_f.split(os.sep)[-1].replace('.fa', '')
    print(species)
    blast_res, count = run_blast(mirg_f)
    total_res += count
    blast_dict, rej_count = check_best_hit(blast_res, 90, 80)
    total_rejected += rej_count
    for key in df_dict:
        df_dict[key].extend(blast_dict[key])
print(total_res)
print(total_rejected)

with open(outf, 'w') as of:
    json.dump(df_dict, of)
# df = pd.DataFrame.from_dict(df_dict)
# rejected = df[~df['confirmed']]
# print(df)
# print(rejected)
