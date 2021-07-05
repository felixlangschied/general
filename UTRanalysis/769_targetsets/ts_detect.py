import json
import subprocess as sp
import os
import pandas as pd

# TODO: check for NDRG1 (https://pubmed.ncbi.nlm.nih.gov/25081069/)

human_path = '/home/felixl/PycharmProjects/general/UTRanalysis/769_targetsets/target_utrs/GCF_000001405.39_GRCh38.p13_UTRs.json'
mouse_path = '/home/felixl/PycharmProjects/general/UTRanalysis/769_targetsets/target_utrs/GCF_000001635.27_GRCm39_UTRs.json'
hmm_path = '/home/felixl/PycharmProjects/general/UTRanalysis/769_targetsets/hsa-miR-769-5p.hmm'
meme_path = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/target_sites/ts_motifs/hsa-miR-769-5p.xml'
outdir = '/home/felixl/PycharmProjects/general/UTRanalysis/769_targetsets/data'
target_path = '/home/felixl/PycharmProjects/general/UTRanalysis/769_targetsets/potential_targets.tsv'

pot_ts = []
with open(target_path, 'r') as th:
    for line in th:
        pot_ts.append(line.strip().split('\t')[3])

if not os.path.isdir(outdir):
    os.mkdir(outdir)
# options: 'hmmer', 'meme', 'string'
analysis = 'string'


def string_analysis(utr_dict, strings):
    df_dict = {'gene': [], 'score': [], 'formula': []}
    for protname in utr_dict:
        # print(f'### human: {protname} ###')
        counter = 0
        for iso_count, seq in enumerate(utr_dict[protname].values(), 1):
            # print('# {}'.format(iso_id))
            for motif in strings:
                if motif in seq:
                   counter += 1
        if iso_count == 0:
            score = 'NaN'
            formula = 'NaN'
            continue
        else:
            formula = '{} / {}'.format(counter, iso_count)
            score = counter / iso_count
        df_dict['gene'].append(protname)
        df_dict['score'].append(score)
        df_dict['formula'].append(formula)
    df = pd.DataFrame.from_dict(df_dict)
    return df

# collect utr sequences and save to file
utr_file = f'{outdir}/all_human_utrs.fa'
if not os.path.isfile(utr_file):
    with open(human_path) as hp:
        hu_d = json.load(hp)
    for protname in hu_d:
        if protname not in pot_ts:
            continue
        print(f'### human: {protname} ###')
        # write UTRs to temporary fasta file
        with open(utr_file, 'a') as tmpf:
            for protid, seq in hu_d[protname].items():
                outstr = f'>{protid}\n{seq}\n'
                tmpf.write(outstr)

if analysis == 'hmmer':
    with open(human_path) as hp:
        hu_d = json.load(hp)
    # run HMMER search
    hmm_cmd = 'hmmsearch {} {}'.format(hmm_path, utr_file)
    rawres = sp.run(hmm_cmd, shell=True, capture_output=True)
    res = rawres.stdout.decode('utf-8')
    for line in res.split('\n'):
        if not line.startswith('#'):
            print(line)
elif analysis == 'meme':
    with open(human_path) as hp:
        hu_d = json.load(hp)
    fimo_cmd = 'fimo {} {}'.format(meme_path, utr_file)
    rawres = sp.run(fimo_cmd, shell=True, capture_output=True)
    res = rawres.stdout.decode('utf-8')
    err = rawres.stderr.decode('utf-8')
    if err:
        print('ERROR:')
        print(err)
    else:
        print(res)
# do simple string search
elif analysis == 'string':
    str_list = ['AGGTCTC', 'AGGTTTC']
    with open(human_path, 'r') as hh:
        raw_d = json.load(hh)
        human_d = {k: v for k, v in raw_d.items() if k in pot_ts}

    with open(mouse_path, 'r') as mh:
        raw_d = json.load(mh)
        mouse_d = {k.upper(): v for k, v in raw_d.items() if k.upper() in pot_ts}

    hu_df = string_analysis(human_d, str_list)
    mu_df = string_analysis(mouse_d, str_list)
    print(hu_df.head())
    print(hu_df.size)

    print(mu_df.head())
    print(mu_df.size)

    df = hu_df.join(mu_df.set_index('gene'), on="gene", lsuffix='_hsa', rsuffix='_mmu', how='inner')
    # df = df.dropna(how='any')
    # print(df.size)
    df = df.sort_values(by=['score_hsa', 'score_mmu'], ascending=[False, True])
    # print(df.head(30))

    above_1 = df['score_hsa'] >= 1
    fil_df = df[above_1].sort_values(by='score_mmu', ascending=False)
    print(fil_df)
    fil_df.to_csv(f'{outdir}/filtered_string_result.tsv', sep='\t')
    # print(f'{outdir}/string_result.tsv')






