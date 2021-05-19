import matplotlib.pyplot as plt
import seaborn as sns
import subprocess as sp
from Bio.Seq import Seq
import os

# inpath = r'C:\Users\felix\project\UTRanalysis\target_analysis\target_sites.tsv'
inpath = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/target_sites/target_sites.tsv'
outdir = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/target_sites/align_data'
# mat_path = r'C:\Users\felix\project\UTRanalysis\target_analysis\all_mats.fa'
mat_path = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/seq_align/all_mats.fa'

out_dict = {}
mirna_set = set()
with open(inpath, 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')
        mat_id = data[1]
        gen_id = data[3]
        site_seq = data[6]
        target = data[3]
        if mat_id not in out_dict:
            tmp_dict = {'seq': site_seq, 'target': target}
            out_dict[mat_id] = [tmp_dict]
        else:
            tmp_dict = {'seq': site_seq, 'target': target}
            out_dict[mat_id].append(tmp_dict)
print(out_dict)

mat_dict = {}
# parse mature sequences
with open(mat_path, 'r') as fh:
    for line in fh:
        if line.startswith('>'):
            mat_id = line.strip().split()[0].replace('>', '')
        else:
            seq = line.strip()
            mat_dict[mat_id] = Seq(seq)
# print(mat_dict)

fa_files = []
for mat_mir in out_dict:
    mat_seq = mat_dict[mat_mir]
    if '-5p' in str(mat_mir):
        typ = '-5p'
        # other_matseq = mat_dict[mat_mir.replace('-5p', '-3p')]
    else:
        typ = '-3p'
        # other_matseq = mat_dict[mat_mir.replace('-3p', '-5p')]
    rc_mat = str(mat_seq.reverse_complement()).replace('U', 'T')
    # other_rc_mat = str(other_matseq.reverse_complement()).replace('U', 'T')
    outname = '{}/{}.fa'.format(outdir, mat_mir)
    fa_files.append(outname)
    with open(outname, 'w') as of:
        for c, t_dict in enumerate(out_dict[mat_mir]):
            target = out_dict[mat_mir][c]['target']
            tseq = out_dict[mat_mir][c]['seq']
            of.write(f'>ts{c+1}_{target}_{typ}\n')
            of.write(tseq)
            of.write('\n')
        # add reverse complement of target miRNA
        of.write(f'>rc_{mat_mir}\n')
        of.write(rc_mat)
        of.write('\n')
        # add reverse complement of the other mature miRNA
        # of.write(f'>rc_otherp_{mat_mir}\n')
        # of.write(other_rc_mat)
        # of.write('\n')
print('# FASTA files per mature miRNA created')

pre_mirna = [mat_mir.replace('-3p', '').replace('-5p', '') for mat_mir in out_dict.keys()]
pre_mirna = list(set(pre_mirna))

# create outputfolder for premirna specific alignments
if not os.path.isdir('{}/pre_mirna'.format(outdir)):
    os.mkdir('{}/pre_mirna'.format(outdir))
for pre_mir in pre_mirna:
    outname = '{}/pre_mirna/{}.fa'.format(outdir, pre_mir)
    fa_files.append(outname)
    # extract 5p information
    five_id = f'{pre_mir}-5p'
    five_mat = mat_dict[five_id]
    five_rc = str(five_mat.reverse_complement()).replace('U', 'T')
    # extract 3p information
    three_id = f'{pre_mir}-3p'
    three_mat = mat_dict[three_id]
    three_rc = str(three_mat.reverse_complement()).replace('U', 'T')

    with open(outname, 'w') as of:
        # write 5p
        ts_count = 1
        if five_id in out_dict:
            for c, t_dict in enumerate(out_dict[five_id]):
                target = out_dict[five_id][c]['target']
                tseq = out_dict[five_id][c]['seq']
                of.write(f'>ts{ts_count}_{target}_5p\n')
                of.write(tseq)
                of.write('\n')
                ts_count += 1
            # add reverse complement of target miRNA
            of.write(f'>rc_{five_id}\n')
            of.write(five_rc)
            of.write('\n')

        # write 3p
        if three_id in out_dict:
            for c, t_dict in enumerate(out_dict[three_id]):
                target = out_dict[three_id][c]['target']
                tseq = out_dict[three_id][c]['seq']
                of.write(f'>ts{ts_count}_{target}_3p\n')
                of.write(tseq)
                of.write('\n')
                ts_count += 1
            # add reverse complement of target miRNA
            of.write(f'>rc_{three_id}\n')
            of.write(three_rc)
            of.write('\n')

for file in fa_files:
    aln_name = file.replace('.fa', '.aln')
    aln_cmd = (
        f'muscle -in {file} -out {aln_name}'
    )
    sp.run(aln_cmd, shell=True)
print('# Alignment done')




##### visualization #####

# # plot length distribution
# plt.figure()
# plt.rcParams.update({'font.size': 14})
# target_lenghts = []
# for mat_mir in out_dict:
#     for target in out_dict[mat_mir]:
#         target_lenghts.append(len(target))
# sns.histplot(target_lenghts)
# plt.xlabel('Target site length [nt]')
# plt.show()

# # plot number of target sites per mat-mir
# plt.figure()
# plt.rcParams.update({'font.size': 14})
# data = [len(out_dict[mat_id]) for mat_id in out_dict]
# labels = out_dict.keys()
# plt.xticks(range(len(data)), labels, rotation=40, ha='right')
# plt.ylabel('Number of target sites')
# plt.bar(range(len(data)), data, edgecolor='k')
# plt.tight_layout()
# plt.show()