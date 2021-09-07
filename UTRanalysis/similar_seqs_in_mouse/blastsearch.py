import subprocess as sp
import pyfaidx
import tempfile
import pandas as pd
import glob
import os


matpath = '/home/felixl/PycharmProjects/general/UTRanalysis/similar_seqs_in_mouse/mat_seqs.fa'
blastdb = '/home/felixl/project/ncOrtho/benchmark/ncortho/Mus_musculus/data/Mus_musculus.fa'
# blastdb = '/home/felixl/PycharmProjects/general/UTRanalysis/similar_seqs_in_mouse/mmu_mat'
# blastdb = '/home/felixl/PycharmProjects/general/UTRanalysis/similar_seqs_in_mouse/data/mmu_hairpin.fa'
mmu_gff = '/home/felixl/PycharmProjects/general/UTRanalysis/similar_seqs_in_mouse/mmu.gff3'
chrom_map = '/home/felixl/PycharmProjects/general/UTRanalysis/similar_seqs_in_mouse/GCF_000001635.27_GRCm39_assembly_report.txt'


def run_blast(seq):
    blast_command = (
                'blastn -task blastn-short -db {0} '
                '-num_threads {1} -evalue {2} '
                '-outfmt "6 std sstrand sseq"'.format(blastdb, 4, 1)
            )
    blast_call = sp.Popen(
        blast_command, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    seq = seq.replace('U', 'T')
    res, err = blast_call.communicate(seq)
    if err:
        print(f'ERROR: {err}')
    return res

# read chrom map
nc_2_chrom = {}
with open(chrom_map, 'r') as ch:
    for line in ch:
        if not line.startswith('#'):
            mapdata = line.strip().split()
            if mapdata[1] == 'assembled-molecule':
                # print(mapdata)
                nc = mapdata[6]
                chrom = mapdata[0]
                nc_2_chrom[nc] = chrom

# read mouse mirbase gff
mirbase_dict = {}
with open(mmu_gff, 'r') as gff:
    for line in gff:
        if not line.startswith('#'):
            linedata = line.strip().split()
            if linedata[2] == 'miRNA_primary_transcript':
                chrom = linedata[0].replace('chr', '')
                start = linedata[3]
                end = linedata[4]
                strand = linedata[6]
                mirid = linedata[-1].split(';')[-1].replace('Name=', '')
                if chrom not in mirbase_dict:
                    mirbase_dict[chrom] = {}
                mirbase_dict[chrom][mirid] = f'{start}|{end}|{strand}'
print(mirbase_dict)

genes = pyfaidx.Fasta(blastdb)
df_dict = {'mirna': [], 'blasthitseq': [], 'score': [], 'seq': [], 'scheme': []}
with open(matpath, 'r') as fh:
    extraregion = 50
    hit_at_start = False
    hit_at_end = False

    for line in fh:
        if line.startswith('>'):
            header = line.strip()
            mirid = header.replace('>', '').split()[0]
            print(f'# {mirid}')
        else:
            seq = line.strip()
            # print(seq.replace('U', 'T'))
            print(seq)
            # print(f'length: {len(seq)}')
            resblop = run_blast(seq)
            for count, result in enumerate(resblop.split('\n')):
                # print(result)
                # RNAfold
                if result:
                    res_data = result.split()
                    hitseq = res_data[-1]

                    # chrom = nc_2_chrom[res_data[1]]
                    chrom = res_data[1]
                    start = int(res_data[8])
                    end = int(res_data[9])
                    if res_data[12] == 'plus':
                        strand = '+'
                    else:
                        strand = '-'
                    # extract flanking regions
                    if start > extraregion:
                        n_start = start - extraregion
                    else:
                        # hit starts at the beginning of the chromosome
                        n_start = 0
                        hit_at_start = True
                    n_end = end + extraregion
                    if n_end > genes[chrom][-1].end:
                        # hit is at the end of the chromosome
                        n_end = genes[chrom][-1].end
                        hit_at_end = True
                    if strand == '+':
                        ext_seq = genes[chrom][n_start:n_end].seq
                    elif strand == '-':
                        ext_seq = genes[chrom][n_start:n_end].reverse.complement.seq
                    with tempfile.NamedTemporaryFile(mode='w+') as fp:
                        header = f'{mirid}_res{count}'
                        fp.write(f'>{header}\n{ext_seq}\n')
                        fp.seek(0)
                        fold_cmd = f'RNAfold -i {fp.name}'
                        res = sp.run(fold_cmd, shell=True, capture_output=True)
                        out = res.stdout.decode('utf-8').strip()
                        folddata = out.split('\n')
                        schemscore = folddata[-1]
                        scheme, score = schemscore.split(' (')
                        score = score.replace(')', '')
                        df_dict['mirna'].append(mirid)
                        df_dict['blasthitseq'].append(hitseq)
                        df_dict['score'].append(score)
                        df_dict['seq'].append(ext_seq)
                        df_dict['scheme'].append(scheme)
psfiles = glob.glob(f'{os.getcwd()}/*.ps')
for file in psfiles:
    os.remove(file)

df = pd.DataFrame.from_dict(df_dict)
print(df)
df.to_csv('/home/felixl/PycharmProjects/general/UTRanalysis/similar_seqs_in_mouse/rnafold.txt', sep='\t', index=False)





                # # check in gff file
                #     for mirid, location in mirbase_dict[chrom].items():
                #         db_start, db_end, db_strand = location.split('|')
                #         db_start = int(db_start)
                #         db_end = int(db_end)
                #         if db_strand != strand:
                #             continue
                #         # first within second
                #         if (
                #                 (start <= db_start <= end)
                #                 or (start <= db_end <= end)
                #         ):
                #             print(mirid)
                #         # second within first
                #         elif (
                #                 (db_start <= start <= db_end)
                #                 or (db_start <= end <= db_end)
                #         ):
                #             print(mirid)
                #         elif (
                #                 (start >= db_start)
                #                 and (end <= db_end)
                #         ):
                #             print(mirid)

# PROBLEM: mirbase uses old assembly, cannot map location
# extract flanking regions of hit and detect folding energies, search for hairpin
# blastsearch in mmu-mirbase pre-miRNAs