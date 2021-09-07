

NCmap = '/home/felixl/project/ncOrtho/data/GCF_000001405.39_GRCh38.p13_assembly_report.txt'
gff_path = '/home/felixl/project/ncOrtho/data/hsa.gff'
pre_path = '/home/felixl/project/ncOrtho/data/hsa-pre.fas'
mat_path = '/home/felixl/project/ncOrtho/data/hsa-mat.fas'

outfile = '/home/felixl/project/ncOrtho/debugging/mat_mirgenedb.tsv'


def read_seqs(path):
    seq_dict = {}
    with open(path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                header = line.strip().replace('>', '').split('_')[0]
            else:
                seq = line.strip()
                seq_dict[header] = seq
    return seq_dict


pre_dict = read_seqs(pre_path)
mat_dict = read_seqs(mat_path)

chrom_map = {}
with open(NCmap, 'r') as fh:
    for line in fh:
        if not line.startswith('#'):
            data = line.strip().split()
            nc = data[6]
            chr = data[-1]
            chrom_map[chr] = nc



with open(gff_path, 'r') as fh, open(outfile, 'w') as of:
    for line in fh:
        if not line.startswith('#'):
            data = line.strip().split()
            if data[2] == 'pre_miRNA':
                start = data[3]
                end = data[4]
                strand = data[6]
                mirid = data[-1].split('ID=')[1].split(';')[0].replace('_pre', '')
                chrom = chrom_map[data[0]]
                pre = pre_dict[mirid]
                mat = mat_dict[mirid]
                outstring = f'{mirid}\t{chrom}\t{start}\t{end}\t{strand}\t{pre}\t{mat}'
                print(outstring)
                of.write(outstring)
                of.write('\n')

