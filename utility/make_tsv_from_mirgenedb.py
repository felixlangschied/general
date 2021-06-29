

NCmap = '/share/project/felixl/ncOrtho/data/hsa_ref/refseq/data/GCF_000001405.39_GRCh38.p13_assembly_report.txt'
gff_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/human.gff'
pre_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/human_precursor.fa'

outfile = '/share/project/felixl/ncOrtho/data/hsa_ref/mirnas/mirgenedb.tsv'

chrom_map = {}
with open(NCmap, 'r') as fh:
    for line in fh:
        if not line.startswith('#'):
            data = line.strip().split()
            nc = data[6]
            chr = data[-1]
            chrom_map[chr] = nc

seq_dict = {}
with open(pre_path, 'r') as fh:
    for line in fh:
        if line.startswith('>'):
            header = line.strip().replace('>', '')
        else:
            seq = line.strip()
            seq_dict[header] = seq

with open(gff_path, 'r') as fh, open(outfile, 'w') as of:
    for line in fh:
        if not line.startswith('#'):
            data = line.strip().split()
            if data[2] == 'pre_miRNA':
                start = data[3]
                end = data[4]
                strand = data[6]
                mirid = data[-1].split('ID=')[1].split(';')[0]
                chrom = chrom_map[data[0]]
                seq = seq_dict[mirid]
                outstring = f'{mirid}\t{chrom}\t{start}\t{end}\t{strand}\t{seq}\tNA'
                print(outstring)
                of.write(outstring)
                of.write('\n')

