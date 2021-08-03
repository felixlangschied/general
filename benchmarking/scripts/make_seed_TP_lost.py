import os

table_path = "/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/MammalianVertebratesRefSeq_taxid2GCF.txt"
flatdir_path = "/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active"
spec_path = '/home/felixl/PycharmProjects/general/benchmarking/data/spec_missing_TP.tsv'
seed_out = "/share/project2/felix/ncOrtho/mirgenedb/TPbenchmark/seed_file.txt"


ncortho_call = (
    "ncSearchcd / "
    "-m /share/project/felixl/ncOrtho/data/mirGeneDB/coreset/output/CMs "
    "-n /share/project/felixl/ncOrtho/data/mirGeneDB/data/missing_TPs.tsv "
    "-o /share/project2/felix/ncOrtho/mirgenedb/TPbenchmark "
    "-q {0} "
    "-r /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/genomic.fna "
    "--refblast /share/project2/felix/ncOrtho/mirgenedb/05cutoff/Homo_sapiens/data/Homo_sapiens.fa"
    "--cpu 4"
    "--queryname {1}"
    "--queryblast /share/project2/felix/ncOrtho/mirgenedb/05cutoff/{1}/data/{1}.fa\n"
)

specs = []
with open(spec_path, 'r') as fh:
    next(fh)
    for line in fh:
        specs.append(line.split('\t'))


with open(table_path, 'r') as fh, open(seed_out, 'w') as of:
    for line in fh:
        data = line.strip().split('|')
        if data[-1] != 'active':
            continue
        name = data[1].replace(' ', '_')
        if name not in specs:
            continue
        accession = data[2]
        filepath = f'{flatdir_path}/{accession}/genomic.fna'
        if not os.path.isfile(filepath):
            print(f'WARNING: No genome file found for {name}')
            continue
        else:
            of.write(ncortho_call.format(filepath, name))
print('# Done')