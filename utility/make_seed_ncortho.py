import os

table_path = "/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/MammalianVertebratesRefSeq_taxid2GCF.txt"
flatdir_path = "/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active"

ncortho_call = (
    "/home/felixl/PycharmProjects/ncOrtho/ncOrtho/ncortho.py "
    "-m /share/project/felixl/ncOrtho/data/hsa_ref/refseq/coreset/output/CMs "
    "-n /share/project/felixl/ncOrtho/data/hsa_ref/mirnas/hsa_highc_505.tsv "
    "-o /share/project/felixl/ncOrtho/data/hsa_ref/refseq/MammalianVertebrates "
    "-q {} "
    "-r /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/genomic.fna "
    "--cpu 4"
    "--queryname {}\n"
)

seed_out = "/share/project/felixl/ncOrtho/data/hsa_ref/refseq/seed_file.txt"


with open(table_path, 'r') as fh, open(seed_out, 'w') as of:
    for line in fh:
        data = line.strip().split('|')
        if data[-1] != 'active':
            continue
        name = data[1].replace(' ', '_')
        accession = data[2]
        filepath = f'{flatdir_path}/{accession}/genomic.fna'
        if not os.path.isfile(filepath):
            print(f'WARNING: No genome file found for {name}')
            continue
        else:
            of.write(ncortho_call.format(filepath, name))
print('# Done')



