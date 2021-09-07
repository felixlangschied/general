
# rodent_path = r'C:\Users\felix\PycharmProjects\general\rodent_missing\rodent_taxids.txt'
# mammal_path = r'C:\Users\felix\PycharmProjects\general\rodent_missing\MammalianVertebratesRefSeq_taxid2GCF.txt'

rodent_path = '/home/felixl/project/ncOrtho/analyses/rodent_missing/rodent_taxids.txt'
mammal_path = '/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/MammalianVertebratesRefSeq_taxid2GCF.txt'
out = '/home/felixl/project/ncOrtho/analyses/rodent_missing/seed.txt'

ncortho_call = (
    'ncSearch '
    '-m /home/felixl/project/ncOrtho/benchmark/CMs '
    '-n /home/felixl/project/ncOrtho/analyses/rodent_missing/candidats_wControl.tsv '
    '-o /home/felixl/project/ncOrtho/analyses/rodent_missing/results '
    '-q {} '
    '-r  /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/genomic.fna '
    '--refblast /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/GCF_000001405.39 '
    '--queryname {} '
    '--queryblast {} '
    '--cpu 4 '
    '\n'
    )

with open(rodent_path, 'r') as rh:
    rodents = rh.read().split('\n')
    # add human as reference
    rodents.append('9606')

with open(mammal_path) as mh, open(out, 'w') as of:
    for line in mh:
        taxid, name, assembly, asname, status = line.strip().split('|')
        if status == 'active' and taxid in rodents:
            qpath = f'/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/{assembly}/genomic.fna'
            qname = name.replace(' ', '_')
            qblast = f'/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/{assembly}/{assembly}'
            seed_str = ncortho_call.format(qpath, qname, qblast)
            of.write(seed_str)
            print(seed_str)