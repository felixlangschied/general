import os

spec_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/available_specs.txt'
mammal_path = '/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/{}/genomic.fna'
nonmammal_path = '/share/gluster/GeneSets/NCBI-Genomes/nonMammalianVertebratesRefSeq/raw_dir/active_June21/{}/genomic.fna'

seed_out = '/share/project2/felix/ncOrtho/06seed.txt'

call = ('time ncSearch '
        '-m /share/project/felixl/ncOrtho/data/mirGeneDB/coreset/output/CMs '
        '-n /share/project/felixl/ncOrtho/data/mirGeneDB/data/mirgenedb.tsv '
        '-o /share/project2/felix/ncOrtho/mirgenedb/06cutoff '
        '-q {0} '
        '-r /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/genomic.fna '
        '--cpu=4 '
        '--queryname={1} '
        '--queryblast=/share/project2/felix/ncOrtho/mirGeneDB/05cutoff/{1}/data/{1}.fa '
        '--cm_cutoff=0.6 '
        '--refblast=/share/project/felixl/ncOrtho/data/hsa_ref/refseq/test/data/genomic.fna')

with open(spec_path, 'r') as fh, open(seed_out, 'w') as of:
    for line in fh:
        taxid, rawname, assembly, asid, status = line.strip().split('|')
        name = rawname.replace(' ', '_')

        mammal = mammal_path.format(assembly)
        nomammal = nonmammal_path.format(assembly)
        # try to find query file
        if os.path.isfile(mammal):
            outcall = call.format(mammal, name)
        elif os.path.isfile(nomammal):
            outcall = call.format(nomammal, name)
        else:
            print('no genome file found for {}'.format(name))
        of.write(outcall)
        of.write('\n')
