import subprocess as sp
import sys

mampath = '/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/MammalianVertebratesRefSeq_taxid2GCF.txt'
nonpath = '/share/gluster/GeneSets/NCBI-Genomes/nonMammalianVertebratesRefSeq/raw_dir/non-mammalianVertebrates-RefSeq_taxid2GCF.txt'
inpath = '/share/gluster/GeneSets/NCBI-Genomes/InvertebratesRefSeq/raw_dir/InvertebratesRefSeq_taxid2GCF.txt'
spec_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb_specs_abb.txt'

out_dir = '/home/felixl/project/ncOrtho/benchmark'

spec_map = {}
spec_list = []
with open(spec_path) as fh:
    for line in fh:
        data = line.strip().split('\t')
        if data[0]:
            species = data[1]
            abb = data[2]
            spec_map[abb] = species
            spec_list.append(species)

# print(spec_list)
remove_dict = {
    'Canis familiaris': 'Canis lupus familiaris',
    'Sarcophilus harrissii': 'Sarcophilus harrisii',
    'Strongylocentrotus purpuratus ': 'Strongylocentrotus purpuratus'
}
for key, value in remove_dict.items():
    spec_list.remove(key)
    spec_list.append(value)

specsfound = []
seed_dict = {}
mammals = '/share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/{}/genomic.fna'
mapfile = f'{out_dir}/mirgenedb_specs.txt'
with open(mapfile, 'w') as of:
    with open(mampath, 'r') as fh:
        for line in fh:
            data = line.strip().split('|')
            if data[-1] == 'active':
                # specname = ' '.join(data[1].split(' ')[:2])
                # print(specname)
                specname = data[1]
                if specname in spec_list:
                    specsfound.append(specname)
                    access = data[2]
                    specname = specname.replace(' ', '_')
                    seed_dict[specname] = mammals.format(access)
                    taxid = data[0]
                    of.write(f'{taxid}\t{specname}\n')
    nonmams = '/share/gluster/GeneSets/NCBI-Genomes/nonMammalianVertebratesRefSeq/raw_dir/active_June21/{}/genomic.fna'
    with open(nonpath, 'r') as fh:
        for line in fh:
            data = line.strip().split('|')
            if data[-1] == 'active':
                # specname = ' '.join(data[1].split(' ')[:2])
                specname = data[1]
                # print(specname)
                if specname in spec_list:
                    specsfound.append(specname)
                    access = data[2]
                    specname = specname.replace(' ', '_')
                    seed_dict[specname] = nonmams.format(access)
                    taxid = data[0]
                    of.write(f'{taxid}\t{specname}\n')
    inverts = '/share/gluster/GeneSets/NCBI-Genomes/InvertebratesRefSeq/raw_dir/active/{}/genomic.fna'
    with open(inpath, 'r') as fh:
        for line in fh:
            data = line.strip().split('|')
            if data[-1] == 'active':
                # specname = ' '.join(data[1].split(' ')[:2])
                specname = data[1]
                # print(specname)
                if specname in spec_list:
                    specsfound.append(specname)
                    access = data[2]
                    specname = specname.replace(' ', '_')
                    seed_dict[specname] = inverts.format(access)
                    taxid = data[0]
                    of.write(f'{taxid}\t{specname}\n')

print(len(specsfound))
print(len(spec_list))
print(set(spec_list) - set(specsfound))
# print(spec_list)





# mirbh_call = (
#     'python /home/felixl/PycharmProjects/mirbh/mirbh/mirna_rbh.py '
#     '-n /share/project2/felix/ncOrtho/mirgenedb/data/mirgenedb.tsv '
#     '-o /home/felixl/project/ncOrtho/benchmark/mirbh '
#     '-q /home/felixl/project/ncOrtho/benchmark/ncortho/{0}/data/{0}.fa '
#     '-r /home/felixl/project/ncOrtho/benchmark/ncortho/Homo_sapiens/data/Homo_sapiens.fa '
#     '--queryname {0} '
#     '--cpu 4'
#     '\n'
# )
# with open(f'/home/felixl/project/ncOrtho/benchmark/mirbh/seed.txt', 'w') as sh:
#     for qname in seed_dict:
#         outstr = mirbh_call.format(qname)
#         sh.write(outstr)


# ncortho_call = (
#     'ncSearch '
#     '-m /home/felixl/project/ncOrtho/benchmark/CMs '
#     '-n /home/felixl/project/ncOrtho/benchmark/mirgenedb.tsv '
#     '-o /home/felixl/project/ncOrtho/benchmark/ncortho '
#     '-q {} '
#     '-r  /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/genomic.fna '
#     '--refblast /home/felixl/project/ncOrtho/benchmark/refblast/Homo_sapiens.fa '
#     '--queryname {} '
#     '--cpu 4 '
#     '\n'
#     )
# with open(f'/home/felixl/project/ncOrtho/benchmark/ncortho/seed.txt', 'w') as sh:
#     for qname in seed_dict:
#         outstr = ncortho_call.format(seed_dict[qname], qname)
#         sh.write(outstr)



# test_call = (
#     'python /home/felixl/PycharmProjects/ncOrtho/ncOrtho/ncortho.py '
#     '-m /home/felixl/project/ncOrtho/benchmark/CMs '
#     # '-n /home/felixl/project/ncOrtho/benchmark/mirgenedb.tsv '
#     '-n /home/felixl/project/ncOrtho/benchmark/test.tsv '
#     '-o /home/felixl/project/ncOrtho/tmp_test '
#     '-q {} '
#     '-r  /share/gluster/GeneSets/NCBI-Genomes/MammalianVertebratesRefSeq/raw_dir/active/GCF_000001405.39/genomic.fna '
#     '--refblast /home/felixl/project/ncOrtho/benchmark/refblast/Homo_sapiens.fa '
#     '--queryname {} '
#     '--cpu 3 '
#     '\n'
#     )
#
#
# for qname in seed_dict:
#     call = test_call.format(seed_dict[qname], qname)
#     # sh.write(outstr)
#     if qname == 'Oryctolagus_cuniculus':
#         res = sp.run(call, shell=True, stdout=sys.stdout)
#     break