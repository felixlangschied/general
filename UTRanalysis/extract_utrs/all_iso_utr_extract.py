from Bio import SeqIO
import json
import os

# genbank_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\human_data\GCF_000001405.39_GRCh38.p13_rna.gbff'
# prot = ['CDK2']

prot_path = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/hsa_rodentmis_genenames.txt'
tax_path = '/share/project/felixl/ncOrtho/data/UTR/taxid_name_assembly.tsv'
genbank_files_dir = '/share/project/felixl/ncOrtho/data/UTR/genbank_files'
out_dir = '/share/project/felixl/ncOrtho/data/UTR/all_iso_all_genes'

with open(prot_path, 'r') as fh:
    prot = fh.read().split('\n')

with open(tax_path, 'r') as fh:
    for line in fh:
        taxid, name, assembly = line.strip().split('\t')
        print('# Starting UTR extraction for {}'.format(name))
        outpath = '{}/{}_UTRs.json'.format(out_dir, assembly)
        if assembly == 'None':
            continue
        elif os.path.isfile(outpath):
            continue
        genbank_loc = '{}/{}_rna.gbff'.format(genbank_files_dir, assembly)
        out_dict = {}
        with open(genbank_loc) as handle:
            for record in SeqIO.parse(handle, "genbank"):
                nuc_seq = record.seq
                for feature in record.features:
                    if feature.type == 'CDS':
                        gene_name = feature.qualifiers['gene'][0]
                        loc = feature.location
                        utrseq = str(nuc_seq[loc.end:])
                        prot_id = feature.qualifiers['protein_id'][0]
                        # test if protein has alread an entry in the output
                        if gene_name not in out_dict:
                            out_dict[gene_name] = {}
                        out_dict[gene_name][prot_id] = utrseq
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        with open(outpath, 'w') as of:
            json.dump(out_dict, of)
