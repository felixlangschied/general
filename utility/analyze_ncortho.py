import os
import subprocess as sp
import sys
import argparse

from analyze_tools import create_overview
from analyze_tools import make_phyloprofile
from analyze_tools import extract_representative
from analyze_tools import align_seqs
from analyze_tools import make_supermatrix


# Define input
analysis_name = 'mirgenedb'

# directory that contains results in the format dir/species/mirna/files.fa
result_path = '/share/project2/felix/ncOrtho/mirGeneDB/05cutoff'

# path to mapping file for phyloprofile
map_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/available_specs.txt'

# output
out_path = '/home/felixl/project/ncOrtho/analyses/diff_sets/{}'.format(analysis_name)
over_file = '{}/overview.txt'.format(out_path)
tree_out = '/home/felixl/project/ncOrtho/analyses/diff_sets/{}/base_tree'.format(analysis_name)
if not os.path.isdir(tree_out):
    os.makedirs(tree_out)
# call to iqtree
iqtree_cmd = '/share/applications/bin/iqtree -s {} -bb 1000 -alrt 1000 -nt AUTO -redo -pre {}/{}'

# species that should not be included in the supermatrix tree analysis
spec_to_skip = ['Petromyzon_marinus']

# tools: mafft or muscle
tool = 'mafft'

# input mirnas as "tsv" or list of "ids"
mode = 'tsv'
mirna_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/mirgenedb.tsv'


def main():
    # Print header
    print('\n'+'#'*39)
    print('###'+' '*33+'###')
    print('###   Analysis of ncOrtho results   ###')
    print('###'+' '*33+'###')
    print('#'*39+'\n')

    # # Parse command-line arguments
    # # Define global variables
    # parser = argparse.ArgumentParser(
    #     description='Analzye ncOrtho results. Gives: PhyloProfile input, supermatrix alignment, species tree'
    # )
    # # covariance models folder
    # parser.add_argument(
    #     '-i', '--input', metavar='<path>', type=str,
    #     help='Path to ncOrtho output directory'
    # )
    # parser.add_argument(
    #     '-n', '--name', metavar='str', type=str,
    #     help='Name of the analysis'
    # )
    # parser.add_argument(
    #     '-i', '--input', metavar='<path>', type=str,
    #     help='Path to ncOrtho output directory'
    # )
    # parser.add_argument(
    #     '-o', '--output', metavar='<path>', type=str,
    #     help='Path to location where output of analysis should be stored'
    # )

    # if overview exists read it
    if os.path.isfile(over_file):
        ortholog_files = {}
        with open(over_file, 'r') as inf:
            for line in inf:
                key, value = line.strip().split()
                ortholog_files[key] = value
    else:
        # create overview from result directory
        ortholog_files = create_overview(result_path, out_path)

    print('# Starting filtering for miRNAs and species')
    # determine miRNAs to analyse
    mirna_ids = []
    if mode == 'tsv':
        with open(mirna_path, 'r') as inf:
            for line in inf:
                mirna_ids.append(line.strip().split()[0])
    elif mode == 'ids':
        with open(mirna_path, 'r') as infile:
            mirna_ids = infile.read().split()
    # determine species to analyse
    all_specs = [path.split('/')[-2] for path in ortholog_files]

    all_specs = set(all_specs)
    specs_to_analyse = all_specs - set(spec_to_skip)
    # filter overview dictionary
    filtered_orthologs = {}
    mirnas_found = set()
    # print(mirna_ids)

    for path in ortholog_files:
        outmirna = path.split('/')[-1].replace('_orthologs.fa', '')
        if outmirna in mirna_ids and path.split('/')[-2] in specs_to_analyse:
            mirnas_found.add(path.split('/')[-1].replace('_orthologs.fa', ''))
            filtered_orthologs[path] = ortholog_files[path]
    print('# Done')

    # double check
    if not filtered_orthologs:
        print('# Nothing left after filtering. Exiting..')
        sys.exit()

    # make phyloprofile input
    # 106582|Maylandia zebra|GCF_000238955.4|MAYZE|active
    name_2_id = {}
    with open(map_path, 'r') as mapfile:
        for line in mapfile:
            taxid, name, assem, bla, status = line.strip().split('|')
            name = name.replace(' ', '_')
            name_2_id[name] = taxid

    make_phyloprofile(filtered_orthologs, name_2_id, analysis_name, out_path)

    # writing fasta files of representatives to directory in out_path
    multi_paths = extract_representative(filtered_orthologs, mirnas_found, result_path, out_path)

    # align representative multifasta files using MAFFT-linsi
    align_seqs(multi_paths, out_path, tool)

    # call Ingo's concat alignment and degap scripts
    superm_path = make_supermatrix(out_path)

    # calculate tree using iqtree
    print('# Starting tree calculation')
    tree_cmd = iqtree_cmd.format(superm_path, tree_out, analysis_name)
    res = sp.run(tree_cmd, shell=True, capture_output=True)
    print(res.stdout.decode('utf-8'))
    if res.returncode != 0:
        print(res.stderr.decode('utf-8'))


if __name__ == '__main__':
    main()