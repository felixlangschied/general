import os

# all_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb\all_mirgenedb.fa'
# abb_path = r'C:\Users\felix\Desktop\mirgene_db.CSV'
# out_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb\{}_mirgenedb.fa'
# spec_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\spec_missing_TP.tsv'

all_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/all_mirgenedb.fa'
abb_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/supp_tab.tsv'
out_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/spec_specific/{}_mirgenedb.fa'
spec_path = '/home/felixl/PycharmProjects/general/benchmarking/data/spec_missing_TP.tsv'

map_dict = {}
with open(abb_path, 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')[1:3]
        if data:
            name, abb = data
            name = name.replace(' ', '_')
            abb = abb.capitalize()
            map_dict[abb] = name

specs = []
with open(spec_path, 'r') as fh:
    next(fh)
    for line in fh:
        specs.append(line.split('\t')[0])


with open(all_path, 'r') as fh:
    data = fh.read().split('>')
    data.pop(0)
    for entry in data:
        print(entry)
        db_abb = entry[0:3]
        name = map_dict[db_abb]
        out_file = out_path.format(name)
        if name not in specs:
            continue
        parsed_entry = '>' + entry.replace(db_abb, name)
        print(parsed_entry.strip())
        with open(out_file, 'a') as of:
            of.write(parsed_entry)

# blastn -task megablast -query benchmark_mirgenedb.fa -db all_results -outfmt "6 std qcovhsp" -perc_identity 90 -qcov_hsp_perc 80 -out blastresults.txt