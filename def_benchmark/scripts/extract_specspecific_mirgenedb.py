import os


all_path = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/ALL-pre.fas'
abb_path = '/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/mirgenedb/supp_tab.tsv'
out_path = '/home/felixl/PycharmProjects/general/def_benchmark/data/spec_specific/{}.fa'


map_dict = {}
with open(abb_path, 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')[1:3]
        if data:
            name, abb = data
            name = name.replace(' ', '_')
            abb = abb.capitalize()
            map_dict[abb] = name


with open(all_path, 'r') as fh:
    data = fh.read().split('>')
    data.pop(0)
    for entry in data:
        print(entry)
        db_abb = entry[0:3]
        name = map_dict[db_abb]
        name = name.replace('Sarcophilus_harrissii', 'Sarcophilus_harrisii')
        name = name.replace('Canis_familiaris', 'Canis_lupus_familiaris')
        name = name.replace('Strongylocentrotus_purpuratus_', 'Strongylocentrotus_purpuratus')
        out_file = out_path.format(name)
        parsed_entry = '>' + entry.replace(db_abb, name)
        print(parsed_entry.strip())
        with open(out_file, 'a') as of:
            of.write(parsed_entry)
