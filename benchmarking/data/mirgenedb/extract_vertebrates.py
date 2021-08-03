

all_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb\all_mirgenedb.fa'
abb_path = r'C:\Users\felix\Desktop\mirgene_db.CSV'
out_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb\benchmark_mirgenedb.fa'
spec_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\spec_missing_TP.tsv'

map_dict = {}
with open(abb_path, 'r') as fh:
    next(fh)
    next(fh)
    for line in fh:
        data = line.strip().split(';')[3:5]
        if data:
            name, abb = data
            name = name.replace(' ', '_')
            abb = abb.capitalize()
            map_dict[abb] = name
# print(map_dict)

specs = []
with open(spec_path, 'r') as fh:
    next(fh)
    for line in fh:
        specs.append(line.split('\t')[0])


with open(all_path, 'r') as fh, open(out_path, 'w') as of:
    data = fh.read().split('>')
    for entry in data:
        db_abb = entry[0:3]
        name = map_dict[db_abb]
        if name not in specs:
            continue
        parsed_entry = '>' + entry.replace(db_abb, name)
        print(parsed_entry.strip())
        of.write(parsed_entry)
