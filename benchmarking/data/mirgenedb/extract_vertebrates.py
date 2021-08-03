

all_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb\all_mirgenedb.fa'
abb_path = r'C:\Users\felix\Desktop\mirgene_db.CSV'
out_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb\named_all_mirgenedb.fa'

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


with open(all_path, 'r') as fh, open(out_path, 'w') as of:
    data = fh.read().split('>')
    for entry in data:
        db_abb = entry[0:3]
        parsed_entry = '>' + entry.replace(db_abb, map_dict[db_abb])
        print(parsed_entry.strip())
        of.write(parsed_entry)
