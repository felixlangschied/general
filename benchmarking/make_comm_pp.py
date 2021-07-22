

five_pp = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\05phyloprofile.long'
six_pp = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\06phyloprofile.long'

pp_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\mirgenedb.long'

map_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\available_specs.txt'
pres_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\present.txt'
# tops = r'C:\Users\felix\PycharmProjects\general\benchmarking\top_30_diff.txt'
# with open(tops, 'r') as fh:
#     mirnas = fh.read().strip().split(',')

with open(pres_path, 'r') as fh:
    present = fh.read().split('\n')
    present.pop()


taxid_2_name = {}
name_2_taxid = {}
# keep_ids = []
with open(map_path, 'r') as fh:
    for line in fh:
        data = line.strip().split('|')
        taxid = data[0]
        name = data[1].replace(' ', '_')
        taxid_2_name[taxid] = name
        name_2_taxid[name] = taxid
        # keep_ids.append(taxid)

keep_ids = [name_2_taxid[name] for name in present]
print(keep_ids)


out_list = []
check_list = []
int_mirnas = []
int_fams = []
with open(six_pp, 'r') as fh:
    old_header = next(fh).strip()
    header = f'{old_header}\twith06\n'
    out_list.append(header)
    for line in fh:
        check_list.append(line)
        linedata = line.strip().split('\t')
        if linedata[1].replace('ncbi', '') not in keep_ids:
            continue
        fam = '-'.join(linedata[0].split('-')[1:3])
        # if fam not in mirnas:
        #     continue
        linedata.append('1')
        linedata[0] = fam
        linestr = '\t'.join(linedata) + '\n'
        out_list.append(linestr)
with open(five_pp, 'r') as fh:
    old_header = next(fh).strip()
    for line in fh:
        if line in check_list:
            continue
        linedata = line.strip().split('\t')
        if linedata[1].replace('ncbi', '') not in keep_ids:
            continue
        int_mirnas.append(linedata[0])
        fam = '-'.join(linedata[0].split('-')[1:3])
        # if fam not in mirnas:
        #     continue
        linedata.append('0')
        linedata[0] = fam
        int_fams.append(fam)
        linestr = '\t'.join(linedata) + '\n'
        out_list.append(linestr)

int_mirnas = set(int_mirnas)
int_fams = set(int_fams)


with open('phylo_input/family_results_19specs.long', 'w') as of:
    # of.write(out_list[0])
    for line in out_list:
        # if line.strip().split('\t')[0] in int_mirnas:
        # if line.strip().split('\t')[0] in int_fams:
            of.write(line)


out_list = []
check_list = []
in_human = []
with open(pp_path, 'r') as fh:
    header = next(fh)
    out_list.append(header)
    for line in fh:
        if line in check_list:
            continue
        check_list.append(line)
        linedata = line.strip().split('\t')
        if linedata[1] == 'ncbi9606':
            in_human.append(linedata[0])
        if linedata[1].replace('ncbi', '') not in keep_ids:
            continue
        linestr = '\t'.join(linedata) + '\n'
        out_list.append(linestr)


with open('phylo_input/family_mirgenedb_19specs.long', 'w') as of:
    of.write(out_list[0])
    for line in out_list:
        if line.strip().split('\t')[0] in in_human:
            of.write(line)