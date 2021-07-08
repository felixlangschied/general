

five_pp = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\05phyloprofile.long'
six_pp = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\06phyloprofile.long'

tops = r'C:\Users\felix\PycharmProjects\general\benchmarking\top_30_diff.txt'

with open(tops, 'r') as fh:
    mirnas = fh.read().strip().split(',')

out_list = []
check_list = []
int_mirnas = []
with open(six_pp, 'r') as fh:
    old_header = next(fh).strip()
    header = f'{old_header}\twith06\n'
    out_list.append(header)
    for line in fh:
        check_list.append(line)
        linedata = line.strip().split('\t')
        fam = '-'.join(linedata[0].split('-')[1:3])
        # if fam not in mirnas:
        #     continue
        linedata.append('1')
        linestr = '\t'.join(linedata) + '\n'
        out_list.append(linestr)
with open(five_pp, 'r') as fh:
    old_header = next(fh).strip()
    for line in fh:
        if line in check_list:
            continue
        linedata = line.strip().split('\t')
        int_mirnas.append(linedata[0])
        fam = '-'.join(linedata[0].split('-')[1:3])
        # if fam not in mirnas:
        #     continue
        linedata.append('0')
        linestr = '\t'.join(linedata) + '\n'
        out_list.append(linestr)

int_mirnas = set(int_mirnas)

with open('comm_pp.long', 'w') as of:
    of.write(out_list[0])
    for line in out_list:
        if line.strip().split('\t')[0] in int_mirnas:
            of.write(line)