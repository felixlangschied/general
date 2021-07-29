import os
import subprocess as sp
import os
import re
from ftplib import FTP
import wget

out_path = '/home/felixl/PycharmProjects/general/benchmarking/data/chrom_map'
avail_path = '/home/felixl/PycharmProjects/general/benchmarking/data/available_specs.txt'
present_path = '/home/felixl/PycharmProjects/general/benchmarking/data/present.txt'

def extract_GCF(ids, output, flavor, naming, dir_mode):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd("/genomes/all/")
    count = -1

    for id in ids:
        count += 1
        ftp.cwd('/genomes/all/')
        id = id.split('.')[0]
        typ = id.split('_')[0]
        id = id.split('_')[1]
        part1 = id[0:3]
        part2 = id[3:6]
        part3 = id[6:9]

        url = '/'.join([typ, part1, part2, part3])
        ftp.cwd(url)
        assemblies = ftp.nlst()
        # will download the top file
        target = assemblies[0]
        if len(assemblies) > 1:
            print('# Found multiple assemblies for {}.\n# Downloading most recent: {}'
                  .format(id, assemblies[0])
                  )
        ftp.cwd(target)
        files = ftp.nlst()

        # find file types that can be downloaded (flavors)
        # ftypes = [file.split('.')[-3:] for file in files if '.gz' in file]
        # ftypes = ['.'.join(parts) for parts in ftypes]
        # ftypes = [ftype.split('_')[2:] for ftype in ftypes]
        # ftypes = ['_'.join(parts) for parts in ftypes]
        # print(ftypes)
        # print('\n'.join(ftypes))

        r = re.compile('.*' + flavor)
        target_file = list(filter(r.match, files))[0]
        download_url = (
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}'
                .format(url, target, target_file)
        )
        os.chdir(output)
        try:
            print('\n# Starting download of {} for: {}'.format(flavor, target))
            wget.download(download_url)
            if naming:
                name = naming[count]
                print(name)
                if dir_mode:
                    if not os.path.isdir('{}/{}'.format(output, name)):
                        cmd = 'mkdir {}'.format(name)
                        sp.run(cmd, shell=True)
                    file_name = '{}_{}'.format(name, flavor)
                    cmd = 'mv {}_{} {}/{}'.format(target, flavor, name, file_name)
                    sp.run(cmd, shell=True)
                else:
                    file_name = '{}_{}'.format(name, flavor)
                    cmd = 'mv {}_{} {}'.format(target, flavor, file_name)
                    sp.run(cmd, shell=True)
            elif not naming and dir_mode:
                print('\n# No column with names detected in the input file. Unable to run in dir_mode')

        except:
            print('\n# {} not found'.format(target_file))


with open(present_path, 'r') as fh:
    present = fh.read().split('\n')
    present.pop()

gcf_list = []
with open(avail_path, 'r') as fh:
    for line in fh:
        data = line.strip().split('|')
        if data[1].replace(' ', '_') not in present:
            continue
        gcf_list.append(data[2])


extract_GCF(gcf_list, out_path, 'assembly_report.txt', False, False)