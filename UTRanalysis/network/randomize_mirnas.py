import random

mirna_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\network\mirnas_in_mirtarbase.txt'

with open(mirna_path, 'r') as fh:
    mirna_list = fh.read().split('\n')

print(mirna_list)
x = random.sample(range(len(mirna_list)), 5)
sample_mirnas = [mirna_list[ind] for ind in x]

for mirna in sample_mirnas:
    print(mirna)