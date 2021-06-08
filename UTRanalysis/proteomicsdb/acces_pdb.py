import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# mapping file I/O
# protname = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\mapping_data\uniprot_genename.txt'
protname = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\mapping_data\genename_uniprot.tab'
omapath = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\mapping_data\hsa_2_mmu_oma.txt'

# target proeints
# target_name_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'
target_name_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\hsa_rodentmis_genenames.txt'

def download_expression_info(prot_id):
    payload = (
        f"PROTEINFILTER='{prot_id}',"
        "MS_LEVEL=1,"
        "TISSUE_ID_SELECTION='',"
        # "TISSUE_CATEGORY_SELECTION='tissue;fluid',"
        "TISSUE_CATEGORY_SELECTION='cell line',"
        "SCOPE_SELECTION=1,"
        "GROUP_BY_TISSUE=0,"
        # 0 - iBAQ, 1 - TOP3 (Integer)
        "CALCULATION_METHOD=0,"
        "EXP_ID=-1"
    )

    res_select = (
        "Results?$select="
        "TISSUE_NAME,NORMALIZED_INTENSITY,SAMPLES"
        "&$format=json"
    )

    r_url = (
        "https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/InputParams({})/{}"
        .format(payload, res_select)
    )
    r = requests.get(r_url)
    return r


# map uniprot 2 genename
name_2_uniprot = {}
uniprot_2_name = {}
with open(protname, 'r') as fh:
    for line in fh:
        if line.startswith('yourlist'):
            continue
        name, uni = line.strip().split()
        name_2_uniprot[name] = uni
        uniprot_2_name[uni] = name
# map uniprot id between mouse and human
hsa_2_mmu = {}
mmu_2_hsa = {}
with open(omapath, 'r') as fh:
    for line in fh:
        data = line.strip().split()
        hsa = data[0]
        mmu = data[1]
        hsa_2_mmu[hsa] = mmu
        mmu_2_hsa[mmu] = hsa

target_protnames = []
with open(target_name_path, 'r') as fh:
    for line in fh:
        target_protnames.append(line.strip())
print(target_protnames)

# GET ALL AVAILABLE TISSUES TO INITIALIZ DATAFRAME
columns = ['gene_name', 'prot_id', 'organism', 'tissue', 'norm_intensity']

# prot_id = 'P00533'
df = pd.DataFrame(columns=columns)


# target_protnames = ['BMPR1A']

for tprot in target_protnames:
    try:
        hsa_id = name_2_uniprot[tprot]
    except KeyError:
        continue
    # FIRST EXTRACT FOR HUMAN
    print('# Extracting expression data for: {}'.format(tprot))
    r_human = download_expression_info(hsa_id)
    # prot_dict = {}
    for result in r_human.json()['d']['results']:
        tissue = result['TISSUE_NAME']
        intensity = float(result['NORMALIZED_INTENSITY'])
        df = df.append({'gene_name': tprot, 'organism': 'human', 'prot_id': hsa_id, 'tissue': tissue, 'norm_intensity': intensity}, ignore_index=True)


    # THEN EXTRACT FOR MOUSE
    try:
        mmu_id = hsa_2_mmu[hsa_id]
    except KeyError:
        continue
    r_mouse = download_expression_info(mmu_id)
    # prot_dict = {}
    for result in r_mouse.json()['d']['results']:
        tissue = result['TISSUE_NAME']
        intensity = float(result['NORMALIZED_INTENSITY'])
        df = df.append({'gene_name': tprot, 'organism': 'mouse', 'prot_id': mmu_id, 'tissue': tissue, 'norm_intensity': intensity},
                       ignore_index=True)

df.to_csv('cell_line_only_proteomics.csv')


# df = df.sort_values(by='norm_intensity')
# print(df)
# to_remove = ['IGF1R', 'BMPR1A', 'RORA', 'TWIST1']
# to_keep = [entry for entry in target_protnames if entry not in to_remove]
# print(to_keep)
#
# df = df.iloc[df.filter(to_keep, axis=1), :]
print(df)
# sns.barplot(x='gene_name', y='norm_intensity', hue='organism', data=df)
# sns.catplot(x='tissue', y='norm_intensity', hue='organism', data=df)
# plt.xticks(rotation=30, ha='right')
# plt.tight_layout()
sns.catplot(x='gene_name', y='norm_intensity', hue='organism', data=df, kind='violin')
plt.tight_layout()
plt.show()

# # MORE RESULT OPTIONS
# res_select = (
#     "/Results?$select="
#     "UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,SAMPLE_ID,SAMPLE_NAME,AFFINITY_PURIFICATION,EXPERIMENT_ID,EXPERIMENT_NAME,EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,PROJECT_NAME,PROJECT_STATUS,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,SAMPLES"
#     "&$format=json"
# )

# prot_dict = {}
# r = requests.get(r_url)
# for result in r.json()['d']['results']:
#     prot_dict[result['TISSUE_NAME']] = result['NORMALIZED_INTENSITY']
# if df.empty:
#     df = pd.DataFrame.from_dict(prot_dict, orient='index', columns=[])
# else:
#     tmp_df = pd.DataFrame.from_dict(prot_dict, orient='index', columns=[prot_id])
#     df = df.join(tmp_df)


#     prot_dict[result['TISSUE_NAME']] = [result['NORMALIZED_INTENSITY']]
# tmp_df = pd.DataFrame.from_dict(prot_dict)
# # tmp_df['organism'] = 'Mouse'
# tmp_df['protid'] = mmu_id
#
# # JOIN DATAFRAMES
# df = df.set_index('protid')
# tmp_df = tmp_df.set_index('protid')
# com_df = pd.concat([df, tmp_df], axis=0)
# print(com_df.head())