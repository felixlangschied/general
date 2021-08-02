import pandas as pd

tp_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\tp_results.csv'

res_df = pd.read_csv(tp_path)
# print(res_df)

has_tps = set(res_df['species'][res_df['True_positive'] == 'True'].values)
print(has_tps)

worked_df = res_df[res_df['species'].isin(has_tps)]
print(worked_df)

all_count = worked_df['True_positive'].size
true_count = worked_df['True_positive'][worked_df['True_positive'] == 'True'].size

print(f'{true_count} / {all_count} = {true_count / all_count}')