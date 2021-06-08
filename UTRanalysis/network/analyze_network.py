import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

csv_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\network\target_network_table.csv'

df = pd.read_csv(csv_path)


df_f = df.filter(items=['name', 'Degree', "NeighborhoodConnectivity","BetweennessCentrality","ClosenessCentrality"])
df_f = df_f.sort_values(by='Degree', ascending=False)
print(df_f.head())
