import json
import matplotlib.pyplot as plt
import seaborn as sns

json_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\all05_RNAfold_results.txt'

all_scores = []
with open(json_path, 'r') as fh:
    for line in fh:
        if line.strip().endswith(')'):
            score = float(line.strip().split(' (')[1].replace(')', ''))
            all_scores.append(score)

print(sorted(all_scores, reverse=True))

sns.histplot(all_scores)
plt.show() m