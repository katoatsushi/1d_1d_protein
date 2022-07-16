# csvモジュールを使ってCSVファイルから1行ずつ読み込む
import csv
import numpy as np
import pandas as pd
import seaborn as sns


full_match_score = []
protein_count = []

filename = 'output_E_0.csv'
with open(filename, encoding='utf8', newline='') as f:
  res = []
  csvreader = csv.reader(f)
  for row in csvreader:
      path_name = "./v3.17/" + row[0] + "/" + row[1] + ".FASTA"

      f = open(path_name, 'r')
      lines = f.readlines()
      count = 0
      for line in lines:
        line = line.replace("\n", "")
        if('pdb' in line):
          count += 1
      f.close()

      percentages = row[2:]
      # 全部一致の配列がいくつあるか
      all_in_percentage = percentages.count('100') / len(percentages)
      # res.append({
      #   "family": row[0],
      #   "subfamily": row[1],
      #   "hundred_ratio": all_in_percentage,
      #   "count": count,
      # })
      full_match_score.append(all_in_percentage)
      protein_count.append(count)


df_full_match = pd.DataFrame(full_match_score, columns=["full_match_score"])
df_protein = pd.DataFrame(protein_count, columns=["protein_count"])

df = pd.concat([df_full_match, df_protein], axis=1)


corr = df.corr()
print(corr)
sns.heatmap(df.corr())
