
# 親 - 子 - 孫
# v3.17/ 以下にあるフォルダの# "./v3.17"以下のフォルダを全て読み込み、該当するFASTAファイルのパスを抽出
import os
import tqdm
import numpy as np
import matplotlib.pyplot as plt

ROOT_PATH = "./v3.17"

no_pdb_file = []
def analysis_each_file(file_name):
  f = open(file_name, 'r')
  lines = f.readlines()
  protein_counts = 0
  protein_in_pdb_counts = 0
  for line in lines:
    line = line.replace("\n", "")
    if('>' in line):
      protein_counts += 1
      if('pdb' in line):
        protein_in_pdb_counts += 1
  f.close()
  if(protein_in_pdb_counts == 0):
    no_pdb_file.append(file_name)
  return protein_counts, protein_in_pdb_counts

def search_all_file():
  files = os.listdir(ROOT_PATH)
  
  family_counts = 0
  childrens = []
  msa_protein_count = []
  msa_protein_in_pdb_count = []

  msa_protein_count_large = []
  msa_protein_count_normal = []
  msa_protein_count_small = []

  msa_protein_in_pdb_count_large = []
  msa_protein_in_pdb_count_normal = []
  msa_protein_in_pdb_count_small = []

  for this_file in tqdm.tqdm(files):
    if ('cd' in this_file):
      family_counts += 1
      NCBI_ID_PATH = ROOT_PATH + "/" + this_file + "/"
      index_file = NCBI_ID_PATH + this_file + ".fam"
      f = open(index_file, 'r')
      lines = f.readlines()
      # 親フォルダ内部のファイルの検査

      childrens.append(len(lines))

      families = []
      for line in lines:
        ncbi_id = line.replace("\n", "")
        families.append(NCBI_ID_PATH + ncbi_id + ".FASTA")
      f.close()

      for family in families:
        protein_counts, protein_in_pdb_counts = analysis_each_file(family)
        msa_protein_count.append(protein_counts)
        msa_protein_in_pdb_count.append(protein_in_pdb_counts)

        if(protein_counts >= 250):
          msa_protein_count_large.append(protein_counts)
        elif((protein_counts < 250) and (100 < protein_counts)):
          msa_protein_count_normal.append(protein_counts)
        else:
          msa_protein_count_small.append(protein_counts)

        if(protein_in_pdb_counts >= 50):
          msa_protein_in_pdb_count_large.append(protein_in_pdb_counts)
        elif((protein_in_pdb_counts < 50) and (10 < protein_in_pdb_counts)):
          msa_protein_in_pdb_count_normal.append(protein_in_pdb_counts)
        else:
          msa_protein_in_pdb_count_small.append(protein_in_pdb_counts)

  # #figure()でグラフを表示する領域をつくり，figというオブジェクトにする．
  # fig = plt.figure()

  # #add_subplot()でグラフを描画する領域を追加する．引数は行，列，場所
  # ax1 = fig.add_subplot(2, 3, 1) #msa large
  # ax2 = fig.add_subplot(2, 3, 2) #msa normal
  # ax3 = fig.add_subplot(2, 3, 3) #msa small
  # ax4 = fig.add_subplot(2, 3, 4) #pdb large
  # ax5 = fig.add_subplot(2, 3, 5) #pdb normal
  # ax6 = fig.add_subplot(2, 3, 6) #pdb small

  # # plt.ylim(0, max(msa_protein_count)+1)
  # index = [ n for n in range(len(msa_protein_count)) ]
  # # msa_protein_count = np.array(msa_protein_count)
  # msa_protein_in_pdb_count = np.array(msa_protein_in_pdb_count)
  # # plt.bar(index, msa_protein_count, color="green")

  # ax1.bar([ n for n in range(len(msa_protein_count_large)) ], msa_protein_count_large, color="orange")
  # ax2.bar([ n for n in range(len(msa_protein_count_normal)) ], msa_protein_count_normal, color="orange")
  # ax3.bar([ n for n in range(len(msa_protein_count_small)) ], msa_protein_count_small, color="orange")

  # ax4.bar([ n for n in range(len(msa_protein_in_pdb_count_large)) ], msa_protein_in_pdb_count_large, color="green")
  # ax5.bar([ n for n in range(len(msa_protein_in_pdb_count_normal)) ], msa_protein_in_pdb_count_normal, color="green")
  # ax6.bar([ n for n in range(len(msa_protein_in_pdb_count_small)) ], msa_protein_in_pdb_count_small, color="green")
  # fig.tight_layout() 
  # plt.show()

  print("ファミリー数: " + str(family_counts))
  print("子ノード数: " + str(sum(childrens)))
  print("全てのMSAに使われたタンパク質の数: " + str(sum(msa_protein_count)))
  print("全てのMSAに使われたタンパク質のうちpdbに該当する数: " + str(sum(msa_protein_in_pdb_count)))

search_all_file()
