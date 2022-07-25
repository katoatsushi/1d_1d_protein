# # python analysis_csv.py
import csv
# データ前処理
import numpy as np
import pandas as pd

# データ可視化
import seaborn
import matplotlib.pyplot as plt


ROOT_PATH = './2022_07_20_SAVED/'

def plot_datas(path, title):
  path = ROOT_PATH + path

  csv_file = open(path, "r", encoding="ms932", errors="", newline="" )
  #リスト形式
  f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
  plot_data = []
  rate_avarage_plot_data = []
  
  # 相関関係のためのリスト
  msa_counts_list = []
  hundred_counts_list = []
  rate_avarage_list = []

  for row in f:
    # 3番目以降のデータを読み取る
    rate_datas = row[3:]
    msa_counts = int(row[0])
    row_length = len(rate_datas)
    hundred_counts = rate_datas.count('100')
    plot_data.append([ msa_counts, int(hundred_counts*100 / row_length)/100])

    total_rate = 0
    for data in rate_datas:
      total_rate += int(data)
    rate_avarage_plot_data.append([ msa_counts, int(total_rate*100 / row_length)/100])
    
    # 相関関係のためのもの
    msa_counts_list.append(msa_counts)
    hundred_counts_list.append(hundred_counts)
    rate_avarage_list.append( int(total_rate*100 / row_length)/100 )

  df_msa_counts_list = pd.DataFrame(msa_counts_list, columns=["msa_counts_list"])
  df_hundred_counts_list = pd.DataFrame(hundred_counts_list, columns=["hundred_counts_list"])
  df_rate_avarage_list = pd.DataFrame(rate_avarage_list, columns=["rate_avarage_list"])

  df_handred_rate = pd.concat([df_msa_counts_list, df_hundred_counts_list], axis=1)
  corr_handred_rate = df_handred_rate.corr()
  print(title)
  print("correration with handred_rate")
  print(str(corr_handred_rate))
  print("="*100)
  df_rate_avarage = pd.concat([df_msa_counts_list, df_rate_avarage_list], axis=1)
  corr_avarage_rate = df_rate_avarage.corr()
  print("correration with avarage rate")
  print(str(corr_avarage_rate))

  return plot_data, rate_avarage_plot_data


# MSAされた本数が少ないと、100%の確率が高いか図る
def match_rate():
  #df_simple = pd.DataFrame(np.array([[2, 0.3],[6, 0.2],[10,0.2], [2, 0.8], [2, 0.9]]))
  #data = [[2, 0.6],[6, 0.2],[10,0.2], [2, 0.8], [2, 0.9],[13,0.2], [10,0.5]]

  datas = [
    {'url':'appearance_rates_without_asa_only_secondary.csv','title': "no ASA secondary structure", 'index': 1 },
    {'url':'appearance_rates_without_asa.csv', 'title': "no ASA", 'index': 2 },
    {'url':'appearance_rates_without_asa_only_secondary.csv', 'title': "ASA", 'index': 3 }
  ]

  #figure()でグラフを表示する領域をつくり，figというオブジェクトにする．
  fig = plt.figure()

  #add_subplot()でグラフを描画する領域を追加する．引数は行，列，場所
  ax1 = fig.add_subplot(3, 2, 1)
  ax2 = fig.add_subplot(3, 2, 2)
  ax3 = fig.add_subplot(3, 2, 3)
  ax4 = fig.add_subplot(3, 2, 4)
  ax5 = fig.add_subplot(3, 2, 5)
  ax6 = fig.add_subplot(3, 2, 6)

  axs = [ax1,ax2,ax3,ax4,ax5,ax6]

  index = 0
  for data in datas:
    plot_data, rate_avarage_plot_data = plot_datas(data['url'], data['title'])
    rate_avarage_plot_data = [ [ n[0], n[1]/100 ] for n in rate_avarage_plot_data]

    x, y = zip(*plot_data)
    x_ava, y_ava = zip(*rate_avarage_plot_data)

    axs[index].scatter(x, y)
    axs[index+1].scatter(x_ava, y_ava)

    axs[index].set_title(data['title'])
    axs[index+1].set_title(data['title'])
    
    index += 2

  ax1.legend(loc = 'upper right') #凡例
  ax2.legend(loc = 'upper right') #凡例
  ax3.legend(loc = 'upper right') #凡例
  ax4.legend(loc = 'upper right') #凡例
  ax5.legend(loc = 'upper right') #凡例
  ax6.legend(loc = 'upper right') #凡例
  fig.tight_layout()              #レイアウトの設定
  plt.show()

match_rate()
