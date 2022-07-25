import csv
import os

COLUMNS = { 
  'SB': 0, 'SE': 1, 'SN': 2, 
  "HB": 3, "HE": 4, "HN": 5,
  "#B": 6, "#E": 7, "#N": 8, 
  "S":  9, "H":  10,"#":  11,
  "--": 12, 'unknown': 13
}
LENGTH = len(COLUMNS.keys())
ROOT_PATH = './appearance_rates'


def save_appearance_rates_in_csv(family_cd_id, cd_id, count_results):
  path = ROOT_PATH + "/" + family_cd_id + "/" +  cd_id + ".csv"
  dir_path = ROOT_PATH + "/" + family_cd_id

  if not (os.path.exists(dir_path)):
    os.makedirs(dir_path)

  # ファイル初期化
  f = open(path,"w")
  f.close()

  for count_result in count_results:
    # 出現頻度をカウントしていく
    init_rates = [0 for n in range(LENGTH)]
    for key, appearance_times in dict(count_result).items():

      if(key in COLUMNS.keys()):
        init_rates[COLUMNS[key]] += int(appearance_times)
      else:
        init_rates[COLUMNS['unknown']] += int(appearance_times)

    # 結果を書き込む
    with open(path, 'a') as f:
        writer = csv.writer( f )
        writer.writerow( init_rates )
  

