import requests
import json
from statistics import mean

BASE_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/"
HELIX = 'helix'
SHEET = 'sheet'


# base_aminos = [ 'リシン','アルギニン','ヒスチジン' ]
# acied_amino = [ 'アスパラギン酸', 'グルタミン酸']

def get_secondary_structure_residues(pdb_id, chain_index):
  if (chain_index == " "):
    chain_index = "A"
  if (chain_index == "X"):
    chain_index = "A"
  print("PDB_ID: " + pdb_id + "  chain: " + chain_index)
  url = BASE_URL + pdb_id + "/" + chain_index
  res = requests.get(url)
  data = res.json()

  #　ChainIDが古い表記でAPIに乗っていない場合、Aに変換してみる
  if (( 'status' in data) and (data['status'] == 404)):
    print("ChainIDをAに変更します！")
    chain_index = "A"
    url = BASE_URL + pdb_id + "/" + chain_index
    res = requests.get(url)
    data = res.json()

  helix_residues = []
  sheet_residues = []
  asa_list = []

  # APIを叩いた結果を見たい場合は以下コメントアウトを外す
  # PDBでrcsb_polymer_instance_featureがないもの
  if ('rcsb_polymer_instance_feature' in data):
    for feature in data['rcsb_polymer_instance_feature']:
      if(('name' in feature) and (feature['name'] == HELIX)):
        for feature_positions in feature['feature_positions']:
          helix_residues += list(range(feature_positions['beg_seq_id'],feature_positions['end_seq_id'] + 1))
      elif(('name' in feature) and (feature['name'] == SHEET)):
        for feature_positions in feature['feature_positions']:
          sheet_residues += list(range(feature_positions['beg_seq_id'],feature_positions['end_seq_id'] + 1))
      elif(('type' in feature) and(feature['type'] == "ASA")):
        # 埋もれ度の値を抽出
        end_positions = []
        for p in feature['feature_positions']:
          end_positions.append(p['end_seq_id'])
        the_end_position = max(end_positions)

        # 埋もれ度の基準の平均値を獲得する
        # 数値が高い = 露出度が高い / 埋もれていない
        values  = []
        for pos in feature['feature_positions']:
            values += pos["values"]
        avarage = mean(values)
        # ここまで

        asa_list = [ "N" for i in range(the_end_position) ]

        for pos in feature['feature_positions']:
          index = pos["beg_seq_id"]
          for val in pos["values"]:
            if(index > pos["end_seq_id"]):
              break
            # 埋もれている: E　埋もれていない: Bにする
            # 1:sortして上、下で判別
            # 2:0なら埋もれている,　? > 0 なら埋もれていない
            if(int(val) == 0):
              asa_list[index - 1] = "E"
            else:
              asa_list[index - 1] = "B"
            # if(val > avarage):
            #   asa_list[index - 1] = "B"
            # else:
            #   asa_list[index - 1] = "E"
            
            # asa_list[index - 1] = int(val)
            index += 1
  return {
            "helix": helix_residues, 
            "sheet": sheet_residues,
            "asa": asa_list
          }


# あるPDBの情報をAPIで取得する
# get_secondary_structure_residues("2HXA", 'A')
# "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/2HXA/A"

get_secondary_structure_residues("1CX0", "A")
