import os
import fetch_pdb
import collections
import csv
import pprint
import tqdm

ROOT_PATH = "./v3.17"

# 塩基と酸性のアミノ酸を判別
ACIDIC_AMINOS = [ 'K', 'R', 'H' ]
BASIC_AMINOS = [ 'D', 'E' ]

def check_amino_acidic_basic(residue):
  if (residue in ACIDIC_AMINOS):
    # 酸性
    return "A"
  elif (residue in BASIC_AMINOS):
    # 塩基性
    return "B"
  else:
    # 中性
    return "N"
  
# FASTAファイルのそれぞれのMSAされた配列情報に二次構造情報を付加する
def attend_secondary_info(secondary_info, residues_str):
  secondary_str = ''
  residue_index = 0
  new_asa_list = []
  acidic_basic_str = ""
  for residue in list(residues_str):
    if(residue == '-'): 
      secondary_str += '-'
      new_asa_list.append("-")
      acidic_basic_str += '-'
    else:
      # 酸性か塩基性か調べる
      acidic_basic_str += check_amino_acidic_basic(residue)
      # 二次構造情報を調べる
      if(residue_index in secondary_info['helix']):
        secondary_str += 'H'
      elif(residue_index in secondary_info['sheet']):
        secondary_str += 'S'
      else:
        secondary_str += '#'
      # residue_indexがASAの配列数より多ければcontinue
      if( residue_index <  len(secondary_info['asa'])):
        new_asa_list.append(secondary_info['asa'][residue_index])
      residue_index += 1

  # 2次構造情報を表示する
  print(residues_str)
  print(secondary_str)
  print(acidic_basic_str)
  print("".join(new_asa_list))

  # EとBの割合を調べる(埋もれているものか埋もれていないものか)
  new_asa_res = {
    'e_count': int((new_asa_list.count("E"))*10/len(new_asa_list))/10,
    'b_count': int((new_asa_list.count("B"))*10/len(new_asa_list))/10,
    'n_count': int((new_asa_list.count("N"))*10/len(new_asa_list))/10,
    'total_length': len(new_asa_list)
  }

  return {
    'residues_str': residues_str,
    'secondary_str': secondary_str,
    'asa': new_asa_list,
    'acidic_basic_str': acidic_basic_str,
    'new_asa_res': new_asa_res,
  }

# それぞれのFASTAのファイルを解析して、pdbに関するものを抽出し、PDBのAPIを叩いて二次構造がどの残基なのかを判別
def analysis_each_file(file_name):
  msa_result = []
  f = open(file_name, 'r')
  lines = f.readlines()
  residues_str = ''
  for line in lines:
    line = line.replace("\n", "")
    if('>' in line):
      if('pdb' in line):
        ids = line.split('|')
        pdb_num = ids.index('pdb')
        # fetch PDB and get secondary structure info
        res = fetch_pdb.get_secondary_structure_residues(ids[pdb_num + 1], ids[pdb_num + 2])
        infos = attend_secondary_info( res, residues_str )
        msa_result.append( infos )

      residues_str = ''
    else:
      residues_str += line
  f.close()
  return msa_result

def output_result(msa_result):
  print('msa["new_asa_res"]')
  for msa in msa_result:
    print(msa["new_asa_res"])

  if( msa_result == []):
    return []
  
  protein_analzed_without_asa = []
  protein_analzed = []
  for res in msa_result:
    protein = []
    protein_without_asa = []
    
    # ASAと２次構造情報の2種類を組み合わせた場合(ASAの情報が配列分ない場合は不足分を"N"で保管する)
    if ( len(res['secondary_str']) > len( res['asa'])):
      null_count = len(res['secondary_str']) - len(res['asa'])
      for x in range(null_count):
          res['asa'].append("N")
    for secondary, asa in zip(list(res['secondary_str']), list( res['asa']) ):
      protein.append(secondary + asa)
    protein_analzed.append(protein)

    # ASAの埋もれ度を加味しない場合
    protein_analzed_without_asa.append(list(res['secondary_str']))

  num_of_sequence = len( protein_analzed )
  sequence_length = len( protein_analzed[0] )

  # 各カラムの出現率をカウントしている
  # ASAと２次構造情報の2種類を活用
  appearance_rates = []
  for index in range(sequence_length):
    analzed_this_column = [each_protein[index] for each_protein in protein_analzed]
    c = collections.Counter(analzed_this_column)
    # ギャップの部分は0にする
    if(c.most_common()[0][0] == "--"):
      appearance_rates.append( 0 )
    else: 
      appearance_rates.append( int( (c.most_common()[0][1]/num_of_sequence) * 100) )

  # ASA(埋もれ度なしの場合)
  appearance_rates_without_asa = []
  for index in range(sequence_length):
    analzed_this_column = [ each_protein[index] for each_protein in protein_analzed_without_asa ]
    c = collections.Counter(analzed_this_column)
    # ギャップの部分は0にする
    if(c.most_common()[0][0] == "-"):
      appearance_rates_without_asa.append( 0 )
    else: 
      appearance_rates_without_asa.append( int( (c.most_common()[0][1]/num_of_sequence) * 100) )

  # ASA(埋もれ度なしの場合)
  appearance_rates_without_asa_only_secondary = []
  for index in range(sequence_length):
    analzed_this_column = [ each_protein[index] for each_protein in protein_analzed_without_asa ]
    c = collections.Counter(analzed_this_column)
    # ギャップの部分は0にする
    if((c.most_common()[0][0] == "S") or (c.most_common()[0][0] == "H")):
      appearance_rates_without_asa_only_secondary.append( int( (c.most_common()[0][1]/num_of_sequence) * 100) )
    else:
      appearance_rates_without_asa_only_secondary.append( 0 )

  print(appearance_rates)
  print(appearance_rates_without_asa)
  print(appearance_rates_without_asa_only_secondary)
  return appearance_rates

  return {
    'appearance_rates': appearance_rates,
    'appearance_rates_without_asa': appearance_rates_without_asa,
    'appearance_rates_without_asa_only_secondary': appearance_rates_without_asa_only_secondary
  }

# "./v3.17"以下のフォルダを全て読み込み、該当するFASTAファイルのパスを抽出
def search_all_file():
  # # ファイルの中身消去
  f = open(".2022_07_19/appearance_rates.csv","w")
  f.close()
  f = open(".2022_07_19/appearance_rates_without_asa.csv","w")
  f.close()
  f = open(".2022_07_19/appearance_rates_without_asa_only_secondary.csv","w")
  f.close()

  files = os.listdir(ROOT_PATH)
  # print(files)
  for this_file in tqdm.tqdm(files):
    NCBI_ID_PATH = ROOT_PATH + "/" + this_file + "/"
    index_file = NCBI_ID_PATH + this_file + ".fam"
    file_paths = []
    f = open(index_file, 'r')
    lines = f.readlines()
    for line in lines:
      ncbi_id = line.replace("\n", "")
      file_paths.append(NCBI_ID_PATH + ncbi_id + ".FASTA")
    f.close()

    for path in file_paths:
      print(path)
      msa_result = analysis_each_file(path)
      if( len(msa_result) > 1):
        res = output_result( msa_result ) # 出現度合いを表した度合い
        appearance_rates = res['appearance_rates']
        appearance_rates_without_asa = res['appearance_rates_without_asa']
        appearance_rates_without_asa_only_secondary = res['appearance_rates_without_asa_only_secondary']

        # 配列の先頭にNCBIのIDを入れる
        ncbi_fam_id = path.split('/')[2]
        ncbi_id = path.split('/')[3].replace(".FASTA", "")
        appearance_rates.insert(0, ncbi_id)
        appearance_rates.insert(0, ncbi_fam_id)
        appearance_rates_without_asa.insert(0, ncbi_id)
        appearance_rates_without_asa.insert(0, ncbi_fam_id)
        appearance_rates_without_asa_only_secondary.insert(0, ncbi_id)
        appearance_rates_without_asa_only_secondary.insert(0, ncbi_fam_id)

        with open('./2022_07_19/appearance_rates.csv', 'a') as f:
            writer = csv.writer( f )
            writer.writerow( appearance_rates )
        with open('./2022_07_19/appearance_rates_without_asa.csv', 'a') as f:
            writer = csv.writer( f )
            writer.writerow( appearance_rates_without_asa )
        with open('./2022_07_19/appearance_rates_without_asa_only_secondary.csv', 'a') as f:
            writer = csv.writer( f )
            writer.writerow( appearance_rates_without_asa_only_secondary )

# 以下のコメントアウトを外す
search_all_file()

# test_path = './v3.17/cd12120/cd12199.FASTA'
# msa_result = analysis_each_file(test_path)
# if(len(msa_result)== 0):
#   print("あああ")
# res = output_result( msa_result )

