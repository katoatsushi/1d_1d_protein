import read_msa_file
import os
import fetch_pdb
import collections
import csv
import pprint
import tqdm

ROOT_PATH = "./v3.17"

# FASTAファイルのそれぞれのMSAされた配列情報に二次構造情報を付加する
def attend_secondary_info(secondary_info, residues_str):
  secondary_str = ''
  residue_index = 0
  new_asa_list = []

  for residue in list(residues_str):
    if(residue == '-'): 
      secondary_str += '-'
      new_asa_list.append("-")
    else:
      if(residue_index in secondary_info['helix']):
        secondary_str += 'H'
      elif(residue_index in secondary_info['sheet']):
        secondary_str += 'S'
      else:
        secondary_str += '#'
      # residue_indexがASAの配列数よりもおお消えればcontinue
      if( residue_index <  len(secondary_info['asa'])):
        new_asa_list.append(secondary_info['asa'][residue_index])
      residue_index += 1

  # 2次構造情報を表示する
  print(residues_str)
  print(secondary_str)
  print("".join(new_asa_list))
  print("\n")

  return {
    'residues_str': residues_str,
    'secondary_str': secondary_str,
    'asa': new_asa_list
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
  if( msa_result == []):
    return []

  protein_analzed = []
  for res in msa_result:
    protein = []
    # for residues, secondary in zip(list(res['residues_str']), list(res['secondary_str']) ):
    # ASAの情報が配列分ない場合は不足分を"N"で保管する
    if ( len(res['secondary_str']) > len( res['asa'])):
      null_count = len(res['secondary_str']) - len(res['asa'])
      for x in range(null_count):
          res['asa'].append("N")
    for residues, secondary in zip(list(res['secondary_str']), list( res['asa']) ):
      protein.append(residues + secondary)
    protein_analzed.append(protein)
    # print(protein_analzed)
  num_of_sequence = len( protein_analzed )
  sequence_length = len( protein_analzed[0] )

  appearance_rates = []
  # 各カラムの出現率をカウントしている
  for index in range(sequence_length):
    analzed_this_column = [each_protein[index] for each_protein in protein_analzed]
    c = collections.Counter(analzed_this_column)
    # ギャップの部分は0にする
    if(c.most_common()[0][0] == "--"):
      appearance_rates.append( 0 )
    else: 
      appearance_rates.append( int( (c.most_common()[0][1]/num_of_sequence) * 100) )


  print(appearance_rates)
  return appearance_rates

# "./v3.17"以下のフォルダを全て読み込み、該当するFASTAファイルのパスを抽出
def search_all_file():
  # ファイルの中身消去
  f = open("output_without_gap.csv","w")
  f.close()

  files = os.listdir(ROOT_PATH)
  print(files)
  # files = ['cd00004']
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
        # 配列の先頭にNCBIのIDを入れる
        ncbi_fam_id = path.split('/')[2]
        ncbi_id = path.split('/')[3].replace(".FASTA", "")
        res.insert(0, ncbi_id)
        res.insert(0, ncbi_fam_id)

        with open('output_without_gap.csv', 'a') as f:
            writer = csv.writer( f )
            writer.writerow( res )

search_all_file()
# test_path = './v3.17/cd12120/cd12199.FASTA'
# msa_result = analysis_each_file(test_path)
# if(len(msa_result)== 0):
#   print("あああ")
# res = output_result( msa_result )

