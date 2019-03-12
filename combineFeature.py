# -*- coding:utf-8 -*-
from utils import combine_feature_v1,svm2csv
import os

input_pos_files = []
input_neg_files = []
pos_file_dir = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_pos_top10/'
neg_file_dir = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_top10/'
for root, dir, files in os.walk(pos_file_dir):
    for file_name in files:
        input_pos_files.append(root+file_name)
for root, dir, files in os.walk(neg_file_dir):
    for file_name in files:
        input_neg_files.append(root+file_name)
# print(input_pos_files)
output_pos_file = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_pos_top10_svm.txt'
output_neg_file = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_top10_svm.txt'

combine_feature_v1(input_pos_files,output_pos_file)
combine_feature_v1(input_neg_files,output_neg_file)

output_pos_file_csv = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_pos_top10_csv.txt'
output_neg_file_csv = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_top10_csv.txt'

svm2csv(output_pos_file,output_pos_file_csv)
svm2csv(output_neg_file,output_neg_file_csv)