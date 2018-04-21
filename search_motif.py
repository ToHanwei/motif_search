#/usr/bin/python
#coding:utf-8

"""
              **********************
              *                    *
              *Sequence Conservasim*
              *                    *
              **********************
This code to find the motif of the proteins family
"""

__author__ = "Wei"
__time__ = "2018-04-18"
__email__ = "hanwei@shanghaitech.edu.cn"


import argparse
import numpy as np
from sys import argv
from pandas import DataFrame
from itertools import groupby
from collections import defaultdict

def analyze_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--inputfile', action='store',
                    help='input the sequence aligment file')
	parser.add_argument('-c', '--conservatism', action='store', default=60,
                    help='input the comservatism as you want,default 60')
	parser.add_argument('-l', '--long', action='store', default=3,
                    help='motif long you need, default 3')
	args = parser.parse_args()
	input_file = args.inputfile
	conservatism = args.conservatism
	motif_long = args.long
	return input_file, conservatism, motif_long

def open_file(fname):
	fopen = open(fname, 'r')
	fread = fopen.read().strip()
	blocks = fread.split('\n\n')[1:]
	return blocks
  
def transformat(blocks):
	blocks_dict = defaultdict(list)
	for block in blocks:
		block = block.strip().split('\n')
		for line in block:
			line = line.split(' ')
			while '' in line:
				line.remove('')
			if line[-1] in ['.', ':', '*']:
				continue
			blocks_dict[line[0]] += list(line[1])
	blocks_df = DataFrame(blocks_dict)
	return blocks_df

def getresult(blocks_df, number):
	leng = len(blocks_df)
	result_dict = defaultdict(list)
	row_dict = defaultdict(list)
	for index, row in blocks_df.iterrows():
		aad_dict = defaultdict(int)
		for aad in row:
			aad_dict[aad] += 1
		result = sorted(aad_dict.items(), key=lambda x:x[1])[-1]
		result = list(result)
		result[1] = 100*round(result[1]/leng, 4)
		row_dict[index] = result
		if result[0] != '-' and result[1] >= number:
			result_dict[index] = result
	return row_dict, result_dict

def search_motif(df, continu_num):
	a_index = np.delete(df.index, 0, axis=0)
	b_index = np.delete(df.index, -1, axis=0)
	c_index = a_index - b_index
	group_list = [(key, len(list(num))) for key, num in groupby(c_index)]
	count = 0; motif_list = []
	for tmp in group_list:
		if tmp[0] == 1 and tmp[1]+1>= continu_num:
			start = df.index[count]
			motif_list.append((start, tmp[1]+1))
		count += tmp[1]
	return motif_list		
	
def construc_motif_matrix(blocks_df, motif_list):
	index_list = list(blocks_df.index)
	for tmp in motif_list:
		start = tmp[0]
		continu = tmp[1]
		for i in range(start, start+continu):
			index_list[i] = '# '+str(i)+' #'
	blocks_df.index = index_list
	motif_matrix = blocks_df.T
	return motif_matrix

def print_motif(result, motif_list):
	count = 1
	for tmp in motif_list:
		str_out = ''
		start = tmp[0]
		continu = tmp[1]
		print("Motif %2d is: " % (count), end=" "*2)
		for i in range(start, start+continu):
			print(result.loc[i][0], end="")
			str_out += (str(i) + '#')
		print()
		print("Motif numbering on aligment is %s" % str_out)
		count += 1
	print("There are %d motif" % (count-1))

if __name__ == "__main__":
	args = analyze_parser()
	fname = args[0]
	conservatism = int(args[1])
	continu_num = int(args[2])
	blocks = open_file(fname)
	blocks_df = transformat(blocks)
	row, result = getresult(blocks_df, conservatism)
	row = DataFrame(row).T
	result = DataFrame(result).T
	if len(result):
		motif_list = search_motif(result, continu_num)
		print_motif(result, motif_list)
		motif_matrix = construc_motif_matrix(blocks_df, motif_list)
		row.to_csv("sequence_all_conservatism.csv", header=True, index=True)
		result.to_csv("sequence_conservatism.csv", header=True, index=True)
		motif_matrix.to_csv("motif_matrix.csv", header=True, index=True)
	else:
		print("Use your parser there are no motif")
		print("You can try decrease the --conservatism or --long parsers")
	
	
	
	