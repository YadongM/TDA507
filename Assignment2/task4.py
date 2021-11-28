# Author: Yadong Mao

# To run this code, use command:
# python ./task4.py ./1HZH_H.pdb
# python ./task4.py ./4GAF_B.pdb

# Ditance threshold is 10
# Find all domains by setting min_residue and min_split_score

import re
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_file(pdb_path):
    lines_atom = []
    with open(pdb_path, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.split()
            if (line[0] == 'ATOM' and line[2] == 'CA'):
                lines_atom.append(line)

    lines_position = [line[6:9] for line in lines_atom]
    lines_position = np.array(lines_position).astype(float)

    lines_index = [re.sub('[A-Z]', '', line[5]) for line in lines_atom]
    lines_index = np.array(lines_index).astype(int)

    return lines_index, lines_position

def cal_distance(lines_index, lines_position, threshold=10):
    pair = []
    all_distance = np.zeros([len(lines_position), len(lines_position)])
    for i in range(len(lines_position)):
        for j in range(len(lines_position)):
            _distance = np.sqrt(np.sum( (lines_position[i]-lines_position[j])**2 ))
            all_distance[i][j] = _distance
            if (i!=j and _distance<=threshold):
                pair.append([lines_index[i], lines_index[j]])
    return np.array(pair)



def cal_domak(lines_index, pair):

    max_split_value = 0
    _counter = 0

    for current in lines_index:
        int_A = np.sum((pair[:,0]<current) & (pair[:,1]<current))
        int_B = np.sum((pair[:,0]>current) & (pair[:,1]>current))
        ext_AB = np.sum((pair[:,0]<=current) & (pair[:,1]>=current))

        if ext_AB != 0:
            _score = (int_A/ext_AB) * (int_B/ext_AB)

            if max_split_value < _score:
                max_split_value = _score
                split_index = current
                index = _counter
        _counter+=1
    
    return max_split_value, split_index, index


def get_all_split(lines_index, lines_position, split_indexs):

    min_residue = 20
    min_split_score = 80
    
    pair = cal_distance(lines_index, lines_position)
    max_split_value, split_index, index = cal_domak(lines_index, pair)

    lines_index1 = lines_index[index:]
    lines_index2 = lines_index[:index]

    lines_position1 = lines_position[index:]
    lines_position2 = lines_position[:index]

    if max_split_value < min_split_score: return 
    if len(lines_position1) < min_residue: return
    if len(lines_position2) < min_residue: return

    split_indexs.append(split_index)


    get_all_split(lines_index1, lines_position1, split_indexs)
    get_all_split(lines_index2, lines_position2, split_indexs)
    

def make_fingure(pair, splits):
    plt.figure(figsize=(8, 8))
    plt.plot(pair[:, 0], pair[:, 1], 'o', markersize=1)

    colors = 'rgbycmkw'
    for inx, split in enumerate(splits):
        plt.axvline(split, color = colors[inx])
        plt.axhline(split, color = colors[inx])
    plt.show()



def main():
    
    pdb_path = sys.argv[1]
    lines_index, lines_position = read_file(pdb_path)
    split_indexs = []
    get_all_split(lines_index, lines_position, split_indexs)
    pair = cal_distance(lines_index, lines_position)
    print("The number of domains is {}".format(len(split_indexs)+1))
    make_fingure(pair, split_indexs)


if __name__ == "__main__":
    main()
