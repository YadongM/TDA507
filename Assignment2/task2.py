# Author: Yadong Mao

# To run this code, use command:
# python ./task2.py ./2csn.pdb

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
    
    lines_index = [line[5] for line in lines_atom]
    lines_index = np.array(lines_index).astype(int)

    return lines_index, lines_position



def cal_distance(lines_index, lines_position, threshold=10):
    pair = []
    # all_distance = np.zeros([len(lines_position), len(lines_position)])
    for i in range(len(lines_position)):
        for j in range(len(lines_position)):

            _distance = np.sqrt(np.sum( (lines_position[i]-lines_position[j])**2 ))
            # all_distance[i][j] = _distance

            if (_distance<=threshold):
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
    threshold=10
    lines_index, lines_position = read_file(pdb_path)
    pair = cal_distance(lines_index, lines_position, threshold)
    max_split_value, split_index, index = cal_domak(lines_index, pair)
    print("When the threshold = {}".format(threshold) )
    print("The max DOMAK score is {:.2f}".format(max_split_value) )
    print("The s_iCode of this position is {}".format(index) )
    print("The index of this position is {}".format(split_index) )
    
    make_fingure(pair, [index])


if __name__ == "__main__":
    main()