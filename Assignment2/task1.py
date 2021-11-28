# Author: Yadong Mao

# To run this code, use command:
# python ./task1.py ./2csn.pdb

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
    all_distance = np.zeros([len(lines_position), len(lines_position)])
    for i in range(len(lines_position)):
        for j in range(len(lines_position)):
            _distance = np.sqrt(np.sum( (lines_position[i]-lines_position[j])**2 ))
            all_distance[i][j] = _distance
            if (i!=j and _distance<=threshold):
                pair.append([lines_index[i], lines_index[j]])
    return np.array(pair)



def make_fingure(pair):
    plt.figure(figsize=(8, 8))
    plt.plot(pair[:, 0], pair[:, 1], 'o', markersize=1)
    plt.show()

    # colors = 'rgbycmkw'
    # for inx, split in enumerate(splits):
    #     plt.axvline(split, color = colors[inx])
    #     plt.axhline(split, color = colors[inx])

    
def main():
    pdb_path = sys.argv[1]
    lines_index, lines_position = read_file(pdb_path)
    pair = cal_distance(lines_index, lines_position, threshold=10)
    names = str(pdb_path).replace("/",'').split('.')

    np.savetxt(names[1]+'.pair', pair, fmt="%d")
    print("The pair file save as {}".format(names[1]+'.pair'))
    make_fingure(pair)


if __name__ == "__main__":
    main()