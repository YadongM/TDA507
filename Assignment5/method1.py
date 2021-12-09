# Author: Yadong Mao

# To run this code, use command:

# 1CDH_2CSN
# python ./method1.py ./1cdh.pdb ./2csn.pdb

# 2CSN_1CDH
# python ./method1.py ./2csn.pdb ./1cdh.pdb


import numpy as np
import sys

ATOM_RADIUS = 2

def read_file(pdb_path):

    with open(pdb_path, 'r', encoding="utf-8") as f:
        other_info = []
        position = []
        for line in f:
            line = line.split()
            if line[0] == 'ATOM' or line[0] == 'HETATM':
                position.append([float(line[6]), float(line[7]), float(line[8])])
                other_info.append([line[1], line[3], line[5], line[11]])

    return position, other_info

def cal_distance(atom1, atom2):
    return np.sqrt(np.sum( (atom1 - atom2)**2 ))


def overla_test(atom1, atom2, ATOM_RADIUS):
    if cal_distance(atom1, atom2) < ATOM_RADIUS*2:
        return 1


def main():

    pdb_path1 = sys.argv[1]
    pdb_path2 = sys.argv[2]

    position1, atom_info1 = read_file(pdb_path1)
    position2, atom_info2 = read_file(pdb_path2)

    cal_times = 0
    overlap_atom = []

    for index, atom2 in enumerate(position2):
        for atom1 in position1:
            cal_times += 1

            if overla_test(np.array(atom1), np.array(atom2), ATOM_RADIUS):
                tmp = [atom_info2[index][0], atom_info2[index][1], atom_info2[index][2], atom_info2[index][3]]
                overlap_atom.append(tmp)
                print(atom_info2[index][0], atom_info2[index][1], atom_info2[index][2], atom_info2[index][3])
                break

    print('Number of clashing atoms: {}'.format(len(overlap_atom)))
    print('Number of comparisons: {}'.format(cal_times))

if __name__ == '__main__':
    main()