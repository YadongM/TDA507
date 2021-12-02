# Author: Yadong Mao

# To run this code, use command:
# python task2_method2.py

import sys
import copy
import numpy as np
import math

def read_file(text_path):
    _tmp = []
    with open(text_path, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.split()
            _tmp.append(line)
    return np.array(_tmp).astype(float)

def find_relation_of_atom(lines):
    atom_relation_dic = {}
    beginnings = {}

    # calculate the distance, 
    # and fit the atom index whoes distance between 3.7 and 3.9
    for i in range(len(lines)):
        _tmp = []
        for j in range(len(lines)):
            _distance = np.sqrt(np.sum( (lines[i][1:]-lines[j][1:])**2 ))
            if (_distance<3.9) and (_distance>3.7):
                _tmp.append(lines[j][0])
        if len(_tmp) == 1:
            beginnings[lines[i][0]] = _tmp
        if len(_tmp) > 0:
            atom_relation_dic[lines[i][0]] = _tmp

    return atom_relation_dic, beginnings


def list2dic(lines):
    atom_dic = {}
    for line in lines:
        atom_dic[int(line[0])] = line[1:]
    return atom_dic

def cal_center_distance(atom_dic, atom_relation_dic):
    center_distance = {}
    for index in atom_relation_dic.keys():
        _tmp = [0,0,0]
        _tmp_len = len(atom_relation_dic[index])
        _tmp_distance = []

        for relate_atom in atom_relation_dic[index]:
            _tmp += atom_dic[relate_atom]
        center_posi = np.array(_tmp)/_tmp_len

        for relate_atom in atom_relation_dic[index]:
            _distance = np.sqrt(np.sum( ( atom_dic[relate_atom] - center_posi )**2 ))
            _tmp_distance.append(_distance)
        
        center_distance[index] = _tmp_distance
        
        # closest_atom[index] = atom_relation_dic[index][np.array(_tmp_distance).argmin()]
    return center_distance


def cal_angle(position1, position2, position3):

    a_x, b_x, c_x = position1[0], position2[0], position3[0]
    a_y, b_y, c_y = position1[1], position2[1], position3[1]
    a_z, b_z, c_z = position1[2], position2[2], position3[2]
    x1,y1,z1 = (a_x-b_x),(a_y-b_y),(a_z-b_z)
    x2,y2,z2 = (c_x-b_x),(c_y-b_y),(c_z-b_z)

    _tmp = (x1*x2 + y1*y2 + z1*z2) / (math.sqrt(x1**2 + y1**2 + z1**2) *(math.sqrt(x2**2 + y2**2 + z2**2)))

    return math.degrees(math.acos(_tmp))


def find_path_2(root, child1, atom_relation_dic, chain, chains, atom_dic):
    _tmp_related_atom1 = atom_relation_dic[child1]

    chain = copy.deepcopy(chain)
    chain.append(child1)

    for _child in _tmp_related_atom1:
        if _child in chain:
            _tmp_related_atom1.remove(_child)
        
        if _child not in list(atom_relation_dic.keys()):
            _tmp_related_atom1.remove(_child)
    
    if len(_tmp_related_atom1) == 0:
        chains.append(chain)
        return

    for child2 in _tmp_related_atom1:
        if 75 < cal_angle(atom_dic[root],atom_dic[child1],atom_dic[child2]) < 160:
            root = child1
            child1 = child2
            find_path_2(root, child1, atom_relation_dic, chain, chains, atom_dic)


text_path = './data_q2.txt'
lines = read_file(text_path)
atom_dic = list2dic(lines)
atom_relation_dic, beginnings = find_relation_of_atom(lines)


chains = []

for child1 in beginnings.keys():
    chain = []
    chain.append(child1)

    for child2 in atom_relation_dic[child1]:
        if (child2 not in chain) and (child2 in atom_relation_dic.keys()):
            find_path_2(child1, child2, atom_relation_dic, chain, chains, atom_dic)

len_chain = []

for chain in chains:
    len_chain.append(len(chain))

for atom in chains[np.argmax(len_chain)]:
    print(atom)
print("The total number of alpha-carbon atoms in the chain {}".format(np.max(len_chain)))

