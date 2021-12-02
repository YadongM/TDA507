import sys
import copy
import numpy as np

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



def dfs(chain, child, atom_relation_dic):
    chain = copy.deepcopy(chain)
    
    if child not in chain:
        chain.append(child) 

        if child not in list(atom_relation_dic.keys()):
            return

        elif set(atom_relation_dic[child]) < set(chain):
            return

        else:
            root = child
            for child in atom_relation_dic[root]:
                dfs(chain, child, atom_relation_dic)
    else:return


text_path = './data_q2.txt'
lines = read_file(text_path)
atom_dic = list2dic(lines)
atom_relation_dic, beginnings = find_relation_of_atom(lines)

for start in beginnings.keys():
    print(start)
    sorted = []
    chain = []
    sorted.append(start)

    dfs(chain, start, atom_relation_dic)