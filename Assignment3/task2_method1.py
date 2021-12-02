# Author: Yadong Mao

# To run this code, use command:
# python task2_method1.py



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

def find_most_relation_atom(center_distance, atom_relation_dic):
    most_relation_atom = {}

    for index in center_distance.keys():
        _tmp = []

        min_distance = np.min(np.array(center_distance[index]))
        for _inx, _dist in enumerate(center_distance[index]):
            if abs(_dist-min_distance)<1.2:
                _tmp.append(int(atom_relation_dic[index][_inx]))

        most_relation_atom[int(index)] = _tmp
    
    return most_relation_atom

def find_path(child, most_relation_atom,chain,chains):
    # explored[child] = most_relation_atom[child]
    _tmp_related_atom = most_relation_atom[child]
    chain = copy.deepcopy(chain)
    chain.append(child)

    for _child in _tmp_related_atom:
        if _child in chain:
            _tmp_related_atom.remove(_child)
        
        if _child not in list(most_relation_atom.keys()):
            _tmp_related_atom.remove(_child)
    
    if len(_tmp_related_atom) == 0:
        chains.append(chain)
        return

    for child in _tmp_related_atom:
        find_path(child, most_relation_atom,chain,chains)
    

def main():
    text_path = './data_q2.txt'
    lines = read_file(text_path)
    atom_dic = list2dic(lines)
    atom_relation_dic, beginnings = find_relation_of_atom(lines)
    all_center_distance = cal_center_distance(atom_dic, atom_relation_dic)
    most_relation_atom = find_most_relation_atom(all_center_distance, atom_relation_dic)


    chains = []
    for start in beginnings.keys():
        find_path(int(start), most_relation_atom,[],chains)
    len_chain = []

    for chain in chains:
        len_chain.append(len(chain))

    for atom in chains[np.argmax(len_chain)]:
        print(atom)
    print("The total number of alpha-carbon atoms in the chain {}".format(np.max(len_chain)))

if __name__ == "__main__":
    main()