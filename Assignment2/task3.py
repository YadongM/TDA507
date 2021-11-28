# Author: Yadong Mao

# To run this code, use command:
# python ./task3.py ./2csn.pdb

import sys
import numpy as np


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


def split_group(lines_index, lines_position):

    def get_swap(group_1, group_2):
        max_count = 0
        for i in group_1.keys():
            _simple_count = 0
            for j in group_2.keys():
                _distance = np.sqrt(np.sum( (group_1[i]-group_2[j])**2 ))
                if (_distance<=threshold):
                    _simple_count += 1

            if max_count < _simple_count:
                max_count = _simple_count
                _swap = i
        return max_count, _swap


    dic_position = dict(zip(lines_index, lines_position))
    _split = int(len(lines_position)/2)
    threshold = 10

    sorted = []

    group_U = dict(list(dic_position.items())[:_split] )
    group_V = dict(list(dic_position.items())[_split:] )

    for _ in range(int(len(lines_position)/2)):
        print('\r' + 'Processing {:.2f}%'.format(_/len(lines_position)*100), end='', flush=True)
        max_count1, _swap1 = get_swap(group_U,group_V)
        max_count2, _swap2 = get_swap(group_V,group_U)
        
        if max_count1 > max_count2:
            group_U.pop(_swap1)
            group_V[_swap1] = dic_position[_swap1]

        else:
            group_V.pop(_swap2)
            group_U[_swap2] = dic_position[_swap2]

        sorted.append(_swap1)
        sorted.append(_swap2)

    best_combine1 = list(set(list(group_U.keys())))
    best_combine2 = list(set(list(group_V.keys())))
            
    return best_combine1, best_combine2



def main():
    
    pdb_path = sys.argv[1]
    lines_index, lines_position = read_file(pdb_path)
    best_combine1, best_combine2 = split_group(lines_index, lines_position)

    print()
    print("Group 1")
    print( best_combine1 )
    print()
    print("Group 2")
    print( best_combine2 )

if __name__ == "__main__":
    main()







