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

    # lines_series = [line[1] for line in line_atom]
    # lines_series = np.array(lines_series).astype(int)

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

# pdb_path = sys.argv[1]
pdb_path = './2csn.pdb'
lines_index, lines_position = read_file(pdb_path)
pair = cal_distance(lines_index, lines_position)
max_split_value, split_index, index = cal_domak(lines_index, pair)





dic_position = dict(zip(lines_index, lines_position))
sorted_keys = []

_split = int(len(lines_position)/2)

group_U = dict(list(dic_position.items())[:_split])
group_V = dict(list(dic_position.items())[_split:])
group_U_keys = list(group_U.keys())
group_V_keys = list(group_V.keys())


threshold = 10
max_count = 0

while True:
    _simple_count = 0
    for i in group_U.keys():
        for j in group_V.keys():
            _distance = np.sqrt(np.sum( (group_U[i]-group_V[j])**2 ))
            if (_distance<=threshold):
                _simple_count += 1
    
    if max_count<_simple_count:
        best_combine1 = list(group_U.keys())
        best_combine2 = list(group_V.keys())

    if (len(group_U_keys)<2) or (len(group_V_keys)<2):
        break

    a = np.random.choice(group_U_keys, 2)
    b = np.random.choice(group_V_keys, 2)

    for i in a:
        group_U.pop(i)
        group_V[i] = dic_position[i]
        group_U_keys.remove(i)
        
    for i in b:
        group_V.pop(i)
        group_U[i] = dic_position[i]
        group_V_keys.remove(i)
