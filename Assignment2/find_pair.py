import sys
import numpy as np


def read_file(pdb_path):
    lines_atom = []
    with open(pdb_path, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.split()
            if (line[0] == 'ATOM' and line[2] == 'CA'):
                lines_atom.append(line)
    return lines_atom

def get_positions(line_atom):
    lines_position = [line[6:9] for line in line_atom]
    lines_position = np.array(lines_position).astype(float)

    lines_index = [line[5] for line in line_atom]
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

def main():
    pdb_path = sys.argv[1]
    lines_atom = read_file(pdb_path)
    lines_index, lines_position = get_positions(lines_atom)
    pair = cal_distance(lines_index, lines_position, threshold=10)
    np.savetxt('result.pair',pair, fmt="%d")

    print('Done, Please check result.pair')


if __name__ == "__main__":
    main()