import numpy as np
import sys

ATOM_RADIUS = 2

def read_file(pdb_path):
    atom = []

    with open(pdb_path, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.split()
            if line[0] == 'ATOM' or line[0] == 'HETATM':
                atom.append(line)
    return np.array(atom)

def cal_distance(atom1, atom2):
    return np.sqrt(np.sum( (atom1 - atom2)**2 ))


def overla_test(atom1, atom2, ATOM_RADIUS):
    if cal_distance(atom1, atom2) < ATOM_RADIUS*2:
        return 1


def main():
    # pdb_path1 = './1cdh.pdb'
    # pdb_path2 = './2csn.pdb'
    pdb_path1 = sys.argv[1]
    pdb_path1 = sys.argv[2]

    all_atom1 = read_file(pdb_path1)
    all_atom2 = read_file(pdb_path2)
    position = all_atom2[:,[6,7,8]].astype(np.float64)

    cal_times = 0
    overlap_atom = []


    for atom1 in all_atom1:
        x = float(atom1[6])
        y = float(atom1[7])
        z = float(atom1[8])

        x_index = np.intersect1d( np.where(position[:, 0]>(x-ATOM_RADIUS*2))[0], np.where(position[:, 0]<(x+ATOM_RADIUS*2))[0] ) 
        y_index = np.intersect1d( np.where(position[:, 1]>(y-ATOM_RADIUS*2))[0], np.where(position[:, 1]<(y+ATOM_RADIUS*2))[0] ) 
        z_index = np.intersect1d( np.where(position[:, 2]>(z-ATOM_RADIUS*2))[0], np.where(position[:, 2]<(z+ATOM_RADIUS*2))[0] )
        
        _index = np.intersect1d(x_index, np.intersect1d(y_index,z_index))

        for atom2 in all_atom2[_index]:
            cal_times += 1
            if overla_test(atom1[6:9].astype(np.float64), atom2[6:9].astype(np.float64), ATOM_RADIUS):
                tmp = [atom1[1], atom1[3], atom1[5], atom1[11]]
                overlap_atom.append(tmp)
                break

    print('\nNumber of clashing atoms: {}'.format(len(overlap_atom)))
    print('Number of comparisons: {}'.format(cal_times))