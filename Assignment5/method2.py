# Author: Yadong Mao

# To run this code, use command:

# 1CDH_2CSN
# python ./method2.py ./1cdh.pdb ./2csn.pdb

# 2CSN_1CDH
# python ./method2.py ./2csn.pdb ./1cdh.pdb


import numpy as np
import sys

ATOM_RADIUS = 2

def read_file(pdb_path):
    '''
    Used to read pdb files.
    Input: The path of pdb file
    Return: A np.array include line contain ATOM and HETATM
    '''
    atom = []
    with open(pdb_path, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.split()
            if line[0] == 'ATOM' or line[0] == 'HETATM':
                tmp = [line[1], line[3], line[5], line[11], line[6], line[7], line[8]]
                atom.append(tmp)
    return np.array(atom)

def cal_distance(atom1, atom2):
    '''
    Calculate the distance between two atoms.
    Input: Coordinates of two atoms
    Output: The distance between them(float)
    '''
    return np.sqrt(np.sum( (atom1 - atom2)**2 ))


def overla_test(atom1, atom2, ATOM_RADIUS):
    '''
    Check if two atoms overlap
    Input: Coordinates of two atoms, Atomic radius
    Output: Return 1 if overlap
    '''
    if cal_distance(atom1, atom2) < ATOM_RADIUS*2:
        return 1


def main():
    # pdb_path1 = './1cdh.pdb'
    # pdb_path2 = './2csn.pdb'
    pdb_path1 = sys.argv[1]
    pdb_path2 = sys.argv[2]

    all_atom1 = read_file(pdb_path1)
    all_atom2 = read_file(pdb_path2)

    '''
    Use np.array to store the location information of the first file, 
    the data format is float64. 
    Use np.array to easily find the index of the atom that meets the condition.
    '''
    position = all_atom1[:,[4,5,6]].astype(np.float64)


    cal_times = 0
    overlap_atom = []


    for atom2 in all_atom2:
        x = float(atom2[4])
        y = float(atom2[5])
        z = float(atom2[6])

        '''
        For each atom in file two, when we want to traversal the atoms in file one, we first filter the atoms in file one.
        Suppose an atom coordinate in the current file 2 is (x, y, z)
        The atomic position in the filter file 1 satisfies:

        x-ATOM_RADIUS*2 < x_coordinate < x+ATOM_RADIUS*2 -------- Set 1 of atomic indexes in file one
        y-ATOM_RADIUS*2 < y_coordinate < y+ATOM_RADIUS*2 -------- Set 2 of atomic indexes in file one
        z-ATOM_RADIUS*2 < z_coordinate < z+ATOM_RADIUS*2 -------- Set 2 of atomic indexes in file one

        Take the intersection of three sets.
        We only need to traverse the atoms of file one whose index is in this intersection.
        '''

        x_index = np.intersect1d( np.where(position[:, 0]>(x-ATOM_RADIUS*2))[0], np.where(position[:, 0]<(x+ATOM_RADIUS*2))[0] ) 
        y_index = np.intersect1d( np.where(position[:, 1]>(y-ATOM_RADIUS*2))[0], np.where(position[:, 1]<(y+ATOM_RADIUS*2))[0] ) 
        z_index = np.intersect1d( np.where(position[:, 2]>(z-ATOM_RADIUS*2))[0], np.where(position[:, 2]<(z+ATOM_RADIUS*2))[0] )
        
        _index = np.intersect1d(x_index, np.intersect1d(y_index,z_index))

        for atom1 in all_atom1[_index]:
            cal_times += 1
            if overla_test(atom1[4:].astype(np.float64), atom2[4:].astype(np.float64), ATOM_RADIUS):
                tmp = [atom2[0], atom2[1], atom2[2], atom2[3]]
                overlap_atom.append(tmp)
                print(atom2[0], atom2[1], atom2[2], atom2[3])
                break

    print('Number of clashing atoms: {}'.format(len(overlap_atom)))
    print('Number of comparisons: {}'.format(cal_times))

if __name__ == '__main__':
    main()