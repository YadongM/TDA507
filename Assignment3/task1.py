# Author: Yadong Mao

# To run this code, use command:
# python .\task1.py .\test_q1.txt

# python .\task1.py .\data_q1.txt



import sys
import numpy as np

def read_file(text_path):
    _tmp = []
    with open(text_path, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.split()
            _tmp.append(line)
    return np.array(_tmp).astype(float)

def find_order_of_chain(lines):
    atom_dic = {}
    beginning = []

    # calculate the distance, 
    # and fit the atom index whoes distance between 3.7 and 3.8
    for i in range(len(lines)):
        _tmp = []
        for j in range(len(lines)):
            _distance = np.sqrt(np.sum( (lines[i][1:]-lines[j][1:])**2 ))
            if (_distance<3.9) and (_distance>3.7):
                _tmp.append(lines[j][0])
        if len(_tmp) == 1:
            beginning.append(lines[i][0])
        atom_dic[lines[i][0]] = _tmp

    # Random choose a atom as start
    _sorted = []
    current_atom = beginning[1]
    _sorted.append(current_atom)
    print(beginning[1])

    # find order by dictory
    while True:
        if len(_sorted) == len(atom_dic):
            break
        
        for j in atom_dic[current_atom]:
            if j not in _sorted:
                _sorted.append(j)
                print(j)
                current_atom = j
    return _sorted
def main():
    # text_path = './test_q1.txt'
    text_path = sys.argv[1]
    lines = read_file(text_path)
    chain = find_order_of_chain(lines)
    print("The total number of alpha-carbon atoms in the chain {}".format(len(chain)))

if __name__ == '__main__':
    main()