#!/usr/bin/env python

def parse_mol(file):
    in_block = 1
    out_block = 0
    f = open(file, 'r')
    state = 1
    counter = 0
    blocks = []
    tokens = []

    for i in f:
        if len(i.strip()) == 0:
            state = 0
            counter = 0
            tokens = []
        else:
            state = 1
        if state == in_block:
            if counter == 0:
                tokens.append(i.split()[0][4:])
            else:
                tokens.append(i.strip())
            counter += 1
        if len(tokens) > 0:
            blocks.append(tokens)
    
    return blocks

def stat_all( block ):
    mol_in_simu = {}
    for i in block:
        for j in range(len(i)):
            if j == 0:
                steps = i[j]
            else:
                mol = i[j].split()[-1]
                number = i[j].split()[0]
                if mol not in mol_in_simu.keys():
                    mol_in_simu[mol] = []
                    mol_in_simu[mol].append([])
                    mol_in_simu[mol].append([])
                mol_in_simu[mol][0].append(steps)
                mol_in_simu[mol][1].append(number)
    return mol_in_simu
            
def test():
    blocks = parse_mol("water.mol")
    test = stat_all(blocks)

    counter = 0
    skip = 1

    for i in test.keys():
        o = open("%s.csv"%i,"w")
        skip = len(test[i][0])/ 1000 + 1
        for j in range(len(test[i][0])):
            if counter % skip == 0:
                o.write(test[i][0][j]+','+ test[i][1][j]+'\n')
            counter += 1
        o.close()

def stat_number():
    blocks = parse_mol("water.mol")
    test = stat_all(blocks)
    counter = 0
        
    for i in test.keys():
        sum = 0
        counter = 0
        for j in range(len(test[i][0])):
            if int(test[i][0][j]) > 8300000:
                sum += int(test[i][1][j])
                counter += 1
        if counter > 0:
            print i, sum


if __name__ == "__main__":
    test()
