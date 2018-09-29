#!/usr/bin/env python

"""generate each fragment as a funtion of time"""

def parse_mol(file):
    """parse .mol file (reaxc)"""
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
            
def output_csv():
    """output the fragment number as a funtion of time.
    The lines of the file can be controled by maxlines"""

    names = []
    blocks = parse_mol("water.mol")
    test = stat_all(blocks)

    counter = 0
    seq = 0
    skip = 1
    maxlines = 1000
    for i in test.keys():
        name = "%010d_%s"%(len(test[i][0]),i)
        names.append(name)
    names.sort()

    for tokens in names:
        i = tokens.split('_')[-1]
        o = open("m_%04d_%s.csv"%((len(names) - seq),i),"w")
        skip = len(test[i][0])/ maxlines + 1
        for j in range(len(test[i][0])):
            if counter % skip == 0:
                o.write(test[i][0][j]+','+ test[i][1][j]+'\n')
            counter += 1
        seq += 1
        o.close()

def stat_number(ntime):
    """counter the number of each fragment in the simulation
    after ntime"""

    blocks = parse_mol("water.mol")
    test = stat_all(blocks)
    counter = 0
        
    for i in test.keys():
        sum = 0
        counter = 0
        for j in range(len(test[i][0])):
            if int(test[i][0][j]) > ntime:
                sum += int(test[i][1][j])
                counter += 1
        if counter > 0:
            print i, sum


if __name__ == "__main__":
    output_csv()
