import os

def get_mol_info():
    molname = []
    molcom = {}
    f = open('../lastFrame.atom', 'r')
    for i in f:
        if i.strip().startswith('#'):
            pass
        else:
            tokens = i.strip().split()
            if len(tokens) == 3:
                if tokens[2] in molcom.keys():
                    molcom[tokens[2]].append(tokens[0])
                else:
                    molcom[tokens[2]] = []
                    molname.append(tokens[1])
                    molcom[tokens[2]].append(tokens[0])
    f.close()
    return molname, molcom

def parse_pdb():
    header = ''
    coords = []
    f = open('output.pdb', 'r')
    for i in f:
        if i.strip().startswith('CRYST1'):
            header = i
        elif i.strip().startswith('ATOM'):
            coords.append(i)
        else:
            pass
    return header, coords
 
def get_fragment(molname, molcom, header, coords):
    for (i,j) in molcom.items():
        pdbfile = '%02d_%s.pdb'%(int(i), molname[(int(i)-1)])
        o = open(pdbfile, 'w')
        o.write(header)
        for m in j:
            for n in coords:
                if m == n.strip().split()[1]:
                    o.write(n)
        o.close()
        os.system('editconf -f %s -o %s -box 2.5 2.5 2.5'%(pdbfile, pdbfile))


if __name__ == '__main__':
    molname, molcom = get_mol_info()
    header, coords = parse_pdb()
    get_fragment(molname, molcom, header, coords)



            
