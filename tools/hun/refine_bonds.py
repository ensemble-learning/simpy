"""
Check the number of bonds for each atom. Delete the redundant bonds. 

@note: Only hard coded for H.
@todo: read input
@date: Mon Apr 28 10:19:03 PDT 2014
@log:
2014-11-22: introduce ingore list
"""
import shutil
import os, sys

class bond():
    """
    bond information
    """
    def __init__(self,):
        self.id = 0
        self.type = 0
        self.nbond = 0
        self.bonded_atoms = []
        self.bond_orders = []
        self.totalbo = 0.0
        self.nlp = 0.0
        self.q = 0.0

def get_nparticles():
    """
    Get the number of particles(or atoms) from input. 
    """
    f = open("reaxbonds.out", "r")
    for i in f:
        tokens = i.strip().split()
        if "Number of particles" in i:
            natom = int(tokens[-1])
    f.close()
    return natom

def get_ncomment():
    """
    Get the number of comments lines from input. 
    """
    f = open("reaxbonds.out", "r")
    n = 0
    for i in f:
        if i.strip().startswith("#"):
            n += 1
        else:
            break
    f.close()
    return n

def get_nblock():
    """
    Get number of blocks.
    """
    f = open("reaxbonds.out", "r")
    n = 0
    for i in f:
        if "Timestep" in i:
            n += 1
    f.close()
    return n

def parse_bond(atom, tokens):
    """
    parse the bond from plain text.
    """
    atom.id = int(tokens[0])
    atom.type = int(tokens[1])
    atom.nbond = int(tokens[2])
    nstart = 3
    nend = nstart + atom.nbond
    for i in range(nstart, nend):
        atom.bonded_atoms.append(int(tokens[i]))
    nstart = nend + 1
    nend = nstart + atom.nbond
    for i in range(nstart, nend):
        atom.bond_orders.append(float(tokens[i]))
    nstart = nend
    atom.totalbo = float(tokens[nstart])
    nstart += 1
    atom.nlp = float(tokens[nstart])
    nstart += 1
    atom.q = float(tokens[nstart])

def check_nbonds(atoms):
    """
    check the number of bonds of each atom
    """
    #@todo: read the atom types and max bond from input.
    atomindex = []
    counter = 0
    for i in atoms:
        # hard coded for H
        if i.type == 2:
            if i.nbond > 1:
                atomindex.append(counter)
        counter += 1
    return atomindex

def build_index(atoms, index, ignore=[]):
    """
    build the bond index to delete.
    """
    #@todo: read the atom types and max bond from input, consistent with check_nbonds()
    index1 = [] # bond list for redudant bonds
    index2 = [] # bond list for tmp use
    for i in index:
        bo = []
        counter = 0
        for j in range(len(atoms[i].bond_orders)):
            tmp = atoms[i].bond_orders[j]*100000
            if ignore:
                if atoms[i].bonded_atoms[j] in ignore:
                    tmp = tmp * 0.0001
            bo.append("%08d_%04d"%(tmp, counter))
            bo.sort()
            counter += 1

        for k in bo[:atoms[i].nbond-1]:
            a1 = atoms[i].id
            a2 = atoms[i].bonded_atoms[int(k.split("_")[-1])]
            index2.append("%08d_%08d"%(a1, a2))
            index2.append("%08d_%08d"%(a2, a1))

    #@note: Here, we use set to delete the repeating terms
    index2 = set(index2)
    index2 = [i for i in index2]
    index2.sort()
    for i in index2:
        tokens = i.split("_")
        a1 = int(tokens[0])
        a2 = int(tokens[1])
        index1.append([a1, a2])
    return index1

def refine_bonds(atoms, index):
    """
    refine the connection by deleting bonds.
    """
    n = 0
    for i in atoms:
        flag = 1
        while(flag and (n < len(index))):
            flag = 0
            if index[n][0] == i.id:
                #print i.id, i.bonded_atoms
                id = i.bonded_atoms.index(index[n][1])
                i.totalbo = i.totalbo - i.bond_orders[id]
                i.nbond = i.nbond - 1
                del i.bonded_atoms[id]
                del i.bond_orders[id]
                flag = 1
                n += 1

def output_out(atoms, comments, outfile):
    """
    Output the data
    """
    if len(comments) < 8:
        outfile.write("#\n")
    for i in comments:
        outfile.write(i)
    for i in atoms:
        outfile.write("%6d"%i.id)
        outfile.write("%4d"%i.type)
        outfile.write("%4d"%i.nbond)
        for j in i.bonded_atoms:
            outfile.write("%6d"%j)
        outfile.write("%4d"%0)
        for j in i.bond_orders:
            outfile.write("%8.4f"%j)
        outfile.write("%8.4f"%i.totalbo)
        outfile.write("%8.4f"%i.nlp)
        outfile.write("%8.4f"%i.q)
        outfile.write("\n")

def read_index(fname):
    ndx = []
    f = open(fname, "r")
    for i in f:
        tokens = i.strip()
        if len(tokens) > 0:
            ndx.append(int(tokens))
    f.close()
    return ndx
                
def main():
    o = open("tmp.out", "w")
    natom = get_nparticles()        
    nblock = get_nblock()
    
    log = open("refine.log", "w")
    log.write("The number of atoms is %d\n"%natom)
    log.write("Total number of snapshots is %d\n"%nblock)
    log.flush()
    
    f = open("reaxbonds.out", "r")

    prev = ''
    flag = 1
    
    ignore = []
    if os.path.exists("ignore.ndx"):
        sys.stdout.write("Reading ignore list.\n")
        sys.stdout.flush()
        ignore = read_index("ignore.ndx")

    # Parse the block
    for n in range(nblock):
        comments = []
        atoms = []

        # Parse the comments
        if len(prev.strip()) > 0 :
            comments.append(prev)

        for i in f:
            tokens = i.strip().split()
            if i.strip().startswith("#"):
                comments.append(i)
                if "Timestep" in i:
                    timestep = int(tokens[-1])
            else:
                break

        # Parse the atoms
        atom = bond()
        parse_bond(atom, tokens)
        atoms.append(atom)
    
        for i in f:
            tokens = i.strip().split()
            if i.strip().startswith("#"):
                break
            else:
                atom = bond()
                parse_bond(atom, tokens)
                atoms.append(atom)
        # Do the analysis
        index1 = check_nbonds(atoms)
        if ignore:
            index2 = build_index(atoms, index1, ignore)
        else:
            index2 = build_index(atoms, index1)
        refine_bonds(atoms, index2)
        output_out(atoms, comments, o)
        
    # clear up
    log.close()
    o.close()
    f.close()
    # copy the files
    #shutil.copy("reaxbonds.out", "reaxbonds.out.bak")
    #shutil.copy("tmp.out", "reaxbonds.out")
    #os.remove("tmp.out")

if __name__ == "__main__":
    main()
