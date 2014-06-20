"""
Make a index file for atoms. Output to gromacs format
@log: 
"""
import re
ELE = re.compile("(\w+)(\d+.*)")

class Atom():
    def __init__(self,):
        self.id = 0
        self.x = [0.0, 0.0, 0.0]
        
class System():
    def __init__(self,):
        self.pbc = []
        self.atoms = []
        self.connect = []

def read_pdb(fname):
    """
    read pdb file
    """
    pattern = re.compile(r'(\D+)(\d*)')
    system = System()
    f = open(fname, "r")
    for i in f:
        tokens = i.strip().split()
        if i.startswith("CRYST1"):
            system.pbc = [float(j) for j in tokens[1:7]]
        elif i.startswith("ATOM"):
            atom = Atom()
            atom.id = int(tokens[1])
            atom.name = tokens[2]
            atom.x[0] = float(tokens[5])
            atom.x[1] = float(tokens[6])
            atom.x[2] = float(tokens[7])
            atom.element = atom.name
            if len(tokens) > 10:
                match = pattern.match(tokens[10])
                if match:
                    atom.element = match.group(1)
            system.atoms.append(atom)
        elif i.startswith("CONECT"):
            system.connect.append([int(j) for j in tokens[1:]])
    return system
    
def output_ndx(ndx, ndxname):
    o = open("index.ndx", "a")
    o.write("[ %s ]\n"%ndxname)
    counter = 0
    for i in ndx:
        if counter > 0 and counter % 15 == 0:
            o.write("\n")
        o.write("%5d"%i)
        counter += 1
    o.write("\n")
    o.close()

def find_site(system, ndx, site):
    ndx_site = []
    for i in system.connect:
        if i[0] in ndx:
            for j in i[1:]:
                if system.atoms[j-1].element == site:
                    ndx_site.append(i[0])
    return ndx_site
            
def main():
    fname = "mmo.pdb"
    system =  read_pdb(fname)
    ndx_sur = find_surface_atoms(system)
    output_ndx(ndx_sur, "surface")
    
    ndx_site =  find_site(system, ndx_sur, "Te")
    output_ndx(ndx_site, "TeSite")
    ndx_site =  find_site(system, ndx_sur, "Mo")
    output_ndx(ndx_site, "MoSite")
    ndx_site =  find_site(system, ndx_sur, "V")
    output_ndx(ndx_site, "VSite")
    ndx_site =  find_site(system, ndx_sur, "Nb")
    output_ndx(ndx_site, "NbSite")



def find_surface_atoms(system):
    ndx = []
    for i in system.atoms:
        if i.x[2] > 30.5 or i.x[2] < 15.7:
            ndx.append(i.id)
    return ndx

if __name__ == "__main__":
    main()
    
    
    