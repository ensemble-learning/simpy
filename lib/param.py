from ffield import Ffield

class Param():
    def __init__(self, filename="param", ntype=0):
        self.params = []
        self.comments = []
        self.ntype = ntype
        self.read(filename)
    def read(self, filename):
        f = open(filename, "r")
        for i in f:
            if "!" in i:
                tokens = i.strip().split("!")
                data = tokens[0].strip().split()
                comment = tokens[1]
            else:
                data = i.strip().split()
                comment = ""
            if len(data) == 6:
                self.params.append(data)
                self.comments.append(comment)
    def update(self,):
        ff = Ffield("ffield", self.ntype)
        for i in self.params:
            nsec = int(i[0])
            if nsec == 1:
                pass
            elif nsec == 2:
                np = int(i[1]) - 1
                nval = int(i[2])
                val = ff.atom[np][nval]
            elif nsec == 3:
                np = int(i[1]) - 1
                nval = int(i[2]) + 1
                val = ff.bond[np][nval]
            elif nsec == 4:
                np = int(i[1]) - 1
                nval = int(i[2]) + 1
                val = ff.off[np][nval]
    def filter(self, a1, a2, eq):
        for i in self.comments:
            print i
          
def test():
    a = Param("params_total")
    a.filter("Li", "Li", "")

if __name__ == "__main__":
    test()
