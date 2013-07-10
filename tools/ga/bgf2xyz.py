import crystal

i = open ("geo", 'r')
fout = open ("geo.xyz", 'w')
fin = i.readlines()

xyz = ""
count = 0
for line in fin:

  token = line.split()

  if (line.find ("DESCRP") >= 0):
    header = "XYZ %s\n"%(token[1])
    count = 0

  elif (line.find ("CRYSTX") >= 0):
    a = token[1]
    b = token[2]
    c = token[3]
    alpha = token[4]
    beta = token[5]
    gamma = token[6]
    
    cry = crystal.Crystal(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    lattice = "Lattice=\"%f %f %f %f %f %f %f %f %f\""%(cry.veca[0],cry.veca[1],cry.veca[2],cry.vecb[0],cry.vecb[1],cry.vecb[2],cry.vecc[0],cry.vecc[1],cry.vecc[2])

  elif (line.find ("HETATM") >= 0):
    xyz = xyz+"%s %s %s %s\n"%(token[2],token[3],token[4],token[5])
    count+=1
  
  elif len(token) == 0:
    fout.write (header)
    fout.write ("%d\n"%count)
    fout.write (lattice+" Properties=species:S:1:pos:R:3\n")
    fout.write (xyz+"\n")
    xyz = ""

