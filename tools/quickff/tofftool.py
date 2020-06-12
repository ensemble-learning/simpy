import datetime
ff_file = 'pars_yaff_mbisvasp.txt'
# read charges
f = open(ff_file, 'r')
lines = f.readlines()
f.close()

charge = []
for i in lines:
    if i.strip().startswith('FIXQ:ATOM'):
        charge.append(i)

# read the bond terms
f = open(ff_file, 'r')
lines = f.readlines()
f.close()

sections = []
for i in lines:
    if i.strip().startswith("# "):
        tokens = i.strip().split("#")[1].strip()
        sections.append(tokens)
        
bond, angle, torsion, oop = [], [], [], []
for i in lines:
    if i.startswith('BONDHARM:PARS'):
        bond.append(i.strip())
    elif i.startswith('BENDAHARM:PARS'):
        angle.append(i.strip())
    elif i.startswith('TORSION:PARS'):
        torsion.append(i.strip())
    elif i.startswith('OOPDIST:PARS'):
        oop.append(i.strip())

o = open('quickff.ff', 'w')
o.write('# %s\n'%datetime.datetime.now())

o.write('\nATOMS\n')
for i in charge:  
    tokens = i.split()
    a1 = tokens[1]
    q  = float(tokens[2])
    o.write('%-4s%4s  '%(a1, a1))
    o.write('%8.3f'%0)
    o.write('%7.3f'%q)
    o.write('  lj  %8.3f%10.5f'%(3.0, 0.5))
    o.write('\n')

o.write('\nBONDS\n')
for i in bond:
    tokens = i.split()
    a1, a2 = tokens[1], tokens[2]
    re, kr = float(tokens[4]), float(tokens[3])
    o.write('%-4s%4s  '%(a1, a2))
    o.write('  harm  ')
    o.write('%8.3f%8.1f'%(re, kr))
    o.write('\n')

o.write('\nANGLES\n')
for i in angle:
    tokens = i.split()
    a1, a2, a3 = tokens[1], tokens[2], tokens[3]
    th, ka = float(tokens[5]), float(tokens[4])
    o.write('%-4s%4s%4s  '%(a1, a2, a3))
    o.write('  harm  ')
    o.write('%8.3f%8.1f'%(th, ka))
    o.write('\n')

o.write('\nDIHEDRALS\n')
for i in torsion:
    tokens = i.split()
    a1, a2, a3, a4 = tokens[1], tokens[2], tokens[3], tokens[4]
    vn = [0.0, 0.0, 0.0, 0.0]
    parms = tokens[5:]
    for n in range(int(len(parms)/3)):
        n_index = int(tokens[n*3+5]) - 1
        vn[n_index] = float(tokens[n*3+6])
    o.write('%-4s%4s%4s%4s  '%(a1, a2, a3, a4))
    o.write('  opls  ')
    for v in vn:
        o.write('%8.4f'%(v))
    o.write('\n')

"""
o.write('\nIMPROPER\n')
"""
o.close()

