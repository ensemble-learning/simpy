import math

from parse_gro import parse_gro
from parse_top import parse_top

DEBYE = 160.2177/3.3564

coords, box = parse_gro('run.gro')
top = parse_top('out.top')

dip_x = 0
dip_y = 0
dip_z = 0

if len(top) == len(coords):
    for i in range(len(top)):
        dip_x += coords[i][-3] * float(top[i].split()[-2])
        dip_y += coords[i][-2] * float(top[i].split()[-2])
        dip_z += coords[i][-1] * float(top[i].split()[-2])

print math.sqrt(dip_x*dip_x + dip_y*dip_y + dip_z*dip_z)*DEBYE
