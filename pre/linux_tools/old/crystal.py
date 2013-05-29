import math
import os

os.system("tail -1 confout.gro > test.log ")

f = open ( "test.log" ,'r')
o = open ( "out.log " , 'w')

a = []

for i in f:
    temp = i.split()
    for m in range(9):
        a.append(float(temp[m]))
boxa = math.sqrt(a[0]*a[0]+a[3]*a[3]+a[4]*a[4])
boxb = math.sqrt(a[5]*a[5]+a[1]*a[1]+a[6]*a[6])
boxc = math.sqrt(a[7]*a[7]+a[8]*a[8]+a[2]*a[2])

alpha = math.acos((a[5]*a[7] + a[1]*a[8] + a[6]*a[2]) / boxb / boxc) / math.pi * 180
beta = math.acos((a[0]*a[7] + a[3]*a[8] + a[4]*a[2]) / boxa /boxc) /  math.pi * 180
gama = math.acos((a[0]*a[5] + a[3]*a[1] + a[4]*a[6] ) /boxa/ boxb) / math.pi * 180

print boxa
print boxb
print boxc
print alpha
print beta
print gama
f.close()
o.close()
