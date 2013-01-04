#!/usr/bin/env python

f = open("lattice.log" ,'r')

a = 1
o1 = open("a.log" ,'w')
o2 = open("b.log" ,'w')
o3 = open("c.log" ,'w')
o4 = open("alpha.log", 'w')
o5 = open("beta.log", 'w')
o6 =open("gama.log", 'w')

for i in f:
   b = a % 6 
   if b == 1:
       o1.write(i)
   elif b == 2:
       o2.write(i)
   elif b == 3:
       o3.write(i)
   elif b == 4:
       o4.write(i) 
   elif b == 5:
       o5.write(i)
   elif b == 0:
       o6.write(i)
   a += 1
o1.close()
o2.close()
o3.close()
o4.close()
o5.close()
o6.close()
f.close()
