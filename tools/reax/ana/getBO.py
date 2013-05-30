"""get specified bond order from fort.8
"""

f = open("fort.8", "r")

n = 1

for i in range(n):
    f.readline()

tokens = f.readline().strip().split()

o = open("bo", "w")
o.write(tokens[13] + '\n')
o.close()

