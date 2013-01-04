f = open("PDMS.top", 'r')

output = []
delta = 0.15
for i in f:
    if len(i.split()) == 8:
        if i.split()[4] == "Si":
            charge = float(i.split()[6]) - delta
            lines = i[:65] + "%7.4f"%charge + i[72:]
        elif i.split()[4] == "O":
            charge = float(i.split()[6]) + delta
            lines = i[:65] + "%7.4f"%charge + i[72:]
        else:
            lines = i
    else:
        lines = i
    output.append(lines)

f.close()

o = open("test.top" ,'w')
for i in output:
    o.write(i)
o.close()

#     1       SiC2O2      1          UNK           Si      0       0.9116      28.0855
#     2        CH3Si      2          UNK            C      1      -0.2889      12.0110
#     3           HC      3          UNK            H      2       0.0430       1.0079
#     4           HC      4          UNK            H      3       0.0430       1.0079
#     5           HC      5          UNK            H      4       0.0430       1.0079
#     6        CH3Si      6          UNK            C      5      -0.2889      12.0110
#     7           HC      7          UNK            H      6       0.0430       1.0079
#     8           HC      8          UNK            H      7       0.0430       1.0079
#     9           HC      9          UNK            H      8       0.0430       1.0079
#    10         OHSi     10          UNK            O      9      -0.6347      15.9994

