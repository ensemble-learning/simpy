def write_usr(fname, ignore=[]):
    natom = 240
    o = open(fname, "w")
    for i in range(natom):
        if i in ignore:
            pass
        else:
            o.write("%d 1 1 1 0.2\n"%i)
    o.close()


#mol1 = "2 106 6 41 77 42 5 109 110 105 78 114 113 1 161 165 162 166" 
# lose two O
#mol2 = "36 7 12 40 3 43 71 68 79 75 107 112 116 120 147 164 168 151"
# lose two O and NO2
#mol3 = "35 8 11 39 4 44 72 67 80 76 108 111 115 119 148 163 167 152"
# No reaction
#mol4 = "34 33 10 37 9 65 38 145 66 117 70 69 146 74 118 73 150 149"
# lose NO2 and break a ring
#mol5 = "18 126 22 57 97 58 21 129 130 125 98 134 133 17 169 173 170 174"
# re-arrangement 10400 NO2 -> No + CO
# N2O + CO
# NO2
# NO

#mol6 = "52 23 28 56 19 59 91 88 99 95 127 132 136 140 155 172 176 159"
# release NO2 O
# release O
# relase NO + N+CN2

#mol7 = "51 24 27 55 20 60 92 87 100 96 128 131 135 139 156 171 175 160"
# break the ring very complex

#mol8 = "50 49 26 53 25 85 54 153 86 137 90 89 154 94 138 93 158 157"
# ring break (N-N)
# release NO2

line = mol8
ignore = [(int(i) - 1)  for i in line.strip().split()]
for i in range(15994):
    write_usr("dump.%d.usr"%i, ignore)



