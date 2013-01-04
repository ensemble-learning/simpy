exp = {
'diisopropyl_ether': ['0.674', '7.07'], 
'methyl_tert-butyl_ether': ['0.702', '6.74'],
'ethyl_methyl_ether':['0.711', '6.31'],
'diethyl_ether':['0.697', '6.47'],
'dimethyl_ether':['0.729', '5.15'],
'dimethyl_ether':['0.729', '5.15'],
'dimethyl_ether':['0.729', '5.15'],
'dimethyl_ether':['0.729', '5.15'],
'dimethyl_ether':['0.729', '5.15'],
'dimethyl_ether':['0.729', '5.15'],
}
f = open("result.csv", 'r')
cal = []
line = []
for i in f:
    terms = i.split(',')
    if len(terms[0].split('#'))<3:
        pass
    else:
        line = [terms[0].split('#')[1],terms[0].split('#')[2],terms[1],terms[2]]
        cal.append(line)

o = open("scan.csv", 'w')
for i in cal:
    o.write("%s,"%i[0])
    o.write("%s,"%(','.join([j for j in i[1].split('_')])))
    o.write("%s,%s,%s,%s\n"%(i[2], exp[i[0]][0], i[3].strip(), exp[i[0]][1]))

o.close()

