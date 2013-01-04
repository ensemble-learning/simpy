f = open("crystal.top", 'r')
o = open("crystal_1.top", 'w')

for i in f:
    i = i.replace("-0.3146", "-0.4520")
    i = i.replace(" 0.9429", " 0.3411")
    i = i.replace("-0.3137", " 0.2424")
    i = i.replace(" 0.3063", " 0.4039")
    i = i.replace("-0.9967", "-0.8109")
    i = i.replace(" 0.3841", " 0.3236")
    o.write(i)

f.close()
o.close()
 
