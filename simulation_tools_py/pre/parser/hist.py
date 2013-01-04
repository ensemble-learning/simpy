f = open('po.log' , 'r')

min = 500
max = 0
ener_list = []
sum = 0

for i in f:
    ener = float(i.strip())
    if min > ener:
        min = ener
    if max < ener:
        max = ener
    ener_list.append(ener)

step = (max-min)/50.0

hist = [0]*101

for i in ener_list:
    sum += i
    hist[int((i-min)/step)] +=1

for i in hist:
    print i

print '------------average--------------'
print sum / len(ener_list)



