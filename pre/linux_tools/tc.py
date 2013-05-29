#!/usr/bin/env python

def readTempgrad(file):
    f = open(file, 'r')

    temp_slice = []

    for i in f:
        for j in range(20):
            temp_slice.append([])
            temp_slice[j].append(i.split()[j+2])
    f.close()
            
    return temp_slice

def averageTempgrad(a2Dlist):
    sum = []
    for i in range(20):
        temp = 0
        for j in a2Dlist[i]:
            temp = temp + float(j)
        sum.append(temp/len(a2Dlist[i]))

    return sum

def write_list(list, file = 'out.log'):
    o = open(file, 'w')
    for i in list:
        o.write(str(i) + '\n')
    o.close()
    
def readEnergy(file):
    f = open(file, 'r')

    sum = 0
    
    for i in f:
        sum = sum + float(i.split()[2])

    return sum
        
    
tempgrad = readTempgrad('tempgradC.energy')
sumall = averageTempgrad(tempgrad)
write_list(sumall, 'tempgrad.log')
energy = readEnergy('exchangeC.energy')
print energy

