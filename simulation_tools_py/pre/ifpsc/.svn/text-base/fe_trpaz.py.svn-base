#parseFile
f = open("freeEnergy.log", 'r')
dmm = {}
water = {}
purewater = [[],[]]
puredmm = [[],[]]
for i in f:
    if "dmm" in i:
        a = i.split()[0]
        ct = a[3:6]  #concentration
        step = float(a[-3:])/100.0
        if a[3:6] not in dmm.keys():
            dmm[ct] = [[], []]
        dmm[ct][0].append(step)
        dmm[ct][1].append(float(i.split()[2]))
    elif "water" in i:
        a = i.split()[0]
        ct = a[3:6]  #concentration
        step = float(a[-3:])/100.0
        if a[3:6] not in water.keys():
            water[ct] = [[], []]
        water[ct][0].append(step)
        water[ct][1].append(float(i.split()[2]))
    elif "a000/a" in i:
        a = i.split()[0]
        step = float(a[-3:])/100.0
        purewater[0].append(step)
        purewater[1].append(float(i.split()[2]))
    elif "a100/a" in i:
        a = i.split()[0]
        step = float(a[-3:])/100.0
        puredmm[0].append(step)
        puredmm[1].append(float(i.split()[2]))

# Integrate
import scipy.integrate
#pure water and pure DMM
if len(puredmm) >1:
    print "pure DMM", scipy.integrate.trapz(puredmm[1], puredmm[0])
if len(purewater) >1:
    print "pure water", scipy.integrate.trapz(purewater[1], purewater[0])

#DMM in water
keys = dmm.keys()
keys.sort()
print "DMM in water"
for i in keys:
    print i, scipy.integrate.trapz(dmm[i][1], dmm[i][0])

#DMM in water
keys = dmm.keys()
keys.sort()
print "water in DMM"
for i in keys:
    print i, scipy.integrate.trapz(water[i][1], water[i][0])

#specified data
#print dmm['090'][0]
#print dmm['090'][1]
#print scipy.integrate.trapz(dmm['090'][1], dmm['090'][0])
