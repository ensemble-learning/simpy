import csv
import numpy as np
import matplotlib.pyplot as plt

x = []
y = []

fig = plt.figure()
with open('Bond_p1.csv', 'rb') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    n = 0
    for row in data:
        if n > 0:
            xv = float(row[1])
            yv = float(row[2])
            x.append(xv)
            y.append(yv)
        else:
            pass
        n += 1
xp = np.array(x)
yp = np.array(y)
#min = np.min([np.min(xp), np.min(yp)])
#max = np.max([np.max(xp), np.max(yp)])
ax = fig.add_subplot(221)
ax.plot([0.5, 1.9], [0.5, 1.9], ls="--")
ax.plot(xp, yp, "o")
ax.set_title("Bond length")
ax.set_xlabel("QM ($\AA$)")
ax.set_ylabel("FF ($\AA$)")

x = []
y = []
with open('Angle_p1.csv', 'rb') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    n = 0
    for row in data:
        if n > 0:
            xv = float(row[1])
            yv = float(row[2])
            x.append(xv)
            y.append(yv)
        else:
            pass
        n += 1
xp = np.array(x)
yp = np.array(y)
min = np.min([np.min(xp), np.min(yp)])
max = np.max([np.max(xp), np.max(yp)])
ax = fig.add_subplot(222)
ax.plot([min, max], [min, max], ls="--")
ax.plot(xp, yp, "o")
ax.set_title("Angle")
ax.set_xlabel("QM (degree)")
ax.set_ylabel("FF (degree)")

x = []
y = []
with open('Torsion_p1.csv', 'rb') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    n = 0
    for row in data:
        if n > 0:
            xv = float(row[1])
            yv = float(row[2])
            x.append(xv)
            y.append(yv)
        else:
            pass
        n += 1
xp = np.array(x)
yp = np.array(y)
min = np.min([np.min(xp), np.min(yp)])
max = np.max([np.max(xp), np.max(yp)])
print min, max
ax = fig.add_subplot(223)
ax.plot([min, max], [min, max], ls="--")
ax.plot(xp, yp, "o")
ax.set_title("Torsion Angle")
ax.set_xlabel("QM (degree)")
ax.set_ylabel("FF (degree)")

fig.subplots_adjust(wspace=0.4, hspace=0.4)
plt.show()
        