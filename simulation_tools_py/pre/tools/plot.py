#!/usr/bin/env python

"""USAGE: plot the csv file
FORMAT: plot.py [-type][csvfile][...]"""

import sys
import random
import matplotlib.pyplot as pl

COLOR = { 0: 'b', 1: 'g', 2: 'r', 3: 'c', 4: 'm', 5: 'y', 6: 'k', 7: 'w'}
FORM = {0: '.', 1: '.', 2: '.', 3: '.', 4: '.', 5: '.', 6: '.'}

def parse_out(filename):
    x = []
    y = []
    f = open(filename, 'r')
    counter = 0
    for i in f:
        tokens = i.strip().split()
        if counter == 0:
            pass
        else:
            x.append(int(tokens[0]))
            y.append(float(tokens[1]))
        counter += 1
    return x, y
            
        

def parse_csv(filename):
    x = []
    y = []
    f = open(filename, "r")
    for i in f:
        tokens = i.strip().split(",")
        x.append(float(tokens[0])/4000)
        y.append(int(tokens[1]))
    return x, y

def plot_data(type):
    pl_x = []
    pl_y = []
    label_pl = []
    if len(sys.argv) < 3:
        print "warning"
    else:
        for i in sys.argv[2:]:
            x = []
            y = []
            if type == 1:
                x, y = parse_csv(i)
            elif type == 2:
                x, y = parse_out(i)
            label_pl.append(i)
            pl_x.append(x)
            pl_y.append(y)

    for i in range(len(pl_x)):
        c = COLOR[i%7]
        f = FORM[i%6]
        pl.plot(pl_x[i], pl_y[i], c+f , label=label_pl[i], lw=3)
    #axis = pl.gca().xaxis
    #for label in axis.get_ticklabels():
    #    label.set_rotation(45)
    pl.legend(loc = 0)
    pl.xlabel('Tims ps')
    pl.ylabel('N')
    pl.grid()
    pl.show()

def parse_argv():
    type = 0
    if len(sys.argv) < 3:
        print "warning"
    else:
        if sys.argv[1] == '-csv':
            type = 1
        elif sys.argv[1] == '-out':
            type = 2
        plot_data(type)

if __name__ == "__main__":
    parse_argv()



    
