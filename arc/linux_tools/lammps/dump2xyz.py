# -*- coding: utf-8 -*-
# !/usr/bin/env python

import os, sys
from sys import argv
"""
Usage: python dump2xyz.py dump.xyz xyzfile.xyz.
Note: data_lammps must be in the same directory
by Lianchi Liu
Mar 20th. 2011
"""

def read_dumpfile(dumpfile):
	'''
	read the dumpfile.lmp to dumplist
	dumplist structure:	["Timestep xxx", "line 1", "line2", ..., "linen"]
	'''
	dumplist = []
	for i in dumpfile:
		if i.count("TIMESTEP") != 0:
			temp_timestep = dumpfile.next()
			dumpfile.next()
			dumplist.append(dumpfile.next())
			dumplist.append("Timestep "+temp_timestep)
			dumpfile.next()
			boxinf = [dumpfile.next(), dumpfile.next(), dumpfile.next()]
			dumpfile.next()
		else:
			dumplist.append(i)
	return dumplist

def get_typedic(datafile):
	typedic = {}
	for i in datafile:
		if len(i.split()) > 2 and i.split()[2] == "#":
			typedic[i.split()[0]] = i.split()[3]
	return typedic

if __name__ == '__main__':
	
	if len(argv) < 2:
		print "Usage: python dump2xyz.py dump.xyz xyzfile.xyz."
		print "Note: data_lammps(exactly the same name) must be in the same directory."
		sys.exit()

	datafile = open("data_lammps")
	typedic = get_typedic(datafile)
	dumpfile = open(argv[1])
	newfile = open(argv[2], 'w')
	#print read_dumpfile(dumpfile)
	for i in read_dumpfile(dumpfile):
		if len(i.split()) <= 2:
			newfile.write(i)
		elif len(i.split()) >= 4:
			newfile.write(typedic[i.split()[1]] + " " + \
							i.split()[2] + " " + \
							i.split()[3] + " " + \
							i.split()[4] + "\n")
	datafile.close()
	dumpfile.close()
	newfile.close()
