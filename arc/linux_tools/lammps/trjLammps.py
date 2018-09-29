# -*- coding: utf-8 -*-
# !/usr/bin/env python

import os
from sys import argv

"""
This is a script to transfer LAMMPS trajectory file to fort.7 file of F-Reax.
Usage: python trj2fort7.py lammps.trj fort.7
coded by Lianchi Liu
Mar 21th 2011
"""

class TrjLammps:
	def read_trj(trjfile):
		'''
		This funtion return the list of all trajectories, result is a LIST with this structure:
		[{typedic}, [frame1], [frame2], ...] in which:
		{typedic} is the type dictionary along with atom number
		[frame1] = [[nstep, natoms, nbondpairs, total_e], [bond_pair1], [bond_pair2], ...]
		'''
		block = None; result = []; typedic = {}
		for i in trjfile:
			if i.find("number_of_atoms:") != -1:
				natoms = i.split()[-1]
			if i.find("molecular_analysis_frequency:") != -1:
				trjfile.next()
				for ii in range(int(natoms)):
					temp = trjfile.next().split()
					typedic[temp[0]] = temp[1]
				result.append(typedic)
			if i.find("step:") != -1:
				nstep = i.split()[-1]
				#Skip and grep useful information 
				trjfile.next();trjfile.next();trjfile.next();trjfile.next()
				trjfile.next();trjfile.next()
				total_energy = trjfile.next().split()[-1]
				total_kinetic = trjfile.next().split()[-1]
				total_potential = trjfile.next().split()[-1]
				bond_energy = trjfile.next().split()[-1]
				atom_energy = trjfile.next().split()[-1]
				lone_pair_energy = trjfile.next().split()[-1]
				valence_angle_energy = trjfile.next().split()[-1]
				threebody_conj = trjfile.next().split()[-1]
				hb_energy = trjfile.next().split()[-1]
				torsion_energy = trjfile.next().split()[-1]
				fourbody_conj = trjfile.next().split()[-1]
				vdw_energy = trjfile.next().split()[-1]
				elec_energy = trjfile.next().split()[-1]
				polar_energy = trjfile.next().split()[-1]
				n_bondpairs = int(trjfile.next().split(',')[-1])
				block = [[nstep, natoms, n_bondpairs, total_energy]]
			elif block is not None and i.startswith("chars_to_skip_frame_header"):
				result.append(block)
				block = None
			elif block is not None:
				block.append(i.split())
				#block.append([int(i.split()[0]), int(i.split()[0]), float(i.split()[0]), float(i.split()[3])])
		result.append(block)
		return result

class DealList:
	def link_list(bondlist, n, maxbonds): 
		'''
		bondlist is a List of bondorders in each list, which is started from the second data of frameList
		n is the atom number
		maxbonds is the length of the list
		'''
		atomlist = []
		bondorderlist = []
		#bondlengthlist = []
		for i in bondlist:
			if i[0:2].count(n) > 0:
				if i[0] == n and atomlist.count(i[1]) == 0:
					atomlist.append(i[1])
					#bondlengthlist.append(i[2])
					bondorderlist.append(i[3])
				elif i[1] == n and atomlist.count(i[0]) == 0:
					atomlist.append(i[0])
					#bondlengthlist.append(i[2])
					bondorderlist.append(i[3])
		while len(atomlist) < maxbonds:
			atomlist.append('0')
			bondorderlist.append('0.000')
		#print len(atomlist)
		atomlist.append('1')
		totalbondorder = sum(float(i) for i in bondorderlist)
		bondorderlist.append(totalbondorder)
		bondorderlist.append('0.000')
		bondorderlist.append('0.000')
		return [atomlist, bondorderlist]

if __name__ == '__main__':
	trjfile = open(argv[1])
	trjlist = TrjLammps.read_trj(trjfile)
	typedic = trjlist[0]
	natoms = len(trjlist[0])
	maxbonds = 12 #Max num of bonds for each atom
	fort7file = open(argv[2], 'w')
	#print trjlist
	for i in trjlist[1:]:
		fort7file.write("%5s SimulationName    Iteration:%10s  #Bonds:%5s\n" %(i[0][1], i[0][0], maxbonds))
		for j in range(natoms):
			templist = DealList.link_list(i[1:], str(j+1), maxbonds)
			fort7file.write("%5s%5s" %(str(j+1), int(typedic[str(j+1)])+1))
			for k in templist[0]:
				fort7file.write("%5s" %(k))
			for k in templist[1]:
				fort7file.write("%7s" %(k))
			fort7file.write("\n")
		fort7file.write("   1234.1234     1234.1234      1234.1234      1234.1234\n")
		
