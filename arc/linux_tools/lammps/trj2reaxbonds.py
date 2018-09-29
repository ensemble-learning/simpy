# !/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
from sys import argv

"""
This is a script to transfer LAMMPS trajectory file to fort.7 file of F-Reax.
Usage: python trj2reaxbonds.py lammps.trj fort.7
coded by Lianchi Liu
Mar 21th 2011
"""

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
			#Check if we need to skip the xyz sections
			check_npairs = int(trjfile.next().split(',')[-1])
			if check_npairs == int(natoms):
				for k in range(int(natoms)):
					trjfile.next()
				n_bondpairs = int(trjfile.next().split(',')[-1])
				block = [[nstep, natoms, n_bondpairs, total_energy]]
			elif check_npairs != int(natoms):
				n_bondpairs = check_npairs
				block = [[nstep, natoms, n_bondpairs, total_energy]]
		elif block is not None and i.startswith("chars_to_skip_frame_header"):
			result.append(block)
			block = None
		elif block is not None:
			block.append(i.split())
			#block.append([int(i.split()[0]), int(i.split()[0]), float(i.split()[0]), float(i.split()[3])])
	result.append(block)
	return result

def link_list(bondlist, natoms):
	'''
	link all the atom pairs to one list
	bondlist is a List of bondorders in each list, which starts from the second item of frameList
	natoms is the atom numbers
	'''
	atomlist = []
	bondorderlist = []
	for i in range(natoms):
		atomlist.append([])
		bondorderlist.append([])
	for j in bondlist:
		atomlist[int(j[0])-1].append(int(j[1]))
		bondorderlist[int(j[0])-1].append(float(j[3]))
		atomlist[int(j[1])-1].append(int(j[0]))
		bondorderlist[int(j[1])-1].append(float(j[3]))
	return [atomlist, bondorderlist]


def convert2graspbonds(): 
	#Convert the bonds trj of lammps to reax.bonds file for grasp
	trjfile = open(argv[1])
	trjlist = read_trj(trjfile)
	typedic = trjlist[0]
	natoms = len(trjlist[0])
	maxbonds = 12 #Max num of bonds for each atom
	fort7file = open(argv[2], 'w')
	for i in trjlist[1:]:
		#fort7file.write("%5s SimulationName    Iteration:%10s  #Bonds:%5s\n" %(i[0][1], i[0][0], maxbonds))
		fort7file.write("# Timestep %-10s\n" %(i[0][0]))
		fort7file.write("#\n")
		fort7file.write("# Number of particles %-10s\n" %(i[0][1]))
		fort7file.write("#\n")
		fort7file.write("# Max.number of bonds per atom 5 with coarse bond order cutoff 0.300\n")
		fort7file.write("# Particle connection table and bond orders\n")
		fort7file.write("# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q\n")
		atomlist = link_list(i[1:], natoms)
		for j in range(natoms):
			fort7file.write("%5s %5s " %(str(j+1), int(typedic[str(j+1)])+1))
			fort7file.write("%5s " %(len(atomlist[0][j])))
			for k in atomlist[0][j]:
				fort7file.write("%5s " %(k))
			fort7file.write("%5s " %(0)) # molecular number
			for k in atomlist[1][j]:
				fort7file.write("%7s " %(k))
			fort7file.write("%7s " %sum(atomlist[1][j]))
			fort7file.write("%7s %7s" %(0.000, 0.000)) # nlp and q
			fort7file.write("\n")
		fort7file.write("#\n")
		
def convert2fort7(): 
	#Convert the bonds trj of lammps to fort.7 
	trjfile = open(argv[1])
	trjlist = read_trj(trjfile)
	typedic = trjlist[0]
	natoms = len(trjlist[0])
	maxbonds = 12 #Max num of bonds for each atom
	fort7file = open(argv[2], 'w')
	for i in trjlist[1:]:
		fort7file.write("%5s SimulationName    Iteration:%10s  #Bonds:%5s\n" %(i[0][1], i[0][0], maxbonds))
		atomlist = link_list(i[1:], natoms)
		for j in range(natoms):
			fort7file.write("%5s %5s " %(str(j+1), int(typedic[str(j+1)])+1))
			while len(atomlist[0][j]) < maxbonds:
				atomlist[0][j].append(0)
				atomlist[1][j].append(0.000)
			for k in atomlist[0][j]:
				fort7file.write("%5s " %(k))
			fort7file.write("%5s " %(0)) # molecular number
			for k in atomlist[1][j]:
				fort7file.write("%7s " %(k))
			fort7file.write("%7s " %sum(atomlist[1][j]))
			fort7file.write("%5s %5s" %(0.0, 0.0)) # nlp and q
			fort7file.write("\n")
		fort7file.write("   1234.1234     1234.1234      1234.1234      1234.1234\n")

if __name__ == '__main__':

	if len(argv) != 3:
		print "################################################"
		print "Usage: python trj2reaxbonds.py lammps.trj output"
		print "coded by Lianchi"
		print "################################################"
		sys.exit()
		
	choice = raw_input("Input the file type you want to convert to: GraspBonds(1) fort.7(2)\n")
	if int(choice) == 1:
		convert2graspbonds()
	elif int(choice) == 2:
		convert2fort7()
	
