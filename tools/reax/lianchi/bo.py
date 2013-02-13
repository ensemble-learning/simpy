#!/home/liulc/src/epd-7.1/epd-7.1-2-rh5-x86/bin/python
#!-*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import os, sys, ffield
from math import exp, pow, sqrt

class FF_PARA:

	def __init__(self, atom_a, atom_b, ff):
		self.cal_bo_pair(atom_a, atom_b, ff)
		self.cal_inner_pair(atom_a, atom_b, ff)
	
	def __com_r(self, a, b):
		return sqrt(a * b)
	
	def cal_bo_pair(self, atom_a, atom_b, ff):
		#Get r_sigma, r_pi, r_pipi
		if ff.offpara.keys().count(atom_a+"-"+atom_b) == 1:
			para = ff.offpara[atom_a+"-"+atom_b]
			r_sigma, r_pi, r_pipi = para[3:6]
		elif ff.offpara.keys().count(atom_b+"-"+atom_a) == 1:
			para = ff.offpara[atom_b+"-"+atom_a]
			r_sigma, r_pi, r_pipi = para[3:6]
		else:
			r_sigma = sqrt(ff.atompara[atom_a][0] * ff.atompara[atom_a][0])
			r_pi = sqrt(ff.atompara[atom_a][6] * ff.atompara[atom_a][6])
			r_pipi = sqrt(ff.atompara[atom_a][16] * ff.atompara[atom_a][16])
		self.r_sigma, self.r_pi, self.r_pipi = r_sigma, r_pi, r_pipi
		#Get bo para
		if ff.bondpara.keys().count(atom_a+"-"+atom_b) == 1:
			para = ff.bondpara[atom_a+"-"+atom_b]
			de_sigma, de_pi, de_pipi = para[0], para[1], para[2]
			p_bo1, p_bo2 = para[12], para[13]
			p_bo3, p_bo4 = para[9], para[10]
			p_bo5, p_bo6 = para[4], para[6]
		elif ff.bondpara.keys().count(atom_b+"-"+atom_a) == 1:
			para = ff.bondpara[atom_b+"-"+atom_a]
			de_sigma, de_pi, de_pipi = para[0], para[1], para[2]
			p_bo1, p_bo2 = para[12], para[13]
			p_bo3, p_bo4 = para[9], para[10]
			p_bo5, p_bo6 = para[4], para[6]
		self.p_bo1, self.p_bo2 = p_bo1, p_bo2
		self.p_bo3, self.p_bo4 = p_bo3, p_bo4
		self.p_bo5, self.p_bo6 = p_bo5, p_bo6
	
	def cal_inner_pair(self, atom_a, atom_b, ff):
		#Get inner-wall
		r_inner = sqrt(ff.atompara[atom_a][29] * ff.atompara[atom_b][29])
		d_inner = sqrt(ff.atompara[atom_a][30] * ff.atompara[atom_b][30])
		alpha_inner = sqrt(ff.atompara[atom_a][31] * ff.atompara[atom_b][31])
		self.r_inner = r_inner
		self.d_inner = d_inner
		self.alpha_inner = alpha_inner

if __name__ == '__main__':
	ff = ffield.FFIELD("ffield")
	data = FF_PARA('C', 'C', ff)
	print data.p_bo3
	
"""
def test():
	#----    cut     ----
	# Get Bond oder list
	this_r = [float(i)*0.01 for i in range(50, 400)]
	def bo_sigma(r):
		return exp(p_bo1 * pow(r/r_sigma, p_bo2))

	def bo_pi(r):
		return exp(p_bo3 * pow(r/r_pi, p_bo4))

	def bo_pipi(r):
		return exp(p_bo5 * pow(r/r_pipi, p_bo6))

	bo1, bo2, bo3 = [], [], []
	allbo = []
	for i in this_r:
		bo1.append(bo_sigma(i))
		temp_bo1 = bo_sigma(i)
		if r_pi > 0:
			bo2.append(bo_pi(i))
			temp_bo2 = bo_pi(i)
		else:
			bo2.append(0)
			temp_bo2 = 0
		if r_pipi > 0:
			bo3.append(bo_pipi(i))
			temp_bo3 = bo_pipi(i)
		else:
			bo3.append(0)
			temp_bo3 = 0
		allbo.append(temp_bo1 + temp_bo2 + temp_bo3)

	#----    cut    ----
	# Get Bond Energy
	def get_ebond(bo, de=0.0):
		return -de*bo*exp(p_be1*(1-pow(bo, p_be2)))

	ebond_sigma, ebond_pi, ebond_pipi = [], [], []
	ebond_total = []
	for i in bo1:
		ebond_sigma.append(get_ebond(i, p_De_sigma))

	for i in bo2:
		ebond_pi.append(get_ebond(i, p_De_pi))

	for i in bo3:
		ebond_pipi.append(get_ebond(i, p_De_pipi))

	for i in range(len(ebond_sigma)):
		ebond_total.append(ebond_sigma[i] + ebond_pi[i] + ebond_pipi[i])

	#----    cut    ----
	# Get VdW Energy
	evdw = []
	def f13(r):
		return pow(r+pow(1.0/10.78, 1.5591),1/1.5591)

	def get_evdw(r):
		return p_De_vdw*(exp(p_alpha_vdw*(1-f13(r)/r_vdw))-2.0*exp(0.5*p_alpha_vdw*(1-f13(r)/r_vdw)))

	for i in this_r:
		evdw.append(get_evdw(i))

	#----    cut    ----
	# Get Total Energy 
	etotal = []
	for i in range(len(ebond_sigma)):
		etotal.append(evdw[i] + ebond_total[i])

	#----    cut    ----
	# Plot figures
	plt.figure(1)
	plt.plot(this_r, bo1, lw = 2, label="bo_sigma")
	plt.plot(this_r, bo2, lw = 2, label="bo_pi")
	plt.plot(this_r, bo3, lw = 2, label="bo_pipi")
	plt.plot(this_r, allbo, lw = 2, label="Total bo")
	plt.grid()
	plt.legend()
	plt.figure(2)
	plt.plot(this_r, ebond_sigma, lw=2, label="be_sigma")
	plt.plot(this_r, ebond_pi, lw=2, label="be_pi")
	plt.plot(this_r, ebond_pipi, lw=2, label="be_pipi")
	plt.plot(this_r, ebond_total, lw=2, label="be_total")
	plt.plot(this_r, evdw, lw=2, label="e_vdw")
	plt.plot(this_r, etotal, lw=2, label="e_total")
	plt.ylim([-300, 300])
	plt.grid()
	plt.legend()
	plt.show()
"""
