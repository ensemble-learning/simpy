#!/usr/bin/env python
#-*- coding: utf-8 -*-

def read_lmplog(lmplog):
	'''
	Read the log.lammps and dump all the data to the list of [result]
	Structure of [result]:   [[column 1], [column 2], ..., [column n]]
	Structure of [column n]: ['Name', 'Value 1', ..., 'Value n']
	'''
	block = None
	for i in lmplog:
		if i.startswith("Loop"):
			block = None
		elif block is not None:
			for ncol in range(len(result)):
				result[ncol].append(i.split()[ncol])
		elif i.startswith("Step"):
			block = []
			result = [[j] for j in i.split()]
	return result

