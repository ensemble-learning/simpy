#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, sys, cmd, math
import matplotlib.pyplot as plt

class LAMMPS_LOG_CDC(cmd.Cmd):
	"Deal with LAMMPS output"
	def __init__(self):
		cmd.Cmd.__init__(self)     #initialize the base class
		self.prompt = "(LAMMPS_LOG)>"
		self.intro = '''
#-------------------------------------------------------------------------
# LAMMPS_LOG Contains:             
# readlmp "name of logfile"                 | Print the column contained
# pt1line "indexnum from readlmp output"    | Plot the figure for each line
# aveall                                    | Calculated all average value
# EOF                                       | Exit from LAMMPS_LOG_CDC
# by liulc AT wag.caltech.edu
#-------------------------------------------------------------------------
'''

	def __average(self, data):
		n = len(data)
		ave = float(sum(data)) / n
		err = 0.0
		for i in data:
			err += i*i - ave*ave
		err = math.sqrt(err / (n-1))
		return (ave, err,)

	def help_EOF(self):
		print "End of the program!"
	def do_EOF(self, line):
		sys.exit()

	def help_ls(self):
		print "same with 'ls' in linux"
	def do_ls(self, line):
		os.system("ls")

	def help_cd(self):
		print "same with 'cd' in linux"
	def do_cd(self, folder):
		try:
			os.chdir(folder)
		except:
			print "folder doesn't exist? ", sys.exc_info()

	def help_pwd(self):
		print "same with 'pwd' in linux"
	def do_pwd(self, folder):
		print os.getcwd()

	def help_readlmp(self):
		print "Read the log.lammps and print the columns index"
		print "Usage: readlmp filename"
	def do_readlmp(self, lmplog):
		try:
			lmpfile = open(lmplog, 'r')
			block = None
			for i in lmpfile:
				if i.startswith("Loop"):
					block = None
				elif block is not None:
					for ncol in range(len(result)):
						result[ncol].append(i.split()[ncol])
				elif i.startswith("Step"):
					block = []
					result = [[] for j in i.split()]
					self.COLNAME = [j for j in i.split()]
			print "%10s %12s" %("Index", "ColName")
			for i in self.COLNAME:
				print "%10s %12s" %(self.COLNAME.index(i), i)
			self.DATA = result
			lmpfile.close()
		except:
			print "Some error in your log file. Empty file?", sys.exc_info()[0]

	def help_pt1line(self):
		print "Plot the figure on the x_index and y_index"
		print "Usage: pt1line 0 10"
	def do_pt1line(self, indextube):
		try:
			index_x = indextube.split()[0]
			index_y = indextube.split()[1]
			plt.figure(int(index_y))
			plt.plot(self.DATA[int(index_x)], self.DATA[int(index_y)])
			plt.show()
		except:
			print sys.exc_info()

	def help_aveall(self):
		print "Calculate all of the data's average value"
		print "Usage: aveall"
	def do_aveall(self, nothing):
		try:
			print "%10s %10s %14s %14s" %("Index", "ColName", "Average", "Std" )
			for i in self.COLNAME:
				temp_i = self.COLNAME.index(i)
				half_n = int(len(self.DATA[temp_i]) / 2)
				this_data = [float(j) for j in self.DATA[temp_i][(half_n + 1):]]
				ave_v = self.__average(this_data)[0]    #Average
				std_v = self.__average(this_data)[1]    #Std
				print "%10s %10s %14.4f %14.4f" %(temp_i, i, ave_v, std_v)
		except:
			print sys.exc_info()
	
if __name__ == '__main__':
	cdc = LAMMPS_LOG_CDC()
	cdc.cmdloop()
