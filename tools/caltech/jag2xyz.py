#!/usr/bin/env python
#
"""Converts a Jaguar input or output file to a xyz file
   Usage: run xyz2xyz.py inputFile and follow instructions interactivly.
   06/23/2008: Andres Jaramillo-Botero 
"""
import os,sys,string
import datetime, time

def splitNum(string):
    str = string
    nums = ['0','1','2','3','4','5','6','7','8','9']
    for num in nums:
        chars = str.split(num)
        char = chars[0]
        if not char:
            break
        elif char != str:
            str = char
            continue
    return char

def jag2xyz(filename):
    """
        Convert Jaguar input or output (converged optimization) file to xyz file (*.xyz).
    """

    basename = ''
    if filename.find('.') != -1:
        parts = filename.split('.')
        for part in parts[:-1]:
            basename += part + '.' 
        basename = basename[:-1]
    else:
        basename = filename

    jagFile = open(filename, 'r')
    input = False
    output = False
    zmatLine = False
    elems = []
    coords = []
    count = 0
    cont=0 	# geometry counter
    geometry_coords={}
    for line in jagFile.readlines():
        if line.find('&zmat') == 0:
            input = True
            zmatLine = True
            continue
        if line.find('  final geometry') == 0:
            output = True
            continue
        if input:
            if zmatLine and line.find('&') != 0:
                if len(line) < 4:
                    sys.exit('zmat line not complete. Abort')
                else:
                    fields = line.split()
                    elems.append(fields[0])
                    for i in range(1,4):
                        if fields[i].find('#') != -1:
                            # Remove the '#'.
                            fields[i] = fields[i].split('#')[0]
                        else:
                            continue
                    coords.append([float(fields[1]),float(fields[2]),float(fields[3])])
                    continue
            if line.find('&') == 0:
                zmatLine = False
                continue
        elif output:
            if count < 2:
                count += 1
                continue
            elif count >= 2:
                fields = line.split()
                if len(fields) != 4:
                    geometry_coords[cont]=coords
                    count=0
                    output=False
                    coords=[]
                    cont += 1
                    continue            
                else:
                    if cont==0: elems.append(fields[0]) # only need elem list for first structure
                    coords.append([float(fields[1]),float(fields[2]),float(fields[3])])
                    continue
        else:
            continue

    if not elems or not geometry_coords:
        sys.exit('Error in reading Jaguar file. Abort')

    numAtoms = len(elems)
    fftype = []
    atomtype = []
    res=[]
    group=[]
    group1=[]
    group2=[]
    for elem in elems:
        char = splitNum(elem) 
        atomtype.append(char)
        fftype.append(char)

    print atomtype
    newlines_start = []
    newlines_end = []
    charge={}
    
    # Need to change DESCRP to short file name instead of tmp
    chunk = '%d\n' % (numAtoms)
    newlines_start.append(chunk)
    chunk = '%s generated from Jaguar output file ...\n' % filename[:-4] 
    newlines_start.append(chunk)

    for geometry in range(len(geometry_coords)):
        newlines = []
        newlines.append(string.join(newlines_start,''))
        for i in range(numAtoms):
            atom_data = (fftype[i], geometry_coords[geometry][i][0], geometry_coords[geometry][i][1], geometry_coords[geometry][i][2])
            format_string = '%-5s %10.5f%10.5f%10.5f\n'
            chunk = format_string % atom_data
            newlines.append(chunk)

        current_dir = os.getcwd()
        xyzname = basename+'.'+str(geometry)+'.xyz'
        full_xyzname = current_dir+'/'+xyzname
        f = open(full_xyzname,'w')
        f.writelines(newlines)
        f.close()
        
        print (' Generated %s !' % (xyzname))

    print "Done ...\n"
    
    return

if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        inputFile = raw_input("Please provide an input file name: ")
    else:
        try:
            inputFile = sys.argv[1]
        except:
            print """Usage: jag2xyz.py <inputfile>
                    - <inputfile> can be either Jaguar input or output file
                  """
            sys.exit(1)
    jag2xyz(inputFile)

