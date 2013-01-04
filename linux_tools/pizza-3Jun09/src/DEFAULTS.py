# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# --------------
# --------------

# 3 variables used by Pizza.py
# to use TOOLS or SCRIPTS, edit and uncomment the line
# to use EXCLUDE, add to the existing list

# TOOLS = list of extra directories that contain Pizza.py tools
#   Pizza.py will load all *.py files as tools from TOOLS dirs, then pizza/src
#   this ordering means your tool can override a Pizza.py tool of the same name
# SCRIPTS = list of extra directories that contain Pizza.py scripts
#   Pizza.py will look in working dir, SCRIPTS dirs, then pizza/scripts
#   this ordering means your script can override a Pizza.py script of same name
# EXCLUDE = Python files to NOT load as tools when Pizza.py starts
#   typically done for auxiliary Python files that are not tools
#   any non-tool Python files from your TOOLS dirs should be added to list

#PIZZA_TOOLS = ["~/mystuff/new_pizza_tools"]
#PIZZA_SCRIPTS = ["~/mystuff/new_pizza_scripts"]
PIZZA_EXCLUDE = ["pizza", "DEFAULTS", "vizinfo"]

# --------------
# --------------

# Pathname for programs executed by various Pizza.py tools

# if you don't use a tool, it's settings can be ignored
# the default values are program names with no path
# to use a default value, the executable must therefore be in your path
# to change a default, edit and uncomment the PIZZA variable line

# --------------

# ImageMagick programs to manipulate image files
# DISPLAY = program to view GIF, PNG, SVG files
# tools that use it: rasmol, raster, svg
# CONVERT = program to convert one image format to another
# MONTAGE = program to stitch 2 images together
# tools that use it: image
# default = display

#PIZZA_DISPLAY = "/usr/bin/display"
#PIZZA_CONVERT = "/usr/bin/convert"
#PIZZA_MONTAGE = "/usr/bin/montage"

# --------------

# GNUPLOT = the GnuPlot plotting package
# GNUTERM = terminal setting used by GnuPlot
# tools that use it: gnu
# default = gnuplot
# default = x11

#PIZZA_GNUPLOT = "gnuplot"
#PIZZA_GNUTERM = "x11"
#PIZZA_GNUTERM = "aqua"                   # for Macs with Aquaterm installed

# --------------

# GUNZIP = program to uncompress gzipped files
# tools that use it: data dump log
# default = gunzip

#PIZZA_GUNZIP = "gunzip"

# --------------

# LABEL3D = program to put a label on a Raster3D image
# tools that use it: raster
# default = label3d

#PIZZA_LABEL3D = "label3d"

# --------------

# MATLAB = the MatLab numerical analysis and plotting package
# tools that use it: matlab
# default = matlab -nosplash -nodesktop -nojvm

#PIZZA_MATLAB = "matlab -nosplash -nodesktop -nojvm"

# --------------

# RASMOL = the RasMol visualization package
# tools that use it: rasmol
# default = rasmol

#PIZZA_RASMOL = "rasmol"

# --------------

# RENDER = the Raster3D visualization rendering engine
# tools that use it: raster
# default = render

#PIZZA_RENDER = "render"
