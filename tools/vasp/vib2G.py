#!/usr/bin/env python

###############################################################################
# Python script to calculate free energy contribution (ZPE included)
# of vibrational modes
# Written by Hai Xiao <haixiao@wag.caltech.edu>
###############################################################################

import sys
from math import exp,log
#import os

if len(sys.argv) == 4:
  T = float(sys.argv[1])
  outcar = open(sys.argv[2])
  state = int(sys.argv[3])

  # Get vibrational eigenvalues from OUTCAR
  vibmatch  = "f  ="
  ivibmatch = "f/i="
  Evib = []
  Nivib = 0
  Eivib = []
  for line in outcar:
    if vibmatch in line:
      # Read in vibrational mode in unit of meV (hv, NOT 1/2*hv)
      Evib.append(float(line.split()[9]))
    if ivibmatch in line:
      Eivib.append(float(line.split()[8]))
      Nivib = Nivib + 1
  outcar.close()
  ###for i in Evib:
  ###  print "%10.6f"%(i)

  # Calculate characteristic vibrational temperatures
  # conversion factor from http://physics.nist.gov/cgi-bin/cuu/Value?evk
  # using the source: 2014 CODATA recommended values
  # T = hv/k
  meV2K = 11.6045221
  Tvib = []
  for i in Evib:
    Tvib.append(i*meV2K)
    print i
  ###for i in Tvib:
  ###  print "%12.6f"%(i)

  # Set the threshold for low frequency modes, which might lead to unphysically large entropy contribution
  # Currently use the peak frequency 60 cm-1 of ice Ih spectra, corresponding to
  # the acoustic translational mode of the six-member rings
  # Planck constant in eV*s from http://physics.nist.gov/cgi-bin/cuu/Value?hev
  # using the source: 2014 CODATA recommended values
  h = float("4.135667662E-15")
  # speed of light in vacuum in unit of cm/s (exact)
  c = float("2.99792458E10")
  # conversion factor for wavenumber 1/lambda (cm-1): h*c
  cm2eV = h*c
  ###print "%10.6E"%(cm2eV)
  threshold = 60.0*cm2eV*1000.0*meV2K
  ###print "%12.6f"%(threshold)
  for i in range(len(Tvib)):
    if Tvib[i] < threshold:
      print "Warning: very small frequency found, which might lead to unphysically large entropy contribution"
      print "very small frequency %3i: %12.6f K"%(i+1,Tvib[i])
      Tvib[i] = threshold
      print "         treated as real: %12.6f K"%(Tvib[i])
  ###for i in Tvib:
  ###  print "%12.6f"%(i)

  # Important check of small imaginary frequencies
  if state == 0:
    if Nivib > 0:
      print "Warning: imaginary frequencies found for supposed local minimum !!!"
      print "Warning: all imaginary frequencies treated as threshold ones for now"
      for i in range(Nivib):
        print "imaginary frequency %i: %12.6f K"%(i+1,Eivib[i]*meV2K)
	Tvib.append(threshold)
        print "       treated as real: %12.6f K"%(threshold)
  elif state == 1:
    if Nivib == 0:
      print "Error: no imaginary frequency found for supposed transition state !!!"
      sys.exit(0)
    elif Nivib > 1:
      print "Warning: more than 1 imaginary frequency found for supposed transition state !!!"
      print "Warning: all imaginary frequencies, except the lowest, treated as threshold ones for now"
      for i in range(Nivib-1):
        print "imaginary frequency %i: %12.6f K"%(i+1,Eivib[i]*meV2K)
	Tvib.append(threshold)
        print "       treated as real: %12.6f K"%(threshold)
      print "imaginary frequency %i: %12.6f K"%(Nivib,Eivib[Nivib-1]*meV2K)
      print "treated as the only 1 imaginary frequency for supposed transition state"
  else:
    print "Error: only accept state 0 or 1 for now"
    sys.exit(0)

  # Free energy contribution, using formulas from
  # http://www.gaussian.com/g_whitepap/thermo.htm
  # R = N_A * k
  # use k here, in order to get final energies in the unit of eV
  # Boltzmann constant from http://physics.nist.gov/cgi-bin/cuu/Value?tkev
  # using the source: 2014 CODATA recommended values
  k = float("8.6173303E-5")
  # ZPE
  H_zpe = 0.0
  for i in Tvib:
    H_zpe = H_zpe + k * i * 0.5
  ###print "%10.6f"%(H_zpe)
  ###H_zpe_test = 0.0
  ###for i in Evib:
  ###  H_zpe_test = H_zpe_test + 0.5 * i / 1000.0
  ###print "%10.6f"%(H_zpe_test)
  # enthalpy and entropy contributions
  H_vib = 0.0
  S_vib = 0.0
  if T == 0.0:
    H_vib = 0.0
    S_vib = 0.0
  else:
    for i in Tvib:
      q = i/T
      H_vib = H_vib + k * i / (exp(q) - 1.0)
      S_vib = S_vib + k * (q / (exp(q) - 1.0) - log(1.0 - exp(-q)))
      print i, q, k * (q / (exp(q) - 1.0) - log(1.0 - exp(-q)))
  ###print "%13.6E"%(H_vib)
  ###print "%13.6E"%(S_vib)
  G_vib = H_vib - T * S_vib
  total_vib = H_zpe + G_vib
  print "ZPE contribution is %10.6f eV"%(H_zpe)
  print "Enthalpy contribution H at %.1f K is %10.6f eV"%(T,H_vib)
  print "Entropy S at %.1f K is %10.6f eV/K"%(T,S_vib)
  print "H - TS at %.1f K is %10.6f eV"%(T,G_vib)
  print "Total contribution at %.1f K with ZPE is %10.6f eV"%(T,total_vib)
else: 
  # Print usage
  print ""
  print "Usage:"
  print "vib2G.py [T(K)] [OUTCAR] [0 or 1 for LM or TS]"
  print ""
  sys.exit(0)

