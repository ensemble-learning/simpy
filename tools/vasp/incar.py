INCAR = """
 System = vasp
##################################
#  Start parameter for this run  #
##################################
 ISTART   = 0
 ICHARG   = 2

##################################
#     Electronic Relaxation      #
##################################
 PREC     = A
 ENCUT    = 400
 NELMIN   = 5
 EDIFF    = 1.0E-6
# IALGO = 48
# ALGO = Very Fast

##################################
# Ionic Relaxation
##################################
 EDIFFG   = -1.0E-3
 NSW      = 500
 IBRION   = 2
 ISIF     = 3

##################################
#General Params
##################################
 ISMEAR   = 0
 SIGMA    = 0.1

##################################
#   Output control               #
##################################
 LWAVE  = .FALSE.
 LCHARG = .FALSE.
 NWRITE = 0

##################################
# Other parameters
##################################
 LVDW = 12      # PBE_ULG
 NPAR = 2

"""

o = open("INCAR", "w")

lines = INCAR
o.write(lines)
o.close()
